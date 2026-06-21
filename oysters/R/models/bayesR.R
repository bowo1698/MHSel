# BayesR Model with Mixture Prior

run_bayesR <- function(matrices, split, gblup_varcomp, config, fold = 0,
                       map = NULL, windsize = NULL) {

  W_train  <- matrices$W_train
  W_test   <- matrices$W_test

  pheno_train_df <- split$train$phenotypes
  pheno_test_df  <- split$test$phenotypes
  has_test       <- nrow(pheno_test_df) > 0

  y <- pheno_train_df$phenotype

  fe_col <- config$fixed_effect
  if (!is.null(fe_col) && fe_col %in% names(pheno_train_df)) {
    fe_levels  <- levels(as.factor(pheno_train_df[[fe_col]]))
    X_cov      <- as.matrix(as.numeric(factor(pheno_train_df[[fe_col]],
                                              levels = fe_levels)))
    storage.mode(X_cov) <- "double"
    if (has_test && fe_col %in% names(pheno_test_df)) {
      X_test_cov <- as.matrix(as.numeric(factor(pheno_test_df[[fe_col]],
                                                levels = fe_levels)))
      storage.mode(X_test_cov) <- "double"
    } else {
      X_test_cov <- NULL
    }
  } else {
    X_cov      <- NULL
    X_test_cov <- NULL
  }

  pi_vec        <- config$bayesR$pi
  sigma2_ah     <- gblup_varcomp$sigma2_g
  sigma2_e_init <- gblup_varcomp$sigma2_e

  prior_params <- list(
    a0_e           = config$bayesR$prior_df$residual,
    a0_g           = config$bayesR$prior_df$genetic,
    variance_class = config$bayesR$var_class
  )

  cat("BayesR Settings:\n")
  cat("π =", pi_vec, "\n")
  cat("σ²_αh =", round(sigma2_ah, 6), "| σ²_e =", round(sigma2_e_init, 6), "\n")
  cat("a0_e =", prior_params$a0_e, "| a0_g =", prior_params$a0_g, "\n")
  cat("var_class =", prior_params$variance_class, "\n")
  cat("fixed effect:", if (is.null(X_cov)) "none" else
                      sprintf("%s (q=%d)", fe_col, ncol(X_cov)), "\n")
  cat("fold:", fold, "\n\n")

  cat("Running BayesR with MCMC algorithm...\n\n")

  res <- masbayes::run_bayesr(
    w             = W_train,
    y             = y,
    X             = X_cov,
    marker_type   = config$marker_type %||% "auto",
    pi_vec        = pi_vec,
    sigma2_e_init = sigma2_e_init,
    sigma2_ah     = sigma2_ah,
    prior_params  = prior_params,
    mcmc_params   = list(
      n_iter = as.integer(config$n_iter),
      n_burn = as.integer(config$n_burn),
      n_thin = as.integer(config$n_thin),
      seed   = as.integer(config$seed)
    ),
    method        = "mcmc",
    response_type = if (!is.null(config$response_type)) config$response_type else "gaussian",
    fold_id       = as.integer(fold),
    save_rds      = FALSE,
    # GWAS surface (masbayes v0.5.0): NULL keeps fit bitwise-identical to v1.4.0.
    # When supplied, masbayes attaches pip / pip_block / gwas / gwas_meta to res.
    map           = map,
    windsize      = windsize
  )

  # Posterior point estimates (raw fields preserved for marker_info etc.)
  beta_hat       <- res$beta_hat
  sigma2_e_hat   <- res$sigma2_e_hat
  alpha_hat      <- res$alpha_hat
  sigma2_g_hat   <- res$sigma2_g
  h2_hat         <- res$h2

  # Predictions via masbayes::predict() — same formula, single source of truth
  pred_train <- as.vector(predict(res)$GEBV)

  if (has_test) {
    pred_test_obj <- predict(res, W_test, pheno_test_df$phenotype,
                             X_new = X_test_cov)
    pred_test <- as.vector(pred_test_obj$GEBV)
  } else {
    pred_test_obj <- NULL
    pred_test     <- numeric(0)
  }

  # Variance components from masbayes (posterior means + π)
  vc <- res$variance_components
  sigma2_small_hat  <- vc$small$mean
  sigma2_medium_hat <- vc$medium$mean
  sigma2_large_hat  <- vc$large$mean
  pi_hat <- c(vc$pi$zero, vc$pi$small, vc$pi$medium, vc$pi$large)

  # Mode gamma per marker — keep manual (not exposed by masbayes)
  gamma_mode <- apply(res$gamma_samples, 2, function(x) {
    tab <- table(x)
    if (length(tab) == 0) return(0)
    as.numeric(names(sort(tab, decreasing = TRUE)[1]))
  })

  # Per-marker 4-class PIP from MCMC class labels (0/1/2/3 = null/small/medium/
  # large). Eq 3 MAP = which.max(row); Eq 4 non-null = 1 - null > 0.5.
  pip_class <- vapply(0:3, function(k) colMeans(res$gamma_samples == k),
                      numeric(ncol(res$gamma_samples)))
  colnames(pip_class) <- c("null", "small", "medium", "large")

  # ESS / Geweke from masbayes diagnostics
  d <- res$diagnostics
  ess <- list(
    sigma2_e      = d$sigma2_e$ess,
    sigma2_small  = d$sigma2_small$ess,
    sigma2_medium = d$sigma2_medium$ess,
    sigma2_large  = d$sigma2_large$ess
  )
  geweke <- list(
    sigma2_e      = d$sigma2_e$geweke,
    sigma2_small  = d$sigma2_small$geweke,
    sigma2_medium = d$sigma2_medium$geweke,
    sigma2_large  = d$sigma2_large$geweke
  )

  cat("Posterior Means:\n")
  cat("σ²_e:", round(sigma2_e_hat, 6), "\n")
  cat("σ²_small:", round(sigma2_small_hat, 8), "\n")
  cat("σ²_medium:", round(sigma2_medium_hat, 8), "\n")
  cat("σ²_large:", round(sigma2_large_hat, 8), "\n")
  cat("π:", round(pi_hat, 3), "\n")
  cat("σ²_g (genetic):", round(sigma2_g_hat, 6), "\n")
  cat("h² (heritability):", round(h2_hat, 4), "\n")
  if (!is.null(alpha_hat))
    cat("α (FE):", round(as.numeric(alpha_hat), 4), "\n")
  cat("\n")

  # Train / test metrics from masbayes — always observed scale for binary
  # (masbayes::compute_metrics_binary applies pnorm before scoring so accuracy
  # is Pearson r on probability scale, R2 = r^2, RMSE^2 ≈ Brier).
  # AUC stored separately so it does not collide with the accuracy field.
  is_binary <- !is.null(config$response_type) && config$response_type == "binary"
  tm <- res$training_metrics
  accuracy_train_pheno <- tm$accuracy
  rmse_train_pheno     <- tm$RMSE
  r2_train_pheno       <- tm$R2
  b_train_pheno        <- tm$bias
  auc_train            <- if (is_binary) tm$AUC else NA_real_

  if (has_test) {
    em <- pred_test_obj$metrics
    accuracy_test_pheno <- em$accuracy
    rmse_test_pheno     <- em$RMSE
    r2_test_pheno       <- em$R2
    b_test_pheno        <- em$bias
    auc_test            <- if (is_binary) em$AUC else NA_real_
  } else {
    accuracy_test_pheno <- NA_real_
    rmse_test_pheno     <- NA_real_
    r2_test_pheno       <- NA_real_
    b_test_pheno        <- NA_real_
    auc_test            <- NA_real_
  }

  # Test heritability (binary uses liability scale; gaussian uses observed)
  if (is_binary) {
    if (has_test) {
      P <- mean(pheno_train_df$phenotype, na.rm = TRUE)
      z <- dnorm(qnorm(P))
      h2_obs <- var(pred_test) / (var(pred_test) + sigma2_e_hat)
      h2_test_bayesr <- h2_obs * P * (1 - P) / z^2
      accuracy_test_gebv <-
        if (!is.na(accuracy_test_pheno) && !is.na(h2_test_bayesr) &&
            h2_test_bayesr > 0)
          accuracy_test_pheno / sqrt(h2_test_bayesr) else NA
    } else {
      h2_test_bayesr     <- NA
      accuracy_test_gebv <- NA
    }
  } else {
    h2_test_bayesr <- if (has_test)
      var(pred_test) / (var(pred_test) + sigma2_e_hat) else NA
    accuracy_test_gebv <-
      if (has_test && !is.na(accuracy_test_pheno) &&
          !is.na(h2_test_bayesr) && h2_test_bayesr > 0.01)
        accuracy_test_pheno / sqrt(h2_test_bayesr) else NA
  }

  # Marker-level CI / PIP / significance — keep custom (project-specific).
  # PIP prefers masbayes's pre-computed res$pip when GWAS is enabled (the
  # full-fit path supplies map → res$pip is populated by masbayes using
  # colMeans(gamma_samples != 0), identical formula). CV per-fold runs
  # don't pass map, so res$pip is NULL and we fall back to the manual
  # computation — bit-identical to pre-v0.5.0 behaviour.
  beta_ci_lower <- apply(res$beta_samples, 2, quantile, probs = 0.025, na.rm = TRUE)
  beta_ci_upper <- apply(res$beta_samples, 2, quantile, probs = 0.975, na.rm = TRUE)
  PIP <- if (!is.null(res$pip)) res$pip else colMeans(res$gamma_samples != 0)
  significant <- (beta_ci_lower > 0 | beta_ci_upper < 0)

  marker_info <- tryCatch({
    m_names <- if (!is.null(matrices$allele_info) &&
                   length(matrices$allele_info$allele_id) == ncol(W_train)) {
      matrices$allele_info$allele_id
    } else {
      colnames(W_train) %||% paste0("M", seq_len(ncol(W_train)))
    }

    df <- data.frame(
      allele_id   = m_names,
      beta_hat    = as.vector(beta_hat),
      beta_abs    = abs(as.vector(beta_hat)),
      ci_lower    = as.vector(beta_ci_lower),
      ci_upper    = as.vector(beta_ci_upper),
      component   = as.integer(gamma_mode),
      component_label = factor(as.integer(gamma_mode),
                               levels = 0:3,
                               labels = c("zero", "small", "medium", "large")),
      PIP         = as.vector(PIP),
      significant = as.logical(significant),
      stringsAsFactors = FALSE
    )
    if (!is.null(matrices$allele_info) &&
        length(matrices$allele_info$allele_id) == ncol(W_train)) {
      df$freq   <- matrices$allele_info$freq
      df$block  <- matrices$allele_info$block
      df$allele <- matrices$allele_info$allele
    }
    df
  }, error = function(e) {
    warning("Marker info for BayesR failed: ", e$message)
    return(NULL)
  })

  list(
    predictions = list(
      train = pred_train,
      test  = pred_test
    ),
    accuracy = list(
      train_pheno = accuracy_train_pheno,
      test_pheno  = accuracy_test_pheno,
      test_gebv   = accuracy_test_gebv,
      auc_train   = auc_train,
      auc_test    = auc_test
    ),
    regression_slopes = list(
      train_pheno = b_train_pheno,
      test_pheno  = b_test_pheno
    ),
    rmse = list(
      train_pheno = rmse_train_pheno,
      test_pheno  = rmse_test_pheno
    ),
    r2 = list(
      train_pheno = r2_train_pheno,
      test_pheno  = r2_test_pheno
    ),
    varcomp = list(
      sigma2_e      = sigma2_e_hat,
      sigma2_small  = sigma2_small_hat,
      sigma2_medium = sigma2_medium_hat,
      sigma2_large  = sigma2_large_hat,
      pi            = pi_hat,
      sigma2_g      = sigma2_g_hat,
      h2            = h2_hat,
      h2_test       = h2_test_bayesr,
      alpha_hat     = if (!is.null(alpha_hat)) as.numeric(alpha_hat) else NULL
    ),
    markers = marker_info,
    samples = if (config$save_mcmc_samples) {
      list(beta = res$beta_samples, gamma = res$gamma_samples)
    } else NULL,
    diagnostics = list(
      ess      = ess,
      geweke_z = geweke,
      mean_ess = mean(unlist(ess)),
      min_ess  = min(unlist(ess))
    ),
    # GWAS pass-through (masbayes v0.5.0). NULL when map = NULL — caller writes
    # these to a separate gwas_full.rds; full_fit.rds consumers can ignore them.
    # pip_class is always populated (derived from gamma_samples, not map).
    pip       = res$pip,
    pip_block = res$pip_block,
    pip_class = pip_class,
    gwas      = res$gwas,
    gwas_meta = res$gwas_meta,
    n_cores = 1,
    runtime = NA
  )
}
