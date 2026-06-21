# BayesA Model with Marker-Specific Variance

run_bayesA <- function(matrices, split, gblup_varcomp, config, fold = 0) {

  W_train  <- matrices$W_train
  W_test   <- matrices$W_test

  pheno_train_df <- split$train$phenotypes
  pheno_test_df  <- split$test$phenotypes
  has_test       <- nrow(pheno_test_df) > 0

  y         <- pheno_train_df$phenotype
  n_alleles <- ncol(W_train)

  fe_col <- config$fixed_effect
  if (!is.null(fe_col) && fe_col %in% names(pheno_train_df)) {
    fe_levels <- levels(as.factor(pheno_train_df[[fe_col]]))
    X_cov     <- as.matrix(as.numeric(factor(pheno_train_df[[fe_col]],
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

  nu            <- config$bayesA$nu
  sigma2_g      <- gblup_varcomp$sigma2_g
  sigma2_e_init <- gblup_varcomp$sigma2_e
  a0_e          <- config$bayesA$prior_df_residual

  prior_params <- list(a0_e = a0_e)

  cat("BayesA Settings:\n")
  cat("ν =", nu, "| σ²_g =", round(sigma2_g, 6), "\n")
  cat("σ²_e_init =", round(sigma2_e_init, 6), "\n")
  cat("fixed effect:", if (is.null(X_cov)) "none" else
                      sprintf("%s (q=%d)", fe_col, ncol(X_cov)), "\n")
  cat("fold:", fold, "\n\n")

  cat("Running BayesA with MCMC algorithm...\n\n")

  res <- masbayes::run_bayesa(
    w             = W_train,
    y             = y,
    X             = X_cov,
    marker_type   = config$marker_type %||% "auto",
    nu            = nu,
    sigma2_g      = sigma2_g,
    sigma2_e_init = sigma2_e_init,
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
    save_rds      = FALSE
  )

  beta_A_hat   <- res$beta_hat
  sigma2_j_hat <- res$sigma2_j_hat
  sigma2_e_hat <- res$sigma2_e_hat
  alpha_hat    <- res$alpha_hat
  sigma2_g_hat <- res$sigma2_g
  h2_hat       <- res$h2

  pred_train <- as.vector(predict(res)$GEBV)

  if (has_test) {
    pred_test_obj <- predict(res, W_test, pheno_test_df$phenotype,
                             X_new = X_test_cov)
    pred_test <- as.vector(pred_test_obj$GEBV)
  } else {
    pred_test_obj <- NULL
    pred_test     <- numeric(0)
  }

  # ESS / Geweke from masbayes diagnostics
  d <- res$diagnostics
  ess_A <- list(
    sigma2_e      = d$sigma2_e$ess,
    sigma2_j_mean = d$sigma2_j$ess_median
  )
  geweke_A <- list(
    sigma2_e      = d$sigma2_e$geweke,
    sigma2_j_mean = d$sigma2_j$geweke_median
  )

  cat("\nBayesA Posterior Means:\n")
  cat("σ²_e:", round(sigma2_e_hat, 6), "\n")
  cat("Mean σ²_j:", round(mean(sigma2_j_hat), 8), "\n")
  cat("Range σ²_j:", round(range(sigma2_j_hat), 8), "\n")
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

  if (is_binary) {
    if (has_test) {
      P <- mean(pheno_train_df$phenotype, na.rm = TRUE)
      z <- dnorm(qnorm(P))
      h2_obs <- var(pred_test) / (var(pred_test) + sigma2_e_hat)
      h2_liability <- h2_obs * P * (1 - P) / z^2
      h2_test_bayesa <- h2_liability
      accuracy_test_gebv <-
        if (!is.na(accuracy_test_pheno) && !is.na(h2_liability) &&
            h2_liability > 0)
          accuracy_test_pheno / sqrt(h2_liability) else NA
    } else {
      h2_test_bayesa     <- NA
      accuracy_test_gebv <- NA
    }
  } else {
    h2_test_bayesa <- if (has_test)
      var(pred_test) / (var(pred_test) + sigma2_e_hat) else NA
    accuracy_test_gebv <-
      if (has_test && !is.na(accuracy_test_pheno) &&
          !is.na(h2_test_bayesa) && h2_test_bayesa > 0.01)
        accuracy_test_pheno / sqrt(h2_test_bayesa) else NA
  }

  # Marker-level CI / significance — keep custom
  beta_ci_lower <- apply(res$beta_samples, 2, quantile, probs = 0.025, na.rm = TRUE)
  beta_ci_upper <- apply(res$beta_samples, 2, quantile, probs = 0.975, na.rm = TRUE)
  significant <- (beta_ci_lower > 0 | beta_ci_upper < 0)

  marker_info <- tryCatch({
    m_names <- colnames(W_train)
    if (is.null(m_names)) m_names <- paste0("M", seq_len(n_alleles))

    data.frame(
      allele_id   = m_names,
      beta_hat    = as.vector(beta_A_hat),
      beta_abs    = abs(as.vector(beta_A_hat)),
      ci_lower    = as.vector(beta_ci_lower),
      ci_upper    = as.vector(beta_ci_upper),
      sigma2_j    = as.vector(sigma2_j_hat),
      significant = as.logical(significant),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("Marker info for BayesA failed: ", e$message)
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
      sigma2_e        = sigma2_e_hat,
      mean_sigma2_j   = mean(sigma2_j_hat),
      range_sigma2_j  = range(sigma2_j_hat),
      sigma2_g        = sigma2_g_hat,
      h2              = h2_hat,
      h2_test         = h2_test_bayesa,
      alpha_hat       = if (!is.null(alpha_hat)) as.numeric(alpha_hat) else NULL
    ),
    markers = marker_info,
    samples = if (config$save_mcmc_samples) {
      list(beta = res$beta_samples, sigma2_j = res$sigma2_j_samples)
    } else NULL,
    diagnostics = list(
      ess      = ess_A,
      geweke_z = geweke_A,
      mean_ess = mean(unlist(ess_A)),
      min_ess  = min(unlist(ess_A))
    ),
    n_cores = 1,
    runtime = NA
  )
}
