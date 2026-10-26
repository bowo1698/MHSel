# BayesR with mixture prior

run_bayesR <- function(matrices, split, gblup_varcomp, config, fold = 0,
                       map = NULL, windsize = NULL) {

  W_train  <- matrices$W_train
  W_test   <- matrices$W_test

  pheno_train <- split$train$phenotypes
  pheno_test  <- split$test$phenotypes
  has_test    <- nrow(pheno_test) > 0

  y <- pheno_train$phenotype

  # Identical prior and kernel for SNP and MH (common/bayesr_prior.R)
  if (!exists("bayesr_prior")) source(file.path(config$common_dir, "bayesr_prior.R"))
  prior         <- bayesr_prior("gaussian")
  pi_vec        <- prior$pi
  sigma2_ah     <- gblup_varcomp$sigma2_g
  sigma2_e_init <- gblup_varcomp$sigma2_e

  prior_params <- list(
    a0_e           = prior$a0_e,
    a0_g           = prior$a0_g,
    variance_class = prior$variance_class
  )

  cat("BayesR Settings:\n")
  cat("π =", pi_vec, "\n")
  cat("σ²_αh =", round(sigma2_ah, 6), "| σ²_e =", round(sigma2_e_init, 6), "\n")
  cat("a0_e =", prior_params$a0_e, "| a0_g =", prior_params$a0_g, "\n")
  cat("var_class =", prior_params$variance_class, "\n")
  cat("fold:", fold, "\n")

  cat("Running BayesR with MCMC algorithm...\n\n")

  res <- masbayes::run_bayesr(
    w             = W_train,
    y             = y,
    marker_type   = prior$marker_type,
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
    fold_id       = as.integer(fold),
    save_rds      = FALSE,
    # NULL keeps the fit bitwise-identical to v1.4.0; when supplied,
    # masbayes attaches pip / pip_block / gwas / gwas_meta to res.
    map           = map,
    windsize      = windsize
  )

  beta_hat     <- res$beta_hat
  mu_hat       <- res$mu_hat
  sigma2_e_hat <- res$sigma2_e_hat
  sigma2_g_hat <- res$sigma2_g
  h2_hat       <- res$h2

  pred_train <- as.vector(predict(res)$GEBV)

  if (has_test) {
    pred_test_obj <- predict(res, W_test, pheno_test$phenotype)
    pred_test     <- as.vector(pred_test_obj$GEBV)
  } else {
    pred_test_obj <- NULL
    pred_test     <- numeric(0)
  }

  # Variance components from masbayes
  vc <- res$variance_components
  sigma2_small_hat  <- vc$small$mean
  sigma2_medium_hat <- vc$medium$mean
  sigma2_large_hat  <- vc$large$mean
  pi_hat <- c(vc$pi$zero, vc$pi$small, vc$pi$medium, vc$pi$large)

  # Mode gamma — keep manual
  gamma_mode <- apply(res$gamma_samples, 2, function(x) {
    tab <- table(x)
    if (length(tab) == 0) return(0)
    as.numeric(names(sort(tab, decreasing = TRUE)[1]))
  })

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
  cat("h² (heritability):", round(h2_hat, 4), "\n\n")

  # Train metrics from masbayes (phenotype scale)
  tm <- res$training_metrics
  accuracy_train_pheno <- tm$accuracy
  rmse_train_pheno     <- tm$RMSE
  r2_train_pheno       <- tm$R2
  b_train_pheno        <- tm$bias

  # tbv-based metrics — custom, masbayes doesn't know tbv
  accuracy_train_tbv <- cor(pred_train, pheno_train$tbv, use = "complete.obs")
  b_train_tbv        <- as.numeric(coef(lm(pheno_train$tbv ~ pred_train))[2])
  rmse_train_tbv     <- sqrt(mean((pred_train - pheno_train$tbv)^2, na.rm = TRUE))
  r2_train_tbv       <- accuracy_train_tbv^2

  if (has_test) {
    em <- pred_test_obj$metrics
    accuracy_test_pheno <- em$accuracy
    rmse_test_pheno     <- em$RMSE
    r2_test_pheno       <- em$R2
    b_test_pheno        <- em$bias

    accuracy_test_tbv <- cor(pred_test, pheno_test$tbv, use = "complete.obs")
    b_test_tbv        <- as.numeric(coef(lm(pheno_test$tbv ~ pred_test))[2])
    rmse_test_tbv     <- sqrt(mean((pred_test - pheno_test$tbv)^2, na.rm = TRUE))
    r2_test_tbv       <- accuracy_test_tbv^2

    h2_test_bayesr <- var(pred_test) / (var(pred_test) + sigma2_e_hat)
    accuracy_test_gebv <-
      if (!is.na(accuracy_test_pheno) && !is.na(h2_test_bayesr))
        accuracy_test_pheno / sqrt(h2_test_bayesr) else NA
  } else {
    accuracy_test_pheno <- NA
    rmse_test_pheno     <- NA
    r2_test_pheno       <- NA
    b_test_pheno        <- NA
    accuracy_test_tbv   <- NA
    b_test_tbv          <- NA
    rmse_test_tbv       <- NA
    r2_test_tbv         <- NA
    h2_test_bayesr      <- NA
    accuracy_test_gebv  <- NA
  }

  # Marker statistics — custom. Uses res$pip when the full-fit path supplies
  # map; split fits pass map = NULL, so we compute it manually (same formula).
  beta_ci_lower <- apply(res$beta_samples, 2, quantile, probs = 0.025, na.rm = TRUE)
  beta_ci_upper <- apply(res$beta_samples, 2, quantile, probs = 0.975, na.rm = TRUE)
  PIP <- if (!is.null(res$pip)) res$pip else colMeans(res$gamma_samples != 0)
  significant <- (beta_ci_lower > 0 | beta_ci_upper < 0)

  marker_info <- tryCatch({
    if (length(matrices$allele_info$allele_id) != ncol(W_train)) {
      stop("allele_info length (", length(matrices$allele_info$allele_id),
           ") != W_train columns (", ncol(W_train), ")")
    }

    data.frame(
      allele_id       = matrices$allele_info$allele_id,
      beta_hat        = as.vector(beta_hat),
      beta_abs        = abs(as.vector(beta_hat)),
      ci_lower        = as.vector(beta_ci_lower),
      ci_upper        = as.vector(beta_ci_upper),
      component       = as.integer(gamma_mode),
      component_label = factor(as.integer(gamma_mode),
                               levels = 0:3,
                               labels = c("zero", "small", "medium", "large")),
      PIP             = as.vector(PIP),
      significant     = as.logical(significant),
      freq            = matrices$allele_info$freq,
      block           = matrices$allele_info$block,
      allele          = matrices$allele_info$allele,
      chr             = matrices$allele_info$chr,
      start_pos       = matrices$allele_info$start_pos,
      end_pos         = matrices$allele_info$end_pos,
      n_snps          = matrices$allele_info$n_snps,
      global_block_id = matrices$allele_info$global_block_id,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("Marker info for BayesR failed: ", e$message)
    return(NULL)
  })

  # ── MCMC diagnostics (opt-in via config$mcmc_diag$enabled) ──
  mcmc_diag <- build_mcmc_diag(
    res, config,
    meta = list(
      model  = "bayesR",
      n_iter = as.integer(config$n_iter),
      n_burn = as.integer(config$n_burn),
      n_thin = as.integer(config$n_thin),
      n_kept = length(res$sigma2_e_samples),
      seed   = as.integer(config$seed)
    )
  )

  list(
    coef = list(
      beta_hat = as.vector(beta_hat),
      mu_hat   = as.numeric(mu_hat)
    ),
    predictions = list(
      train = pred_train,
      test  = pred_test
    ),
    accuracy = list(
      train_pheno = accuracy_train_pheno,
      train_tbv   = accuracy_train_tbv,
      test_pheno  = accuracy_test_pheno,
      test_tbv    = accuracy_test_tbv,
      test_gebv   = accuracy_test_gebv
    ),
    regression_slopes = list(
      train_pheno = b_train_pheno,
      train_tbv   = b_train_tbv,
      test_pheno  = b_test_pheno,
      test_tbv    = b_test_tbv
    ),
    rmse = list(
      train_pheno = rmse_train_pheno,
      train_tbv   = rmse_train_tbv,
      test_pheno  = rmse_test_pheno,
      test_tbv    = rmse_test_tbv
    ),
    r2 = list(
      train_pheno = r2_train_pheno,
      train_tbv   = r2_train_tbv,
      test_pheno  = r2_test_pheno,
      test_tbv    = r2_test_tbv
    ),
    varcomp = list(
      sigma2_e      = sigma2_e_hat,
      sigma2_small  = sigma2_small_hat,
      sigma2_medium = sigma2_medium_hat,
      sigma2_large  = sigma2_large_hat,
      pi            = pi_hat,
      sigma2_g      = sigma2_g_hat,
      h2            = h2_hat,
      h2_test       = h2_test_bayesr
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
    # NULL when map = NULL — caller writes these to a separate gwas_full.rds.
    pip       = res$pip,
    pip_block = res$pip_block,
    gwas      = res$gwas,
    gwas_meta = res$gwas_meta,
    n_cores   = 1,
    runtime   = NA,
    mcmc_diag = mcmc_diag
  )
}
