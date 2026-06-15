# BayesA Model with Marker-Specific Variance

run_bayesA <- function(matrices, split, gblup_varcomp, config, fold = 0) {

  W_train  <- matrices$W_train
  W_test   <- matrices$W_test
  WtW_diag <- matrices$WtW_diag

  pheno_train <- split$train$phenotypes
  pheno_test  <- split$test$phenotypes
  has_test    <- nrow(pheno_test) > 0

  y <- pheno_train$phenotype

  nu            <- config$bayesA$nu
  is_snp        <- identical(config$marker_type, "snp")
  vary          <- var(y)
  sigma2_g      <- if (is_snp) vary / 4 else gblup_varcomp$sigma2_g
  sigma2_e_init <- if (is_snp) vary / 2 else gblup_varcomp$sigma2_e
  a0_e          <- config$bayesA$prior_df_residual

  prior_params <- list(a0_e = a0_e)

  cat("BayesA Settings:\n")
  cat("ν =", nu, "| σ²_g =", round(sigma2_g, 6), "\n")
  cat("σ²_e_init =", round(sigma2_e_init, 6), "\n")
  cat("fold:", fold, "\n\n")

  cat("Running BayesA with MCMC algorithm...\n\n")

  res <- masbayes::run_bayesa(
    w             = W_train,
    y             = y,
    wtw_diag      = WtW_diag,
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
    fold_id       = as.integer(fold),
    save_rds      = FALSE
  )

  beta_A_hat   <- res$beta_hat
  mu_hat       <- res$mu_hat
  sigma2_j_hat <- res$sigma2_j_hat
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
  cat("h² (heritability):", round(h2_hat, 4), "\n\n")

  # Train metrics from masbayes (phenotype scale)
  tm <- res$training_metrics
  accuracy_train_pheno <- tm$accuracy
  rmse_train_pheno     <- tm$RMSE
  r2_train_pheno       <- tm$R2
  b_train_pheno        <- tm$bias

  # tbv-based metrics — keep custom
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

    h2_test_bayesa <- var(pred_test) / (var(pred_test) + sigma2_e_hat)
    accuracy_test_gebv <-
      if (!is.na(accuracy_test_pheno) && !is.na(h2_test_bayesa))
        accuracy_test_pheno / sqrt(h2_test_bayesa) else NA
  } else {
    accuracy_test_pheno <- NA
    rmse_test_pheno     <- NA
    r2_test_pheno       <- NA
    b_test_pheno        <- NA
    accuracy_test_tbv   <- NA
    b_test_tbv          <- NA
    rmse_test_tbv       <- NA
    r2_test_tbv         <- NA
    h2_test_bayesa      <- NA
    accuracy_test_gebv  <- NA
  }

  # Marker statistics — keep custom
  beta_ci_lower <- apply(res$beta_samples, 2, quantile, probs = 0.025, na.rm = TRUE)
  beta_ci_upper <- apply(res$beta_samples, 2, quantile, probs = 0.975, na.rm = TRUE)
  significant <- (beta_ci_lower > 0 | beta_ci_upper < 0)

  marker_info <- tryCatch({
    data.frame(
      allele_id       = matrices$allele_info$allele_id,
      beta_hat        = as.vector(beta_A_hat),
      beta_abs        = abs(as.vector(beta_A_hat)),
      ci_lower        = as.vector(beta_ci_lower),
      ci_upper        = as.vector(beta_ci_upper),
      sigma2_j        = as.vector(sigma2_j_hat),
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
    warning("Marker info for BayesA failed: ", e$message)
    return(NULL)
  })

  # ── MCMC diagnostics (opt-in via config$mcmc_diag$enabled) ──
  mcmc_diag <- build_mcmc_diag(
    res, config,
    meta = list(
      model  = "bayesA",
      n_iter = as.integer(config$n_iter),
      n_burn = as.integer(config$n_burn),
      n_thin = as.integer(config$n_thin),
      n_kept = length(res$sigma2_e_samples),
      seed   = as.integer(config$seed)
    )
  )

  list(
    coef = list(
      beta_hat = as.vector(beta_A_hat),
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
      sigma2_e        = sigma2_e_hat,
      mean_sigma2_j   = mean(sigma2_j_hat),
      range_sigma2_j  = range(sigma2_j_hat),
      sigma2_g        = sigma2_g_hat,
      h2              = h2_hat,
      h2_test         = h2_test_bayesa
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
    n_cores   = 1,
    runtime   = NA,
    mcmc_diag = mcmc_diag
  )
}
