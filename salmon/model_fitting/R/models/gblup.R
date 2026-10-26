# GBLUP via masreml (REML for VC, Cholesky/PCG for BLUP).
# G comes pre-built in matrices$A_train / A_combined.

run_gblup <- function(matrices, split, config) {

  pheno_train <- split$train$phenotypes
  pheno_test  <- split$test$phenotypes
  has_test    <- nrow(pheno_test) > 0

  pheno_train <- pheno_train %>%
    dplyr::select(individual_id, phenotype, tbv) %>%
    dplyr::rename(y = phenotype)
  pheno_train$individual_id <- as.character(pheno_train$individual_id)

  pheno_test <- pheno_test %>%
    dplyr::select(individual_id, phenotype, tbv) %>%
    dplyr::rename(y = phenotype)
  pheno_test$individual_id <- as.character(pheno_test$individual_id)

  g_label   <- matrices$g_label %||% "g"
  train_ids <- matrices$train_ids %||% pheno_train$individual_id
  test_ids  <- matrices$test_ids  %||% pheno_test$individual_id

  G_train <- matrices$A_train
  G_full  <- matrices$A_combined
  rownames(G_train) <- colnames(G_train) <- train_ids
  rownames(G_full)  <- colnames(G_full)  <- c(train_ids, test_ids)

  y_train_named <- setNames(pheno_train$y, pheno_train$individual_id)

  gblup_fit <- masreml(
    y      = y_train_named,
    G      = setNames(list(G_train), g_label),
    method = "AI",
    tol    = 1e-6
  )

  sigma2_g <- as.numeric(gblup_fit$varcomp$sigma2[g_label])
  sigma2_e <- as.numeric(gblup_fit$varcomp$sigma2["residual"])
  h2       <- as.numeric(gblup_fit$varcomp$h2[g_label])
  mu_hat   <- as.numeric(gblup_fit$fixed_effects[1])

  pred_train <- as.numeric(gblup_fit$total_gebv) + mu_hat

  # Candidate GEBV = BLUP G_te,tr (G_tr + lambda I)^-1 (y - mu). masreml's
  # predict(G_full = ...) solves against the training BLUP instead of (y - mu),
  # shrinking candidates twice (review C1); computed here directly.
  # Candidate phenotypes are never used (no h2 fit on the test set).
  h2_test_gblup <- NA_real_
  if (has_test) {
    lambda    <- sigma2_e / sigma2_g
    r_train   <- pheno_train$y[match(train_ids, pheno_train$individual_id)] - mu_hat
    pred_test <- as.vector(G_full[test_ids, train_ids, drop = FALSE] %*%
                             solve(G_train + diag(lambda, nrow(G_train)), r_train)) + mu_hat
  } else {
    pred_test <- numeric(0)
  }

  calculate_regression_slope <- function(predicted, observed) {
    fit <- lm(observed ~ predicted)
    coef(fit)[2]
  }

  accuracy_train_pheno <- cor(pred_train, pheno_train$y)
  accuracy_train_tbv   <- cor(pred_train, pheno_train$tbv)
  b_train_pheno        <- calculate_regression_slope(pred_train, pheno_train$y)
  b_train_tbv          <- calculate_regression_slope(pred_train, pheno_train$tbv)
  rmse_train_pheno     <- sqrt(mean((pred_train - pheno_train$y)^2, na.rm = TRUE))
  rmse_train_tbv       <- sqrt(mean((pred_train - pheno_train$tbv)^2, na.rm = TRUE))
  r2_train_pheno       <- accuracy_train_pheno^2
  r2_train_tbv         <- accuracy_train_tbv^2

  if (has_test) {
    accuracy_test_pheno <- cor(pred_test, pheno_test$y)
    accuracy_test_tbv   <- cor(pred_test, pheno_test$tbv)
    accuracy_test_gebv  <- if (!is.na(h2_test_gblup) && h2_test_gblup > 0)
      accuracy_test_pheno / sqrt(h2_test_gblup) else NA_real_
    b_test_pheno    <- calculate_regression_slope(pred_test, pheno_test$y)
    b_test_tbv      <- calculate_regression_slope(pred_test, pheno_test$tbv)
    rmse_test_pheno <- sqrt(mean((pred_test - pheno_test$y)^2, na.rm = TRUE))
    rmse_test_tbv   <- sqrt(mean((pred_test - pheno_test$tbv)^2, na.rm = TRUE))
    r2_test_pheno   <- accuracy_test_pheno^2
    r2_test_tbv     <- accuracy_test_tbv^2
  } else {
    accuracy_test_pheno <- NA_real_
    accuracy_test_tbv   <- NA_real_
    accuracy_test_gebv  <- NA_real_
    b_test_pheno        <- NA_real_
    b_test_tbv          <- NA_real_
    rmse_test_pheno     <- NA_real_
    rmse_test_tbv       <- NA_real_
    r2_test_pheno       <- NA_real_
    r2_test_tbv         <- NA_real_
  }

  list(
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
      sigma2_g           = sigma2_g,
      sigma2_e           = sigma2_e,
      h2                 = h2,
      h2_test            = h2_test_gblup,
      lambda             = sigma2_e / sigma2_g,
      mean_diag_train    = mean(diag(G_train)),
      mean_diag_combined = mean(diag(G_full))
    ),
    model   = gblup_fit,
    n_cores = 1,
    runtime = NA
  )
}
