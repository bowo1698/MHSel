# GBLUP via masreml (REML for VC, Cholesky/PCG for BLUP).
# G comes pre-built in matrices$A_train / A_combined.

run_gblup <- function(matrices, split, config) {

  pheno_train <- split$train$phenotypes
  pheno_test  <- split$test$phenotypes
  has_test    <- nrow(pheno_test) > 0

  fe_col <- config$fixed_effect

  pheno_train <- pheno_train %>%
    dplyr::select(ID, phenotype, any_of(fe_col)) %>%
    dplyr::rename(y = phenotype)
  pheno_train$ID <- as.character(pheno_train$ID)

  pheno_test <- pheno_test %>%
    dplyr::select(ID, phenotype, any_of(fe_col)) %>%
    dplyr::rename(y = phenotype)
  pheno_test$ID <- as.character(pheno_test$ID)

  g_label   <- matrices$g_label %||% "g"
  train_ids <- matrices$train_ids %||% as.character(pheno_train$ID)
  test_ids  <- matrices$test_ids  %||% as.character(pheno_test$ID)

  # G matrices (pre-built upstream)
  G_train <- matrices$A_train
  G_full  <- matrices$A_combined
  rownames(G_train) <- colnames(G_train) <- train_ids
  rownames(G_full)  <- colnames(G_full)  <- c(train_ids, test_ids)

  trait_type <- if (!is.null(config$response_type) &&
                    config$response_type == "binary") "binary" else "continuous"
  link_type  <- if (!is.null(config$link)) config$link else "logit"

  # Fixed effects: intercept + factor(fe_col) if set
  X_train <- if (!is.null(fe_col) && fe_col %in% names(pheno_train)) {
    cbind(1, as.numeric(as.factor(pheno_train[[fe_col]])))
  } else {
    matrix(1, nrow(pheno_train), 1)
  }
  rownames(X_train) <- pheno_train$ID

  y_train_named <- setNames(pheno_train$y, pheno_train$ID)

  gblup_fit <- masreml(
    y      = y_train_named,
    X      = X_train,
    G      = setNames(list(G_train), g_label),
    method = "AI",
    trait  = trait_type,
    link   = link_type,
    tol    = 1e-6
  )

  sigma2_g <- as.numeric(gblup_fit$varcomp$sigma2[g_label])
  sigma2_e <- as.numeric(gblup_fit$varcomp$sigma2["residual"])
  h2       <- as.numeric(gblup_fit$varcomp$h2[g_label])
  mu_hat   <- as.numeric(gblup_fit$fixed_effects[1])

  pred_train <- as.numeric(gblup_fit$total_gebv) + mu_hat

  # Test-set h2 from an independent fit on test individuals
  h2_test_gblup <- NA_real_
  if (has_test && nrow(pheno_test) > 1) {
    h2_test_gblup <- tryCatch({
      G_test_block <- G_full[test_ids, test_ids, drop = FALSE]
      X_test_fe <- if (!is.null(fe_col) && fe_col %in% names(pheno_test)) {
        m <- cbind(1, as.numeric(as.factor(pheno_test[[fe_col]])))
        rownames(m) <- pheno_test$ID
        m
      } else {
        m <- matrix(1, nrow(pheno_test), 1)
        rownames(m) <- pheno_test$ID
        m
      }
      y_test_named <- setNames(pheno_test$y, pheno_test$ID)
      fit_test <- masreml(
        y      = y_test_named,
        X      = X_test_fe,
        G      = setNames(list(G_test_block), g_label),
        method = "AI",
        trait  = trait_type,
        link   = link_type,
        tol    = 1e-6
      )
      h2_raw <- as.numeric(fit_test$varcomp$h2[g_label])
      if (is.na(h2_raw) || h2_raw < 0 || h2_raw > 1) NA_real_ else h2_raw
    }, error = function(e) NA_real_)
  }

  # Test prediction via masreml S3 predict (Cholesky on G_train + lambda*I)
  if (has_test) {
    gblup_pred <- predict(
      gblup_fit,
      G_full    = setNames(list(G_full), g_label),
      train_ids = train_ids,
      test_ids  = test_ids
    )
    pred_test <- as.numeric(gblup_pred$GEBV) + mu_hat

    auc_test <- if (trait_type == "binary") {
      tryCatch(
        as.numeric(pROC::auc(pROC::roc(pheno_test$y, pred_test, quiet = TRUE))),
        error = function(e) NA_real_
      )
    } else NA_real_
  } else {
    pred_test <- numeric(0)
    auc_test  <- NA_real_
  }

  calculate_regression_slope <- function(predicted, observed) {
    fit <- lm(observed ~ predicted)
    coef(fit)[2]
  }

  # Binary: score on the probability scale, to match training_metrics and AUC.
  # Continuous: predictions are already on the phenotype scale.
  if (trait_type == "binary") {
    score_train <- as.numeric(gblup_fit$binary$fitted)
    if (has_test) {
      link_fn     <- if (link_type == "probit") stats::pnorm else stats::plogis
      score_test  <- pmax(pmin(link_fn(pred_test), 1 - 1e-6), 1e-6)
    }
  } else {
    score_train <- pred_train
    if (has_test) score_test <- pred_test
  }

  accuracy_train_pheno <- cor(score_train, pheno_train$y, use = "complete.obs")
  b_train_pheno        <- calculate_regression_slope(score_train, pheno_train$y)
  rmse_train_pheno     <- sqrt(mean((score_train - pheno_train$y)^2, na.rm = TRUE))
  r2_train_pheno       <- accuracy_train_pheno^2

  if (has_test) {
    accuracy_test_pheno <- cor(score_test, pheno_test$y, use = "complete.obs")
    accuracy_test_gebv  <- if (!is.na(h2_test_gblup) && h2_test_gblup > 0)
      accuracy_test_pheno / sqrt(h2_test_gblup) else NA_real_
    b_test_pheno    <- calculate_regression_slope(score_test, pheno_test$y)
    rmse_test_pheno <- sqrt(mean((score_test - pheno_test$y)^2, na.rm = TRUE))
    r2_test_pheno   <- accuracy_test_pheno^2
  } else {
    accuracy_test_pheno <- NA_real_
    accuracy_test_gebv  <- NA_real_
    b_test_pheno        <- NA_real_
    rmse_test_pheno     <- NA_real_
    r2_test_pheno       <- NA_real_
  }

  list(
    predictions = list(
      train = pred_train,
      test  = pred_test
    ),
    accuracy = list(
      train_pheno = accuracy_train_pheno,
      test_pheno  = accuracy_test_pheno,
      test_gebv   = accuracy_test_gebv,
      auc_train   = if (trait_type == "binary") gblup_fit$binary$auc else NA_real_,
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
