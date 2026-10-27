# Shared by main-mh-fit-eval.R and main-snp-fit-eval.R: build train/test
# inputs, fit GBLUP + BayesR on reference individuals, score candidates.
# Requires helper.R, models/gblup.R, models/bayesR.R and common/bayesr_prior.R.

grm_summary <- function(A) {
  off <- A[lower.tri(A)]
  list(mean_diag = mean(diag(A)), sd_diag = sd(diag(A)),
       mean_offdiag = mean(off), sd_offdiag = sd(off))
}

# SNP: BayesR W = centred + globally scaled (as W_ah) or per-column z-score
# (encoding = "zscore"); GRM = VanRaden 1 either way.
# Frequencies, scale and sd come from the training (reference) set only.
build_snp_inputs <- function(X_train, X_test, allele_info = NULL, encoding = c("global", "zscore")) {
  encoding <- match.arg(encoding)
  ctr <- center_snp(X_train); cte <- center_snp(X_test, ctr$freq)
  k   <- 2 * sum(ctr$freq * (1 - ctr$freq))
  if (encoding == "global") {
    str_ <- scale_marker_w(ctr$W); ste <- scale_marker_w(cte$W, str_$scale)
  } else {
    ztr <- masbayes::construct_snp_matrix(X_train, encoding = "zscore")
    zte <- masbayes::construct_snp_matrix(X_test, encoding = "zscore", ref_freq = ztr$freq, ref_sd = ztr$sd)
    str_ <- list(W = ztr$W, scale = NA_real_); ste <- list(W = zte$W)
    dimnames(str_$W) <- dimnames(X_train); dimnames(ste$W) <- dimnames(X_test)
  }
  A_train    <- tcrossprod(ctr$W) / k
  A_combined <- tcrossprod(rbind(ctr$W, cte$W)) / k
  if (is.null(allele_info))
    allele_info <- data.frame(allele_id = colnames(X_train), freq = ctr$freq,
                              block = colnames(X_train), allele = 1L, chr = NA_integer_,
                              start_pos = NA_integer_, end_pos = NA_integer_, n_snps = 1L,
                              global_block_id = colnames(X_train))
  list(W_train = str_$W, W_test = ste$W, A_train = A_train, A_combined = A_combined,
       k_ah = k, w_scale = str_$scale, allele_info = allele_info, dropped_alleles = NULL,
       train_ids = rownames(X_train), test_ids = rownames(X_test), g_label = "snp_add",
       grm_stats = grm_summary(A_train))
}

# MH: W_ah from training allele frequencies; test reuses the training structure.
build_mh_inputs <- function(hap_train, hap_test, config, coords_file) {
  hcols_tr <- hap_train %>% select(-individual_id)
  hcols_te <- hap_test  %>% select(-individual_id)
  af   <- compute_allele_frequencies(hcols_tr, config)
  W_tr <- construct_W_matrix(hcols_tr, af, drop_baseline = config$drop_baseline)
  W_te <- construct_W_matrix(hcols_te, allele_freq_filtered = NULL,
                             reference_structure = W_tr, drop_baseline = config$drop_baseline)
  block_id <- attr(W_tr$W_ah, "block_id")
  keep <- apply(W_tr$W_ah, 2, var) > 0
  s_tr <- scale_marker_w(W_tr$W_ah[, keep, drop = FALSE])
  s_te <- scale_marker_w(W_te$W_ah[, keep, drop = FALSE], s_tr$scale)
  W_train <- s_tr$W; W_test <- s_te$W
  rownames(W_train) <- hap_train$individual_id; rownames(W_test) <- hap_test$individual_id
  if (!is.null(block_id)) attr(W_train, "block_id") <- block_id[keep]
  G_tr <- tcrossprod(W_train); k_ah <- sum(diag(G_tr)) / nrow(W_train)
  A_train <- G_tr / k_ah
  list(W_train = W_train, W_test = W_test, A_train = A_train,
       A_combined = tcrossprod(rbind(W_train, W_test)) / k_ah, k_ah = k_ah,
       w_scale = s_tr$scale,
       allele_info = enrich_allele_info(data.frame(W_tr$allele_info)[keep, ], coords_file),
       dropped_alleles = W_tr$dropped_alleles,
       train_ids = hap_train$individual_id, test_ids = hap_test$individual_id,
       g_label = "mh_add", grm_stats = grm_summary(A_train))
}

# Phenotypes of test (candidate) individuals are passed only for scoring.
fit_one_scenario <- function(pheno_all, inputs, config, map = NULL, windsize = NULL) {
  pick <- function(ids) pheno_all[match(ids, pheno_all$individual_id), ]
  pheno_train <- pick(inputs$train_ids); pheno_test <- pick(inputs$test_ids)
  stopifnot(!anyNA(pheno_train$individual_id), !anyNA(pheno_test$individual_id))
  split <- list(train = list(phenotypes = pheno_train), test = list(phenotypes = pheno_test))
  t0 <- proc.time()[["elapsed"]]
  gblup  <- run_gblup(inputs, split, config)
  t1 <- proc.time()[["elapsed"]]
  bayesR <- run_bayesR(inputs, split, gblup$varcomp, config, fold = 0, map = map, windsize = windsize)
  t2 <- proc.time()[["elapsed"]]
  list(gblup = gblup, bayesR = bayesR, runtime = c(gblup = t1 - t0, bayesR = t2 - t1),
       pheno_train = pheno_train, pheno_test = pheno_test,
       train_ids = inputs$train_ids, test_ids = inputs$test_ids, grm_stats = inputs$grm_stats)
}

# One row per candidate; accuracy/bias = candidates, *_train = reference (in-sample).
candidate_results <- function(fit, model, meta) {
  m <- fit[[model]]; g <- m$predictions$test; tbv <- fit$pheno_test$tbv
  gs <- fit$grm_stats %||% list(mean_diag = NA, sd_diag = NA, mean_offdiag = NA, sd_offdiag = NA)
  tibble(iter = meta$iter, individual_id = fit$test_ids, split = "candidate",
         tbv = tbv, GEBV = g, scenario = meta$scenario, model = model,
         marker = meta$marker, variant = meta$variant,
         accuracy = cor(g, tbv), bias = unname(coef(lm(tbv ~ g))[2]),
         rmse = sqrt(mean((g - tbv)^2)), h2 = m$varcomp$h2,
         n_train = length(fit$train_ids), n_test = length(fit$test_ids),
         acc_train = cor(m$predictions$train, fit$pheno_train$tbv),
         bias_train = unname(coef(lm(fit$pheno_train$tbv ~ m$predictions$train))[2]),
         grm_mean_diag = gs$mean_diag, grm_sd_diag = gs$sd_diag,
         grm_mean_offdiag = gs$mean_offdiag, grm_sd_offdiag = gs$sd_offdiag)
}
