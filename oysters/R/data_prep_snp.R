# SNP data loading and preprocessing

suppressPackageStartupMessages({ library(tidyverse) })

load_genomic_data <- function(config) {

  pheno <- read_csv(file.path(config$data_dir, config$phenotype_file), show_col_types = FALSE) %>%
    rename(phenotype = all_of(config$phenotype_trait)) %>%
    { df <- .
      cols_keep <- c("ID", "phenotype",
                     intersect(c(config$fixed_effect, config$extra_cols), names(df)))
      df[, cols_keep] }
  pheno <- pheno[!is.na(pheno$phenotype), ]

  geno <- read_csv(file.path(config$data_dir, config$genotype_file), show_col_types = FALSE)

  common_ids <- intersect(pheno$ID, geno$ID)
  cat("Phenotype N:", nrow(pheno), "| Genotype N:", nrow(geno), "| Common:", length(common_ids), "\n")

  pheno <- pheno %>% filter(ID %in% common_ids) %>% arrange(match(ID, common_ids))
  geno  <- geno  %>% filter(ID %in% common_ids) %>% arrange(match(ID, common_ids))

  geno_matrix <- as.matrix(geno %>% select(-ID))
  storage.mode(geno_matrix) <- "numeric"

  qc  <- quality_control(geno_matrix, config)
  imp <- impute_only(qc$filtered)

  cat("Final: N =", length(common_ids), "| SNPs =", qc$n_retained, "\n\n")

  af  <- colMeans(imp$imputed) / 2
  Z   <- sweep(imp$imputed, 2, 2 * af, "-")
  k   <- 2 * sum(af * (1 - af))
  GRM <- tcrossprod(Z) / k

  list(
    phenotypes = pheno,
    genotypes  = data.frame(ID = common_ids, imp$imputed),
    n_snps     = qc$n_retained,
    GRM        = GRM
  )
}

quality_control <- function(geno_matrix, config) {
  allele_freq <- colMeans(geno_matrix, na.rm = TRUE) / 2
  maf         <- pmin(allele_freq, 1 - allele_freq)
  callrate    <- colMeans(!is.na(geno_matrix))
  keep        <- (maf >= config$min_allele_freq) &
                 (callrate >= config$callrate_threshold)
  cat("QC: retained", sum(keep), "/", ncol(geno_matrix), "SNPs\n")
  list(filtered = geno_matrix[, keep], keep_snps = keep, n_retained = sum(keep))
}

impute_only <- function(geno_filtered) {
  geno_imp <- apply(geno_filtered, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE); x
  })
  zero_var <- apply(geno_imp, 2, var) < 1e-10
  if (sum(zero_var) > 0) {
    cat("Removing", sum(zero_var), "zero-variance SNPs after imputation\n")
    geno_imp <- geno_imp[, !zero_var]
  }
  list(imputed = geno_imp)
}

setup_cv <- function(data, config) {
  set.seed(config$cv_seed)
  k      <- config$k_folds
  n_clus <- config$n_grm_clusters %||% 10

  if (!is.null(config$cv_method) && config$cv_method == "cluster") {
    if (is.null(data$GRM)) stop("GRM is required for CV clusters")
    eig <- eigen(data$GRM, symmetric = TRUE)
    pcs <- eig$vectors[, 1:5]
    km  <- kmeans(pcs, centers = k, nstart = 20, iter.max = 100)
    cat("Cluster-based CV | k =", k, "\n")
    cat("Cluster sizes:\n"); print(table(km$cluster))
    return(list(folds = km$cluster, method = "cluster", reason = "grm_kmeans_pc5"))
  }

  if (!is.null(config$cv_method) && config$cv_method == "grm_stratified") {
    if (is.null(data$GRM)) stop("GRM is required for grm_stratified CV")
    eig <- eigen(data$GRM, symmetric = TRUE)
    pcs <- eig$vectors[, 1:5]
    km  <- kmeans(pcs, centers = n_clus, nstart = 20, iter.max = 100)
    fold_assignments <- integer(nrow(data$phenotypes))
    for (cl in 1:n_clus) {
      idx <- which(km$cluster == cl)
      fold_assignments[idx] <- sample(rep(1:k, length.out = length(idx)))
    }
    cat("GRM-stratified CV | k_folds =", k, "| n_clusters =", n_clus, "\n")
    cat("Cluster sizes:\n"); print(table(km$cluster))
    cat("Fold sizes:\n"); print(table(fold_assignments))
    return(list(folds = fold_assignments, method = "grm_stratified",
                reason = "grm_kmeans_stratified", n_clusters = n_clus))
  }

  folds <- sample(rep(1:k, length.out = nrow(data$phenotypes)))
  list(folds = folds, method = "random", reason = "random_10fold")
}

get_fold_split_info <- function(data, cv_setup, fold_id) {
  test_idx  <- which(cv_setup$folds == fold_id)
  train_idx <- which(cv_setup$folds != fold_id)
  list(
    fold       = fold_id,
    n_train    = length(train_idx),
    n_test     = length(test_idx),
    prop_train = length(train_idx) / nrow(data$phenotypes),
    prop_test  = length(test_idx)  / nrow(data$phenotypes)
  )
}

split_fold_data <- function(data, fold_assignments, fold_id) {
  test_idx  <- which(fold_assignments == fold_id)
  train_idx <- which(fold_assignments != fold_id)
  list(
    train = list(
      phenotypes = data$phenotypes[train_idx, ],
      genotypes  = data$genotypes[train_idx, ]
    ),
    test = list(
      phenotypes = data$phenotypes[test_idx, ],
      genotypes  = data$genotypes[test_idx, ]
    )
  )
}

# Build matrices for one (train, test) split.
# G GBLUP always from VanRaden W = X - 2p (classical GBLUP semantics).
# config$snp_encoding = "zscore" z-scores the Bayesian W only; A stays VanRaden.
# Missing snp_encoding defaults to "vanRaden": one W shared by both.
construct_matrices <- function(train_data, test_data, config) {

  geno_train <- as.matrix(train_data$genotypes %>% select(-ID))
  geno_test  <- as.matrix(test_data$genotypes  %>% select(-ID))
  storage.mode(geno_train) <- "double"
  storage.mode(geno_test)  <- "double"

  train_ids <- as.character(train_data$phenotypes$ID)
  test_ids  <- as.character(test_data$phenotypes$ID)
  rownames(geno_train) <- train_ids
  rownames(geno_test)  <- test_ids

  # Filter zero-variance SNPs based on training set
  col_var <- apply(geno_train, 2, var)
  keep    <- col_var > 1e-10
  if (sum(!keep) > 0) cat("Removing", sum(!keep), "zero-variance SNPs\n")
  geno_train <- geno_train[, keep, drop = FALSE]
  geno_test  <- geno_test[, keep, drop = FALSE]

  encoding <- config$snp_encoding %||% "vanRaden"

  # Always built for G GBLUP; also the Bayesian W when encoding == "vanRaden".
  vR_train <- masbayes::construct_snp_matrix(geno_train, encoding = "vanRaden")
  vR_test  <- masbayes::construct_snp_matrix(geno_test, encoding = "vanRaden",
                                             ref_freq = vR_train$freq)
  W_vR_train <- vR_train$W
  W_vR_test  <- vR_test$W
  rownames(W_vR_train) <- train_ids
  rownames(W_vR_test)  <- test_ids

  freq    <- vR_train$freq
  k_grm   <- 2 * sum(freq * (1 - freq))
  G_train <- tcrossprod(W_vR_train) / k_grm
  W_vR_full <- rbind(W_vR_train, W_vR_test)
  G_full    <- tcrossprod(W_vR_full) / k_grm

  rownames(G_train) <- colnames(G_train) <- train_ids
  rownames(G_full)  <- colnames(G_full)  <- c(train_ids, test_ids)

  if (encoding == "zscore") {
    zs_train <- masbayes::construct_snp_matrix(geno_train,
                                               encoding = "zscore")
    zs_test  <- masbayes::construct_snp_matrix(geno_test,
                                               encoding = "zscore",
                                               ref_freq = zs_train$freq,
                                               ref_sd   = zs_train$sd)
    W_train <- zs_train$W
    W_test  <- zs_test$W
    rownames(W_train) <- train_ids
    rownames(W_test)  <- test_ids
  } else {
    W_train <- W_vR_train
    W_test  <- W_vR_test
  }

  cat("W_train:", dim(W_train), "| encoding:", encoding,
      "| k_grm:", round(k_grm, 4), "\n")
  cat("G_train diagonal mean:", round(mean(diag(G_train)), 3), "\n\n")

  list(
    W_train    = W_train,
    W_test     = W_test,
    A_train    = G_train,
    A_combined = G_full,
    k_ah       = k_grm,
    train_ids  = train_ids,
    test_ids   = test_ids,
    g_label    = "snp_add",
    n_snps     = ncol(W_train)
  )
}
