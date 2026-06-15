# Cross-Validation Setup Functions

check_family_structure <- function(phenotypes, min_offspring = 5) {
  
  # Check family structure for stratification feasibility
  family_stats <- phenotypes %>%
    count(father_id) %>%
    summarise(
      n_sires = n(),
      min_offspring = min(n),
      median_offspring = median(n),
      max_offspring = max(n),
      sires_with_enough = sum(n >= min_offspring),
      prop_adequate = sires_with_enough / n_sires
    )
  
  list(
    stats = family_stats,
    feasible = family_stats$n_sires > 1 && 
               family_stats$median_offspring >= min_offspring &&
               family_stats$prop_adequate >= 0.8
  )
}

setup_cluster_cv <- function(data, config, G = NULL) {
  set.seed(config$cv_seed)
  k <- config$k_folds

  if (is.null(G)) stop("GRM is required for cluster-based CV")

  eig <- eigen(G, symmetric = TRUE)
  pcs <- eig$vectors[, 1:5]

  km <- kmeans(pcs, centers = k, nstart = 20, iter.max = 100)
  fold_assignments <- km$cluster

  cat("Cluster-based CV | k =", k, "\n")
  cat("Cluster sizes:\n")
  print(table(fold_assignments))

  list(
    folds  = fold_assignments,
    method = "cluster",
    reason = "grm_kmeans_pc5"
  )
}

setup_grm_stratified_cv <- function(data, config, G = NULL) {
  set.seed(config$cv_seed)
  k      <- config$k_folds
  n_clus <- config$n_grm_clusters %||% 10  # default 10 cluster

  if (is.null(G)) stop("GRM is required for GRM-stratified CV")

  # PCA dari GRM
  eig <- eigen(G, symmetric = TRUE)
  pcs <- eig$vectors[, 1:5]

  # K-means dengan n_clus cluster
  km <- kmeans(pcs, centers = n_clus, nstart = 20, iter.max = 100)

  # Stratified assignment: dari tiap cluster, bagi merata ke k fold
  n_ind <- nrow(data$phenotypes)
  fold_assignments <- integer(n_ind)
  for (cl in 1:n_clus) {
    idx <- which(km$cluster == cl)
    fold_assignments[idx] <- sample(rep(1:k, length.out = length(idx)))
  }

  cat("GRM-stratified CV | k_folds =", k, "| n_clusters =", n_clus, "\n")
  cat("Cluster sizes:\n"); print(table(km$cluster))
  cat("Fold sizes:\n"); print(table(fold_assignments))

  list(
    folds       = fold_assignments,
    method      = "grm_stratified",
    reason      = "grm_kmeans_stratified",
    km_clusters = km$cluster,
    n_clusters  = n_clus
  )
}

setup_stratified_cv <- function(data, config) {
  
  set.seed(config$cv_seed)
  
  pheno <- data$phenotypes
  k <- config$k_folds
  
  if(config$stratify_by_sire && "father_id" %in% colnames(pheno)) {
    family_check <- check_family_structure(pheno, min_offspring = config$min_offspring_per_sire)
  } else {
    family_check <- list(feasible = FALSE, stats = list(n_sires = 0))
  }

  if(config$stratify_by_sire && family_check$feasible) {
    
    fold_assignments <- createFolds(
      y = pheno$father_id,
      k = k,
      list = FALSE,
      returnTrain = FALSE
    )
    
    stratification_check <- table(pheno$father_id, fold_assignments)
    
    list(
      folds = fold_assignments,
      method = "stratified",
      family_stats = family_check$stats,
      sires_per_fold = colSums(stratification_check > 0),
      avg_folds_per_sire = mean(rowSums(stratification_check > 0))
    )
    
  } else if (!is.null(config$cv_method) && config$cv_method == "cluster") {
    return(setup_cluster_cv(data, config, G = data$GRM))
  } else if (!is.null(config$cv_method) && config$cv_method == "grm_stratified") {
    return(setup_grm_stratified_cv(data, config, G = data$GRM))
  } else {
    fold_assignments <- sample(rep(1:k, length.out = nrow(pheno)))
    list(folds = fold_assignments, method = "random",
         family_stats = family_check$stats,
         reason = if(!config$stratify_by_sire) "stratification_disabled" else "no_pedigree_available")
  }
}

get_fold_split_info <- function(data, cv_setup, fold_id) {
  
  test_idx <- which(cv_setup$folds == fold_id)
  train_idx <- which(cv_setup$folds != fold_id)
  
  list(
    fold = fold_id,
    n_train = length(train_idx),
    n_test = length(test_idx),
    prop_train = length(train_idx) / nrow(data$phenotypes),
    prop_test = length(test_idx) / nrow(data$phenotypes)
  )
}