# MH data loading and preparation

load_genomic_data <- function(config) {
  
  # Load phenotypes
  pheno <- read_csv(
    file.path(config$data_dir, config$phenotype_file),
    show_col_types = FALSE
  ) %>%
    rename(phenotype = all_of(config$phenotype_trait)) %>%
    { df <- .
      cols_keep <- c("ID", "phenotype",
                     intersect(c(config$fixed_effect, config$extra_cols), names(df)))
      df[, cols_keep] }
  
  # Load haplotype genotypes
  hap_geno_list <- lapply(1:config$n_chromosomes, function(chr) {
    read_table(
      file.path(config$data_dir, config$block_dir, paste0("hap_geno_", chr)),
      show_col_types = FALSE
    )
  })
  
  # Merge chromosomes
  hap_geno <- hap_geno_list[[1]]
  for(chr in 2:config$n_chromosomes) {
    hap_geno <- bind_cols(
      hap_geno,
      hap_geno_list[[chr]] %>% select(-ID)
    )
  }
  
  # Keep individuals present in both datasets
  common_ids <- intersect(pheno$ID, hap_geno$ID)
  if (length(common_ids) < nrow(pheno)) {
    cat("Individuals in phenotype but not genotype:", nrow(pheno) - length(common_ids), "\n")
  }
  if (length(common_ids) < nrow(hap_geno)) {
    cat("Individuals in genotype but not phenotype:", nrow(hap_geno) - length(common_ids), "\n")
  }
  pheno <- pheno[!is.na(pheno$phenotype), ]
  common_ids <- intersect(pheno$ID, hap_geno$ID)
  pheno_ids <- pheno$ID[pheno$ID %in% common_ids]  # pheno order is the reference
  pheno    <- pheno    %>% filter(ID %in% common_ids) %>% arrange(match(ID, pheno_ids))
  hap_geno <- hap_geno %>% filter(ID %in% common_ids) %>% arrange(match(ID, pheno_ids))
  cat("Common individuals used:", length(common_ids), "\n\n")

  GRM <- NULL
  if (!is.null(config$cv_method) && config$cv_method %in% c("cluster", "grm_stratified")) {
    cat("Constructing GRM for cluster CV...\n")
    hap_cols    <- hap_geno %>% select(-ID)
    allele_freq <- compute_allele_frequencies(hap_cols, config)
    W_result    <- construct_W_matrix(hap_cols, allele_freq, drop_baseline = config$drop_baseline)
    col_var     <- apply(W_result$W_ah, 2, var)
    W_ah        <- W_result$W_ah[, col_var > 0]
    scale_f     <- sqrt(mean(apply(W_ah, 2, var)))
    W_ah        <- W_ah / scale_f
    k_ah        <- sum(diag(tcrossprod(W_ah))) / nrow(W_ah)
    GRM         <- tcrossprod(W_ah) / k_ah
    cat("GRM constructed:", dim(GRM), "\n\n")
  }

  list(
    phenotypes = pheno,
    haplotypes = hap_geno,
    GRM        = GRM
  )
}

split_fold_data <- function(data, fold_assignments, fold_id) {
  
  test_idx <- which(fold_assignments == fold_id)
  train_idx <- which(fold_assignments != fold_id)
  
  pheno_train <- data$phenotypes[train_idx, ]
  pheno_test <- data$phenotypes[test_idx, ]
  
  hap_train <- data$haplotypes %>%
    filter(ID %in% pheno_train$ID) %>%
    arrange(match(ID, pheno_train$ID))

  hap_test <- data$haplotypes %>%
    filter(ID %in% pheno_test$ID) %>%
    arrange(match(ID, pheno_test$ID))
  
  cat("\n=== Genotype Data Preview ===\n")
  print(head(hap_train[, 1:min(10, ncol(hap_train))]))
  cat("\nDimensions:", dim(hap_train), "\n\n")
  
  list(
    train = list(phenotypes = pheno_train, haplotypes = hap_train),
    test = list(phenotypes = pheno_test, haplotypes = hap_test)
  )
}

compute_allele_frequencies <- function(hap_cols, config) {
  
  # Pool both copies per block, stripping the _1 suffix
  col_names <- colnames(hap_cols)
  block_names <- sub("_1$", "", col_names)
  hap_cols_pooled <- hap_cols
  colnames(hap_cols_pooled) <- block_names
  
  allele_freq <- hap_cols_pooled %>%
    pivot_longer(everything(), names_to = "haplotype", values_to = "allele") %>%
    count(haplotype, allele) %>%
    group_by(haplotype) %>%
    mutate(freq = n / sum(n)) %>%
    select(-n)
  
  allele_freq_filtered <- allele_freq %>%
    filter(freq >= config$min_allele_freq & freq <= config$max_allele_freq) %>%
    ungroup()

  cat("Total alleles before filter:", nrow(allele_freq), "\n")
  cat("Total alleles after filter:", nrow(allele_freq_filtered), "\n")
  cat("Alleles removed:", nrow(allele_freq) - nrow(allele_freq_filtered), "\n\n")
  
  allele_freq_filtered
}