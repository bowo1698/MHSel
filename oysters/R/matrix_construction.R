# W matrix and GRM construction (MH path)
#
# W_ah from masbayes::construct_wah_matrix(), shared by Bayesian models
# (direct) and GBLUP (G = tcrossprod(W_ah)/k_ah) — no duplicate encoding.
construct_W_matrix <- function(hap_cols, allele_freq_filtered, reference_structure = NULL, drop_baseline = TRUE) {

  hap_matrix <- as.matrix(hap_cols)
  storage.mode(hap_matrix) <- "integer"
  colnames_vec <- colnames(hap_cols)

  allele_freq_list <- if (is.null(reference_structure)) {
    list(
      haplotype = as.character(allele_freq_filtered$haplotype),
      allele    = as.integer(allele_freq_filtered$allele),
      freq      = as.numeric(allele_freq_filtered$freq)
    )
  } else NULL

  ref_struct <- if (!is.null(reference_structure)) {
    list(
      allele_info = list(
        allele_id = as.character(reference_structure$allele_info$allele_id),
        freq      = as.numeric(reference_structure$allele_info$freq)
      ),
      dropped_alleles = if (!is.null(reference_structure$dropped_alleles) &&
                            nrow(reference_structure$dropped_alleles) > 0) {
        list(
          block  = as.character(reference_structure$dropped_alleles$block),
          allele = as.integer(reference_structure$dropped_alleles$allele),
          freq   = as.numeric(reference_structure$dropped_alleles$freq)
        )
      } else list()
    )
  } else NULL

  result <- masbayes::construct_wah_matrix(
    hap_matrix           = hap_matrix,
    colnames             = colnames_vec,
    allele_freq_filtered = allele_freq_list,
    reference_structure  = ref_struct,
    drop_baseline        = drop_baseline
  )

  if (is.null(colnames(result$W_ah)) &&
      !is.null(result$allele_info) &&
      "allele_id" %in% names(result$allele_info)) {
    colnames(result$W_ah) <- result$allele_info$allele_id
  }

  cat("W_ah final columns:", ncol(result$W_ah), "\n")
  cat("Metadata rows:", length(result$allele_info$allele_id), "\n")
  cat("Dropped alleles:", length(result$dropped_alleles$block), "\n")

  result
}

construct_matrices <- function(train_data, test_data, allele_freq, config) {

  hap_train <- train_data$haplotypes %>% select(-ID)
  hap_test  <- test_data$haplotypes %>% select(-ID)

  W_result_train <- construct_W_matrix(
    hap_train,
    allele_freq,
    drop_baseline = config$drop_baseline
  )

  cat("\n=== W_ah Matrix Preview (Training) ===\n")
  cat("Dimensions:", dim(W_result_train$W_ah), "\n")
  print(head(W_result_train$W_ah[, 1:min(10, ncol(W_result_train$W_ah))]))
  cat("\nFirst 10 column names:\n")
  print(colnames(W_result_train$W_ah)[1:10])
  cat("\n")

  W_result_test <- construct_W_matrix(
    hap_test,
    allele_freq_filtered = NULL,
    reference_structure  = W_result_train,
    drop_baseline        = config$drop_baseline
  )

  # Capture block_id BEFORE subset/arithmetic — both drop matrix attributes,
  # and run_bayesr(map = ...) needs attr(W, "block_id") to detect the MH path.
  block_id_orig <- attr(W_result_train$W_ah, "block_id")

  # Filter zero-variance columns
  col_var <- apply(W_result_train$W_ah, 2, var)
  keep    <- col_var > 0
  cat("Removing", sum(!keep), "zero-variance columns from", length(keep), "\n")

  W_result_train$W_ah <- W_result_train$W_ah[, keep]
  W_result_test$W_ah  <- W_result_test$W_ah[, keep]
  block_id_keep       <- if (!is.null(block_id_orig)) block_id_orig[keep] else NULL

  # Standardize W_ah
  vx_before <- apply(W_result_train$W_ah, 2, var)
  cat("W_ah BEFORE scaling | mean col-var:", round(mean(vx_before), 4),
      "| sum_vx:", round(sum(vx_before), 2), "\n")

  scale_factor <- sqrt(mean(vx_before))
  W_result_train$W_ah <- W_result_train$W_ah / scale_factor
  W_result_test$W_ah  <- W_result_test$W_ah  / scale_factor

  # Re-attach after subset + scale. Train and test share columns by
  # construction, so one vector serves both; test needs it for predict().
  if (!is.null(block_id_keep)) {
    stopifnot(length(block_id_keep) == ncol(W_result_train$W_ah))
    attr(W_result_train$W_ah, "block_id") <- block_id_keep
    attr(W_result_test$W_ah,  "block_id") <- block_id_keep
  }

  vx_after <- apply(W_result_train$W_ah, 2, var)
  cat("W_ah AFTER  scaling | mean col-var:", round(mean(vx_after), 4),
      "| sum_vx:", round(sum(vx_after), 2),
      "| scale_factor:", round(scale_factor, 4), "\n")

  if (!is.null(W_result_train$allele_info)) {
    W_result_train$allele_info <- lapply(W_result_train$allele_info, function(x) x[keep])
  }

  W_train <- W_result_train$W_ah
  W_test  <- W_result_test$W_ah

  # GBLUP relationship matrix from shared W_ah
  n_train <- nrow(W_train)
  G_train <- tcrossprod(W_train)
  k_ah    <- sum(diag(G_train)) / n_train
  G_train <- G_train / k_ah

  W_full <- rbind(W_train, W_test)
  G_full <- tcrossprod(W_full) / k_ah

  train_ids <- as.character(train_data$phenotypes$ID)
  test_ids  <- as.character(test_data$phenotypes$ID)
  rownames(G_train) <- colnames(G_train) <- train_ids
  rownames(G_full)  <- colnames(G_full)  <- c(train_ids, test_ids)

  cat("G_full dimension:", dim(G_full), "\n")
  cat("Mean diagonal:", mean(diag(G_full)), "\n")
  cat("Mean off-diagonal:", mean(G_full[upper.tri(G_full)]), "\n\n")

  list(
    W_train         = W_train,
    W_test          = W_test,
    allele_info     = W_result_train$allele_info,
    dropped_alleles = W_result_train$dropped_alleles,
    A_train         = G_train,
    A_combined      = G_full,
    k_ah            = k_ah,
    train_ids       = train_ids,
    test_ids        = test_ids,
    g_label         = "mh_add"
  )
}
