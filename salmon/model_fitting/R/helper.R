# Allele freqs
compute_allele_frequencies <- function(hap_cols, config) {

  # Pool both copy per blok with strip suffix _1
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

# W Matrix and Genomic Relationship Matrix Construction
construct_W_matrix <- function(hap_cols, allele_freq_filtered, reference_structure = NULL, drop_baseline = TRUE) {
  
  # Convert to matrix format
  hap_matrix <- as.matrix(hap_cols)
  storage.mode(hap_matrix) <- "integer"
  colnames_vec <- colnames(hap_cols)
  
  # Prepare allele_freq for Rust (training set only)
  allele_freq_list <- if(is.null(reference_structure)) {
    list(
      haplotype = as.character(allele_freq_filtered$haplotype),
      allele = as.integer(allele_freq_filtered$allele),
      freq = as.numeric(allele_freq_filtered$freq)
    )
  } else {
    NULL
  }
  
  # Prepare reference structure for Rust (test set only)
  ref_struct <- if(!is.null(reference_structure)) {
    list(
      allele_info = list(
        allele_id = as.character(reference_structure$allele_info$allele_id),
        freq = as.numeric(reference_structure$allele_info$freq)
      ),
      dropped_alleles = if(!is.null(reference_structure$dropped_alleles) && 
                           nrow(reference_structure$dropped_alleles) > 0) {
        list(
          block = as.character(reference_structure$dropped_alleles$block),
          allele = as.integer(reference_structure$dropped_alleles$allele),
          freq = as.numeric(reference_structure$dropped_alleles$freq)
        )
      } else {
        list()
      }
    )
  } else {
    NULL
  }
  
  # Call Rust implementation
  result <- masbayes::construct_wah_matrix(
    hap_matrix = hap_matrix,
    colnames = colnames_vec,
    allele_freq_filtered = allele_freq_list,
    reference_structure = ref_struct,
    drop_baseline = drop_baseline
  )
  
  # Ensure column names are set
  if(is.null(colnames(result$W_ah)) && 
     !is.null(result$allele_info) && 
     "allele_id" %in% names(result$allele_info)) {
    colnames(result$W_ah) <- result$allele_info$allele_id
  }
  
  # Print summary
  cat("W_ah final columns:", ncol(result$W_ah), "\n")
  cat("Metadata rows:", length(result$allele_info$allele_id), "\n")
  cat("Dropped alleles:", length(result$dropped_alleles$block), "\n")
  
  result
}

suppressPackageStartupMessages({
  library(tidyverse)
})

quality_control <- function(geno_matrix, config) {
  if(!is.numeric(geno_matrix)) {
    stop("Genotype matrix must be numeric. Check if ID column was removed.")
  }
  
  allele_freq <- colMeans(geno_matrix, na.rm = TRUE) / 2
  maf <- pmin(allele_freq, 1 - allele_freq)
  callrate <- colMeans(!is.na(geno_matrix))
  
  keep_snps <- (maf >= config$min_allele_freq) & 
               (callrate >= config$callrate_threshold)
  
  cat("=== Quality Control ===\n")
  cat("SNPs retained:", sum(keep_snps), "out of", ncol(geno_matrix), "\n")
  
  list(
    filtered = geno_matrix[, keep_snps],
    keep_snps = keep_snps,
    n_retained = sum(keep_snps)
  )
}

impute_and_standardize <- function(geno_filtered) {
  geno_imputed <- apply(geno_filtered, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    x
  })

  snp_means <- colMeans(geno_imputed)
  snp_sds <- apply(geno_imputed, 2, sd)
  zero_var <- snp_sds < 1e-10

  if(sum(zero_var) > 0) {
    cat("Removing", sum(zero_var), "zero-variance SNPs\n")
    geno_imputed <- geno_imputed[, !zero_var]
    snp_means <- snp_means[!zero_var]
    snp_sds <- snp_sds[!zero_var]
  }
  
  geno_std <- scale(geno_imputed, center = snp_means, scale = snp_sds)

  if(any(is.na(geno_std))) {
    stop("Standardization produced NA values - check for zero variance SNPs")
  }
  
  cat("Clean matrix:", ncol(geno_std), "SNPs\n\n")
  
  list(
    standardized = geno_std,
    imputed = geno_imputed
  )
}

enrich_allele_info <- function(allele_info, coords_file) {
  coords <- read_csv(coords_file, show_col_types = FALSE) %>%
    mutate(chr = as.integer(chr)) %>%
    mutate(global_block_id = block_id) %>%
    group_by(chr) %>%
    mutate(block_num = row_number()) %>%
    ungroup() %>%
    select(chr, block_num, global_block_id, start_pos, end_pos, n_snps)

  allele_info %>%
    mutate(
      chr       = as.integer(str_extract(block, "(?<=hap_)\\d+")),
      block_num = as.integer(str_extract(block, "\\d+$"))
    ) %>%
    left_join(coords, by = c("chr", "block_num")) %>%
    select(-block_num)
}

# QTL Recovery Evaluation
##
## CHANGE 2026-05-20: positional matching for MH×MH-block path
##   Previously: `top_blocks %in% qtl_mh_blocks` (block_id string match).
##   That breaks when G1 block definition is re-derived fresh (block_id namespace
##   differs G0 ↔ G1). New approach: did fitted MH block contain any of the
##   *causal SNPs* (positions known from qtl_s{3,4,5,6}_effects)? Consistent
##   with the 5 other paths that already use positional matching.
##   New param `qtl_block_causal_snps` — character vector of causal SNP IDs
##   for the current scenario's block QTLs (extracted by caller from
##   qtl_sN_effects). NULL → fallback to old string-match behaviour.
compute_qtl_recovery <- function(marker_info, sc, qtl_snp_ids, qtl_mh_blocks,
                                  qtl_hap_blocks, snp_meta, marker_type,
                                  mh_coords = NULL, hap_coords = NULL,
                                  qtl_block_causal_snps = NULL,
                                  top_k = c(50, 100, 200)) {
  if (is.null(marker_info)) return(NULL)

  # Define QTL 
  qtl_unit <- dplyr::case_when(
    sc %in% c("S1_equal", "S2_unequal") ~ "snp",
    sc %in% c("S3_equal", "S4_unequal") ~ "mh_block",
    sc %in% c("S5_equal", "S6_unequal") ~ "hap_block",
    TRUE ~ "unknown"
  )

  results <- lapply(top_k, function(k) {
    # Rank marker 
    if ("PIP" %in% names(marker_info)) {
      ranked <- marker_info %>% arrange(desc(PIP))
    } else {
      ranked <- marker_info %>% arrange(desc(abs(beta_hat)))
    }
    top_markers <- ranked %>% slice_head(n = k)

    if (qtl_unit == "snp") {
      # SNP marker vs SNP QTL — direct match
      if (marker_type == "snp") {
        n_hit <- sum(top_markers$allele_id %in% qtl_snp_ids)
      } else {
        # MH marker — check
        qtl_snp_pos <- snp_meta %>% filter(snp_id %in% qtl_snp_ids)
        top_blocks  <- top_markers %>% select(block, chr, start_pos, end_pos) %>% distinct()
        n_hit <- nrow(top_blocks %>%
          inner_join(qtl_snp_pos, by = "chr", relationship = "many-to-many") %>%
          filter(pos >= start_pos & pos <= end_pos) %>%
          distinct(block))
      }
      n_qtl <- length(qtl_snp_ids)

    } else if (qtl_unit == "mh_block") {
      if (marker_type == "mh") {
        ## NEW 2026-05-20: positional overlap with causal SNP positions.
        ## Fallback to legacy block-id string match if qtl_block_causal_snps not provided.
        block_scores <- marker_info %>%
          group_by(global_block_id) %>%
          summarise(block_score = sum(beta_hat^2),
                    chr        = dplyr::first(chr),
                    start_pos  = dplyr::first(start_pos),
                    end_pos    = dplyr::first(end_pos),
                    .groups    = "drop") %>%
          arrange(desc(block_score))
        top_blocks_df <- block_scores %>% slice_head(n = k)

        if (!is.null(qtl_block_causal_snps)) {
          causal_pos <- snp_meta %>%
            filter(snp_id %in% qtl_block_causal_snps) %>%
            select(chr, pos) %>%
            mutate(chr = as.integer(chr))
          top_blocks_df <- top_blocks_df %>% mutate(chr = as.integer(chr))
          n_hit <- nrow(top_blocks_df %>%
            inner_join(causal_pos, by = "chr", relationship = "many-to-many") %>%
            filter(pos >= start_pos & pos <= end_pos) %>%
            distinct(global_block_id))
        } else {
          ## Legacy fallback (only valid when block_id namespace shared G0↔G1, e.g. via --reuse-blocks)
          n_hit <- sum(top_blocks_df$global_block_id %in% qtl_mh_blocks)
        }
      } else if (marker_type == "snp" && !is.null(mh_coords)) {
        qtl_block_ranges <- mh_coords %>%
          filter(block_id %in% qtl_mh_blocks) %>%
          select(chr, start_pos, end_pos) %>%
          mutate(chr = as.integer(chr))
        top_snps <- top_markers %>%
          select(allele_id, chr, start_pos) %>%
          rename(snp_pos = start_pos)
        n_hit <- nrow(top_snps %>%
          inner_join(qtl_block_ranges, by = "chr", relationship = "many-to-many") %>%
          filter(snp_pos >= start_pos & snp_pos <= end_pos) %>%
          distinct(allele_id))
      } else {
        n_hit <- NA
      }
      n_qtl <- length(qtl_mh_blocks)

    } else if (qtl_unit == "hap_block") {
      if (!is.null(hap_coords)) {
        qtl_hap_ranges <- hap_coords %>%
          filter(block_id %in% qtl_hap_blocks) %>%
          select(chr, start_pos, end_pos) %>%
          mutate(chr = as.integer(chr))

        if (marker_type == "snp") {
          top_snps <- top_markers %>%
            select(allele_id, chr, start_pos) %>%
            rename(snp_pos = start_pos)
          n_hit <- nrow(top_snps %>%
            inner_join(qtl_hap_ranges, by = "chr", relationship = "many-to-many") %>%
            filter(snp_pos >= start_pos & snp_pos <= end_pos) %>%
            distinct(allele_id))

        } else if (marker_type == "mh") {
          top_blocks <- top_markers %>%
            select(global_block_id, chr, start_pos, end_pos) %>%
            distinct()
          n_hit <- nrow(top_blocks %>%
            inner_join(qtl_hap_ranges, by = "chr", relationship = "many-to-many") %>%
            filter(start_pos.x >= start_pos.y & end_pos.x <= end_pos.y) %>%
            distinct(global_block_id))
        } else {
          n_hit <- NA
        }
      } else {
        n_hit <- NA
      }
      n_qtl <- length(qtl_hap_blocks)

    } else {
      n_hit <- NA
      n_qtl <- NA
    }

    tibble(
      top_k        = k,
      n_hit        = n_hit,
      n_qtl        = n_qtl,
      recovery_rate = if (!is.na(n_hit)) n_hit / n_qtl else NA_real_,
      qtl_unit     = qtl_unit
    )
  })

  bind_rows(results)
}

# ---------------------------------------------------------------------------
# MCMC diagnostics helpers
# ---------------------------------------------------------------------------

# Polyfill for %||% (base since R 4.4.0). Safe to redefine; identical semantics.
`%||%` <- function(a, b) if (is.null(a)) b else a

# MCMC diagnostics -- moving-window Geweke z-score.
# Split each scalar chain into k contiguous windows of equal length.
# For each window, compute Geweke z (10% first vs 50% last within the window).
# Returns NULL when n_kept < k * 10 (insufficient samples).
compute_geweke_moving <- function(scalar_chains, k = 20L) {
  n_kept <- length(scalar_chains[[1L]])
  if (n_kept < k * 10L) {
    warning("n_kept (", n_kept, ") too small for k_window = ", k,
            " -- geweke_window skipped")
    return(NULL)
  }
  window_size <- floor(n_kept / k)
  out <- data.frame(window_idx = seq_len(k))
  for (nm in names(scalar_chains)) {
    z <- numeric(k)
    for (i in seq_len(k)) {
      start_i <- (i - 1L) * window_size + 1L
      end_i   <- min(i * window_size, n_kept)
      sub     <- scalar_chains[[nm]][start_i:end_i]
      z[i]    <- coda::geweke.diag(coda::as.mcmc(sub))$z
    }
    out[[paste0(nm, "_z")]] <- z
  }
  out
}

# MCMC diagnostics -- summary stats per chain.
# scalar_chains: list of numeric vectors (one per scalar VC parameter).
# pi_matrix: matrix n_kept x 4 for BayesR's pi chain, or NULL.
compute_mcmc_summary <- function(scalar_chains, pi_matrix = NULL,
                                 k_window = 20L) {

  ess <- lapply(scalar_chains, function(x) {
    unname(coda::effectiveSize(coda::as.mcmc(x)))
  })
  geweke_z <- lapply(scalar_chains, function(x) {
    unname(coda::geweke.diag(coda::as.mcmc(x))$z)
  })
  hw <- lapply(scalar_chains, function(x) {
    # heidel.diag aborts on degenerate chains (I=Inf -> pcramer(Inf)=NaN); NA-fallback keeps the fit
    h <- tryCatch(coda::heidel.diag(coda::as.mcmc(x)), error = function(e) NULL)
    if (is.null(h)) return(list(pvalue = NA_real_, passed = NA))
    list(pvalue = unname(h[1L, "pvalue"]),
         passed = as.logical(h[1L, "stest"]))
  })

  geweke_window <- compute_geweke_moving(scalar_chains, k = k_window)

  posterior_mean <- lapply(scalar_chains, mean, na.rm = TRUE)
  posterior_q025 <- lapply(scalar_chains, function(x) {
    unname(quantile(x, probs = 0.025, na.rm = TRUE))
  })
  posterior_q975 <- lapply(scalar_chains, function(x) {
    unname(quantile(x, probs = 0.975, na.rm = TRUE))
  })

  if (!is.null(pi_matrix)) {
    ess$pi <- apply(pi_matrix, 2L, function(x) {
      unname(coda::effectiveSize(coda::as.mcmc(x)))
    })
    geweke_z$pi <- apply(pi_matrix, 2L, function(x) {
      unname(coda::geweke.diag(coda::as.mcmc(x))$z)
    })
    hw$pi <- lapply(seq_len(ncol(pi_matrix)), function(j) {
      # heidel.diag aborts on degenerate chains (I=Inf -> pcramer(Inf)=NaN); NA-fallback keeps the fit
      h <- tryCatch(coda::heidel.diag(coda::as.mcmc(pi_matrix[, j])),
                    error = function(e) NULL)
      if (is.null(h)) return(list(pvalue = NA_real_, passed = NA))
      list(pvalue = unname(h[1L, "pvalue"]),
           passed = as.logical(h[1L, "stest"]))
    })
    posterior_mean$pi <- apply(pi_matrix, 2L, mean, na.rm = TRUE)
    posterior_q025$pi <- apply(pi_matrix, 2L, function(x) {
      unname(quantile(x, probs = 0.025, na.rm = TRUE))
    })
    posterior_q975$pi <- apply(pi_matrix, 2L, function(x) {
      unname(quantile(x, probs = 0.975, na.rm = TRUE))
    })
  }

  list(
    ess             = ess,
    geweke_z        = geweke_z,
    geweke_window   = geweke_window,
    hw_stationarity = hw,
    posterior_mean  = posterior_mean,
    posterior_q025  = posterior_q025,
    posterior_q975  = posterior_q975
  )
}

# MCMC diagnostics -- build mcmc_diag list dari masbayes res object.
# Returns NULL kalau config$mcmc_diag$enabled bukan TRUE.
build_mcmc_diag <- function(res, config, meta) {
  if (!isTRUE(config$mcmc_diag$enabled)) return(NULL)
  if (is.null(res$sigma2_e_samples)) {
    warning("res$sigma2_e_samples missing -- mcmc_diag skipped")
    return(NULL)
  }
  if (is.null(res$mu_samples)) {
    warning("res$mu_samples missing -- mcmc_diag skipped")
    return(NULL)
  }

  model <- meta$model
  k_window <- as.integer(config$mcmc_diag$k_window %||% 20L)

  chains <- list(
    sigma2_e = as.numeric(res$sigma2_e_samples),
    mu       = as.numeric(res$mu_samples)
  )

  pi_matrix <- NULL
  if (model == "bayesR") {
    chains$sigma2_small  <- as.numeric(res$sigma2_small_samples)
    chains$sigma2_medium <- as.numeric(res$sigma2_medium_samples)
    chains$sigma2_large  <- as.numeric(res$sigma2_large_samples)
    if (!is.matrix(res$pi_samples) || ncol(res$pi_samples) != 4L) {
      stop("res$pi_samples must be a matrix with 4 columns (got ",
           paste(dim(res$pi_samples), collapse = " x "), ")")
    }
    pi_matrix <- res$pi_samples
    sigma2_g_chain <- pi_matrix[, 2L] * chains$sigma2_small +
      pi_matrix[, 3L] * chains$sigma2_medium +
      pi_matrix[, 4L] * chains$sigma2_large
    chains$h2 <- sigma2_g_chain / (sigma2_g_chain + chains$sigma2_e)
  } else if (model == "bayesA") {
    if (is.null(res$sigma2_j_samples)) {
      warning("res$sigma2_j_samples missing for bayesA -- mcmc_diag skipped")
      return(NULL)
    }
    mean_sigma2_j_chain <- rowMeans(res$sigma2_j_samples)
    chains$h2 <- mean_sigma2_j_chain /
      (mean_sigma2_j_chain + chains$sigma2_e)
  } else {
    stop("Unknown model: ", model)
  }

  summary_out <- compute_mcmc_summary(
    scalar_chains = chains,
    pi_matrix     = pi_matrix,
    k_window      = k_window
  )

  out_chains <- chains
  if (!is.null(pi_matrix)) out_chains$pi <- pi_matrix

  list(
    meta    = meta,
    chains  = out_chains,
    summary = summary_out
  )
}