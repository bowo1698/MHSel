suppressPackageStartupMessages({
  library(tidyverse)
  library(sommer)
  library(coda)
  library(masbayes)
  library(masreml)
  library(RhpcBLASctl)
})

select <- dplyr::select
filter <- dplyr::filter
blas_set_num_threads(1)

config_file <- Sys.getenv("GENOMIC_CONFIG", "config.R")
source(config_file)
setwd(config$base_dir)

source("R/helper.R")
source("R/models/gblup.R")
source("R/models/bayesR.R")

pipeline_start <- Sys.time()
cat(sprintf("=== Full-data Pipeline started: %s ===\n", Sys.time()))

# ── 1. Load genotype ──────────────────────────────────────────
load_hap_geno <- function(data_dir, n_chr, block_dir) {
  geno_list <- lapply(1:n_chr, function(chr) {
    read_table(
      file.path(data_dir, block_dir, paste0("hap_geno_", chr)),
      show_col_types = FALSE
    )
  })
  geno <- geno_list[[1]]
  for (chr in 2:n_chr)
    geno <- bind_cols(geno, geno_list[[chr]] %>% select(-ID))
  geno %>% rename(individual_id = ID)
}

hap_geno_current <- load_hap_geno(config$data_dir, config$n_chromosomes, config$block_dir)

cat(sprintf("Genotype loaded: %d individuals\n", nrow(hap_geno_current)))

# ── 2. All individuals used for fitting ───────────────────────
pedigree_current <- read_csv(
  file.path(config$data_dir, config$pedigree),
  show_col_types = FALSE
) %>% filter(population_type == config$population_type)

train_ids <- pedigree_current$individual_id
test_ids  <- character(0)

cat(sprintf("Full-data fit: train=%d (no test split)\n", length(train_ids)))

# ── 3. Construct W matrix ─────────────────────────────────────────────────────
# Full-data W matrix
hap_all           <- hap_geno_current %>%
                     arrange(match(individual_id, c(train_ids, test_ids)))
allele_freq_full  <- compute_allele_frequencies(hap_all %>% select(-individual_id), config)
W_result_full     <- construct_W_matrix(hap_all %>% select(-individual_id),
                                        allele_freq_full, drop_baseline = config$drop_baseline)
# Capture block_id BEFORE [, keep_full] and division
# attr(W, "block_id"); the attr must be present when run_bayesr(map = ...)
# is called downstream in the scenario loop.
block_id_full_orig <- attr(W_result_full$W_ah, "block_id")
col_var_full      <- apply(W_result_full$W_ah, 2, var)
keep_full         <- col_var_full > 0
W_result_full$W_ah <- W_result_full$W_ah[, keep_full]
block_id_full_keep <- if (!is.null(block_id_full_orig)) block_id_full_orig[keep_full] else NULL
scale_factor_full <- sqrt(mean(apply(W_result_full$W_ah, 2, var)))
W_full            <- W_result_full$W_ah / scale_factor_full
rownames(W_full)  <- hap_all$individual_id
# Re-attach after subset + scale; keep_full mask is shared between W and
# block_id, so alignment with ncol(W_full) is guaranteed by construction.
if (!is.null(block_id_full_keep)) {
  stopifnot(length(block_id_full_keep) == ncol(W_full))
  attr(W_full, "block_id") <- block_id_full_keep
}
n_full            <- nrow(W_full)
A_full            <- tcrossprod(W_full)
k_ah_full         <- sum(diag(A_full)) / n_full
A_full            <- A_full / k_ah_full

diag_vals  <- diag(A_full)
off_vals   <- A_full[lower.tri(A_full)]
grm_stats  <- list(
  mean_diag    = mean(diag_vals),
  sd_diag      = sd(diag_vals),
  mean_offdiag = mean(off_vals),
  sd_offdiag   = sd(off_vals)
)

coords_file       <- file.path(config$data_dir,
                               "mh_info_ld_haploblock/stats/microhaplotype_coordinates.csv")
allele_info_full_keep <- enrich_allele_info(
  data.frame(W_result_full$allele_info)[keep_full, ], coords_file
)

hap_train <- hap_geno_current %>%
  filter(individual_id %in% train_ids) %>%
  arrange(match(individual_id, train_ids))
hap_test <- hap_geno_current %>%
  filter(individual_id %in% test_ids) %>%
  arrange(match(individual_id, test_ids))

allele_freq <- compute_allele_frequencies(hap_train %>% select(-individual_id), config)

W_result_train <- construct_W_matrix(hap_train %>% select(-individual_id),
                                     allele_freq, drop_baseline = config$drop_baseline)
W_result_test  <- construct_W_matrix(hap_test %>% select(-individual_id),
                                     allele_freq_filtered = NULL,
                                     reference_structure  = W_result_train,
                                     drop_baseline        = config$drop_baseline)

col_var <- apply(W_result_train$W_ah, 2, var)
keep    <- col_var > 0
W_result_train$W_ah <- W_result_train$W_ah[, keep]
W_result_test$W_ah  <- W_result_test$W_ah[, keep]

scale_factor        <- sqrt(mean(apply(W_result_train$W_ah, 2, var)))
W_train <- W_result_train$W_ah / scale_factor
W_test  <- W_result_test$W_ah  / scale_factor
rownames(W_train) <- train_ids
rownames(W_test)  <- test_ids
n_train <- nrow(W_train)

G_train  <- tcrossprod(W_train)
k_ah     <- sum(diag(G_train)) / n_train
A_train  <- G_train / k_ah

W_combined <- rbind(W_train, W_test)
G_combined <- tcrossprod(W_combined)
A_combined <- G_combined / k_ah
all_ids    <- c(train_ids, test_ids)
rownames(A_combined) <- colnames(A_combined) <- all_ids

allele_info_keep <- enrich_allele_info(
  data.frame(W_result_train$allele_info)[keep, ], coords_file
)

cat(sprintf("W matrix ready: %d train x %d alleles | k_ah=%.4f | scale_factor=%.4f\n",
            n_train, ncol(W_train), k_ah, scale_factor))

# Load QTL
qtl_snp_ids    <- readRDS(file.path(config$qtl_dir, "qtl_snp_ids.rds"))
qtl_mh_blocks  <- readRDS(file.path(config$qtl_dir, "qtl_mh_blocks.rds"))
qtl_hap_blocks <- readRDS(file.path(config$qtl_dir, "qtl_hap_blocks.rds"))
snp_meta       <- readRDS(file.path(config$qtl_dir, "snp_meta.rds"))
hap_coords     <- read_csv(config$hap_coords_file, show_col_types = FALSE) %>%
                 mutate(chr = as.integer(chr))
## REDESIGN 2026-05-22: block QTL is now HAPLOTYPE-ADDITIVE (not dosage epistasis).
## Each list element is list(chr, snp_ids, snp_sorted_idx, causal_alleles,
##   allele_freqs, beta) — we use $snp_ids (all block SNPs) for positional
## recovery in compute_qtl_recovery (MH×MH-block path). Pre-2026-05-22 schema
## was list(chr, snp_ids=K-causal-SNPs, beta=K-betas, gamma=scalar).
qtl_s3_effects <- readRDS(file.path(config$qtl_dir, "qtl_s3_effects.rds"))
qtl_s4_effects <- readRDS(file.path(config$qtl_dir, "qtl_s4_effects.rds"))
qtl_s5_effects <- readRDS(file.path(config$qtl_dir, "qtl_s5_effects.rds"))
qtl_s6_effects <- readRDS(file.path(config$qtl_dir, "qtl_s6_effects.rds"))
extract_causal_snps <- function(effects_list) {
  unique(unlist(lapply(effects_list, `[[`, "snp_ids")))
}
qtl_block_effects_lookup <- list(
  S3_equal   = qtl_s3_effects,
  S4_unequal = qtl_s4_effects,
  S5_equal   = qtl_s5_effects,
  S6_unequal = qtl_s6_effects
)
cat(sprintf("QTL loaded: %d SNP QTL | %d MH block QTL | %d HAP block QTL\n",
            length(qtl_snp_ids), length(qtl_mh_blocks), length(qtl_hap_blocks)))

# ── GWAS map (MH): per-iter coords from this iter's mh_info_ld_haploblock ────
gwas_enabled <- isTRUE(config$gwas$enabled)
if (gwas_enabled) {
  block_id_W_full <- attr(W_full, "block_id")
  if (is.null(block_id_W_full)) {
    stop("attr(W_full, 'block_id') is NULL — masbayes v0.5.0 required for GWAS.")
  }

  coords_raw <- read_csv(config$mh_coords_file, show_col_types = FALSE) %>%
    select(block_id, chr, start_pos, end_pos) %>%
    mutate(chr = as.integer(chr))

  unique_w_blocks <- unique(block_id_W_full)
  w_block_chr     <- as.integer(sub("^hap_(\\d+)_\\d+$", "\\1", unique_w_blocks))
  if (any(is.na(w_block_chr))) {
    stop("Cannot parse chromosome from W block_id values (expected 'hap_<chr>_<block>').\n  Example values: ",
         paste(head(unique_w_blocks[is.na(w_block_chr)], 3), collapse = ", "))
  }
  translation <- character(length(unique_w_blocks))
  for (chr_i in unique(w_block_chr)) {
    w_in_chr      <- which(w_block_chr == chr_i)
    coords_in_chr <- which(coords_raw$chr == chr_i)
    if (length(w_in_chr) != length(coords_in_chr)) {
      stop(sprintf(
        "Chr %d block count mismatch: W has %d unique blocks but coords CSV has %d rows. Position-based translation requires identical per-chr counts.",
        chr_i, length(w_in_chr), length(coords_in_chr)))
    }
    translation[w_in_chr] <- coords_raw$block_id[coords_in_chr]
  }
  names(translation) <- unique_w_blocks

  attr(W_full, "block_id") <- unname(translation[block_id_W_full])
  block_id_W_full <- attr(W_full, "block_id")

  map_mh_gwas <- coords_raw %>%
    filter(block_id %in% unique(block_id_W_full))
  cat(sprintf("GWAS enabled | W unique blocks: %d | map blocks: %d | windsize: %.0f bp\n",
              length(unique_w_blocks), nrow(map_mh_gwas), config$gwas$windsize))
} else {
  map_mh_gwas      <- NULL
  block_id_W_full  <- NULL
  cat("GWAS disabled (config$gwas$enabled != TRUE)\n")
}

# ── 4. Loop scenario ──────────────────────────────────────────────────────────
for (sc in names(config$scenarios)) {

  out_dir    <- file.path(config$output_dir, sc)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sc_start   <- Sys.time()

  # Load pheno
  pheno_all <- read_csv(
    file.path(config$data_dir, config$scenarios[[sc]]),
    show_col_types = FALSE
  )

  pheno_train <- pheno_all %>%
    filter(individual_id %in% train_ids) %>%
    arrange(match(individual_id, train_ids))

  pheno_test <- pheno_all %>%
    filter(individual_id %in% test_ids) %>%
    arrange(match(individual_id, test_ids))

  cat(sprintf("\n[%s] Scenario: %s | train=%d | test=%d\n",
              sc, Sys.time(), nrow(pheno_train), nrow(pheno_test)))
  
  # ── Full-data fitting ──
  pheno_full     <- pheno_all %>% arrange(match(individual_id, hap_all$individual_id))
  y_full         <- pheno_full$phenotype
  matrices_full  <- list(
    W_train         = W_full,
    W_test          = W_full[0, , drop = FALSE],
    allele_info     = allele_info_full_keep,
    dropped_alleles = W_result_full$dropped_alleles,
    A_train         = A_full,
    A_combined      = A_full,
    k_ah            = k_ah_full,
    train_ids       = as.character(hap_all$individual_id),
    test_ids        = character(0),
    g_label         = "mh_add"
  )
  split_full <- list(
    train = list(phenotypes = pheno_full, haplotypes = hap_all),
    test  = list(phenotypes = pheno_full[0, ], haplotypes = hap_all[0, ])
  )

  cat("  Running full-data fitting...\n")
  t0_gblup  <- proc.time()
  gblup_full <- run_gblup(matrices_full, split_full, config)
  rt_gblup  <- (proc.time() - t0_gblup)[["elapsed"]]
  cat(sprintf("    [full] gblup  | acc=%.4f bias=%.4f h2=%.3f\n",
    cor(gblup_full$predictions$train, pheno_full$tbv, use="complete.obs"),
    coef(lm(pheno_full$tbv ~ gblup_full$predictions$train))[2],
    gblup_full$varcomp$h2))

  t0_bayesR  <- proc.time()
  bayesR_full <- run_bayesR(
    matrices_full, split_full, gblup_full$varcomp, config, fold = 0,
    map      = if (gwas_enabled) map_mh_gwas              else NULL,
    windsize = if (gwas_enabled) config$gwas$windsize     else NULL
  )
  rt_bayesR  <- (proc.time() - t0_bayesR)[["elapsed"]]
  cat(sprintf("    [full] bayesR | acc=%.4f bias=%.4f h2=%.3f\n",
    cor(bayesR_full$predictions$train, pheno_full$tbv, use="complete.obs"),
    coef(lm(pheno_full$tbv ~ bayesR_full$predictions$train))[2],
    bayesR_full$varcomp$h2))

  # ── GWAS Rds (separate from results_all.rds per Rds-contract policy) ──
  # max-per-block PIP: more discriminative than masbayes pip_block on small
  # blocks. pip_block is a union statistic (P(any allele has gamma != 0))
  # which inflates with block size k roughly as 1 - (1-p)^k; max-per-block
  # is bounded by max(PIP) over alleles. See masbayes/examples/04_gwas.R
  # lines 130-138 for the demo discussion.
  if (gwas_enabled && !is.null(bayesR_full$pip)) {
    max_pip_block <- vapply(map_mh_gwas$block_id, function(b) {
      idx <- which(block_id_W_full == b)
      if (length(idx) == 0L) NA_real_ else max(bayesR_full$pip[idx])
    }, numeric(1L))

    iter_label <- basename(config$data_dir)
    gen_label  <- sub("^output_", "", basename(dirname(config$data_dir)))

    saveRDS(
      list(
        map           = map_mh_gwas,
        pip           = bayesR_full$pip,
        pip_block     = bayesR_full$pip_block,
        max_pip_block = max_pip_block,
        pip_class     = bayesR_full$pip_class,
        gwas          = bayesR_full$gwas,
        gwas_meta     = bayesR_full$gwas_meta,
        marker_type   = "mh",
        windsize      = config$gwas$windsize,
        scenario      = sc,
        iter          = iter_label,
        generation    = gen_label
      ),
      file.path(out_dir, "gwas_full.rds")
    )
    cat(sprintf("    [full] gwas   | saved: %s\n",
                file.path(out_dir, "gwas_full.rds")))
  }

  full_results <- list(
    gblup  = list(GEBV = gblup_full$predictions$train,  h2 = gblup_full$varcomp$h2),
    bayesR = list(GEBV = bayesR_full$predictions$train, h2 = bayesR_full$varcomp$h2)
  )

  # ── Save full-data fit per model ──
  for (model_name in c("gblup", "bayesR")) {
    dir.create(file.path(out_dir, model_name),
               recursive = TRUE, showWarnings = FALSE)

    GEBV <- full_results[[model_name]]$GEBV
    h2   <- full_results[[model_name]]$h2
    acc  <- cor(GEBV, pheno_full$tbv, use = "complete.obs")
    bias <- coef(lm(pheno_full$tbv ~ GEBV))[2]
    rmse <- sqrt(mean((GEBV - pheno_full$tbv)^2, na.rm = TRUE))

    results_all <- tibble(
      individual_id    = pheno_full$individual_id,
      tbv              = pheno_full$tbv,
      GEBV             = GEBV,
      scenario         = sc,
      model            = model_name,
      marker           = "mh",
      accuracy         = acc,
      bias             = bias,
      rmse             = rmse,
      h2               = h2,
      n_train          = nrow(pheno_full),
      grm_mean_diag    = grm_stats$mean_diag,
      grm_sd_diag      = grm_stats$sd_diag,
      grm_mean_offdiag = grm_stats$mean_offdiag,
      grm_sd_offdiag   = grm_stats$sd_offdiag
    )

    markers <- switch(model_name,
      gblup  = NULL,
      bayesR = bayesR_full$markers
    )
    causal_snps_sc <- if (sc %in% names(qtl_block_effects_lookup))
                        extract_causal_snps(qtl_block_effects_lookup[[sc]]) else NULL
    qtl_recovery <- compute_qtl_recovery(markers, sc, qtl_snp_ids,
                                         qtl_mh_blocks, qtl_hap_blocks,
                                         snp_meta, "mh",
                                         mh_coords = NULL,
                                         hap_coords = hap_coords,
                                         qtl_block_causal_snps = causal_snps_sc)
    diagnostics_fit <- switch(model_name,
      gblup  = NULL,
      bayesR = bayesR_full$diagnostics
    )
    runtime_fit <- switch(model_name,
      gblup  = rt_gblup,
      bayesR = rt_bayesR
    )

    mcmc_diag_obj <- switch(model_name,
      gblup  = NULL,
      bayesR = bayesR_full$mcmc_diag
    )
    if (!is.null(mcmc_diag_obj)) {
      iter_label_mcmc <- basename(config$data_dir)
      gen_label_mcmc  <- sub("^output_", "",
                             basename(dirname(config$data_dir)))
      mcmc_diag_obj$meta$scenario    <- sc
      mcmc_diag_obj$meta$iter        <- iter_label_mcmc
      mcmc_diag_obj$meta$generation  <- gen_label_mcmc
      mcmc_diag_obj$meta$marker_type <- "mh"
      saveRDS(mcmc_diag_obj,
              file.path(out_dir, model_name, "mcmc_diag.rds"))
      cat(sprintf("    [%s] mcmc_diag | saved\n", model_name))
    }

    saveRDS(
      list(
        predictions       = results_all,
        markers           = markers,
        markers_full      = markers,
        qtl_recovery      = qtl_recovery,
        qtl_recovery_full = qtl_recovery,
        diagnostics_fit   = diagnostics_fit,
        runtime_fit       = runtime_fit
      ),
      file.path(out_dir, model_name, "results_all.rds")
    )

    cat(sprintf("  %s | acc=%.4f bias=%.4f rmse=%.4f h2=%.3f\n",
                model_name, acc, bias, rmse, h2))
  }

  cat(sprintf("Scenario %s done in %.1f mins\n", sc,
              as.numeric(difftime(Sys.time(), sc_start, units = "mins"))))
}

cat(sprintf("\nAll scenarios complete | Total: %.1f mins\n",
            as.numeric(difftime(Sys.time(), pipeline_start, units = "mins"))))