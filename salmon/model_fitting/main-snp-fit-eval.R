suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
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
source("R/models/bayesA.R")

pipeline_start <- Sys.time()
cat(sprintf("=== CV Pipeline started: %s ===\n", Sys.time()))

# ── 1. Load genotypes ──────────────────────────────────────────
cat("Loading genotype...\n")
genotype_raw <- {
  g <- as_tibble(fread(file.path(config$data_dir, config$genotype)))
  if ("ID" %in% names(g) && !"individual_id" %in% names(g))
    g <- g %>% rename(individual_id = ID)
  g %>% select(-any_of(c("generation", "population_type")))
}
cat(sprintf("Genotype loaded: %d individuals\n", nrow(genotype_raw)))

# ── 2. Split generasi terbaru per family ──────────────────────────────────────
pedigree_current <- read_csv(
  file.path(config$data_dir, config$pedigree),
  show_col_types = FALSE
) %>% filter(population_type == config$population_type)

set.seed(config$seed_split)
split_ids <- pedigree_current %>%
  group_by(family_id) %>%
  summarise(
    train_ids = list(sample(individual_id, round(n() * config$train_prop))),
    test_ids  = list(setdiff(individual_id, train_ids[[1]])),
    .groups   = "drop"
  )

current_train_ids <- unlist(split_ids$train_ids)
test_ids          <- unlist(split_ids$test_ids)
train_ids         <- current_train_ids

cat(sprintf("Split: train=%d | test=%d\n", length(train_ids), length(test_ids)))

# ── 3. QC + standardize dari training set ────────────────────────────────────
ind_ids_raw <- genotype_raw$individual_id
geno_matrix <- genotype_raw %>%
  select(-individual_id) %>%
  as.matrix()

# Full-data GRM via shared masbayes W
all_idx_raw      <- which(ind_ids_raw %in% c(train_ids, test_ids))
ind_ids_all      <- ind_ids_raw[all_idx_raw]
geno_all_raw     <- geno_matrix[all_idx_raw, ]
qc_full          <- quality_control(geno_all_raw, config)
geno_all_proc    <- impute_and_standardize(qc_full$filtered)
geno_all_imputed <- geno_all_proc$imputed
storage.mode(geno_all_imputed) <- "double"
rownames(geno_all_imputed) <- as.character(ind_ids_all)

# SNP encoding dispatch: G GBLUP always derives from VanRaden W; Bayesian
# W may additionally be z-scored when `config$snp_encoding = "zscore"`.
# Backward compat: missing config field → vanRaden everywhere, identical
# to the pre-PR-2 path.
snp_encoding <- config$snp_encoding %||% "vanRaden"

vR_full           <- masbayes::construct_snp_matrix(geno_all_imputed,
                                                    encoding = "vanRaden")
W_vR_full         <- vR_full$W
rownames(W_vR_full) <- as.character(ind_ids_all)
allele_freq_full  <- vR_full$freq
k_grm_full        <- 2 * sum(allele_freq_full * (1 - allele_freq_full))
A_full            <- tcrossprod(W_vR_full) / k_grm_full
rownames(A_full)  <- colnames(A_full) <- as.character(ind_ids_all)

if (snp_encoding == "zscore") {
  zs_full <- masbayes::construct_snp_matrix(geno_all_imputed,
                                            encoding = "zscore")
  W_full  <- zs_full$W
  rownames(W_full) <- as.character(ind_ids_all)
} else {
  W_full  <- W_vR_full
}

diag_vals  <- diag(A_full)
off_vals   <- A_full[lower.tri(A_full)]
grm_stats  <- list(
  mean_diag    = mean(diag_vals),
  sd_diag      = sd(diag_vals),
  mean_offdiag = mean(off_vals),
  sd_offdiag   = sd(off_vals)
)

train_idx_raw  <- which(ind_ids_raw %in% train_ids)
test_idx_raw   <- which(ind_ids_raw %in% test_ids)
geno_train_raw <- geno_matrix[train_idx_raw, ]
geno_test_raw  <- geno_matrix[test_idx_raw, ]
ind_ids_train  <- ind_ids_raw[train_idx_raw]
ind_ids_test   <- ind_ids_raw[test_idx_raw]

# QC dari training saja
qc_result           <- quality_control(geno_train_raw, config)
keep_snps           <- qc_result$keep_snps
geno_train_filtered <- geno_train_raw[, keep_snps]
geno_test_filtered  <- geno_test_raw[, keep_snps]

# Impute (training-only) and apply to test
geno_train_proc    <- impute_and_standardize(qc_result$filtered)
snp_means          <- colMeans(geno_train_proc$imputed)
geno_train_imputed <- geno_train_proc$imputed

geno_test_imputed <- geno_test_filtered
for (j in seq_len(ncol(geno_test_imputed))) {
  na_idx <- is.na(geno_test_imputed[, j])
  if (any(na_idx)) geno_test_imputed[na_idx, j] <- snp_means[j]
}
storage.mode(geno_train_imputed) <- "double"
storage.mode(geno_test_imputed)  <- "double"
rownames(geno_train_imputed) <- as.character(ind_ids_train)
rownames(geno_test_imputed)  <- as.character(ind_ids_test)

# VanRaden W: always built — supplies G GBLUP. Bayesian W may be a
# separately z-scored matrix when `config$snp_encoding = "zscore"`.
vR_train <- masbayes::construct_snp_matrix(geno_train_imputed,
                                           encoding = "vanRaden")
vR_test  <- masbayes::construct_snp_matrix(geno_test_imputed,
                                           encoding = "vanRaden",
                                           ref_freq = vR_train$freq)
W_vR_train <- vR_train$W
W_vR_test  <- vR_test$W
rownames(W_vR_train) <- as.character(ind_ids_train)
rownames(W_vR_test)  <- as.character(ind_ids_test)

allele_freq_train <- vR_train$freq
k_grm             <- 2 * sum(allele_freq_train * (1 - allele_freq_train))
GRM_train         <- tcrossprod(W_vR_train) / k_grm
W_vR_combined     <- rbind(W_vR_train, W_vR_test)
A_combined        <- tcrossprod(W_vR_combined) / k_grm

if (snp_encoding == "zscore") {
  zs_train <- masbayes::construct_snp_matrix(geno_train_imputed,
                                             encoding = "zscore")
  zs_test  <- masbayes::construct_snp_matrix(geno_test_imputed,
                                             encoding = "zscore",
                                             ref_freq = zs_train$freq,
                                             ref_sd   = zs_train$sd)
  W_train  <- zs_train$W
  W_test   <- zs_test$W
  rownames(W_train) <- as.character(ind_ids_train)
  rownames(W_test)  <- as.character(ind_ids_test)
} else {
  W_train  <- W_vR_train
  W_test   <- W_vR_test
}

rownames(GRM_train) <- colnames(GRM_train) <- as.character(ind_ids_train)
rownames(A_combined) <- colnames(A_combined) <-
  c(as.character(ind_ids_train), as.character(ind_ids_test))

cat(sprintf("Genotype ready: %d train | %d test | %d SNPs\n",
            length(ind_ids_train), length(ind_ids_test), sum(keep_snps)))
            
snp_map <- read_tsv(
  config$snp_map_file,
  show_col_types = FALSE
) %>%
  rename(snp_id = SNPID, chr = Chr, position = Position) %>%
  mutate(chr = as.integer(chr))

# Buat allele_info untuk SNP — match ke kolom W_train
snp_names_train <- colnames(geno_train_imputed)  # nama SNP setelah QC
snp_allele_info <- snp_map %>%
  filter(snp_id %in% snp_names_train) %>%
  arrange(match(snp_id, snp_names_train)) %>%
  rename(allele_id = snp_id, start_pos = position) %>%
  mutate(end_pos = start_pos, block = allele_id, allele = 1L,
         freq = NA_real_, n_snps = 1L, global_block_id = allele_id)

snp_allele_info_full <- snp_map %>%
  filter(snp_id %in% colnames(W_full)) %>%
  arrange(match(snp_id, colnames(W_full))) %>%
  rename(allele_id = snp_id, start_pos = position) %>%
  mutate(end_pos = start_pos, block = allele_id, allele = 1L,
         freq = NA_real_, n_snps = 1L, global_block_id = allele_id)

# ── GWAS map (SNP): masbayes schema, aligned to W_full columns ───────────────
# Reuses snp_allele_info_full ordering (already match()-ed to colnames(W_full))
# so masbayes's "one row per W column" contract holds without a second match.
gwas_enabled <- isTRUE(config$gwas$enabled)
if (gwas_enabled) {
  map_snp_gwas <- snp_allele_info_full %>%
    transmute(SNP   = allele_id,
              CHROM = as.integer(chr),
              POS   = as.integer(start_pos))
  if (nrow(map_snp_gwas) != ncol(W_full)) {
    stop(sprintf("map_snp_gwas rows (%d) != ncol(W_full) (%d) — SNP map and W_full out of sync.",
                 nrow(map_snp_gwas), ncol(W_full)))
  }
  cat(sprintf("GWAS enabled | map SNPs: %d | windsize: %.0f bp\n",
              nrow(map_snp_gwas), config$gwas$windsize))
} else {
  map_snp_gwas <- NULL
  cat("GWAS disabled (config$gwas$enabled != TRUE)\n")
}


qtl_snp_ids    <- readRDS(file.path(config$qtl_dir, "qtl_snp_ids.rds"))
qtl_mh_blocks  <- readRDS(file.path(config$qtl_dir, "qtl_mh_blocks.rds"))
qtl_hap_blocks <- readRDS(file.path(config$qtl_dir, "qtl_hap_blocks.rds"))
snp_meta_qtl   <- readRDS(file.path(config$qtl_dir, "snp_meta.rds"))
mh_coords      <- read_csv(config$mh_coords_file, show_col_types = FALSE) %>%
                  mutate(chr = as.integer(chr))
hap_coords     <- read_csv(config$hap_coords_file, show_col_types = FALSE) %>%
                  mutate(chr = as.integer(chr))
cat(sprintf("QTL loaded: %d SNP QTL | %d MH block QTL | %d HAP block QTL\n",
            length(qtl_snp_ids), length(qtl_mh_blocks), length(qtl_hap_blocks)))  

# ── 4. Loop skenario ──────────────────────────────────────────────────────────
for (sc in names(config$scenarios)) {

  out_dir  <- file.path(config$output_dir, sc)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  sc_start <- Sys.time()

  pheno_all <- read_csv(
    file.path(config$data_dir, config$scenarios[[sc]]),
    show_col_types = FALSE
  )

  pheno_train <- pheno_all %>%
    filter(individual_id %in% ind_ids_train) %>%
    arrange(match(individual_id, ind_ids_train))

  pheno_test <- pheno_all %>%
    filter(individual_id %in% ind_ids_test) %>%
    arrange(match(individual_id, ind_ids_test))

  y_train <- pheno_train$phenotype

  cat(sprintf("\n[%s] Scenario: %s | train=%d | test=%d\n",
              sc, Sys.time(), nrow(pheno_train), nrow(pheno_test)))

  # ── Full-data fitting ──
  pheno_full    <- pheno_all %>% arrange(match(individual_id, ind_ids_all))
  y_full        <- pheno_full$phenotype
  matrices_full <- list(
    W_train     = W_full,
    W_test      = W_full[0, , drop = FALSE],
    A_train     = A_full,
    A_combined  = A_full,
    k_ah        = k_grm_full,
    WtW_diag    = colSums(W_full^2),
    allele_info = snp_allele_info_full,
    train_ids   = as.character(ind_ids_all),
    test_ids    = character(0),
    g_label     = "snp_add"
  )
  split_full <- list(
    train = list(phenotypes = pheno_full),
    test  = list(phenotypes = pheno_full[0, ])
  )

  cat("  Running full-data fitting...\n")
  t0_gblup  <- proc.time()
  gblup_full  <- run_gblup(matrices_full, split_full, config)
  rt_gblup  <- (proc.time() - t0_gblup)[["elapsed"]]
  cat(sprintf("    [full] gblup  | acc=%.4f bias=%.4f h2=%.3f\n",
    cor(gblup_full$predictions$train, pheno_full$tbv, use="complete.obs"),
    coef(lm(pheno_full$tbv ~ gblup_full$predictions$train))[2],
    gblup_full$varcomp$h2))

  t0_bayesR  <- proc.time()
  bayesR_full <- run_bayesR(
    matrices_full, split_full, gblup_full$varcomp, config, fold = 0,
    map      = if (gwas_enabled) map_snp_gwas          else NULL,
    windsize = if (gwas_enabled) config$gwas$windsize  else NULL
  )
  rt_bayesR  <- (proc.time() - t0_bayesR)[["elapsed"]]
  cat(sprintf("    [full] bayesR | acc=%.4f bias=%.4f h2=%.3f\n",
    cor(bayesR_full$predictions$train, pheno_full$tbv, use="complete.obs"),
    coef(lm(pheno_full$tbv ~ bayesR_full$predictions$train))[2],
    bayesR_full$varcomp$h2))

  # ── GWAS Rds (separate from results_all.rds per Rds-contract policy) ──
  # For SNP path, pip_block == pip by masbayes convention, so no max-per-block
  # derivation is needed; the per-SNP PIP is already the finest resolution.
  if (gwas_enabled && !is.null(bayesR_full$pip)) {
    iter_label <- basename(config$data_dir)
    gen_label  <- sub("^output_", "", basename(dirname(config$data_dir)))

    saveRDS(
      list(
        map         = map_snp_gwas,
        pip         = bayesR_full$pip,
        pip_block   = bayesR_full$pip_block,
        gwas        = bayesR_full$gwas,
        gwas_meta   = bayesR_full$gwas_meta,
        marker_type = "snp",
        windsize    = config$gwas$windsize,
        scenario    = sc,
        iter        = iter_label,
        generation  = gen_label
      ),
      file.path(out_dir, "gwas_full.rds")
    )
    cat(sprintf("    [full] gwas   | saved: %s\n",
                file.path(out_dir, "gwas_full.rds")))
  }

  t0_bayesA  <- proc.time()
  bayesA_full <- run_bayesA(matrices_full, split_full, gblup_full$varcomp, config, fold = 0)
  rt_bayesA  <- (proc.time() - t0_bayesA)[["elapsed"]]
  cat(sprintf("    [full] bayesA | acc=%.4f bias=%.4f h2=%.3f\n",
    cor(bayesA_full$predictions$train, pheno_full$tbv, use="complete.obs"),
    coef(lm(pheno_full$tbv ~ bayesA_full$predictions$train))[2],
    bayesA_full$varcomp$h2))

  full_results <- list(
    gblup  = list(GEBV = gblup_full$predictions$train,  h2 = gblup_full$varcomp$h2),
    bayesR = list(GEBV = bayesR_full$predictions$train, h2 = bayesR_full$varcomp$h2),
    bayesA = list(GEBV = bayesA_full$predictions$train, h2 = bayesA_full$varcomp$h2)
  )

  matrices <- list(
    W_train     = W_train,
    W_test      = W_test,
    A_train     = GRM_train,
    A_combined  = A_combined,
    k_ah        = k_grm,
    WtW_diag    = colSums(W_train^2),
    allele_info = snp_allele_info,
    train_ids   = as.character(ind_ids_train),
    test_ids    = as.character(ind_ids_test),
    g_label     = "snp_add"
  )

  split <- list(
    train = list(phenotypes = pheno_train),
    test  = list(phenotypes = pheno_test)
  )

  # ── GBLUP ──
  cat("  Running GBLUP...\n")
  dir.create(file.path(out_dir, "gblup"), recursive = TRUE, showWarnings = FALSE)
  gblup_res        <- run_gblup(matrices, split, config)
  GEBV_test_gblup  <- gblup_res$predictions$test
  GEBV_train_gblup <- gblup_res$predictions$train
  cat(sprintf("    [cv]   gblup  | train: acc=%.4f bias=%.4f | test: acc=%.4f bias=%.4f | h2=%.3f\n",
    cor(GEBV_train_gblup, pheno_train$tbv, use="complete.obs"),
    coef(lm(pheno_train$tbv ~ GEBV_train_gblup))[2],
    cor(GEBV_test_gblup, pheno_test$tbv, use="complete.obs"),
    coef(lm(pheno_test$tbv ~ GEBV_test_gblup))[2],
    gblup_res$varcomp$h2))

  # ── BayesR ──
  cat("  Running BayesR...\n")
  dir.create(file.path(out_dir, "bayesR"), recursive = TRUE, showWarnings = FALSE)
  bayesR_res        <- run_bayesR(matrices, split, gblup_res$varcomp, config, fold = 0)
  GEBV_test_bayesR  <- bayesR_res$predictions$test
  GEBV_train_bayesR <- bayesR_res$predictions$train
  cat(sprintf("    [cv]   bayesR | train: acc=%.4f bias=%.4f | test: acc=%.4f bias=%.4f | h2=%.3f\n",
    cor(GEBV_train_bayesR, pheno_train$tbv, use="complete.obs"),
    coef(lm(pheno_train$tbv ~ GEBV_train_bayesR))[2],
    cor(GEBV_test_bayesR, pheno_test$tbv, use="complete.obs"),
    coef(lm(pheno_test$tbv ~ GEBV_test_bayesR))[2],
    bayesR_res$varcomp$h2))

  # ── BayesA ──
  cat("  Running BayesA...\n")
  dir.create(file.path(out_dir, "bayesA"), recursive = TRUE, showWarnings = FALSE)
  bayesA_res        <- run_bayesA(matrices, split, gblup_res$varcomp, config, fold = 0)
  GEBV_test_bayesA  <- bayesA_res$predictions$test
  GEBV_train_bayesA <- bayesA_res$predictions$train
  cat(sprintf("    [cv]   bayesA | train: acc=%.4f bias=%.4f | test: acc=%.4f bias=%.4f | h2=%.3f\n",
    cor(GEBV_train_bayesA, pheno_train$tbv, use="complete.obs"),
    coef(lm(pheno_train$tbv ~ GEBV_train_bayesA))[2],
    cor(GEBV_test_bayesA, pheno_test$tbv, use="complete.obs"),
    coef(lm(pheno_test$tbv ~ GEBV_test_bayesA))[2],
    bayesA_res$varcomp$h2))
  
  # ── Save predict_cv per model ──
  for (model_name in c("gblup", "bayesR", "bayesA")) {

    GEBV_test <- switch(model_name,
      gblup  = GEBV_test_gblup,
      bayesR = GEBV_test_bayesR,
      bayesA = GEBV_test_bayesA
    )
    h2 <- switch(model_name,
      gblup  = gblup_res$varcomp$h2,
      bayesR = bayesR_res$varcomp$h2,
      bayesA = bayesA_res$varcomp$h2
    )
    GEBV_train <- switch(model_name,
      gblup  = GEBV_train_gblup,
      bayesR = GEBV_train_bayesR,
      bayesA = GEBV_train_bayesA
    )

    acc  <- cor(GEBV_test,  pheno_test$tbv,  use = "complete.obs")
    bias <- coef(lm(pheno_test$tbv  ~ GEBV_test))[2]
    rmse <- sqrt(mean((GEBV_test  - pheno_test$tbv)^2,  na.rm = TRUE))

    acc_train  <- cor(GEBV_train, pheno_train$tbv, use = "complete.obs")
    bias_train <- coef(lm(pheno_train$tbv ~ GEBV_train))[2]
    rmse_train <- sqrt(mean((GEBV_train - pheno_train$tbv)^2, na.rm = TRUE))

    GEBV_full <- full_results[[model_name]]$GEBV
    h2_full   <- full_results[[model_name]]$h2
    acc_full  <- cor(GEBV_full, pheno_full$tbv, use = "complete.obs")
    bias_full <- coef(lm(pheno_full$tbv ~ GEBV_full))[2]
    rmse_full <- sqrt(mean((GEBV_full - pheno_full$tbv)^2, na.rm = TRUE))

    results_all <- tibble(
      individual_id    = pheno_test$individual_id,
      tbv              = pheno_test$tbv,
      GEBV             = GEBV_test,
      scenario         = sc,
      model            = model_name,
      marker           = "snp",
      split            = "test",
      accuracy         = acc,
      bias             = bias,
      rmse             = rmse,
      h2               = h2,
      n_train          = length(ind_ids_train),
      n_test           = length(ind_ids_test),
      h2_full          = h2_full,
      acc_full         = acc_full,
      bias_full        = bias_full,
      rmse_full        = rmse_full,
      grm_mean_diag    = grm_stats$mean_diag,
      grm_sd_diag      = grm_stats$sd_diag,
      grm_mean_offdiag = grm_stats$mean_offdiag,
      grm_sd_offdiag   = grm_stats$sd_offdiag
    )

    markers_cv <- switch(model_name,
      gblup  = NULL,
      bayesR = bayesR_res$markers,
      bayesA = bayesA_res$markers
    )
    markers_full_mod <- switch(model_name,
      gblup  = NULL,
      bayesR = bayesR_full$markers,
      bayesA = bayesA_full$markers
    )

    qtl_recovery_cv   <- compute_qtl_recovery(markers_cv, sc, qtl_snp_ids, qtl_mh_blocks, qtl_hap_blocks, snp_meta_qtl, "snp", mh_coords = mh_coords, hap_coords = hap_coords)
    qtl_recovery_full <- compute_qtl_recovery(markers_full_mod, sc, qtl_snp_ids, qtl_mh_blocks, qtl_hap_blocks, snp_meta_qtl, "snp", mh_coords = mh_coords, hap_coords = hap_coords)

    diagnostics_fit <- switch(model_name,
      gblup  = NULL,
      bayesR = bayesR_full$diagnostics,
      bayesA = bayesA_full$diagnostics
    )

    runtime_fit <- switch(model_name,
      gblup  = rt_gblup,
      bayesR = rt_bayesR,
      bayesA = rt_bayesA
    )

    # ── MCMC diagnostic Rds (separate from results_all.rds; opt-in) ──
    # Wrapper attaches $mcmc_diag if config$mcmc_diag$enabled = TRUE.
    # Saved separately so the results_all.rds contract stays unchanged.
    mcmc_diag_obj <- switch(model_name,
      gblup  = NULL,
      bayesR = bayesR_full$mcmc_diag,
      bayesA = bayesA_full$mcmc_diag
    )
    if (!is.null(mcmc_diag_obj)) {
      iter_label_mcmc <- basename(config$data_dir)
      gen_label_mcmc  <- sub("^output_", "",
                             basename(dirname(config$data_dir)))
      mcmc_diag_obj$meta$scenario    <- sc
      mcmc_diag_obj$meta$iter        <- iter_label_mcmc
      mcmc_diag_obj$meta$generation  <- gen_label_mcmc
      mcmc_diag_obj$meta$marker_type <- "snp"
      saveRDS(
        mcmc_diag_obj,
        file.path(out_dir, model_name, "mcmc_diag.rds")
      )
      cat(sprintf("    [%s] mcmc_diag | saved: %s\n", model_name,
                  file.path(out_dir, model_name, "mcmc_diag.rds")))
    }

    saveRDS(
      list(
        predictions   = results_all,
        markers       = markers_cv,
        markers_full  = markers_full_mod,
        qtl_recovery  = qtl_recovery_cv,
        qtl_recovery_full = qtl_recovery_full,
        diagnostics_fit = diagnostics_fit,
        runtime_fit = runtime_fit
      ),
      file.path(out_dir, model_name, "results_all.rds")
    )
    
    cat(sprintf("  %s | full: acc=%.4f bias=%.4f h2=%.3f | train: acc=%.4f bias=%.4f rmse=%.4f | test: acc=%.4f bias=%.4f rmse=%.4f\n",
                model_name, acc_full, bias_full, h2_full,
                acc_train, bias_train, rmse_train, acc, bias, rmse))
  }

  cat(sprintf("Scenario %s done in %.1f mins\n", sc,
              as.numeric(difftime(Sys.time(), sc_start, units = "mins"))))
}

cat(sprintf("\nAll scenarios complete | Total: %.1f mins\n",
            as.numeric(difftime(Sys.time(), pipeline_start, units = "mins"))))