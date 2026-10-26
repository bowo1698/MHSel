suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(coda)
  library(masbayes)
  library(masreml)
  library(RhpcBLASctl)
})

select <- dplyr::select
filter <- dplyr::filter
blas_set_num_threads(1)

# One task = one (iteration, scenario). Train on reference (phenotyped)
# offspring, evaluate on unphenotyped candidates (Delomas et al. 2024).
config_file <- Sys.getenv("GENOMIC_CONFIG", "config.R")
source(config_file)
TASK  <- as.integer(commandArgs(trailingOnly = TRUE)[1])
FORCE <- Sys.getenv("FORCE") == "1"

mf <- file.path(config$code_dir, "salmon", "model_fitting", "R")
for (f in c("helper.R", "models/gblup.R", "models/bayesR.R", "fit_scenario.R",
            "causal.R", "task_index.R")) source(file.path(mf, f))
source(file.path(config$common_dir, "bayesr_prior.R"))
source(file.path(config$common_dir, "manifest.R"))

ts       <- task_to_iter_scenario(TASK)
sc       <- ts$scenario
iter_lab <- sprintf("iter%02d", ts$iter)
data_dir <- file.path(config$g1_dir, iter_lab)
out_dir  <- file.path(config$output_dir, iter_lab, sc)
if (should_skip(out_dir, FORCE)) { cat("Done already:", out_dir, "\n"); quit(save = "no") }
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("=== SNP | variant=%s | %s | %s | %s ===\n", config$variant, iter_lab, sc, Sys.time()))

# ── 1. Genotypes and reference/candidate split ───────────────────────────────
geno_file <- file.path(data_dir, config$genotype)
geno <- as_tibble(fread(geno_file)) %>%
  rename(individual_id = ID) %>%
  select(-any_of(c("generation", "population_type")))
pedigree  <- read_csv(file.path(data_dir, "pedigree_G1.csv"), show_col_types = FALSE)
train_ids <- pedigree$individual_id[pedigree$population_type == "reference"]
test_ids  <- pedigree$individual_id[pedigree$population_type == "candidate"]

# ── 2. QTL truth; optional removal of causal SNPs (variant "nocausal") ───────
qtl <- list(snp_ids = readRDS(file.path(config$qtl_dir, "qtl_snp_ids.rds")),
            s3 = readRDS(file.path(config$qtl_dir, "qtl_s3_effects.rds")),
            s5 = readRDS(file.path(config$qtl_dir, "qtl_s5_effects.rds")))
qtl_mh_blocks  <- readRDS(file.path(config$qtl_dir, "qtl_mh_blocks.rds"))
qtl_hap_blocks <- readRDS(file.path(config$qtl_dir, "qtl_hap_blocks.rds"))
snp_meta       <- readRDS(file.path(config$qtl_dir, "snp_meta.rds"))
dropped_snps <- character(0)
if (isTRUE(config$exclude_causal)) {
  dropped_snps <- intersect(causal_snp_ids(sc, qtl), names(geno))
  geno <- drop_causal_snps(geno, dropped_snps)
  cat(sprintf("nocausal: dropped %d causal SNPs\n", length(dropped_snps)))
}

# QC, imputation and frequencies from the reference set only
X <- as.matrix(geno %>% select(-individual_id)); rownames(X) <- geno$individual_id
X_train_raw <- X[train_ids, , drop = FALSE]; X_test_raw <- X[test_ids, , drop = FALSE]
qc <- quality_control(X_train_raw, config)
X_train <- impute_and_standardize(qc$filtered)$imputed
X_test  <- X_test_raw[, colnames(X_train), drop = FALSE]
mu_tr   <- colMeans(X_train)
for (j in which(colSums(is.na(X_test)) > 0)) X_test[is.na(X_test[, j]), j] <- mu_tr[j]
storage.mode(X_train) <- "double"; storage.mode(X_test) <- "double"
rownames(X_train) <- train_ids; rownames(X_test) <- test_ids

snp_map <- read_tsv(config$snp_map_file, show_col_types = FALSE) %>%
  rename(snp_id = SNPID, chr = Chr, position = Position) %>% mutate(chr = as.integer(chr))
allele_info <- snp_map[match(colnames(X_train), snp_map$snp_id), ] %>%
  transmute(allele_id = snp_id, freq = colMeans(X_train) / 2, block = snp_id, allele = 1L,
            chr, start_pos = position, end_pos = position, n_snps = 1L, global_block_id = snp_id)
inputs <- build_snp_inputs(X_train, X_test, as.data.frame(allele_info))
cat(sprintf("W ready: %d train x %d SNPs | %d candidates\n",
            nrow(inputs$W_train), ncol(inputs$W_train), nrow(inputs$W_test)))

map_snp <- if (isTRUE(config$gwas$enabled))
  data.frame(SNP = allele_info$allele_id, CHROM = allele_info$chr, POS = as.integer(allele_info$start_pos))

# ── 3. Fit and save ──────────────────────────────────────────────────────────
pheno_file <- file.path(data_dir, sprintf("pheno_%s_G1.csv", sub("_.*", "", sc)))
pheno_all  <- read_csv(pheno_file, show_col_types = FALSE)
fit <- fit_one_scenario(pheno_all, inputs, config, map = map_snp,
                        windsize = if (isTRUE(config$gwas$enabled)) config$gwas$windsize)
meta <- list(iter = iter_lab, scenario = sc, marker = "snp", variant = config$variant)

if (!is.null(fit$bayesR$pip)) {
  saveRDS(list(map = map_snp, pip = fit$bayesR$pip, pip_block = fit$bayesR$pip_block,
               pip_class = fit$bayesR$pip_class, gwas = fit$bayesR$gwas, gwas_meta = fit$bayesR$gwas_meta,
               marker_type = "snp", windsize = config$gwas$windsize, scenario = sc, iter = iter_lab,
               variant = config$variant),
          file.path(out_dir, "gwas_full.rds"))
}

mh_coords  <- read_csv(config$mh_coords_file,  show_col_types = FALSE) %>% mutate(chr = as.integer(chr))
hap_coords <- read_csv(config$hap_coords_file, show_col_types = FALSE) %>% mutate(chr = as.integer(chr))
for (model in c("gblup", "bayesR")) {
  dir.create(file.path(out_dir, model), showWarnings = FALSE)
  markers <- if (model == "bayesR") fit$bayesR$markers else NULL
  qrec <- compute_qtl_recovery(markers, sc, qtl$snp_ids, qtl_mh_blocks, qtl_hap_blocks,
                               snp_meta, "snp", mh_coords = mh_coords, hap_coords = hap_coords)
  if (model == "bayesR" && !is.null(fit$bayesR$mcmc_diag)) {
    md <- fit$bayesR$mcmc_diag
    md$meta[c("scenario", "iter", "marker_type", "variant")] <- list(sc, iter_lab, "snp", config$variant)
    md$meta$p <- ncol(inputs$W_train)
    saveRDS(md, file.path(out_dir, model, "mcmc_diag.rds"))
  }
  preds <- candidate_results(fit, model, meta)
  saveRDS(list(predictions = preds, markers = markers, markers_full = markers,
               qtl_recovery = qrec, qtl_recovery_full = qrec,
               diagnostics_fit = if (model == "bayesR") fit$bayesR$diagnostics,
               runtime_fit = fit$runtime[[model]], dropped_snps = dropped_snps),
          file.path(out_dir, model, "results_all.rds"))
  cat(sprintf("  %s | candidate acc=%.4f bias=%.4f | reference acc=%.4f | h2=%.3f\n",
              model, preds$accuracy[1], preds$bias[1], preds$acc_train[1], preds$h2[1]))
}

write_manifest(c(geno_file, file.path(data_dir, "pedigree_G1.csv"), pheno_file, config$qtl_dir,
                 config$snp_map_file, config$mh_coords_file, config$hap_coords_file),
               file.path(out_dir, "manifest_inputs.tsv"), "input")
log_run_context(c(config, list(TASK = TASK, iter = iter_lab, scenario = sc)), out_dir)
write_manifest(out_dir, file.path(out_dir, "manifest_outputs.tsv"), "output")
cat(sprintf("=== done %s ===\n", Sys.time()))
