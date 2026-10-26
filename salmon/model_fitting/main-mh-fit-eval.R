suppressPackageStartupMessages({
  library(tidyverse)
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
cat(sprintf("=== MH | variant=%s | %s | %s | %s ===\n", config$variant, iter_lab, sc, Sys.time()))

# ── 1. Genotypes and reference/candidate split ───────────────────────────────
load_hap_geno <- function(data_dir, n_chr, block_dir) {
  geno_list <- lapply(1:n_chr, function(chr)
    read_table(file.path(data_dir, block_dir, paste0("hap_geno_", chr)), show_col_types = FALSE))
  geno <- geno_list[[1]]
  for (chr in 2:n_chr) geno <- bind_cols(geno, geno_list[[chr]] %>% select(-ID))
  geno %>% rename(individual_id = ID)
}
hap_geno  <- load_hap_geno(data_dir, config$n_chromosomes, config$block_dir)
pedigree  <- read_csv(file.path(data_dir, "pedigree_G1.csv"), show_col_types = FALSE)
train_ids <- pedigree$individual_id[pedigree$population_type == "reference"]
test_ids  <- pedigree$individual_id[pedigree$population_type == "candidate"]

# ── 2. QTL truth; optional removal of causal blocks (variant "nocausal") ─────
qtl <- list(snp_ids = readRDS(file.path(config$qtl_dir, "qtl_snp_ids.rds")),
            s3 = readRDS(file.path(config$qtl_dir, "qtl_s3_effects.rds")),
            s5 = readRDS(file.path(config$qtl_dir, "qtl_s5_effects.rds")))
qtl_mh_blocks  <- readRDS(file.path(config$qtl_dir, "qtl_mh_blocks.rds"))
qtl_hap_blocks <- readRDS(file.path(config$qtl_dir, "qtl_hap_blocks.rds"))
snp_meta       <- readRDS(file.path(config$qtl_dir, "snp_meta.rds"))
qtl_effects    <- list(S3_equal = qtl$s3,
                       S4_unequal = readRDS(file.path(config$qtl_dir, "qtl_s4_effects.rds")),
                       S5_equal = qtl$s5,
                       S6_unequal = readRDS(file.path(config$qtl_dir, "qtl_s6_effects.rds")))
coords_file <- file.path(data_dir, "mh_info_ld_haploblock/stats/microhaplotype_coordinates.csv")
dropped_blocks <- character(0)
if (isTRUE(config$exclude_causal)) {
  sel <- read_csv(file.path(data_dir, "mh_info_ld_haploblock/stats/snp_selection_detailed.csv"),
                  show_col_types = FALSE)
  smi <- snp_meta %>% group_by(chr) %>% mutate(snp_index = row_number() - 1L) %>% ungroup()
  hap_geno <- drop_causal_mh_blocks(hap_geno, sel, smi, causal_snp_ids(sc, qtl))
  dropped_blocks <- attr(hap_geno, "dropped_blocks")
  cat(sprintf("nocausal: dropped %d MH blocks containing causal SNPs\n", length(dropped_blocks)))
}

pick_rows <- function(d, ids) d[match(ids, d$individual_id), ]
inputs <- build_mh_inputs(pick_rows(hap_geno, train_ids), pick_rows(hap_geno, test_ids),
                          config, coords_file)
cat(sprintf("W ready: %d train x %d alleles | %d candidates\n",
            nrow(inputs$W_train), ncol(inputs$W_train), nrow(inputs$W_test)))

# ── 3. GWAS map (MH) aligned to W block ids ──────────────────────────────────
map_mh <- NULL
if (isTRUE(config$gwas$enabled)) {
  coords <- read_csv(coords_file, show_col_types = FALSE) %>% mutate(chr = as.integer(chr))
  attr(inputs$W_train, "block_id") <- translate_block_ids(attr(inputs$W_train, "block_id"), coords)
  map_mh <- coords %>% filter(block_id %in% unique(attr(inputs$W_train, "block_id"))) %>%
    select(block_id, chr, start_pos, end_pos)
}

# ── 4. Fit and save ──────────────────────────────────────────────────────────
pheno_file <- file.path(data_dir, sprintf("pheno_%s_G1.csv", sub("_.*", "", sc)))
pheno_all  <- read_csv(pheno_file, show_col_types = FALSE)
fit <- fit_one_scenario(pheno_all, inputs, config, map = map_mh,
                        windsize = if (isTRUE(config$gwas$enabled)) config$gwas$windsize)
meta <- list(iter = iter_lab, scenario = sc, marker = "mh", variant = config$variant)

if (!is.null(fit$bayesR$pip)) {
  bid <- attr(inputs$W_train, "block_id")
  saveRDS(list(map = map_mh, pip = fit$bayesR$pip, pip_block = fit$bayesR$pip_block,
               max_pip_block = vapply(map_mh$block_id, function(b) max(fit$bayesR$pip[bid == b]), 0),
               pip_class = fit$bayesR$pip_class, gwas = fit$bayesR$gwas, gwas_meta = fit$bayesR$gwas_meta,
               marker_type = "mh", windsize = config$gwas$windsize, scenario = sc, iter = iter_lab,
               variant = config$variant),
          file.path(out_dir, "gwas_full.rds"))
}

for (model in c("gblup", "bayesR")) {
  dir.create(file.path(out_dir, model), showWarnings = FALSE)
  markers <- if (model == "bayesR") fit$bayesR$markers else NULL
  qrec <- compute_qtl_recovery(markers, sc, qtl$snp_ids, qtl_mh_blocks, qtl_hap_blocks,
                               snp_meta, "mh", mh_coords = NULL,
                               hap_coords = read_csv(config$hap_coords_file, show_col_types = FALSE) %>%
                                 mutate(chr = as.integer(chr)),
                               qtl_block_causal_snps = if (sc %in% names(qtl_effects))
                                 unique(unlist(lapply(qtl_effects[[sc]], `[[`, "snp_ids"))))
  if (model == "bayesR" && !is.null(fit$bayesR$mcmc_diag)) {
    md <- fit$bayesR$mcmc_diag
    md$meta[c("scenario", "iter", "marker_type", "variant")] <- list(sc, iter_lab, "mh", config$variant)
    md$meta$p <- ncol(inputs$W_train)
    saveRDS(md, file.path(out_dir, model, "mcmc_diag.rds"))
  }
  preds <- candidate_results(fit, model, meta)
  saveRDS(list(predictions = preds, markers = markers, markers_full = markers,
               qtl_recovery = qrec, qtl_recovery_full = qrec,
               diagnostics_fit = if (model == "bayesR") fit$bayesR$diagnostics,
               runtime_fit = fit$runtime[[model]], dropped_blocks = dropped_blocks),
          file.path(out_dir, model, "results_all.rds"))
  cat(sprintf("  %s | candidate acc=%.4f bias=%.4f | reference acc=%.4f | h2=%.3f\n",
              model, preds$accuracy[1], preds$bias[1], preds$acc_train[1], preds$h2[1]))
}

write_manifest(c(file.path(data_dir, config$block_dir), file.path(data_dir, "pedigree_G1.csv"),
                 pheno_file, coords_file, config$qtl_dir, config$hap_coords_file),
               file.path(out_dir, "manifest_inputs.tsv"), "input")
log_run_context(c(config, list(TASK = TASK, iter = iter_lab, scenario = sc)), out_dir)
write_manifest(out_dir, file.path(out_dir, "manifest_outputs.tsv"), "output")
cat(sprintf("=== done %s ===\n", Sys.time()))
