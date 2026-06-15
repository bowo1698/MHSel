suppressPackageStartupMessages({
  library(tidyverse)
  library(coda)
  library(foreach)
  library(doParallel)
  library(masbayes)
  library(masreml)
  library(pROC)
})

select <- dplyr::select
filter <- dplyr::filter

config_file <- Sys.getenv("GENOMIC_CONFIG", "config.R")
source(config_file)

cat("\n=== Configuration Loaded ===\n")
cat("Trait:", config$phenotype_trait, "| Response:", config$response_type %||% "gaussian", "\n")
cat("Genotype:", config$genotype_file, "\n")
cat("Output:", config$output_dir, "| Folds:", config$k_folds, "\n")
cat("============================\n\n")

pipeline_start_time <- Sys.time()
setwd(config$base_dir)

source("R/data_prep_snp.R")
source("R/models/gblup.R")
source("R/models/bayesR.R")
source("R/models/bayesA.R")

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(config$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

is_binary <- !is.null(config$response_type) && config$response_type == "binary"

cat("Loading genomic data...\n")
data <- load_genomic_data(config)
cat("Total individuals:", nrow(data$phenotypes), "\n")
cat("Total SNPs:", data$n_snps, "\n\n")

cat("Setting up cross-validation...\n")
cv_setup <- setup_cv(data, config)
cat("CV method:", cv_setup$method, "| Folds:", config$k_folds, "\n")
cat("Fold distribution:\n"); print(table(cv_setup$folds)); cat("\n")

# ── Full data fit ──
base_output_dir <- config$base_output_dir %||% config$output_dir
full_fit_path   <- file.path(base_output_dir, "full_fit.rds")

if (file.exists(full_fit_path)) {
  cat("Full fit found at:", full_fit_path, "— loading...\n")
  full_fit <- readRDS(full_fit_path)
  h2_full  <- full_fit$h2_full
  cat(sprintf("h2_full: GBLUP=%.4f | BayesR=%.4f | BayesA=%.4f\n\n",
              h2_full$gblup %||% NA, h2_full$bayesR, h2_full$bayesA))
} else {
  cat("\n=== Full Data Fit ===\n")
  full_start <- Sys.time()

  full_split <- list(
    train = list(phenotypes = data$phenotypes, genotypes = data$genotypes),
    test  = list(phenotypes = data$phenotypes[0, ], genotypes = data$genotypes[0, ])
  )
  matrices_full <- construct_matrices(full_split$train, full_split$test, config)

  grm_stats_full <- list(
    mean_diag    = mean(diag(matrices_full$A_train)),
    mean_offdiag = mean(matrices_full$A_train[lower.tri(matrices_full$A_train)]),
    range_diag   = range(diag(matrices_full$A_train)),
    k_ah         = matrices_full$k_ah
  )
  cat("GRM full | mean_diag:", round(grm_stats_full$mean_diag, 4),
      "| mean_offdiag:", round(grm_stats_full$mean_offdiag, 4), "\n\n")

  cat("Running GBLUP full...\n")
  t0_gblup  <- proc.time()
  gblup_full <- run_gblup(matrices_full, full_split, config)
  rt_gblup  <- (proc.time() - t0_gblup)[["elapsed"]]
  full_auc_tag <- function(res) {
    if (!is_binary) return("")
    sprintf(" | AUC_train=%.4f", res$accuracy$auc_train %||% NA_real_)
  }

  cat(sprintf("GBLUP full: h2=%.4f | train_acc=%.4f%s\n\n",
              gblup_full$varcomp$h2, gblup_full$accuracy$train_pheno,
              full_auc_tag(gblup_full)))
  gblup_varcomp_full <- gblup_full$varcomp

  config_full <- config
  config_full$save_mcmc_samples <- TRUE

  # ── GWAS map (SNP): read config$snp_map_file, align to colnames(W_full) ──
  # masbayes contract: nrow(map) == ncol(W) with rows in the same order as W
  # columns. We reorder the map to colnames(W_full) and stop with a clear
  # message if any column is unmatched, so misaligned SNP IDs surface at
  # config load time rather than as silent NULL fields in gwas_full.rds.
  gwas_enabled <- isTRUE(config$gwas$enabled)
  if (gwas_enabled) {
    snp_col_names <- colnames(matrices_full$W_train)
    if (is.null(snp_col_names)) {
      stop("colnames(W_full) is NULL — cannot align SNP map to W columns.")
    }
    map_snp_raw <- read_tsv(config$snp_map_file, show_col_types = FALSE)
    map_snp_gwas <- map_snp_raw %>%
      transmute(SNP   = SNPID,
                CHROM = as.integer(Chr),
                POS   = as.integer(Position)) %>%
      filter(SNP %in% snp_col_names) %>%
      arrange(match(SNP, snp_col_names))
    if (nrow(map_snp_gwas) != length(snp_col_names)) {
      # Surface first few unmatched IDs so the user can tell whether the
      # genotype panel uses a different ID namespace (e.g. custom imputed
      # target SNPs vs. the 12K MH-derived map) and decide between supplying
      # a panel-specific map or disabling GWAS for this panel.
      missing_ids    <- setdiff(snp_col_names, map_snp_gwas$SNP)
      sample_missing <- paste(head(missing_ids, 5), collapse = ", ")
      stop(sprintf(
        paste0("map_snp rows (%d) != ncol(W_full) (%d) — SNP map missing %d W columns.\n",
               "  First missing IDs: %s\n",
               "  Hint: periksa config$snp_map_file (panel-specific map) ",
               "atau set config$gwas$enabled = FALSE untuk panel ini."),
        nrow(map_snp_gwas), length(snp_col_names),
        length(missing_ids), sample_missing))
    }
    cat(sprintf("GWAS enabled | map SNPs: %d | windsize: %.0f bp\n",
                nrow(map_snp_gwas), config$gwas$windsize))
  } else {
    map_snp_gwas <- NULL
    cat("GWAS disabled (config$gwas$enabled != TRUE)\n")
  }

  cat("Running BayesR full...\n")
  t0_bayesR  <- proc.time()
  bayesR_full <- run_bayesR(
    matrices_full, full_split, gblup_varcomp_full, config_full, fold = 0,
    map      = if (gwas_enabled) map_snp_gwas         else NULL,
    windsize = if (gwas_enabled) config$gwas$windsize else NULL
  )
  rt_bayesR  <- (proc.time() - t0_bayesR)[["elapsed"]]
  cat(sprintf("BayesR full: h2=%.4f | train_acc=%.4f%s\n\n",
              bayesR_full$varcomp$h2, bayesR_full$accuracy$train_pheno,
              full_auc_tag(bayesR_full)))

  # ── GWAS Rds (separate from full_fit.rds per Rds-contract policy) ──
  # For SNP, pip_block == pip by masbayes convention; per-SNP PIP is the
  # finest resolution and no max-per-block derivation is needed.
  if (gwas_enabled && !is.null(bayesR_full$pip)) {
    gwas_path <- file.path(base_output_dir, "gwas_full.rds")
    dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(
      list(
        map         = map_snp_gwas,
        pip         = bayesR_full$pip,
        pip_block   = bayesR_full$pip_block,
        pip_class   = bayesR_full$pip_class,
        gwas        = bayesR_full$gwas,
        gwas_meta   = bayesR_full$gwas_meta,
        marker_type = "snp",
        windsize    = config$gwas$windsize,
        trait       = config$phenotype_trait,
        panel       = config$genotype_file,
        rep_id      = config$rep_id %||% 1
      ),
      gwas_path
    )
    cat(sprintf("GWAS RDS saved: %s\n\n", gwas_path))
  }

  cat("Running BayesA full...\n")
  t0_bayesA  <- proc.time()
  bayesA_full <- run_bayesA(matrices_full, full_split, gblup_varcomp_full, config_full, fold = 0)
  rt_bayesA  <- (proc.time() - t0_bayesA)[["elapsed"]]
  cat(sprintf("BayesA full: h2=%.4f | train_acc=%.4f%s\n\n",
              bayesA_full$varcomp$h2, bayesA_full$accuracy$train_pheno,
              full_auc_tag(bayesA_full)))

  h2_full <- list(
    gblup  = gblup_full$varcomp$h2,
    bayesR = bayesR_full$varcomp$h2,
    bayesA = bayesA_full$varcomp$h2
  )

  # All three models fit in a single full_fit. Store diagnostics/runtime as
  # named lists per model (oysters pipeline differs from simulation, where a
  # for-loop selects one model_name per iteration).
  diagnostics_fit <- list(
    gblup  = NULL,
    bayesR = bayesR_full$diagnostics,
    bayesA = bayesA_full$diagnostics
  )

  runtime_fit <- list(
    gblup  = rt_gblup,
    bayesR = rt_bayesR,
    bayesA = rt_bayesA
  )

  full_fit <- list(
    grm_stats = grm_stats_full,
    gblup     = gblup_full,
    bayesR    = bayesR_full,
    bayesA    = bayesA_full,
    h2_full   = h2_full,
    diagnostics_fit = diagnostics_fit,
    runtime_fit = runtime_fit
  )

  dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(full_fit, full_fit_path)
  cat(sprintf("Full data fit complete | saved: %s | %.2fh\n\n",
              full_fit_path,
              as.numeric(difftime(Sys.time(), full_start, units = "hours"))))
}

# Parallel setup
n_cores <- config$k_folds
total_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores()))
n_threads_per_fold <- max(1, floor(total_cpus / (n_cores + 2)))
cat("Parallel:", n_cores, "cores |", n_threads_per_fold, "threads/fold\n\n")

cl <- makeCluster(n_cores)
registerDoParallel(cl)
clusterExport(cl, c("n_threads_per_fold", "h2_full"))
clusterEvalQ(cl, {
  library(RhpcBLASctl); blas_set_num_threads(n_threads_per_fold)
  library(tidyverse); library(coda); library(masreml)
  library(masbayes); library(pROC)
})

# CV loop
all_results <- foreach(
  fold = 1:config$k_folds,
  .packages = c("tidyverse", "coda", "masreml", "masbayes", "pROC"),
  .errorhandling = "pass",
  .verbose = FALSE
) %dopar% {

  source("R/data_prep_snp.R")
  source("R/models/gblup.R")
  source("R/models/bayesR.R")
  source("R/models/bayesA.R")
  select <- dplyr::select; filter <- dplyr::filter
  is_binary <- !is.null(config$response_type) && config$response_type == "binary"

  tryCatch({
    fold_start <- Sys.time()
    log_file <- file.path(config$output_dir, "logs", paste0("fold_", fold, ".log"))
    log_con  <- file(log_file, open = "wt")
    sink(log_con, type = "output"); sink(log_con, type = "message")

    cat(sprintf("\n========== FOLD %d/%d ==========\n", fold, config$k_folds))
    cat("Started:", as.character(Sys.time()), "\n\n")
    cat(sprintf("h2_full: GBLUP=%.4f | BayesR=%.4f | BayesA=%.4f\n\n",
                h2_full$gblup %||% NA, h2_full$bayesR, h2_full$bayesA))

    fold_info <- get_fold_split_info(data, cv_setup, fold)
    cat(sprintf("Train N=%d | Test N=%d\n", fold_info$n_train, fold_info$n_test))

    split    <- split_fold_data(data, cv_setup$folds, fold)
    matrices <- construct_matrices(split$train, split$test, config)
    cat("W_train:", dim(matrices$W_train), "| k_grm:", round(matrices$k_ah, 4), "\n\n")

    # Helper: append AUC to per-fold log line for binary trait
    auc_tag <- function(res) {
      if (!is_binary) return("")
      sprintf(" | AUC_train=%.4f | AUC_test=%.4f",
              res$accuracy$auc_train %||% NA_real_,
              res$accuracy$auc_test  %||% NA_real_)
    }

    # GBLUP
    cat("Running GBLUP...\n")
    t0 <- Sys.time()
    gblup_res <- run_gblup(matrices, split, config)
    gblup_res$runtime <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    r_MG_gblup <- if (!is.na(h2_full$gblup) && h2_full$gblup > 0.01)
      gblup_res$accuracy$test_pheno / sqrt(h2_full$gblup) else NA
    gblup_res$accuracy$test_gebv <- r_MG_gblup
    cat(sprintf("GBLUP train=%.4f | test=%.4f | r_MG=%.4f | h2=%.4f%s | %.2fh\n\n",
                gblup_res$accuracy$train_pheno, gblup_res$accuracy$test_pheno,
                ifelse(is.na(r_MG_gblup), NA, r_MG_gblup),
                gblup_res$varcomp$h2, auc_tag(gblup_res), gblup_res$runtime))
    gblup_varcomp <- gblup_res$varcomp

    # BayesR
    cat("Running BayesR...\n")
    t0 <- Sys.time()
    bayesR_res <- run_bayesR(matrices, split, gblup_varcomp, config, fold = fold)
    bayesR_res$runtime <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    r_MG_bayesR <- if (!is.na(h2_full$bayesR) && h2_full$bayesR > 0.01)
      bayesR_res$accuracy$test_pheno / sqrt(h2_full$bayesR) else NA
    bayesR_res$accuracy$test_gebv <- r_MG_bayesR
    cat(sprintf("BayesR train=%.4f | test=%.4f | r_MG=%.4f | h2=%.4f%s | ESS=%.1f | %.2fh\n\n",
                bayesR_res$accuracy$train_pheno, bayesR_res$accuracy$test_pheno,
                ifelse(is.na(r_MG_bayesR), NA, r_MG_bayesR),
                bayesR_res$varcomp$h2, auc_tag(bayesR_res),
                bayesR_res$diagnostics$mean_ess, bayesR_res$runtime))

    # BayesA
    cat("Running BayesA...\n")
    t0 <- Sys.time()
    bayesA_res <- run_bayesA(matrices, split, gblup_varcomp, config, fold = fold)
    bayesA_res$runtime <- as.numeric(difftime(Sys.time(), t0, units = "hours"))
    r_MG_bayesA <- if (!is.na(h2_full$bayesA) && h2_full$bayesA > 0.01)
      bayesA_res$accuracy$test_pheno / sqrt(h2_full$bayesA) else NA
    bayesA_res$accuracy$test_gebv <- r_MG_bayesA
    cat(sprintf("BayesA train=%.4f | test=%.4f | r_MG=%.4f | h2=%.4f%s | ESS=%.1f | %.2fh\n\n",
                bayesA_res$accuracy$train_pheno, bayesA_res$accuracy$test_pheno,
                ifelse(is.na(r_MG_bayesA), NA, r_MG_bayesA),
                bayesA_res$varcomp$h2, auc_tag(bayesA_res),
                bayesA_res$diagnostics$mean_ess,
                bayesA_res$runtime))

    fold_runtime <- as.numeric(difftime(Sys.time(), fold_start, units = "hours"))
    cat(sprintf("========== FOLD %d DONE | %.2fh ==========\n", fold, fold_runtime))

    sink(type = "message"); sink(type = "output"); close(log_con)
    gc()

    list(
      fold   = fold,
      gblup  = gblup_res,
      bayesR = bayesR_res,
      bayesA = bayesA_res,
      runtime = list(
        fold_total = fold_runtime,
        gblup      = gblup_res$runtime,
        bayesR     = bayesR_res$runtime,
        bayesA     = bayesA_res$runtime
      )
    )

  }, error = function(e) {
    sink(type = "message"); sink(type = "output")
    if (exists("log_con")) close(log_con)
    list(fold = fold, error = e$message)
  })
}

stopCluster(cl)

cat(sprintf("\n=== Parallel Complete: %.2f hours ===\n",
            as.numeric(difftime(Sys.time(), pipeline_start_time, units = "hours"))))

# Check failures
failed <- sapply(all_results, function(x) !is.null(x$error))
if (any(failed)) {
  cat("Failed folds:", which(failed), "\n")
  for (i in which(failed)) cat("  Fold", i, ":", all_results[[i]]$error, "\n")
  all_results <- all_results[!failed]
}
cat("Completed", length(all_results), "/", config$k_folds, "folds\n\n")

# Summary log per model
# - Continuous: train/test = Pearson r on phenotype scale (label: r_MP)
# - Binary: train/test = Pearson r on observed (probability) scale, plus
#   separate AUC_train/AUC_test from rank-based ROC.
summarize_model <- function(results, model_name) {
  is_bin <- !is.null(config$response_type) && config$response_type == "binary"
  vals <- Filter(Negate(is.null), lapply(results, function(x) {
    if (!is.null(x$error)) return(NULL)
    r <- x[[model_name]]
    if (is.null(r) || !is.list(r)) return(NULL)
    if (is.null(r$accuracy) || is.null(r$varcomp)) return(NULL)
    r
  }))
  if (length(vals) == 0) { cat("No valid results for", model_name, "\n"); return(invisible(NULL)) }

  cat(sprintf("\n=== %s | %s | %s ===\n",
              toupper(model_name), config$phenotype_trait, config$genotype_file))
  for (i in seq_along(vals)) {
    v <- vals[[i]]
    base <- sprintf("  Fold %2d: r_train=%.4f | r_test=%.4f | r_MG=%.4f | h2=%.4f",
                    i, v$accuracy$train_pheno, v$accuracy$test_pheno,
                    ifelse(is.na(v$accuracy$test_gebv), NA, v$accuracy$test_gebv),
                    v$varcomp$h2)
    auc_part <- if (is_bin)
      sprintf(" | AUC_train=%.4f | AUC_test=%.4f",
              v$accuracy$auc_train %||% NA_real_,
              v$accuracy$auc_test  %||% NA_real_)
    else ""
    ess_part <- if (!is.null(v$diagnostics$mean_ess))
      sprintf(" | ESS=%.1f", v$diagnostics$mean_ess) else ""
    cat(base, auc_part, ess_part, "\n", sep = "")
  }
  test_vals <- sapply(vals, function(v) v$accuracy$test_pheno)
  rMG_vals  <- sapply(vals, function(v) v$accuracy$test_gebv)
  h2_vals   <- sapply(vals, function(v) v$varcomp$h2)
  cat(sprintf("  Mean+/-SD: r=%.4f+/-%.4f | r_MG=%.4f+/-%.4f | h2=%.4f+/-%.4f",
              mean(test_vals, na.rm = TRUE), sd(test_vals, na.rm = TRUE),
              mean(rMG_vals,  na.rm = TRUE), sd(rMG_vals,  na.rm = TRUE),
              mean(h2_vals,   na.rm = TRUE), sd(h2_vals,   na.rm = TRUE)))
  if (is_bin) {
    auc_test_vals <- sapply(vals, function(v) v$accuracy$auc_test)
    cat(sprintf(" | AUC_test=%.4f+/-%.4f",
                mean(auc_test_vals, na.rm = TRUE),
                sd(auc_test_vals,   na.rm = TRUE)))
  }
  cat("\n")
}

summarize_model(all_results, "gblup")
summarize_model(all_results, "bayesR")
summarize_model(all_results, "bayesA")

# Build summary
build_summary <- function(results, model_name) {
  vals <- Filter(Negate(is.null), lapply(results, function(x) {
    if (!is.null(x$error)) return(NULL)
    r <- x[[model_name]]
    if (is.null(r) || !is.list(r)) return(NULL)
    if (is.null(r$accuracy) || is.null(r$varcomp)) return(NULL)
    r
  }))
  if (length(vals) == 0) return(NULL)
  is_bin <- !is.null(config$response_type) && config$response_type == "binary"
  out <- list(
    train_acc_mean = mean(sapply(vals, function(v) v$accuracy$train_pheno), na.rm=TRUE),
    train_acc_sd   = sd(sapply(vals,   function(v) v$accuracy$train_pheno), na.rm=TRUE),
    test_acc_mean  = mean(sapply(vals, function(v) v$accuracy$test_pheno),  na.rm=TRUE),
    test_acc_sd    = sd(sapply(vals,   function(v) v$accuracy$test_pheno),  na.rm=TRUE),
    r_MG_mean      = mean(sapply(vals, function(v) v$accuracy$test_gebv),   na.rm=TRUE),
    r_MG_sd        = sd(sapply(vals,   function(v) v$accuracy$test_gebv),   na.rm=TRUE),
    h2_mean        = mean(sapply(vals, function(v) v$varcomp$h2),           na.rm=TRUE),
    h2_sd          = sd(sapply(vals,   function(v) v$varcomp$h2),           na.rm=TRUE),
    h2_full        = h2_full[[model_name]],
    # diagnostics$mean_ess can be NULL for folds that skip diagnostics (e.g.
    # GBLUP) or whose Bayesian fit returned without an ESS field; sapply on
    # mixed NULL/numeric collapses to a list and trips mean.default. vapply
    # with NA_real_ default keeps the output type numeric and the log clean.
    ess_mean       = mean(vapply(vals, function(v) {
                            e <- v$diagnostics$mean_ess
                            if (is.null(e) || !is.numeric(e)) NA_real_ else as.numeric(e)[1]
                          }, NA_real_), na.rm = TRUE),
    runtime_mean   = mean(sapply(vals, function(v) v$runtime),              na.rm=TRUE)
  )
  if (is_bin) {
    out$auc_train_mean <- mean(sapply(vals, function(v) v$accuracy$auc_train), na.rm=TRUE)
    out$auc_train_sd   <- sd(sapply(vals,   function(v) v$accuracy$auc_train), na.rm=TRUE)
    out$auc_test_mean  <- mean(sapply(vals, function(v) v$accuracy$auc_test),  na.rm=TRUE)
    out$auc_test_sd    <- sd(sapply(vals,   function(v) v$accuracy$auc_test),  na.rm=TRUE)
  }
  out
}

summary_list <- list(
  gblup  = build_summary(all_results, "gblup"),
  bayesR = build_summary(all_results, "bayesR"),
  bayesA = build_summary(all_results, "bayesA")
)

# Save RDS
rep_id   <- config$rep_id %||% 1
rds_path <- file.path(config$output_dir, sprintf("results_cv_rep%02d.rds", rep_id))

output_rds <- list(
  config       = config,
  cv_setup     = cv_setup,
  rep_id       = rep_id,
  fold_results = setNames(all_results, paste0("fold_", sapply(all_results, `[[`, "fold"))),
  summary      = summary_list,
  runtime      = list(
    total_hours = as.numeric(difftime(Sys.time(), pipeline_start_time, units = "hours")),
    per_fold    = sapply(all_results, function(x) x$runtime$fold_total)
  )
)

saveRDS(output_rds, rds_path)
cat(sprintf("\nResults saved: %s\n", rds_path))

pipeline_end_time <- Sys.time()
cat(sprintf("\n=== Pipeline Complete | Total: %.2f hours ===\n",
            as.numeric(difftime(pipeline_end_time, pipeline_start_time, units = "hours"))))
cat("Started:", as.character(pipeline_start_time), "\n")
cat("Ended:  ", as.character(pipeline_end_time), "\n")