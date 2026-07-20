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
cat("Output:", config$output_dir, "| Folds:", config$k_folds, "\n")
cat("============================\n\n")

pipeline_start_time <- Sys.time()
setwd(config$base_dir)

source("R/data_prep_mh.R")
source("R/matrix_construction.R")
source("R/cv_setup.R")
source("R/models/gblup.R")
source("R/models/bayesR.R")

dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(config$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

is_binary <- !is.null(config$response_type) && config$response_type == "binary"

cat("Loading genomic data...\n")
data <- load_genomic_data(config) # R/data_preparation.R
cat("Total individuals:", nrow(data$phenotypes), "\n")
cat("Total haplotype markers:", ncol(data$haplotypes) - 1, "\n\n")

cat("Setting up cross-validation...\n")
cv_setup <- setup_stratified_cv(data, config) # R/cv_setup.R
cat("CV method:", cv_setup$method, "| Folds:", config$k_folds, "\n")
cat("Fold distribution:\n"); print(table(cv_setup$folds)); cat("\n")

# ── Full data fit ──
base_output_dir <- config$base_output_dir %||% config$output_dir
full_fit_path   <- file.path(base_output_dir, "full_fit.rds")

if (file.exists(full_fit_path)) {
  cat("Full fit found at:", full_fit_path, "— loading...\n")
  full_fit <- readRDS(full_fit_path)
  h2_full  <- full_fit$h2_full
  cat(sprintf("h2_full: GBLUP=%.4f | BayesR=%.4f\n\n",
              h2_full$gblup %||% NA, h2_full$bayesR))
} else {
  cat("\n=== Full Data Fit ===\n")
  full_start <- Sys.time()

  hap_full         <- data$haplotypes %>% select(-ID)
  allele_freq_full <- compute_allele_frequencies(hap_full, config)

  full_split <- list(
    train = list(phenotypes = data$phenotypes, haplotypes = data$haplotypes),
    test  = list(phenotypes = data$phenotypes[0, ], haplotypes = data$haplotypes[0, ])
  )
  matrices_full <- construct_matrices(full_split$train, full_split$test, allele_freq_full, config) #R/matrix_construction.R

  grm_stats_full <- list(
    mean_diag    = mean(diag(matrices_full$A_train)),
    mean_offdiag = mean(matrices_full$A_train[lower.tri(matrices_full$A_train)]),
    range_diag   = range(diag(matrices_full$A_train)),
    k_ah         = matrices_full$k_ah
  )
  cat("GRM full | mean_diag:", round(grm_stats_full$mean_diag, 4),
      "| mean_offdiag:", round(grm_stats_full$mean_offdiag, 4), "\n\n")

  # Append AUC_train to the full-fit log line (binary only)
  full_auc_tag <- function(res) {
    if (!is_binary) return("")
    sprintf(" | AUC_train=%.4f", res$accuracy$auc_train %||% NA_real_)
  }

  cat("Running GBLUP full...\n")
  t0_gblup  <- proc.time()
  gblup_full <- run_gblup(matrices_full, full_split, config)
  rt_gblup  <- (proc.time() - t0_gblup)[["elapsed"]]
  cat(sprintf("GBLUP full: h2=%.4f | train_acc=%.4f%s\n\n",
              gblup_full$varcomp$h2, gblup_full$accuracy$train_pheno,
              full_auc_tag(gblup_full)))
  gblup_varcomp_full <- gblup_full$varcomp

  config_full <- config
  config_full$save_mcmc_samples <- TRUE

  # ── GWAS map (MH): per-trait file at config$coords_mh_file ──
  gwas_enabled <- isTRUE(config$gwas$enabled)
  if (gwas_enabled) {
    block_id_W_full <- attr(matrices_full$W_train, "block_id")
    if (is.null(block_id_W_full)) {
      stop("attr(W_full, 'block_id') is NULL — masbayes v0.5.0 required for GWAS.")
    }

    coords_raw <- read_csv(config$coords_mh_file, show_col_types = FALSE) %>%
      select(block_id, chr, start_pos, end_pos) %>%
      mutate(chr = as.integer(chr))

    # Chr-aware positional translation: W namespace → coords namespace.
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

    attr(matrices_full$W_train, "block_id") <- unname(translation[block_id_W_full])
    block_id_W_full <- attr(matrices_full$W_train, "block_id")

    map_mh_gwas <- coords_raw %>%
      filter(block_id %in% unique(block_id_W_full))
    cat(sprintf("GWAS enabled | W unique blocks: %d | map blocks: %d | windsize: %.0f bp\n",
                length(unique_w_blocks), nrow(map_mh_gwas), config$gwas$windsize))
  } else {
    map_mh_gwas     <- NULL
    block_id_W_full <- NULL
    cat("GWAS disabled (config$gwas$enabled != TRUE)\n")
  }

  cat("Running BayesR full...\n")
  t0_bayesR  <- proc.time()
  bayesR_full <- run_bayesR(
    matrices_full, full_split, gblup_varcomp_full, config_full, fold = 0,
    map      = if (gwas_enabled) map_mh_gwas          else NULL,
    windsize = if (gwas_enabled) config$gwas$windsize else NULL
  )
  rt_bayesR  <- (proc.time() - t0_bayesR)[["elapsed"]]
  cat(sprintf("BayesR full: h2=%.4f | train_acc=%.4f%s\n\n",
              bayesR_full$varcomp$h2, bayesR_full$accuracy$train_pheno,
              full_auc_tag(bayesR_full)))

  # ── GWAS Rds (kept separate from full_fit.rds) ──
  if (gwas_enabled && !is.null(bayesR_full$pip)) {
    max_pip_block <- vapply(map_mh_gwas$block_id, function(b) {
      idx <- which(block_id_W_full == b)
      if (length(idx) == 0L) NA_real_ else max(bayesR_full$pip[idx])
    }, numeric(1L))

    gwas_path <- file.path(base_output_dir, "gwas_full.rds")
    dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
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
        trait         = config$phenotype_trait,
        rep_id        = config$rep_id %||% 1
      ),
      gwas_path
    )
    cat(sprintf("GWAS RDS saved: %s\n\n", gwas_path))
  }

  h2_full <- list(
    gblup  = gblup_full$varcomp$h2,
    bayesR = bayesR_full$varcomp$h2
  )

  # Both models fit in a single full_fit.
  diagnostics_fit <- list(
    gblup  = NULL,
    bayesR = bayesR_full$diagnostics
  )

  runtime_fit <- list(
    gblup  = rt_gblup,
    bayesR = rt_bayesR
  )

  full_fit <- list(
    grm_stats = grm_stats_full,
    gblup     = gblup_full,
    bayesR    = bayesR_full,
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
  .packages = c("tidyverse", "coda", "masbayes", "masreml", "pROC"),
  .errorhandling = "pass",
  .verbose = FALSE
) %dopar% {

  source("R/data_prep_mh.R")
  source("R/matrix_construction.R")
  source("R/cv_setup.R")
  source("R/models/gblup.R")
  source("R/models/bayesR.R")
  select <- dplyr::select; filter <- dplyr::filter
  is_binary <- !is.null(config$response_type) && config$response_type == "binary"

  tryCatch({
    fold_start <- Sys.time()
    log_file <- file.path(config$output_dir, "logs", paste0("fold_", fold, ".log"))
    log_con  <- file(log_file, open = "wt")
    sink(log_con, type = "output"); sink(log_con, type = "message")

    cat(sprintf("\n========== FOLD %d/%d ==========\n", fold, config$k_folds))
    cat("Started:", as.character(Sys.time()), "\n\n")
    cat(sprintf("h2_full: GBLUP=%.4f | BayesR=%.4f\n\n",
                h2_full$gblup %||% NA, h2_full$bayesR))

    fold_info <- get_fold_split_info(data, cv_setup, fold) #R/cv_setup.R
    cat(sprintf("Train N=%d | Test N=%d\n", fold_info$n_train, fold_info$n_test))

    split       <- split_fold_data(data, cv_setup$folds, fold) 
    hap_train   <- split$train$haplotypes %>% select(-ID)
    allele_freq <- compute_allele_frequencies(hap_train, config)

    cat("\nConstructing W matrices...\n")
    matrices <- construct_matrices(split$train, split$test, allele_freq, config) # R/matrix_construction.R
    cat("W_train:", dim(matrices$W_train), "| k_ah =", round(matrices$k_ah, 4), "\n\n")

    # Append AUC to the per-fold log line (binary only)
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

    fold_runtime <- as.numeric(difftime(Sys.time(), fold_start, units = "hours"))
    cat(sprintf("========== FOLD %d DONE | %.2fh ==========\n", fold, fold_runtime))

    sink(type = "message"); sink(type = "output"); close(log_con)
    gc()

    list(
      fold      = fold,
      gblup     = gblup_res,
      bayesR    = bayesR_res,
      runtime   = list(
        fold_total = fold_runtime,
        gblup      = gblup_res$runtime,
        bayesR     = bayesR_res$runtime
      )
    )

  }, error = function(e) {
    sink(type = "message"); sink(type = "output")
    if (exists("log_con")) close(log_con)
    list(fold = fold, error = e$message)
  })
}

stopCluster(cl)

parallel_end_time <- Sys.time()
cat(sprintf("\n=== Parallel Complete: %.2f hours ===\n",
            as.numeric(difftime(parallel_end_time, pipeline_start_time, units = "hours"))))

# Check failures
failed <- sapply(all_results, function(x) !is.null(x$error))
if (any(failed)) {
  cat("Failed folds:", which(failed), "\n")
  for (i in which(failed)) cat("  Fold", i, ":", all_results[[i]]$error, "\n")
  all_results <- all_results[!failed]
}
cat("Completed", length(all_results), "/", config$k_folds, "folds\n\n")

# Summary log per model. Continuous: Pearson r on phenotype scale (r_MP).
# Binary: Pearson r on probability scale, plus separate AUC from rank-based ROC.
summarize_model <- function(results, model_name) {
  is_bin <- !is.null(config$response_type) && config$response_type == "binary"
  vals <- Filter(Negate(is.null), lapply(results, function(x) {
    r <- x[[model_name]]
    if (is.null(r)) return(NULL)
    list(train     = r$accuracy$train_pheno,
         test      = r$accuracy$test_pheno,
         auc_train = r$accuracy$auc_train,
         auc_test  = r$accuracy$auc_test,
         r_MG      = r$accuracy$test_gebv,
         h2        = r$varcomp$h2,
         ess       = r$diagnostics$mean_ess)
  }))

  cat(sprintf("\n=== %s | %s ===\n", toupper(model_name), config$phenotype_trait))
  for (i in seq_along(vals)) {
    v <- vals[[i]]
    base <- sprintf("  Fold %2d: r_train=%.4f | r_test=%.4f | r_MG=%.4f | h2=%.4f",
                    i, v$train, v$test,
                    ifelse(is.na(v$r_MG), NA, v$r_MG), v$h2)
    auc_part <- if (is_bin)
      sprintf(" | AUC_train=%.4f | AUC_test=%.4f",
              v$auc_train %||% NA_real_, v$auc_test %||% NA_real_)
    else ""
    ess_part <- if (!is.null(v$ess)) sprintf(" | ESS=%.1f", v$ess) else ""
    cat(base, auc_part, ess_part, "\n", sep = "")
  }

  test_vals     <- sapply(vals, `[[`, "test")
  rMG_vals      <- sapply(vals, `[[`, "r_MG")
  h2_vals       <- sapply(vals, `[[`, "h2")
  cat(sprintf("  Mean+/-SD: r=%.4f+/-%.4f | r_MG=%.4f+/-%.4f | h2=%.4f+/-%.4f",
              mean(test_vals, na.rm = TRUE), sd(test_vals, na.rm = TRUE),
              mean(rMG_vals,  na.rm = TRUE), sd(rMG_vals,  na.rm = TRUE),
              mean(h2_vals,   na.rm = TRUE), sd(h2_vals,   na.rm = TRUE)))
  if (is_bin) {
    auc_test_vals <- sapply(vals, `[[`, "auc_test")
    cat(sprintf(" | AUC_test=%.4f+/-%.4f",
                mean(auc_test_vals, na.rm = TRUE),
                sd(auc_test_vals,   na.rm = TRUE)))
  }
  cat("\n")
}

summarize_model(all_results, "gblup")
summarize_model(all_results, "bayesR")

# Build summary
build_summary <- function(results, model_name) {
  vals <- Filter(Negate(is.null), lapply(results, function(x) {
    if (!is.null(x$error)) return(NULL)
    r <- x[[model_name]]
    if (is.null(r) || !is.list(r)) return(NULL)
    if (is.null(r$accuracy) || is.null(r$varcomp)) return(NULL)
    r
  }))
  if (length(vals) == 0) { cat("No valid results for", model_name, "\n"); return(invisible(NULL)) }
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
  bayesR = build_summary(all_results, "bayesR")
)

rep_id   <- config$rep_id %||% 1
rds_path <- file.path(config$output_dir, sprintf("results_cv_rep%02d.rds", rep_id))

# Save RDS
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