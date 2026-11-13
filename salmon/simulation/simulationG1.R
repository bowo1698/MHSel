suppressPackageStartupMessages({
  library(AlphaSimR)
  library(tidyverse)
  library(vcfR)
})

# ── Config ────────────────────────────────────────────────────────────────────
CODE        <- Sys.getenv("CODE")
RAW_DIR     <- Sys.getenv("RAW_DIR")
OUT         <- file.path(Sys.getenv("RUN_ROOT"), "salmon", "output_G1")
stopifnot(nzchar(CODE), nzchar(RAW_DIR), nzchar(Sys.getenv("RUN_ROOT")))
source(file.path(CODE, "salmon", "simulation", "sim_helpers.R"))
source(file.path(CODE, "common", "manifest.R"))
VCF_FILE    <- file.path(RAW_DIR, "salmon_60k_pairs125bp.vcf.gz")
MH_STATS    <- file.path(RAW_DIR, "mh_info_ld_haploblock_G0/stats")
HAP_STATS   <- file.path(RAW_DIR, "all_info_haploblock_G0/stats")

args        <- commandArgs(trailingOnly = TRUE)
ITER        <- as.integer(args[1])
N_FAM_PHENO <- 20L   # phenotyped reference offspring per family
N_FAM_CAND  <- 5L    # unphenotyped selection candidates per family
N_OFF_FAM   <- N_FAM_PHENO + N_FAM_CAND
CM_PER_MB   <- 1     # uniform recombination rate
N_QTL       <- 100
H2          <- 0.3

# QTL architecture
K_BLK_QTL          <- 2L     # causal haplotype-alleles per block QTL (S3-S6)
RARE_PROP          <- 0.30
MAF_RARE_MIN       <- 0.01
MAF_RARE_MAX       <- 0.05

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RDS_DIR  <- file.path(OUT, "rds")
ITER_DIR <- file.path(OUT, sprintf("iter%02d", ITER))
dir.create(RDS_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(ITER_DIR, recursive = TRUE, showWarnings = FALSE)

# ── 1. Read VCF → MapPop ──────────────────────────────────────────────────────
vcf       <- read.vcfR(VCF_FILE, verbose = FALSE)

gt        <- extract.gt(vcf, element = "GT")
haplo_mat <- do.call(cbind, lapply(1:ncol(gt), function(i) {
  m <- do.call(rbind, strsplit(gt[, i], "\\|"))
  matrix(as.integer(m), nrow = nrow(m), ncol = 2)
}))

# Genotype as-is: the panel already carries natural rare variants (MAF < 0.05).

gen_map      <- data.frame(chr = as.integer(vcf@fix[, "CHROM"]),
                            pos = as.integer(vcf@fix[, "POS"]))
n_ind_vcf    <- ncol(gt)
chrs_ordered <- sort(unique(gen_map$chr))

# Morgans at CM_PER_MB (previous version used Mb as Morgans, i.e. 100 cM/Mb)
genMap_list <- build_gen_map(gen_map$chr, gen_map$pos, CM_PER_MB)
map_length_M <- sum(vapply(genMap_list, max, numeric(1)))

haplo_list <- lapply(chrs_ordered, function(c) {
  snp_idx <- which(gen_map$chr == c)
  t(haplo_mat[snp_idx, , drop = FALSE])
})

founderPop <- newMapPop(genMap = genMap_list, haplotypes = haplo_list, inbred = FALSE)

vcf_fix  <- as.data.frame(getFIX(vcf))
# global_idx preserves VCF row order through arrange, so haplo_mat can still
# be indexed by snp_id when extracting per-block haplotype patterns.
snp_meta <- tibble(snp_id     = vcf_fix$ID,
                   chr        = as.integer(vcf_fix$CHROM),
                   pos        = as.integer(vcf_fix$POS),
                   global_idx = seq_len(nrow(vcf_fix))) %>%
            arrange(chr, pos)
n_snp    <- nrow(snp_meta)

# G0 SNP allele frequencies + MAF for rare/common stratification
snp_freq_g0 <- setNames(rowMeans(haplo_mat), vcf_fix$ID)[snp_meta$snp_id]
maf_g0      <- pmin(snp_freq_g0, 1 - snp_freq_g0)
is_rare_snp <- maf_g0 >= MAF_RARE_MIN & maf_g0 < MAF_RARE_MAX

# Within-chr SNP index lookup for joining with block_snp_map
snp_meta_idx_df <- snp_meta %>%
  group_by(chr) %>%
  mutate(snp_index = row_number() - 1L) %>%
  ungroup()

# ── 2. Load block info ────────────────────────────────────────────────────────
mh_coords  <- read_csv(file.path(MH_STATS,  "microhaplotype_coordinates.csv"), show_col_types = FALSE)
hap_coords <- read_csv(file.path(HAP_STATS, "microhaplotype_coordinates.csv"), show_col_types = FALSE)
mh_snps    <- read_csv(file.path(MH_STATS,  "snp_selection_detailed.csv"),     show_col_types = FALSE)
hap_snps   <- read_csv(file.path(HAP_STATS, "snp_selection_detailed.csv"),     show_col_types = FALSE)

# ── 3. Define QTL effects from G0 ─────────────────────────────────────────────
set.seed(123)

n_qtl_rare   <- round(N_QTL * RARE_PROP)
n_qtl_common <- N_QTL - n_qtl_rare

# S1/S2: stratified SNP QTL (30 rare + 70 common, deterministic)
rare_snp_pool   <- snp_meta$snp_id[is_rare_snp]
common_snp_pool <- snp_meta$snp_id[!is_rare_snp & maf_g0 >= MAF_RARE_MAX]
if (length(rare_snp_pool)   < n_qtl_rare   ||
    length(common_snp_pool) < n_qtl_common)
  stop("Not enough SNPs in rare or common pool; lower RARE_PROP or relax MAF bounds.")

qtl_snp_ids <- c(sample(rare_snp_pool,   n_qtl_rare,   replace = FALSE),
                 sample(common_snp_pool, n_qtl_common, replace = FALSE))

define_snp_effects <- function(qtl_snp_ids, type = "equal") {
  n   <- length(qtl_snp_ids)
  mag <- if (type == "equal") rep(1, n) else rgamma(n, shape = 0.4, scale = 1)
  sgn <- sample(c(-1, 1), n, replace = TRUE)
  eff <- mag * sgn
  tibble(snp_id = qtl_snp_ids, effect = eff / sqrt(sum(eff^2)))
}
qtl_s1 <- define_snp_effects(qtl_snp_ids, "equal")
qtl_s2 <- define_snp_effects(qtl_snp_ids, "unequal")

# S3-S6: block QTL classification by rare/common haplotype-allele containment
classify_blocks_by_rare <- function(block_snp_map, snp_meta_idx_df, is_rare_named) {
  joined <- block_snp_map %>%
    inner_join(snp_meta_idx_df %>% select(chr, snp_index, snp_id_meta = snp_id),
               by = c("chr", "snp_index")) %>%
    mutate(is_rare = unname(is_rare_named[snp_id_meta]))
  by_blk <- joined %>%
    group_by(block_id) %>%
    summarise(any_rare = any(is_rare, na.rm = TRUE), .groups = "drop")
  list(rare_containing = by_blk$block_id[by_blk$any_rare],
       common_only     = by_blk$block_id[!by_blk$any_rare])
}

is_rare_named <- setNames(is_rare_snp, snp_meta$snp_id)
mh_blk_pools  <- classify_blocks_by_rare(mh_snps,  snp_meta_idx_df, is_rare_named)
hap_blk_pools <- classify_blocks_by_rare(hap_snps, snp_meta_idx_df, is_rare_named)

if (length(mh_blk_pools$rare_containing)  < n_qtl_rare   ||
    length(mh_blk_pools$common_only)      < n_qtl_common ||
    length(hap_blk_pools$rare_containing) < n_qtl_rare   ||
    length(hap_blk_pools$common_only)     < n_qtl_common)
  stop("Insufficient blocks in rare or common pool — lower RARE_PROP or N_QTL")

qtl_mh_blocks  <- c(sample(mh_blk_pools$rare_containing,  n_qtl_rare,   replace = FALSE),
                    sample(mh_blk_pools$common_only,      n_qtl_common, replace = FALSE))
qtl_hap_blocks <- c(sample(hap_blk_pools$rare_containing, n_qtl_rare,   replace = FALSE),
                    sample(hap_blk_pools$common_only,     n_qtl_common, replace = FALSE))

# Encode a K-SNP haplotype as sum(allele_k * 2^(k-1)).
# Stable across generations while snp_meta order is preserved.
encode_block_haplotype <- function(hap_vec) {
  as.integer(sum(as.integer(hap_vec) * 2L^(seq_along(hap_vec) - 1L)))
}

# Haplotype-additive QTL per block: take segregating haplotype-alleles at G0,
# pick K causal ones, assign β. TBV = Σ β_a · count(allele_a in ind).
# This captures cis/trans — AD and BC differ in TBV at identical SNP dosages.
define_block_qtl_haplotype <- function(block_ids, block_snp_map, snp_meta_idx_df,
                                       haplo_mat_g0, k_alleles_per_block,
                                       type = "equal", is_rare_named = NULL) {
  joined <- block_snp_map %>%
    filter(block_id %in% block_ids) %>%
    inner_join(snp_meta_idx_df %>% select(chr, snp_index, snp_id_meta = snp_id),
               by = c("chr", "snp_index"))

  res <- vector("list", length(block_ids))
  names(res) <- block_ids

  for (blk in block_ids) {
    blk_snps <- joined %>% filter(block_id == blk)
    if (nrow(blk_snps) < 1L) {
      res[[blk]] <- list(chr = NA_integer_, snp_ids = character(0),
                         snp_sorted_idx = integer(0), causal_alleles = integer(0),
                         allele_freqs = numeric(0), beta = numeric(0))
      next
    }
    chr_blk        <- blk_snps$chr[1]
    snp_ids_blk    <- blk_snps$snp_id_meta
    snp_sorted_idx <- match(snp_ids_blk, snp_meta_idx_df$snp_id)

    haps_at_block <- haplo_mat_g0[snp_sorted_idx, , drop = FALSE]
    hap_codes_g0  <- apply(haps_at_block, 2, encode_block_haplotype)

    allele_tbl       <- table(hap_codes_g0)
    allele_codes_all <- as.integer(names(allele_tbl))
    allele_freqs_all <- as.numeric(allele_tbl) / length(hap_codes_g0)

    if (length(allele_codes_all) < 2L) {
      res[[blk]] <- list(chr = chr_blk, snp_ids = snp_ids_blk,
                         snp_sorted_idx = snp_sorted_idx,
                         causal_alleles = integer(0),
                         allele_freqs   = numeric(0), beta = numeric(0))
      next
    }

    k_actual       <- min(k_alleles_per_block, length(allele_codes_all))
    is_rare_allele <- allele_freqs_all >= MAF_RARE_MIN &
                      allele_freqs_all < MAF_RARE_MAX

    # Force one rare haplotype-allele causal when the block has any.
    if (!is.null(is_rare_named) && any(is_rare_allele) && k_actual >= 2L) {
      rare_idx    <- which(is_rare_allele)
      forced_rare <- if (length(rare_idx) == 1L) rare_idx
                     else sample(rare_idx, 1L)
      remaining   <- setdiff(seq_along(allele_codes_all), forced_rare)
      rest        <- if (length(remaining) <= (k_actual - 1L)) remaining
                     else sample(remaining, k_actual - 1L)
      chosen_idx  <- c(forced_rare, rest)
    } else {
      chosen_idx  <- if (k_actual == length(allele_codes_all))
                       seq_along(allele_codes_all)
                     else
                       sample(seq_along(allele_codes_all), k_actual)
    }

    beta_mag <- if (type == "equal") rep(1, length(chosen_idx))
                else                 rgamma(length(chosen_idx), shape = 0.4, scale = 1)
    beta_sign <- sample(c(-1, 1), length(chosen_idx), replace = TRUE)

    res[[blk]] <- list(
      chr            = chr_blk,
      snp_ids        = snp_ids_blk,
      snp_sorted_idx = snp_sorted_idx,
      causal_alleles = allele_codes_all[chosen_idx],
      allele_freqs   = allele_freqs_all[chosen_idx],
      beta           = beta_mag * beta_sign
    )
  }

  # Global normalisation across all causal effects so ‖β‖² = 1
  all_eff    <- unlist(lapply(res, `[[`, "beta"))
  total_norm <- sqrt(sum(all_eff^2))
  if (total_norm > 0) {
    for (blk in names(res)) {
      if (length(res[[blk]]$beta) > 0)
        res[[blk]]$beta <- res[[blk]]$beta / total_norm
    }
  }
  res
}

# Reorder haplo_mat rows to match sorted snp_meta, so row i is snp_meta$snp_id[i].
haplo_mat_sorted <- haplo_mat[snp_meta$global_idx, , drop = FALSE]

qtl_s3_effects <- define_block_qtl_haplotype(qtl_mh_blocks,  mh_snps,  snp_meta_idx_df,
                                             haplo_mat_sorted, K_BLK_QTL, "equal",
                                             is_rare_named = is_rare_named)
qtl_s4_effects <- define_block_qtl_haplotype(qtl_mh_blocks,  mh_snps,  snp_meta_idx_df,
                                             haplo_mat_sorted, K_BLK_QTL, "unequal",
                                             is_rare_named = is_rare_named)
qtl_s5_effects <- define_block_qtl_haplotype(qtl_hap_blocks, hap_snps, snp_meta_idx_df,
                                             haplo_mat_sorted, K_BLK_QTL, "equal",
                                             is_rare_named = is_rare_named)
qtl_s6_effects <- define_block_qtl_haplotype(qtl_hap_blocks, hap_snps, snp_meta_idx_df,
                                             haplo_mat_sorted, K_BLK_QTL, "unequal",
                                             is_rare_named = is_rare_named)

# ── 4. TBV calculation functions ──────────────────────────────────────────────
calc_tbv_snp <- function(dosage_mat, qtl_effects) {
  cols <- intersect(qtl_effects$snp_id, colnames(dosage_mat))
  eff  <- qtl_effects$effect[match(cols, qtl_effects$snp_id)]
  tbv  <- as.matrix(dosage_mat[, cols]) %*% eff
  as.vector(tbv - mean(tbv))
}

# Block TBV. haplo_mat is n_snp × (2*n_ind) in sorted snp_meta order. Per
# block: encode the haplotype at causal SNPs, count copies per individual, × β.
calc_tbv_block <- function(haplo_mat, qtl_block_list) {
  n_alleles_total <- ncol(haplo_mat)
  if (n_alleles_total %% 2L != 0L)
    stop("calc_tbv_block: ncol(haplo_mat) must be even (2 haplotypes per ind)")
  n_ind <- n_alleles_total %/% 2L
  tbv   <- numeric(n_ind)

  for (blk in names(qtl_block_list)) {
    q <- qtl_block_list[[blk]]
    if (length(q$causal_alleles) == 0L || length(q$snp_sorted_idx) < 1L) next
    if (max(q$snp_sorted_idx) > nrow(haplo_mat)) next

    haps_at_block <- haplo_mat[q$snp_sorted_idx, , drop = FALSE]
    hap_codes     <- apply(haps_at_block, 2, function(v)
                            sum(as.integer(v) * 2L^(seq_along(v) - 1L)))
    pat_codes <- hap_codes[seq(1L, length(hap_codes), by = 2L)]
    mat_codes <- hap_codes[seq(2L, length(hap_codes), by = 2L)]
    for (i in seq_along(q$causal_alleles)) {
      a      <- q$causal_alleles[i]
      counts <- (pat_codes == a) + (mat_codes == a)
      tbv    <- tbv + q$beta[i] * counts
    }
  }
  as.vector(tbv - mean(tbv))
}

# ── 5. Simulation loop ────────────────────────────────────────────────────────
SP <- SimParam$new(founderPop)
n_snp_per_chr <- snp_meta %>% count(chr) %>% pull(n)
SP$addSnpChip(nSnpPerChr = n_snp_per_chr)

all_geno     <- list()
all_pheno    <- list(S1 = list(), S2 = list(), S3 = list(),
                     S4 = list(), S5 = list(), S6 = list())
all_pedigree <- list()

for (iter in ITER) {
  set.seed(iter * 100)

  G0 <- newPop(founderPop, simParam = SP)

  n_founders <- G0@nInd
  sires <- 1:floor(n_founders / 2)
  dams  <- (floor(n_founders / 2) + 1):n_founders
  crosses <- data.frame(
    sire = sample(sires, min(length(sires), length(dams)), replace = FALSE),
    dam  = sample(dams,  min(length(sires), length(dams)), replace = FALSE)
  )

  G1_list <- lapply(seq_len(nrow(crosses)), function(f) {
    makeCross(G0, crossPlan = matrix(c(crosses$sire[f], crosses$dam[f]), nrow = 1),
              nProgeny = N_OFF_FAM, simParam = SP)
  })
  G1 <- mergePops(G1_list)

  n_G1    <- G1@nInd
  fam_ids  <- rep(seq_len(nrow(crosses)), each = N_OFF_FAM)
  pop_type <- population_type_vector(nrow(crosses), N_FAM_PHENO, N_FAM_CAND)
  ind_ids <- sprintf("iter%d_gen1_fam%03d_ind%02d",
                    iter, fam_ids, sequence(rep(N_OFF_FAM, nrow(crosses))))

  pedigree_g1 <- tibble(
    individual_id   = ind_ids,
    generation      = 1,
    population_type = pop_type,
    mother_id       = rep(sprintf("founder_%d", crosses$dam),  each = N_OFF_FAM),
    father_id       = rep(sprintf("founder_%d", crosses$sire), each = N_OFF_FAM),
    iteration       = iter,
    family_id       = fam_ids
  )
  all_pedigree[[iter]] <- pedigree_g1

  saveRDS(G1, file.path(ITER_DIR, "G1_pop.rds"))

  dosage <- pullSnpGeno(G1, simParam = SP)
  colnames(dosage) <- snp_meta$snp_id
  rownames(dosage) <- ind_ids
  all_geno[[iter]] <- as_tibble(dosage, rownames = "ID") %>%
    mutate(generation = 1, population_type = pop_type, .after = "ID")

  # Phased haplotypes G1: transpose to n_snp × (2*n_ind) for calc_tbv_block
  haplo_mat_g1 <- t(pullSegSiteHaplo(G1, simParam = SP))
  stopifnot(nrow(haplo_mat_g1) == nrow(snp_meta))

  rescale_tbv <- function(tbv) (tbv - mean(tbv)) / sd(tbv)

  tbv_s1 <- rescale_tbv(calc_tbv_snp(dosage, qtl_s1))
  tbv_s2 <- rescale_tbv(calc_tbv_snp(dosage, qtl_s2))
  tbv_s3 <- rescale_tbv(calc_tbv_block(haplo_mat_g1, qtl_s3_effects))
  tbv_s4 <- rescale_tbv(calc_tbv_block(haplo_mat_g1, qtl_s4_effects))
  tbv_s5 <- rescale_tbv(calc_tbv_block(haplo_mat_g1, qtl_s5_effects))
  tbv_s6 <- rescale_tbv(calc_tbv_block(haplo_mat_g1, qtl_s6_effects))

  sim_pheno <- function(tbv, h2 = H2) {
    var_e <- var(tbv) * (1 - h2) / h2
    tbv + rnorm(length(tbv), 0, sqrt(var_e))
  }
  # TBV scaling, residual variance and h2 use all 1,000 G1 individuals;
  # candidate phenotypes are kept for audit only (never used for fitting).
  make_pheno_df <- function(tbv, ind_ids, iter, fam_ids, h2 = H2) {
    pheno   <- sim_pheno(tbv, h2)
    h2_real <- var(tbv) / var(pheno)
    tibble(individual_id   = ind_ids,
           tbv             = tbv,
           phenotype       = pheno,
           generation      = 1,
           population_type = pop_type,
           iteration       = iter,
           family_id       = fam_ids,
           h2_realized     = h2_real)
  }
  all_pheno$S1[[iter]] <- make_pheno_df(tbv_s1, ind_ids, iter, fam_ids)
  all_pheno$S2[[iter]] <- make_pheno_df(tbv_s2, ind_ids, iter, fam_ids)
  all_pheno$S3[[iter]] <- make_pheno_df(tbv_s3, ind_ids, iter, fam_ids)
  all_pheno$S4[[iter]] <- make_pheno_df(tbv_s4, ind_ids, iter, fam_ids)
  all_pheno$S5[[iter]] <- make_pheno_df(tbv_s5, ind_ids, iter, fam_ids)
  all_pheno$S6[[iter]] <- make_pheno_df(tbv_s6, ind_ids, iter, fam_ids)

  cat("Iter", iter, "done | n_G1:", n_G1, "\n")
}

# ── 6. Output ─────────────────────────────────────────────────────────────────
write_csv(all_geno[[ITER]],    file.path(ITER_DIR, "geno_G1.csv"))
walk(names(all_pheno), ~write_csv(all_pheno[[.x]][[ITER]],
       file.path(ITER_DIR, sprintf("pheno_%s_G1.csv", .x))))
write_csv(all_pedigree[[ITER]], file.path(ITER_DIR, "pedigree_G1.csv"))

# QTL design + founderPop + pedigree G0 written only by iter 1
if (ITER == 1) {
  pedigree_g0 <- tibble(
    individual_id   = sprintf("founder_%d", 1:n_ind_vcf),
    generation      = 0,
    population_type = "founder",
    mother_id       = NA, father_id = NA,
    iteration       = NA, family_id = NA
  )
  write_csv(pedigree_g0, file.path(OUT, "pedigree_G0.csv"))
  saveRDS(founderPop,             file.path(RDS_DIR, "founderPop.rds"))
  saveRDS(SP,                     file.path(RDS_DIR, "SP.rds"))
  saveRDS(qtl_s1,                 file.path(RDS_DIR, "qtl_s1.rds"))
  saveRDS(qtl_s2,                 file.path(RDS_DIR, "qtl_s2.rds"))
  saveRDS(qtl_s3_effects,         file.path(RDS_DIR, "qtl_s3_effects.rds"))
  saveRDS(qtl_s4_effects,         file.path(RDS_DIR, "qtl_s4_effects.rds"))
  saveRDS(qtl_s5_effects,         file.path(RDS_DIR, "qtl_s5_effects.rds"))
  saveRDS(qtl_s6_effects,         file.path(RDS_DIR, "qtl_s6_effects.rds"))
  saveRDS(qtl_snp_ids,            file.path(RDS_DIR, "qtl_snp_ids.rds"))
  saveRDS(qtl_mh_blocks,          file.path(RDS_DIR, "qtl_mh_blocks.rds"))
  saveRDS(qtl_hap_blocks,         file.path(RDS_DIR, "qtl_hap_blocks.rds"))
  saveRDS(snp_meta,               file.path(RDS_DIR, "snp_meta.rds"))
  saveRDS(snp_freq_g0,            file.path(RDS_DIR, "snp_freq_g0.rds"))
  saveRDS(maf_g0,                 file.path(RDS_DIR, "maf_g0.rds"))
  saveRDS(list(
    K_BLK_QTL       = K_BLK_QTL,
    RARE_PROP       = RARE_PROP,
    MAF_RARE_MIN    = MAF_RARE_MIN,
    MAF_RARE_MAX    = MAF_RARE_MAX,
    n_qtl_rare      = n_qtl_rare,
    n_qtl_common    = n_qtl_common
  ), file.path(RDS_DIR, "qtl_design_params.rds"))
}

# ── 7. Validation log: confirm block-level haplotype-additive design ─────────
# Per QTL block in S3-S6, fit SNP-additive vs MH-additive lm on tbv_blk.
# R²_MH ≈ 1.0 confirms TBV is haplotype-additive by construction; R²_SNP < 1.0
# shows the block signal is spread across SNPs. The gap grows with the number
# of haplotype-alleles h: ≈ 0 when biallelic, substantial at h ≥ 5.
cat("\n=== Validation iter ", ITER, ": block-level QTL vs SNP-additive bias ===\n", sep = "")

per_block_r2 <- function(qtl_block_list, haplo_mat_g, dosage_mat) {
  n_ind <- ncol(haplo_mat_g) %/% 2L
  rows  <- list()
  for (blk in names(qtl_block_list)) {
    q <- qtl_block_list[[blk]]
    if (length(q$causal_alleles) == 0L || length(q$snp_sorted_idx) < 1L) next

    haps_at_block <- haplo_mat_g[q$snp_sorted_idx, , drop = FALSE]
    hap_codes <- apply(haps_at_block, 2, function(v)
                       sum(as.integer(v) * 2L^(seq_along(v) - 1L)))
    pat <- hap_codes[seq(1L, length(hap_codes), by = 2L)]
    mat <- hap_codes[seq(2L, length(hap_codes), by = 2L)]

    tbv_blk <- numeric(n_ind)
    for (i in seq_along(q$causal_alleles)) {
      counts  <- (pat == q$causal_alleles[i]) + (mat == q$causal_alleles[i])
      tbv_blk <- tbv_blk + q$beta[i] * counts
    }
    if (var(tbv_blk) < 1e-12) next

    snp_cols <- intersect(q$snp_ids, colnames(dosage_mat))
    r2_snp <- NA_real_
    if (length(snp_cols) >= 1L && length(snp_cols) < n_ind - 1L) {
      fit <- tryCatch(lm(tbv_blk ~ as.matrix(dosage_mat[, snp_cols, drop = FALSE])),
                      error = function(e) NULL)
      if (!is.null(fit)) r2_snp <- summary(fit)$r.squared
    }

    all_alleles <- sort(unique(c(pat, mat)))
    r2_mh   <- NA_real_
    h_count <- length(all_alleles)
    if (h_count >= 2L) {
      baseline <- all_alleles[which.max(tabulate(match(c(pat, mat), all_alleles)))]
      kept     <- setdiff(all_alleles, baseline)
      if (length(kept) >= 1L && length(kept) < n_ind - 1L) {
        Z <- vapply(kept, function(a) (pat == a) + (mat == a), numeric(n_ind))
        fit_mh <- tryCatch(lm(tbv_blk ~ Z), error = function(e) NULL)
        if (!is.null(fit_mh)) r2_mh <- summary(fit_mh)$r.squared
      }
    }

    rows[[length(rows) + 1L]] <- data.frame(
      block_id     = blk,
      h_alleles_g1 = h_count,
      var_tbv_blk  = var(tbv_blk),
      r2_snp       = r2_snp,
      r2_mh        = r2_mh,
      gap          = r2_mh - r2_snp,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
}

validation <- list(
  S3 = per_block_r2(qtl_s3_effects, haplo_mat_g1, dosage),
  S4 = per_block_r2(qtl_s4_effects, haplo_mat_g1, dosage),
  S5 = per_block_r2(qtl_s5_effects, haplo_mat_g1, dosage),
  S6 = per_block_r2(qtl_s6_effects, haplo_mat_g1, dosage)
)

# Variance-weighted global R² per scenario
cat(sprintf("%-4s %8s %8s %9s   %s\n", "Scn", "R²_SNP", "R²_MH", "Gap",
            "(weighted by var(TBV_block))"))
for (sc in names(validation)) {
  df <- validation[[sc]]
  if (nrow(df) == 0L) {
    cat(sprintf("%-4s     NA       NA        NA   no blocks\n", sc))
    next
  }
  w      <- df$var_tbv_blk; w <- w / sum(w)
  r2_snp <- sum(w * df$r2_snp, na.rm = TRUE)
  r2_mh  <- sum(w * df$r2_mh,  na.rm = TRUE)
  cat(sprintf("%-4s %8.3f %8.3f %+9.3f\n", sc, r2_snp, r2_mh, r2_mh - r2_snp))
}

# Per-h-bin breakdown (S5 HAP-block — most multi-allelic)
h_bin <- function(h) {
  if (is.na(h))      NA_character_
  else if (h == 2L)  "h=2"
  else if (h <= 4L)  "h=3-4"
  else if (h <= 9L)  "h=5-9"
  else               "h>=10"
}
cat("\nS5 (HAP-block) per-h-bin: gap grows with block multi-allelism\n")
df_s5 <- validation$S5
if (nrow(df_s5) > 0L) {
  df_s5$bin <- vapply(df_s5$h_alleles_g1, h_bin, character(1))
  bin_order <- c("h=2", "h=3-4", "h=5-9", "h>=10")
  for (bin_label in bin_order) {
    d <- df_s5[df_s5$bin == bin_label, , drop = FALSE]
    if (nrow(d) == 0L) next
    cat(sprintf("  %-6s n=%3d  R²_SNP=%.3f  R²_MH=%.3f  gap=%+.3f\n",
                bin_label, nrow(d),
                mean(d$r2_snp, na.rm = TRUE),
                mean(d$r2_mh,  na.rm = TRUE),
                mean(d$gap,    na.rm = TRUE)))
  }
}

cat("\nExpected pattern (= design is haplotype-additive at block level):\n")
cat("  R²_MH ≈ 1.000          (MH-allele model fits block-level TBV by construction)\n")
cat("  R²_SNP < 1.000          (SNP-additive missing variance — signal spread)\n")
cat("  gap grows with h        (biallelic h=2: gap≈0; multi-allelic h≥5: gap large)\n")

saveRDS(validation, file.path(ITER_DIR, "validation_per_block.rds"))

# Summary per iter
pheno_s1 <- all_pheno$S1[[ITER]]
n_chr <- length(unique(snp_meta$chr))
summary_lines <- c(
  "=== Simulation Summary ===",
  sprintf("Date         : %s",   Sys.time()),
  sprintf("ITER         : %d",   ITER),
  sprintf("N_QTL        : %d",   N_QTL),
  sprintf("H2_target    : %.2f", H2),
  sprintf("Map          : %.1f cM/Mb | length %.2f M", CM_PER_MB, map_length_M),
  sprintf("G1 per family: %d reference + %d candidate", N_FAM_PHENO, N_FAM_CAND),
  sprintf("N_founders   : %d",   founderPop@nInd),
  sprintf("N_chr        : %d",   n_chr),
  sprintf("N_SNP_total  : %d",   n_snp),
  sprintf("N_MH_blocks  : %d",   nrow(mh_coords)),
  sprintf("N_HAP_blocks : %d",   nrow(hap_coords)),
  sprintf("n_G1         : %d",   n_G1),
  "",
  "=== h2_realized per scenario ===",
  sprintf("S1: %.4f", unique(all_pheno$S1[[ITER]]$h2_realized)),
  sprintf("S2: %.4f", unique(all_pheno$S2[[ITER]]$h2_realized)),
  sprintf("S3: %.4f", unique(all_pheno$S3[[ITER]]$h2_realized)),
  sprintf("S4: %.4f", unique(all_pheno$S4[[ITER]]$h2_realized)),
  sprintf("S5: %.4f", unique(all_pheno$S5[[ITER]]$h2_realized)),
  sprintf("S6: %.4f", unique(all_pheno$S6[[ITER]]$h2_realized))
)
writeLines(summary_lines, file.path(ITER_DIR, "simulation_summary.txt"))

write_manifest(c(VCF_FILE, MH_STATS, HAP_STATS), file.path(ITER_DIR, "manifest_inputs.tsv"), "input")
write_manifest(ITER_DIR, file.path(ITER_DIR, "manifest_outputs.tsv"), "output")
log_run_context(list(ITER = ITER, N_FAM_PHENO = N_FAM_PHENO, N_FAM_CAND = N_FAM_CAND,
                     CM_PER_MB = CM_PER_MB, N_QTL = N_QTL, H2 = H2, VCF_FILE = VCF_FILE,
                     map_length_M = map_length_M), ITER_DIR)

cat("Output written to", OUT, "\n")
