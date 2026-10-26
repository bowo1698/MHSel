# Helpers for simulationG1.R.

# Uniform recombination: Morgans = (pos - min pos) * cM/Mb / 1e8, one vector per chr.
# Input order is kept (it must match the haplotype columns; VCF is position-sorted).
build_gen_map <- function(chr, pos_bp, cm_per_mb = 1) {
  lapply(sort(unique(chr)), function(c) {
    p <- pos_bp[chr == c]
    (p - min(p)) * cm_per_mb / 1e8
  })
}

# Per family: first n_pheno offspring phenotyped (reference), rest candidates.
population_type_vector <- function(n_fam, n_pheno = 25L, n_cand = 10L) {
  rep(rep(c("reference", "candidate"), c(n_pheno, n_cand)), n_fam)
}
