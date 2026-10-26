# Sensitivity variant "nocausal": remove causal SNPs and MH blocks containing them.

causal_snp_ids <- function(scenario, qtl) {
  blk_snps <- function(eff) unique(unlist(lapply(eff, `[[`, "snp_ids")))
  switch(sub("_.*", "", scenario),
         S1 = , S2 = qtl$snp_ids,
         S3 = , S4 = blk_snps(qtl$s3),
         S5 = , S6 = blk_snps(qtl$s5),
         stop("unknown scenario: ", scenario))
}

drop_causal_snps <- function(geno, causal_ids) {
  geno[, !names(geno) %in% causal_ids, drop = FALSE]
}

# hap_<chr>_<k> is the k-th block of chr in snp_selection order (blkN);
# block -> SNP via (chr, snp_index) as in simulationG1.R.
drop_causal_mh_blocks <- function(hap_geno, snp_selection, snp_meta_idx, causal_ids) {
  sel <- merge(snp_selection[, c("block_id", "chr", "snp_index")],
               snp_meta_idx[, c("chr", "snp_index", "snp_id")], by = c("chr", "snp_index"))
  blk_order <- unique(snp_selection[, c("block_id", "chr")])
  blk_order$col <- paste0("hap_", blk_order$chr, "_", ave(seq_along(blk_order$chr), blk_order$chr, FUN = seq_along))
  hit  <- unique(sel$block_id[sel$snp_id %in% causal_ids])
  drop <- blk_order$col[blk_order$block_id %in% hit]
  # both haplotype copies: "<block>" and readr's deduplicated "<block>_1"
  out  <- hap_geno[, !names(hap_geno) %in% c(drop, paste0(drop, "_1")), drop = FALSE]
  attr(out, "dropped_blocks") <- drop
  out
}

# W block ids "hap_<chr>_<k>" -> coords block_id of the k-th block on chr.
translate_block_ids <- function(block_ids, coords) {
  chr <- as.integer(sub("^hap_(\\d+)_\\d+$", "\\1", block_ids))
  k   <- as.integer(sub("^hap_\\d+_(\\d+)$", "\\1", block_ids))
  key <- paste(as.integer(coords$chr), ave(seq_along(coords$chr), coords$chr, FUN = seq_along))
  out <- coords$block_id[match(paste(chr, k), key)]
  stopifnot(!anyNA(out))
  out
}
