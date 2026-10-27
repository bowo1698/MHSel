# Single source of the BayesR prior and W scaling, shared by SNP and MH fits.
# Both marker types use the same masbayes kernel ("multiallelic"): a SNP is a
# two-allele block with one non-baseline column, so priors are identical.

# set "erbe": one prior and kernel for both markers (classes 1e-4..1e-2, df 10/10,
#   sigma2 start from GBLUP, SNP W centred + globally scaled).
# set "v1": shared pi and classes 0/1e-3/1e-2/1e-1 x sigma2_g; each marker keeps its
#   masbayes kernel defaults: SNP = "snp" kernel, z-score W, df -1/2, sigma2_ah = var(y)/4,
#   sigma2_e = var(y)/2; MH = "auto" (multiallelic) kernel, W_ah, df 10/10, sigma2 from GBLUP.
bayesr_prior <- function(response = c("gaussian", "binary"), set = c("erbe", "v1"), marker = c("mh", "snp")) {
  response <- match.arg(response); set <- match.arg(set); marker <- match.arg(marker)
  p <- list(pi = if (response == "gaussian") c(.90, .05, .03, .02) else c(.995, .003, .001, .001),
            variance_class = c(0, 1e-4, 1e-3, 1e-2),  # x sigma2_g
            a0_e = 10, a0_g = 10, marker_type = "multiallelic",
            sigma2_init = "gblup", snp_encoding = "global")
  if (set == "v1") {
    p$variance_class <- c(0, 1e-3, 1e-2, 1e-1)
    p$marker_type <- if (marker == "snp") "snp" else "auto"
    if (marker == "snp") p[c("a0_e", "a0_g", "sigma2_init", "snp_encoding")] <- list(-1, 2, "vary", "zscore")
  }
  p
}

# Centre SNP dosages on 2p; test sets reuse the training freq.
center_snp <- function(X, freq = NULL) {
  if (is.null(freq)) freq <- colMeans(X) / 2
  list(W = sweep(X, 2, 2 * freq), freq = freq)
}

# Global scale to mean column variance 1 (as for W_ah); test sets reuse it.
scale_marker_w <- function(W, scale = NULL) {
  if (is.null(scale)) scale <- sqrt(mean(apply(W, 2, var)))
  list(W = W / scale, scale = scale)
}
