# Single source of the BayesR prior and W scaling, shared by SNP and MH fits.
# Both marker types use the same masbayes kernel ("multiallelic"): a SNP is a
# two-allele block with one non-baseline column, so priors are identical.

bayesr_prior <- function(response = c("gaussian", "binary")) {
  response <- match.arg(response)
  list(pi = if (response == "gaussian") c(.90, .05, .03, .02) else c(.995, .003, .001, .001),
       variance_class = c(0, 1e-4, 1e-3, 1e-2),  # Erbe et al. 2012, x sigma2_g
       a0_e = 10, a0_g = 10,
       marker_type = "multiallelic")
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
