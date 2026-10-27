# Single source of the BayesR prior and W scaling, shared by SNP and MH fits.
# Both marker types use the same masbayes kernel ("multiallelic"): a SNP is a
# two-allele block with one non-baseline column, so priors are identical.

# set: "erbe" (classes 1e-4..1e-2, df 10/10), "v1" (classes 1e-3..1e-1, df 10/10),
# "v2" (classes 1e-3..1e-1, df -1/2). Each set is applied to SNP and MH alike.
# In the multiallelic kernel masbayes sets b0_e = sigma2_e*(a0_e-1) and
# b0_g = sigma2_g*(a0_g-2)/a0_g/(1-pi0), so v2 gives b0_e < 0 and b0_g = 0.
bayesr_prior <- function(response = c("gaussian", "binary"), set = c("erbe", "v1", "v2")) {
  response <- match.arg(response); set <- match.arg(set)
  list(pi = if (response == "gaussian") c(.90, .05, .03, .02) else c(.995, .003, .001, .001),
       variance_class = if (set == "erbe") c(0, 1e-4, 1e-3, 1e-2) else c(0, 1e-3, 1e-2, 1e-1),  # x sigma2_g
       a0_e = if (set == "v2") -1 else 10, a0_g = if (set == "v2") 2 else 10,
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
