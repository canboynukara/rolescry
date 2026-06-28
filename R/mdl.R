# Two-part Minimum Description Length (MDL) code lengths, in bits -- v0.2.0 phase G2.
#
# Role assignment reframed as compression: for a candidate role MODEL M of a column
# (or column-set), total code length L(M) = L(model) + L(data | model). A role is
# asserted iff its model compresses the data better (fewer bits) than the NULL
# (marginal / independent) model; the bit MARGIN is the confidence and the silence
# criterion (calibrated in G4). This REPLACES the v0.1.x ad-hoc cutoffs (TAU,
# CORR_MIN, GAMMA, the pct>=30 accept, the survival assoc floor) with one derived
# criterion. The ONLY retained free quantity is the bit-margin; the parameter cost
# is the DERIVED Rissanen term and MIN_PRED=2 is a STRUCTURAL (definitional)
# constant, not a tuned threshold.
#
# Everything is computed from value-based, permutation/relabel-invariant features
# (residual sums of squares, value-defined predictor SETS, fingerprint tie-breaks),
# so the RELABEL / S_n gates hold. n-thresholds in the detectors are untouched: MDL
# does NOT override the inviolable small-n protection (a small n simply makes the
# parameter cost (k/2)log2(n) easier to pay only when the data truly compress).
# No exports.

.MDL_MAXPRED <- 10L   # structural cap on selected predictors (p < n); matches v0.1.x
.MDL_NMAX    <- 1024L # row cap for the regression cost (fixed-seed, column-invariant
                      # subsample via .dcor_idx). 1024 leaves the whole dev simulation
                      # corpus (n <= 1000) UNCAPPED and bounds only the largest real
                      # datasets. Used consistently for RSS AND the (k/2)log2(n)
                      # parameter cost so the MDL break-even stays coherent.

# Rissanen two-part parameter cost: k free parameters -> (k/2) log2(n) bits. DERIVED.
.mdl_param_bits <- function(k, n) (k / 2) * log2(n)

# Gaussian residual code length (bits) for residual sum of squares `rss` over n obs:
# (n/2) log2(2*pi*e*sigma^2), sigma^2 = rss/n. Full form kept so absolute bit margins
# are meaningful (not only model differences).
.mdl_gauss_bits <- function(rss, n) {
  s2 <- rss / n
  if (!is.finite(s2) || s2 <= 0) s2 <- .Machine$double.eps
  (n / 2) * log2(2 * pi * exp(1) * s2)
}

# Binary (0/1) outcomes are scored with the SAME Gaussian MDL code length on the 0/1
# vector (a linear-probability model). A logistic cross-entropy code length is
# theoretically purer but quasi-separates on real classification targets (glm.fit
# diverges -> garbage code length -> total miss); the LPM MDL cannot diverge, shares
# the continuous code path, and keeps the n-aware break-even. (See _dev2/G2_mdl.json.)

# OLS residuals via lm.fit (no formula overhead); rank-deficient X handled (NA coefs,
# residuals still defined). Returns residual vector.
.lm_resid <- function(y, X) {
  f <- suppressWarnings(stats::lm.fit(cbind(1, X), y))
  r <- f$residuals
  if (any(!is.finite(r))) r[!is.finite(r)] <- 0
  r
}
.lm_rss <- function(y, X) sum(.lm_resid(y, X)^2)

# Forward MDL selection of predictors for a CONTINUOUS y. Greedy: each step adds the
# predictor most correlated with the current residual IFF the augmented model reduces
# total code length (delta bits < 0). Order-invariant: ties on the screen broken by
# value-fingerprint `fp` (a per-predictor character key). Returns the selected SET and
# bit lengths. bits_saved = bits_null - bits_model is the MDL evidence.
.mdl_forward_gauss <- function(y, X, fp, maxpred = .MDL_MAXPRED) {
  n <- length(y); ybar <- mean(y); rss0 <- sum((y - ybar)^2)
  bits_null <- .mdl_param_bits(1, n) + .mdl_gauss_bits(rss0, n)
  sel <- integer(0); resid <- y - ybar; cur_bits <- bits_null
  remaining <- seq_len(ncol(X))
  repeat {
    if (length(sel) >= maxpred || length(remaining) == 0L) break
    cors <- vapply(remaining, function(j) {
      v <- X[, j]; s <- stats::sd(v)
      if (!is.finite(s) || s < 1e-12) return(0)
      a <- abs(suppressWarnings(stats::cor(resid, v))); if (is.na(a)) 0 else a
    }, numeric(1))
    if (max(cors) <= 0) break
    top <- remaining[cors >= max(cors) - 1e-12]
    j <- top[order(fp[top])][1]                              # canonical (fingerprint) tie-break
    cand <- c(sel, j)
    rss <- .lm_rss(y, X[, cand, drop = FALSE])
    new_bits <- .mdl_param_bits(length(cand) + 1, n) + .mdl_gauss_bits(rss, n)
    if (new_bits < cur_bits - 1e-9) {                        # pays for its parameter bits
      sel <- cand; cur_bits <- new_bits
      resid <- .lm_resid(y, X[, sel, drop = FALSE])
      remaining <- setdiff(remaining, j)
    } else break
  }
  list(sel = sel, bits_null = bits_null, bits_model = cur_bits,
       bits_saved = bits_null - cur_bits, k = length(sel))
}

# MDL mutual-independence test for a predictor SET (columns of XP). The causes of a
# genuine collider are mutually independent; a cluster member's "predictors" are its
# collinear mates (mutually dependent). Threshold-free: predictors are mutually
# independent iff NO column is MDL-compressed by the others (each fails to pay its
# parameter bits). This is the identifiability guard that distinguishes an OUTCOME
# sink from a CLUSTER source -- replacing the v0.1.x tuned GAMMA with the same MDL
# break-even. `fp` = per-column value-fingerprints for order-free tie-breaks.
.mdl_set_independent <- function(XP, fp) {
  m <- ncol(XP)
  if (m < 2L) return(TRUE)
  for (j in seq_len(m)) {
    others <- setdiff(seq_len(m), j)
    r <- .mdl_forward_gauss(XP[, j], XP[, others, drop = FALSE], fp[others])
    if (r$k >= 1L && r$bits_saved > 0) return(FALSE)   # column j is explained by another "cause" -> not independent
  }
  TRUE
}

# (A variance-aware 2-group MDL association `.mdl_assoc_bin_cont` was prototyped in C1
# for non-monotone survival but reverted -- it over-admitted on monotone real survival;
# removed at C4 cleanup. Non-monotone survival remains a documented Tier-2 gap.)

# MDL measurement-pair criterion (description-length form of MI). Reparametrize a
# pair (a, b) into consensus m=(a+b)/2 and difference d=a-b: bits saved coding (m, d)
# independently vs (a, b) independently = (n/2) log2( var(a) var(b) / (var(m) var(d)) ),
# which for equal-variance Gaussian equals (n/2) log2(1/(1-rho^2)) = n * MI(a,b) in
# bits. Positive => the columns share information (compresses). This REPLACES the
# pct>=30 accept, the paired/agreement r-windows, and the dCor 0.30 gate with one
# derived criterion (no tuned constant). Returns the bits; the caller also requires
# the STRUCTURAL measurement-pair condition var(d) < var(m) ("agree more than they
# differ", definitional, not tuned).
.mdl_pair_bits <- function(a, b) {
  ok <- is.finite(a) & is.finite(b); a <- a[ok]; b <- b[ok]; n <- length(a)
  if (n < 10L) return(list(bits = -Inf, agree = FALSE))
  va <- stats::var(a); vb <- stats::var(b)
  if (!is.finite(va) || !is.finite(vb) || va < 1e-12 || vb < 1e-12) return(list(bits = -Inf, agree = FALSE))
  vm <- stats::var((a + b) / 2); vd <- stats::var(a - b)
  if (vd < 1e-12) vd <- 1e-12
  list(bits = (n / 2) * (log2(va) + log2(vb) - log2(vm) - log2(vd)), agree = vd < vm)
}
