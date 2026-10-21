# Two-part Minimum Description Length (MDL) code lengths, in bits -- v0.2.0 phase G2.
#
# Role assignment reframed as compression: for a candidate role MODEL M of a column
# (or column-set), total code length L(M) = L(model) + L(data | model). A role is
# asserted iff its model compresses the data better (fewer bits) than the NULL
# (marginal / independent) model; the bit MARGIN is the confidence and the silence
# criterion (calibrated in G4). For pairs, outcomes, group and ordinal scale items this
# REPLACED the v0.1.x ad-hoc cutoffs (TAU, CORR_MIN, GAMMA, the pct>=30 accept) with one
# derived criterion; survival, wide repeated_measures and the binary scale path remain on
# the ported point scorers (measured derived replacements were rejected under the
# no-regression rule). The ONLY retained free quantity is the bit-margin; the parameter cost
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

# Measured-and-rejected alternatives: see the v0.3.0 development report (rolescry_v030/v030_dev_report.md, not distributed with the package).
# (Non-monotone survival remains a documented Tier-2 gap.)
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

# Bits saved by coding d = a - b with a free mean (systematic shift) over mean zero,
# net of the one shift parameter. > 0 => paired (shift), <= 0 => agreement.
# Sentinels: +Inf = exact constant non-zero offset (purest paired); -Inf = n < 2 or exact duplicate.
.mdl_shift_bits <- function(a, b) {
  ok <- is.finite(a) & is.finite(b); d <- (a - b)[ok]; n <- length(d)
  ssc <- sum((d - mean(d))^2)
  if (n < 2L) return(-Inf)
  if (ssc < 1e-12) return(if (sum(d^2) >= 1e-12) Inf else -Inf)  # constant shift = purest paired; exact duplicate = none
  (n / 2) * log2(sum(d^2) / ssc) - .mdl_param_bits(1, n)
}

# Bernoulli marginal MDL evidence that predictor `xv` compresses a two-valued
# outcome `y01` (0/1): bits saved by a one-slope logistic code over the null
# (intercept-only) Bernoulli code, MINUS the slope's parameter bits (1/2 log2 n)
# and the log2(m) model-index bits that NAME which of the m candidate predictors
# was selected. The naming term is what makes the evidence selection-consistent:
# scanning m predictors on noise buys ~log2(m) bits of chance compression, and
# the code must pay to say which one it used (Rissanen two-part code). No tuned
# constants; positive bits = net compression = admissible cause.
.mdl_bern_bits <- function(y01, xv, m) {
  ok <- !is.na(y01) & is.finite(xv)
  y <- y01[ok]; x <- xv[ok]; n <- length(y)
  if (n < 10L || stats::sd(x) < 1e-12) return(0)
  k1 <- sum(y)
  if (k1 == 0L || k1 == n) return(0)
  p0 <- k1 / n
  dev0 <- -2 * (k1 * log(p0) + (n - k1) * log(1 - p0))
  fit <- tryCatch(suppressWarnings(
    stats::glm.fit(cbind(1, x), y, family = stats::binomial())), error = function(e) NULL)
  b1 <- if (is.null(fit) || !is.finite(fit$deviance)) 0 else
    (dev0 - fit$deviance) / (2 * log(2)) - .mdl_param_bits(1, n) - log2(max(m, 1))
  # Second route, same two-part-code economics: a NONMONOTONE dependence (dev-band
  # nonlinear/nongauss regimes) leaves the logit slope near zero but still moves the
  # conditional distribution of x given y. Code x with one Gaussian (2 params) vs
  # two Gaussians split by y (4 params): location AND scale shifts both compress.
  # The split must pay its 2 extra parameters plus the same log2(m) naming bits.
  b2 <- 0
  x0 <- x[y == 0L]; x1 <- x[y == 1L]
  if (length(x0) >= 5L && length(x1) >= 5L &&
      stats::sd(x0) > 1e-12 && stats::sd(x1) > 1e-12) {
    rss  <- sum((x - mean(x))^2)
    rss0 <- sum((x0 - mean(x0))^2); rss1 <- sum((x1 - mean(x1))^2)
    b2 <- (.mdl_param_bits(2, n) + .mdl_gauss_bits(rss, n)) -
          (.mdl_param_bits(4, n) + .mdl_gauss_bits(rss0, length(x0)) +
           .mdl_gauss_bits(rss1, length(x1))) - log2(max(m, 1))
  }
  max(b1, b2)
}

# One-factor Gaussian block code: Sigma1 = lambda lambda' + diag(psi) vs independence
# (the identity after standardization). The loadings come from ITERATED PRINCIPAL-AXIS
# factoring of the correlation matrix R: start at communalities h2 = 1, repeatedly
# replace diag(R) by h2, take the leading eigenpair of that reduced matrix,
# lambda = sqrt(e1) v1, h2 <- lambda^2, until the communalities stop moving. This is a
# one-factor (common-variance) fit, not a first principal component. bits = code saved
# over independence, net of k loading parameters. same_quantity: common-factor variance
# exceeds mean unique variance (columns measure one quantity, not just co-vary).
# Identification: for k = 2 the one-factor model is exact (lambda1 lambda2 = r) and for
# k = 3 it is just-identified, so the code cannot reject unidimensionality there;
# the unidimensionality constraint bites only for k >= 4.
# The fit is on STANDARDIZED columns (center, divide by sd): a per-column linear map, so the decision
# is invariant to column units (the absolute raw-variance guard v0 < 1e-12 below is the one exception, for
# near-constant columns); S0 = I (ln det S0 = 0) and same_quantity reads mean communality > mean uniqueness.
# bits = (n/2)/ln2 * [ln det S0 - ln det S1 - (tr(S1^-1 S) - k)] - (k/2) log2 n, with ln det S0 = 0.
# The Gaussian code length under S1 is (n/2)[ln det S1 + tr(S1^-1 S)] nats; independence has
# tr = k exactly, so this is the exact two-part comparison (the log-det-only form holds only at the ML optimum).
.mdl_factor_bits <- function(M) {
  M <- as.matrix(M); k <- ncol(M)
  fail <- function(n) list(bits = -Inf, same_quantity = FALSE, k = as.integer(k), n = as.integer(n))
  M <- M[stats::complete.cases(M) & apply(is.finite(M), 1, all), , drop = FALSE]
  n <- nrow(M)
  if (n > .MDL_NMAX) { M <- M[.dcor_idx(n, .MDL_NMAX, seed = 7L), , drop = FALSE]; n <- nrow(M) }
  if (n < 10L || k < 2L) return(fail(n))
  v0 <- apply(M, 2, stats::var)
  if (any(!is.finite(v0)) || any(v0 < 1e-12)) return(fail(n))
  S <- stats::cov(scale(M)); v <- diag(S)   # standardized columns: S is the correlation matrix, S0 = I
  h2 <- rep(1, k)
  for (it in 1:100) {   # 1e-8 / 100 are numeric-convergence cost caps, not decision thresholds
    Sr <- S; diag(Sr) <- h2
    e <- eigen(Sr, symmetric = TRUE)
    lambda <- sqrt(max(e$values[1], 0)) * e$vectors[, 1]
    h2n <- pmin(lambda^2, 1 - 1e-12)
    done <- max(abs(h2n - h2)) < 1e-8
    h2 <- h2n
    if (done) break
  }
  psi <- pmax(v - lambda^2, 1e-12 * v)
  Sigma1 <- tcrossprod(lambda) + diag(psi, k)
  ld <- tryCatch(determinant(Sigma1, logarithm = TRUE),
                 error = function(e) NULL)
  if (is.null(ld) || ld$sign <= 0 || !is.finite(ld$modulus)) return(fail(n))
  tr <- tryCatch(sum(diag(solve(Sigma1, S))), error = function(e) NA_real_)
  if (!is.finite(tr)) return(fail(n))
  list(bits = (n / 2) / log(2) * (sum(log(v)) - as.numeric(ld$modulus) - (tr - k)) - .mdl_param_bits(k, n),
       same_quantity = mean(lambda^2) > mean(psi), k = as.integer(k), n = as.integer(n))
}
