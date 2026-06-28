# Distance-correlation dependency backbone (v0.2.0, phase G1). Replaces linear
# stats::cor as the name-blind measure of inter-column dependence.
#
# Distance correlation (Szekely, Rizzo & Bakirov 2007) is 0 iff the two columns
# are independent (finite second moments) and detects nonlinear, non-monotone
# dependence that Pearson r misses -- the U-shaped, quadratic, and saturating
# relations that real outcomes/pairs exhibit. It reads only the joint value
# distribution, never names ("Data inspice, non nomen").
#
# INVARIANCE: dCor is a function of pairwise VALUE distances only, so it is
# invariant to row order, column order, and column renaming -- the RELABEL
# (TURNUSOL) and S_n gates are preserved by construction (no column index or
# name ever enters).
#
# BOUNDED COST: the estimator is O(n^2) in time and memory. For n above a cap it
# is computed on a deterministic, RNG-state-preserving row subsample (rows never
# permute under the column-relabel / column-permutation gates, so subsampling is
# invariant). Semantics are preserved; only cost is capped. n-thresholds in the
# detectors are untouched. No exports.

.DCOR_NMAX <- 160L   # row cap for the O(n^2) estimator (cost-bound). dCor is stable
                     # for the two threshold-ROBUST decisions it now drives (repeated
                     # block adjacency >=0.5; pair dependence gate >=0.30) by n~160;
                     # subsample is fixed-seed and column-invariant (see .dcor_idx).
                     # Mission-sanctioned cost bound (estimator semantics preserved).

# Deterministic row-subsample indices that do NOT disturb the global RNG (so the
# eval harness's permutation RNG is unaffected). Fixed seed -> same data yields
# the same subsample on every call (determinism), and is identical regardless of
# column order/names (the gates permute columns, never rows).
.dcor_idx <- function(n, cap = .DCOR_NMAX, seed = 42L) {
  if (n <= cap) return(seq_len(n))
  has_old <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old <- if (has_old) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (has_old) assign(".Random.seed", old, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  set.seed(seed)
  sort(sample.int(n, cap))
}

# Value-based numeric coercion: factor/character -> value-sorted integer codes
# (match on sorted-unique levels, name-free), numeric kept as-is.
.dcor_num <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    lv <- sort(unique(x[!is.na(x)]), method = "radix")
    return(as.numeric(match(x, lv)))
  }
  suppressWarnings(as.numeric(x))
}

# Double-centre a (symmetric) distance matrix: D_ij - Dbar_i. - Dbar_.j + Dbar_..
.dcor_center <- function(D) {
  n <- nrow(D)
  D - rowMeans(D) - rep(colMeans(D), each = n) + mean(D)
}

# Distance correlation in [0, 1]. Returns 0 for degenerate input (a constant
# column, <4 complete pairs). Symmetric: .dcor(x, y) == .dcor(y, x). `cap` bounds
# the O(n^2) cost (default .DCOR_NMAX); ANM passes a larger cap for a cleaner
# residual-independence signal.
.dcor <- function(x, y, cap = .DCOR_NMAX) {
  x <- .dcor_num(x); y <- .dcor_num(y)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L) return(0)
  x <- x[ok]; y <- y[ok]
  idx <- .dcor_idx(length(x), cap); x <- x[idx]; y <- y[idx]
  n <- length(x)
  if (stats::sd(x) < 1e-12 || stats::sd(y) < 1e-12) return(0)   # constant -> independent
  A <- .dcor_center(as.matrix(stats::dist(x)))
  B <- .dcor_center(as.matrix(stats::dist(y)))
  dvarx <- sum(A * A); dvary <- sum(B * B)
  if (dvarx < 1e-14 || dvary < 1e-14) return(0)
  v <- sum(A * B) / sqrt(dvarx * dvary)
  sqrt(max(0, min(1, v)))
}

# Bias-corrected (U-centered) distance-correlation independence T-statistic
# (Szekely & Rizzo 2013/2014). The biased .dcor has a positive floor (~0.27 on a
# binary-vs-continuous pair) that would FLOOD a permissive survival-association test;
# the U-centered statistic has mean ~0 under independence, so under H0 (independence)
# T ~ t_{v-1} with v = M(M-3)/2, and a one-sided T > 2 is the 2-sigma reject -- the SAME
# universal as the ANM margin and .gauss_margin (NO new tuned constant). It is n-INVARIANT
# as a decision (a t-threshold), so a genuine NON-MONOTONE dependence is detected without
# non-survival variance-splits flooding as n grows. Row-capped + RNG-preserving like
# .dcor; value-only => RELABEL/S_n invariant. Returns the T-statistic (0 for degenerate).
.ucenter <- function(D) {
  M <- nrow(D); rs <- rowSums(D); tot <- sum(D)
  A <- D - rs / (M - 2) - rep(rs, each = M) / (M - 2) + tot / ((M - 1) * (M - 2))
  diag(A) <- 0
  A
}
.dcor_t <- function(x, y, cap = .DCOR_NMAX) {
  x <- .dcor_num(x); y <- .dcor_num(y)
  ok <- is.finite(x) & is.finite(y); if (sum(ok) < 8L) return(0)
  x <- x[ok]; y <- y[ok]; idx <- .dcor_idx(length(x), cap); x <- x[idx]; y <- y[idx]
  M <- length(x)
  if (M < 8L || stats::sd(x) < 1e-12 || stats::sd(y) < 1e-12) return(0)
  A <- .ucenter(as.matrix(stats::dist(x))); B <- .ucenter(as.matrix(stats::dist(y)))
  XX <- sum(A * A); YY <- sum(B * B); if (XX <= 0 || YY <= 0) return(0)
  R <- sum(A * B) / sqrt(XX * YY)                  # bias-corrected dCor (can be < 0)
  v <- M * (M - 3) / 2; if (v <= 1 || abs(R) >= 1) return(0)
  sqrt(v - 1) * R / sqrt(1 - R^2)                  # T-statistic; reject independence at T > 2 (2-sigma)
}

# Kraskov-Stogbauer-Grassberger (2004) mutual information, k-nearest-neighbour
# estimator (algorithm 1), at a FIXED k=3 -- an information-theoretic convention
# (universal, NOT corpus-tuned), which removes the discretization/binning knob for
# continuous MI. Max-norm joint metric on standardized values (z-scoring is linear,
# hence MI-invariant, and makes the two coordinates scale-comparable). O(n^2),
# row-capped like .dcor; invariant (function of value vectors; rows never permute
# under the column gates). Returns MI in nats (>= 0; estimator noise clamped at 0).
.KSG_K <- 3L
.ksg_mi <- function(x, y, k = .KSG_K, cap = .DCOR_NMAX) {
  x <- .dcor_num(x); y <- .dcor_num(y)
  ok <- is.finite(x) & is.finite(y); if (sum(ok) < k + 2L) return(0)
  x <- x[ok]; y <- y[ok]
  idx <- .dcor_idx(length(x), cap); x <- x[idx]; y <- y[idx]; n <- length(x)
  if (n < k + 2L || stats::sd(x) < 1e-12 || stats::sd(y) < 1e-12) return(0)
  x <- (x - mean(x)) / stats::sd(x); y <- (y - mean(y)) / stats::sd(y)
  Dx <- abs(outer(x, x, "-")); Dy <- abs(outer(y, y, "-"))
  Dj <- pmax(Dx, Dy)                                   # joint max-norm distance
  acc <- 0
  for (i in seq_len(n)) {
    eps <- sort.int(Dj[i, ], partial = k + 1L)[k + 1L] # k-th NN distance (self at 0 is index 1)
    nx <- sum(Dx[i, ] < eps) - 1L                      # marginal neighbours strictly within eps (excl self)
    ny <- sum(Dy[i, ] < eps) - 1L
    acc <- acc + digamma(nx + 1) + digamma(ny + 1)
  }
  max(0, digamma(k) + digamma(n) - acc / n)            # KSG algorithm 1
}

# Low-cardinality / categorical test for the dependency dispatch.
.is_discrete <- function(v) {
  if (is.character(v) || is.factor(v)) return(TRUE)
  u <- unique(v[!is.na(v)])
  length(u) <= 10L && all(abs(u - round(u)) < 1e-9)
}

# Unified name-blind dependency in [0, 1] -- the completed signature's dependency
# component. Discrete/categorical pair -> plug-in normalized mutual information (no
# binning). Continuous/mixed pair -> KSG MI mapped to its CORRELATION-EQUIVALENT
# sqrt(1 - exp(-2*MI)) (which equals |rho| exactly for a bivariate Gaussian, so .dep
# is a drop-in for linear cor on Gaussian data -- no re-tuning of existing cutoffs),
# supplemented by distance correlation for non-monotone dependence. Symmetric;
# invariant; cost-bounded.
.dep <- function(x, y) {
  if (.is_discrete(x) && .is_discrete(y)) return(compute_nmi(x, y))
  mi <- .ksg_mi(x, y)
  ksg_corr <- sqrt(1 - exp(-2 * mi))
  max(ksg_corr, .dcor(x, y))
}

# Gaussian-reduction admissibility guard (D1). The normalized MI used as a
# dependency magnitude is nmi(X, Y) = sqrt(1 - exp(-2 I(X;Y))); for a BIVARIATE
# GAUSSIAN, I = -0.5 ln(1 - rho^2) so nmi = |rho| EXACTLY, and any second-moment
# (variance/covariance) criterion that equals n*I under Gaussianity -- e.g. the
# measurement-pair bits (n/2)log2(1/(1-rho^2)) -- is then the closed-form NMI.
# `.gauss_margin` is the cheap, n-adaptive, DERIVED normality screen on one column:
# under the Gaussian null, sample skewness ~ N(0, 6/n) and excess kurtosis
# ~ N(0, 24/n), so a 2-sigma two-sided band -- |skew| < 2*sqrt(6/n) and
# |exkurt| < 2*sqrt(24/n) -- is the admissibility certificate. The "2" is the SAME
# 2-sigma statistical universal as the ANM decisiveness margin (no new tuned
# constant); 6/n and 24/n are the exact asymptotic Gaussian sampling variances.
# Value-only (skew/kurt are functions of the column VALUES) => permutation/relabel
# invariant. A constant, short (n<20, the inviolable estimator floor), or
# non-numeric column is NOT admissible -> the caller uses the KSG/dCor estimator.
.gauss_margin <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]; n <- length(x)
  if (n < 20L) return(FALSE)                      # below the inviolable n<20 estimator floor
  s <- stats::sd(x); if (!is.finite(s) || s < 1e-12) return(FALSE)
  z <- (x - mean(x)) / s
  skew <- mean(z^3); exkurt <- mean(z^4) - 3
  is.finite(skew) && is.finite(exkurt) &&
    abs(skew) < 2 * sqrt(6 / n) && abs(exkurt) < 2 * sqrt(24 / n)
}
