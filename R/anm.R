# Additive-noise-model (ANM) causal direction -- v0.2.0 phase G3.
#
# For a candidate cause/effect PAIR, the additive-noise model y = f(x) + e is
# identifiable when f is nonlinear OR the noise e is non-Gaussian (Hoyer et al.
# 2009; Shimizu et al. 2006 LiNGAM): in the CAUSAL direction the residual is
# independent of the cause, in the ANTI-causal direction it is not. We fit each
# direction with a flexible cubic model and measure the residual's dependence on
# the putative cause via distance correlation (.dcor). The direction with the more
# INDEPENDENT residual is cause->effect.
#
# The linear-Gaussian case is SYMMETRIC: both residuals are independent of their
# regressor, so the direction is UNIDENTIFIABLE -> honest SILENCE (NA). The aim is
# to reduce silence only where theory permits, never to force a pick.
#
# Invariant: a function of the two value vectors only (no names/positions). Bounded:
# residual independence uses the dCor row cap. n<20 -> silent (literature threshold,
# inviolable: ANM needs enough data to estimate residual independence). No exports.

# Decisiveness floor on the residual-independence asymmetry: set ABOVE the finite-
# sample dCor-difference noise floor at .ANM_CAP rows, so a noise-level asymmetry is
# silent (never a forced pick). This is the ANM analogue of the bit-margin (the lone
# confidence knob, calibrated in G4). RELIABLE for nonlinear additive-noise pairs;
# the non-Gaussian-LINEAR (LiNGAM) signal is below this floor with a dCor-residual
# statistic and is conservatively left SILENT (documented limitation).
.ANM_CAP <- 512L   # cost budget: rows for the residual-independence dCor (row subsample,
                   # data-independent; does not change semantics). NOT a tuned detection knob.

# Residual of effect regressed on a cubic expansion of cause (captures nonlinearity;
# the linear term alone covers the non-Gaussian-noise/LiNGAM case). Reuses mdl.R's
# rank-safe lm.fit residual.
.anm_resid <- function(cause, effect) {
  cc <- scale(cause)[, 1]                          # stabilise the cubic powers
  .lm_resid(effect, cbind(cc, cc^2, cc^3))
}

# Returns list(direction in {"x->y","y->x",NA}, asym, d_forward, d_backward). NA =
# unidentifiable (linear-Gaussian symmetric, or below the decisiveness margin).
# C2: the decisiveness margin is DERIVED, not a tuned 0.08 -- it is the 2-sigma
# finite-sample noise floor of the dCor difference, which under independence decays
# as O(1/sqrt(m)); margin = 2/sqrt(m_eff). The "2" is a 2-sigma statistical universal
# (favour silence), m_eff the effective sample size. Larger n -> sharper (smaller
# margin, more decisive); small n -> larger margin -> silence (Tier-3 honest).
.anm_direction <- function(x, y) {
  x <- suppressWarnings(as.numeric(x)); y <- suppressWarnings(as.numeric(y))
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20L) return(list(direction = NA_character_, asym = 0, reason = "n<20"))
  x <- x[ok]; y <- y[ok]
  if (stats::sd(x) < 1e-12 || stats::sd(y) < 1e-12)
    return(list(direction = NA_character_, asym = 0, reason = "constant"))
  margin <- 2 / sqrt(min(length(x), .ANM_CAP))     # derived 2-sigma dCor-difference noise floor
  d_fwd <- .dcor(.anm_resid(x, y), x, cap = .ANM_CAP)   # x->y: residual of y|x vs cause x
  d_bwd <- .dcor(.anm_resid(y, x), y, cap = .ANM_CAP)   # y->x: residual of x|y vs cause y
  asym <- d_bwd - d_fwd                            # >0 => forward residual more independent => x->y
  if (!is.finite(asym) || abs(asym) < margin)
    return(list(direction = NA_character_, asym = asym, d_forward = d_fwd, d_backward = d_bwd,
                reason = "symmetric/unidentifiable"))
  list(direction = if (asym > 0) "x->y" else "y->x", asym = asym,
       d_forward = d_fwd, d_backward = d_bwd)
}
