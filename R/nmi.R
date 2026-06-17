#' Normalized mutual information
#'
#' Computes the normalized mutual information (NMI) between two discrete
#' variables, a name-blind, information-theoretic measure of association in
#' \eqn{[0, 1]}. NMI is the mutual information divided by the smaller of the two
#' marginal Shannon entropies; it is 0 for independent variables and 1 for a
#' perfect (deterministic) association, and unlike a raw chi-squared it is
#' comparable across variables with different numbers of levels.
#'
#' @param x Either a two-way contingency table / matrix of counts, or a vector
#'   (factor, character, or numeric) of the first variable.
#' @param y Optional. If \code{x} is a vector, the second variable's vector;
#'   a contingency table is formed via \code{table(x, y)} on complete cases.
#'   Ignored when \code{x} is already a table/matrix.
#'
#' @return A single numeric in \eqn{[0, 1]}. Returns 0 for degenerate input
#'   (fewer than two rows/columns, zero total, or near-zero marginal entropy).
#'
#' @examples
#' set.seed(1)
#' g <- sample(c("A", "B", "C"), 200, replace = TRUE)
#' y <- ifelse(g == "A", "yes", sample(c("yes", "no"), 200, replace = TRUE))
#' compute_nmi(g, y)            # > 0: g carries information about y
#' compute_nmi(g, sample(g))    # ~0: shuffled -> independent
#'
#' @export
compute_nmi <- function(x, y = NULL) {
  ct <- if (is.table(x) || is.matrix(x)) {
    x
  } else {
    if (is.null(y)) stop("compute_nmi(): supply a contingency table, or both x and y.")
    cc <- stats::complete.cases(x, y)
    table(x[cc], y[cc])
  }
  if (!is.matrix(ct) && !is.table(ct)) {
    return(0)
  }
  if (nrow(ct) < 2 || ncol(ct) < 2) {
    return(0)
  }
  total <- sum(ct)
  if (total == 0) {
    return(0)
  }
  pxy <- ct / total
  px <- rowSums(pxy)
  py <- colSums(pxy)
  mi <- 0
  for (i in seq_len(nrow(pxy))) {
    for (j in seq_len(ncol(pxy))) {
      if (pxy[i, j] > 0 && px[i] > 0 && py[j] > 0) {
        mi <- mi + pxy[i, j] * log2(pxy[i, j] / (px[i] * py[j]))
      }
    }
  }
  hx <- -sum(px[px > 0] * log2(px[px > 0]))
  hy <- -sum(py[py > 0] * log2(py[py > 0]))
  denom <- min(hx, hy)
  if (denom < 1e-10) {
    return(0)
  }
  max(0, min(1, mi / denom))
}
