# Seeded data-generating processes for the test fixtures. Everything is
# generated in-process (no real MDStatR CSVs, no patient data); each frame is a
# few hundred rows (well under 5 MB).

# Turnusol twin pair: byte-identical numeric data under col_N names and under
# meaningful names. Mirrors the DS40/DS41 DGP (col_7 = col_3 - col_5 + noise).
make_turnusol_twins <- function(seed = 2026L) {
  set.seed(seed)
  n <- 200L
  col_1 <- stats::rnorm(n, 50, 10)
  col_2 <- stats::rnorm(n, 0, 1)
  col_3 <- stats::rnorm(n, 10, 3)
  col_4 <- stats::rnorm(n, 5, 2)
  col_5 <- stats::rnorm(n, 8, 2)
  col_6 <- stats::rnorm(n, 100, 15)
  col_7 <- col_3 - col_5 + stats::rnorm(n, 0, 5)
  m <- data.frame(col_1, col_2, col_3, col_4, col_5, col_6, col_7, check.names = FALSE)
  named <- m
  names(named) <- c("age", "weight", "biomarker_X", "site_score", "biomarker_Y", "lab_panel", "outcome")
  list(col_n = m, named = named)
}

# Clinical-shape twins with a real categorical group and TWO binary columns: a
# demographic (positionally first) and the intended outcome (named "death").
# Used to exercise the name_bonus tie-breaker that the all-continuous DS40/DS41
# fixtures could not (Phase A caveat).
make_namebonus_twins <- function(seed = 7L) {
  set.seed(seed)
  n <- 160L
  treat <- rep(c("A", "B"), each = n / 2L)            # categorical group
  male  <- stats::rbinom(n, 1, 0.5)                    # demographic binary (first)
  pre   <- stats::rnorm(n, 10, 2)
  post  <- pre + stats::rnorm(n, 1, 1)                 # paired with pre
  death <- stats::rbinom(n, 1, 0.35)                   # intended binary outcome
  named <- data.frame(treatment_arm = treat, male = male, pre = pre,
                      post = post, death = death, check.names = FALSE)
  col_n <- named
  names(col_n) <- paste0("col_", seq_along(col_n))
  list(col_n = col_n, named = named)
}

# Map a role's detected columns to a sorted integer index set within a frame.
role_index <- function(res, role, frame) {
  cols <- res$roles[[role]]$columns
  if (length(cols) == 0) integer(0) else sort(match(cols, names(frame)))
}
