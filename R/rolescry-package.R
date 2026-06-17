#' rolescry: Name-blind Variable-Role Detection
#'
#' Deterministic, name-blind detection of variable roles (group, outcome,
#' survival time/event, paired, agreement, repeated measures, scale items,
#' subject id, covariate) in tabular data, using information-theoretic
#' signatures (Shannon entropy, normalized mutual information) and
#' distributional shape rather than column names. The guiding principle is
#' "Data inspice, non nomen" -- inspect the data, not the name.
#'
#' The single public entry point is [detect_roles()]. Header-aware data
#' loading is available via [read_data()]. No LLMs, no external data
#' transmission; detection is >= 90 percent mathematical signature with an
#' optional, capped (<= 10 percent) name bonus (see the \code{name_bonus}
#' argument).
#'
#' Extracted from the MDStatR biostatistics engine.
#'
#' @keywords internal
"_PACKAGE"
