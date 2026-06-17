#' Default name-bonus keyword dictionary
#'
#' Returns a ready-made, ASCII-English keyword dictionary suitable for the
#' \code{name_bonus} argument of [detect_roles()]. It externalizes the
#' hard-coded keyword lists that lived inside the original MDStatR engine
#' (group/treatment, outcome, survival-time and event, subject-id terms) into a
#' plain, inspectable, locale-neutral list.
#'
#' Passing this turns column names into a small, capped tie-breaker only
#' (\eqn{\le} 10 percent of the selection score); the mathematical signature
#' still dominates (>= 90 percent), satisfying the name-blindness contract.
#' Detection without it (\code{name_bonus = NULL}) is purely mathematical.
#'
#' @return A named list of character vectors (regex fragments), with keys
#'   \code{group_var}, \code{outcome_continuous}, \code{outcome_binary},
#'   \code{subject_id}, \code{time_variable}, \code{event_variable}.
#'
#' @examples
#' nb <- rolescry_default_name_bonus()
#' names(nb)
#' set.seed(1)
#' d <- data.frame(
#'   treatment_arm = rep(c("A", "B"), each = 60),
#'   biomarker     = rnorm(120),
#'   death         = rbinom(120, 1, 0.3)
#' )
#' detect_roles(d, name_bonus = nb)$roles$group_var$columns
#'
#' @export
rolescry_default_name_bonus <- function() {
  list(
    group_var = c(
      "group", "treat", "arm", "cohort", "regimen", "protocol", "drug",
      "therapy", "interv", "exposure"
    ),
    outcome_continuous = c(
      "outcome", "response", "score", "change", "delta", "result", "level"
    ),
    outcome_binary = c(
      "outcome", "response", "endpoint", "death", "died", "mortality", "event",
      "status", "relapse", "recur", "progression", "failure", "remission",
      "complication", "readmission", "composite", "primary"
    ),
    subject_id = c(
      "id", "subject", "patient", "participant", "record", "case"
    ),
    time_variable = c(
      "time", "days", "months", "years", "duration", "futime", "surv",
      "followup", "follow_up", "tstart", "tstop"
    ),
    event_variable = c(
      "event", "death", "died", "status", "censor", "fail", "delta",
      "relapse", "recur", "outcome"
    )
  )
}
