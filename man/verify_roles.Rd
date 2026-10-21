% Hand-written to mirror the roxygen block in R/verify_roles.R
\name{verify_roles}
\alias{verify_roles}
\alias{print.role_verification}
\alias{summary.role_verification}
\title{Verify a declared analysis plan against the data signature}
\usage{
verify_roles(data, claim)

\method{print}{role_verification}(x, ...)

\method{summary}{role_verification}(object, ...)
}
\arguments{
\item{data}{A \code{data.frame} with unique, non-empty column names.}

\item{claim}{A named list. Names are role keys: scalar roles
(\code{subject_id}, \code{group_var}, \code{outcome_continuous},
\code{outcome_binary}, \code{time_variable}, \code{event_variable},
\code{time_index}) take one column name; set roles
(\code{repeated_measures}, \code{scale_items}, \code{covariate}) take one or
more; pair roles (\code{coupled_pair}, \code{paired_pairs},
\code{agreement_pairs}) take a length-2 character vector or a list of them.}

\item{x, object}{A \code{role_verification} object.}

\item{...}{Ignored.}
}
\value{
An object of class \code{"role_verification"}: a list with
\code{table} (role, columns, verdict, sub_reason, statistic, reason, blind),
\code{unclaimed} (role, columns, statistic), \code{counts}
(claims, testable, contradicted, consistent, unverifiable), \code{n_obs}
and \code{n_var}.
}
\description{
Checks a user-declared assignment of variable roles against what the data
themselves can support. The claim is mapped to column positions once and the
data are then examined name-blind (the same information-theoretic signature
used by \code{\link[=detect_roles]{detect_roles()}}). Each claimed role gets one verdict:
\code{contradicted} (a falsification backed by a value statistic),
\code{consistent} (the data support it) or \code{unverifiable}
(\code{sub_reason} \code{identifiability} or \code{coverage}). Structure the
blind detector finds that the claim does not mention is listed separately.
}
\details{
\code{contradicted} is issued only for these falsifications, each backed by a value
statistic: \code{subject_id} whose values repeat while a detected long panel has a
different entity cluster; a \code{paired_pairs}/\code{agreement_pairs} pair that
shares no information beyond the independence noise floor, whose members differ
more than they agree, or with a constant member; a \code{coupled_pair} whose members
are independent; an event column that is not a two-level or 0/1/2 coding, or
negative follow-up times; \code{repeated_measures} that are not the measures of a
detected long panel, or a wide set with no common factor; ordinal
\code{scale_items} with no common factor, or binary items forming a one-hot
partition; a \code{time_index} that is not a calendar column; a \code{group_var}
that is all-distinct or calendar; a non-numeric \code{outcome_continuous}; an
\code{outcome_binary} without exactly two observed values. Absence of a time-event
association is not a falsification (constant-hazard censoring makes the event
independent of the observed time) and is reported as
\code{unverifiable(coverage)}. A \code{consistent} pair means the two columns share
information and their difference is smaller than their mean. Each claimed pair is
judged on its own. A survival claim (\code{time_variable} plus
\code{event_variable}) is one claim shown as two rows and is counted once in
\code{counts}.

}
\examples{
set.seed(1)
rct <- data.frame(id = 1:60, arm = rep(c("a", "b"), 30),
                  pre = rnorm(60), post = rnorm(60))
verify_roles(rct, list(subject_id = "id", group_var = "arm",
                       paired_pairs = c("pre", "post")))
}
