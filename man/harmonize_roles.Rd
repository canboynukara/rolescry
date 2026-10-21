% Hand-written to mirror the roxygen block in R/harmonize_roles.R
\name{harmonize_roles}
\alias{harmonize_roles}
\alias{print.role_harmonization}
\alias{summary.role_harmonization}
\title{Align detected variable roles across sites}
\usage{
harmonize_roles(sites)

\method{print}{role_harmonization}(x, ...)

\method{summary}{role_harmonization}(object, ...)
}
\arguments{
\item{sites}{A named list (unique, non-empty names) of at least two data.frames.}

\item{x, object}{A \code{role_harmonization} object.}

\item{...}{Ignored.}
}
\value{
An object of class \code{"role_harmonization"}: a list with \code{map}
(role, status, one column per site), \code{checks} (role, check, one column per
site, flag), \code{site_only} (site, role, columns), \code{n_sites} and
\code{sites}. Flags are \code{type_mismatch}, \code{coding_mismatch} and
\code{size_mismatch}; no numeric threshold decides a flag. \code{summary}
returns \code{map}.
}
\description{
Runs \code{\link[=detect_roles]{detect_roles()}} name-blind on each site and reports, per role, which
column each site holds, plus descriptive consistency checks (value type, event
coding, group-label overlap, range overlap, battery size). Alignment is by role
only; column names appear only in the returned cells. The residual
\code{covariate} role is excluded from \code{map}, \code{checks} and
\code{site_only}. \code{aligned} means the detector found the role at every site; it does not establish that the columns measure the same thing; read \code{checks} (levels, Jaccard, ranges) before pooling. Cell order is alphabetical and carries no pairing direction.
}
\examples{
set.seed(1)
mk <- function(n) data.frame(id = seq_len(n), arm = sample(c("a", "b"), n, TRUE),
                             x = rnorm(n), y = rnorm(n))
harmonize_roles(list(north = mk(60), south = mk(80)))
}
