.ROLE_ARITY <- c(
  subject_id = "scalar", group_var = "scalar", outcome_continuous = "scalar",
  outcome_binary = "scalar", time_variable = "scalar", event_variable = "scalar",
  time_index = "scalar", coupled_pair = "pair", paired_pairs = "pair",
  agreement_pairs = "pair", repeated_measures = "set", scale_items = "set",
  covariate = "set"
)

.vr_stop <- function(...) stop("verify_roles(): ", ..., call. = FALSE)

# Validate the claim against `data`; returns nothing, stops on the first problem.
.validate_claim <- function(data, claim) {
  if (!is.data.frame(data)) .vr_stop("`data` must be a data.frame.")
  nm <- names(data)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) .vr_stop("`data` has empty column names.")
  if (anyDuplicated(nm)) {
    .vr_stop("`data` has duplicate column names (", paste(unique(nm[duplicated(nm)]), collapse = ", "),
             "); the claim needs unambiguous names.")
  }
  if (!is.list(claim) || (length(claim) > 0 && (is.null(names(claim)) || any(!nzchar(names(claim)))))) {
    .vr_stop("`claim` must be a named list.")
  }
  if (anyDuplicated(names(claim))) {
    .vr_stop("duplicate role key: ", paste(unique(names(claim)[duplicated(names(claim))]), collapse = ", "), ".")
  }
  bad <- setdiff(names(claim), names(.ROLE_ARITY))
  if (length(bad)) {
    .vr_stop("unknown role(s) ", paste(bad, collapse = ", "), "; valid roles are ",
             paste(names(.ROLE_ARITY), collapse = ", "), ".")
  }
  for (role in names(claim)) {
    els <- if (.ROLE_ARITY[[role]] == "pair" && is.list(claim[[role]])) claim[[role]] else list(claim[[role]])
    for (el in els) {
      if (!is.character(el)) .vr_stop("claim for `", role, "` must be character.")
      if (anyNA(el)) .vr_stop("claim contains NA (role `", role, "`).")
      if (.ROLE_ARITY[[role]] == "scalar" && length(el) != 1L) {
        .vr_stop("`", role, "` needs a single column name (length 1), got length ", length(el), ".")
      }
      if (.ROLE_ARITY[[role]] == "pair" && length(el) != 2L) {
        .vr_stop("each `", role, "` element needs two column names (length 2), got length ", length(el), ".")
      }
      if (.ROLE_ARITY[[role]] == "set" && length(el) == 0L) {
        .vr_stop("`", role, "` needs at least one column name (length 0 given).")
      }
      if (.ROLE_ARITY[[role]] == "set" && anyDuplicated(el)) {
        .vr_stop("duplicate columns in `", role, "`.")
      }
      if (.ROLE_ARITY[[role]] == "pair" && el[1] == el[2]) {
        .vr_stop("the two columns of a `", role, "` pair are identical (", el[1], ").")
      }
      unk <- setdiff(el, nm)
      if (length(unk)) {
        .vr_stop("unknown column(s) ", paste(unk, collapse = ", "), "; available columns are ",
                 paste(nm, collapse = ", "), ".")
      }
    }
  }
  invisible(TRUE)
}

# One row per scalar/set role, one per pair for pair roles.
.claim_rows <- function(claim) {
  roles <- character(0); cols <- list()
  for (role in names(claim)) {
    x <- claim[[role]]
    els <- if (.ROLE_ARITY[[role]] == "pair" && is.list(x)) x else list(x)
    for (el in els) { roles <- c(roles, role); cols[[length(cols) + 1L]] <- as.character(el) }
  }
  out <- data.frame(role = roles, stringsAsFactors = FALSE)
  out$cols <- I(cols)
  out
}

.verdict_row <- function(role, cols, verdict, sub_reason = NA_character_,
                         statistic = NA_character_, reason = NA_character_, blind = "silent") {
  data.frame(role = role, columns = paste(cols, collapse = ", "), verdict = verdict,
             sub_reason = sub_reason, statistic = statistic, reason = reason, blind = blind,
             stringsAsFactors = FALSE)
}

.vr <- function(verdict, statistic, reason, sub_reason = NA_character_) {
  if (verdict == "unverifiable" && is.na(sub_reason)) sub_reason <- "coverage"
  list(verdict = verdict, sub_reason = sub_reason, statistic = statistic, reason = reason)
}
# Shared empty-evidence guard: a coverage row when a column is all NA or < 10 complete cases.
.vr_thin <- function(...) {
  m <- data.frame(..., check.names = FALSE)
  if (any(colSums(!is.na(m)) == 0L)) return(.vr("unverifiable", NA_character_, "No observed values in a claimed column."))
  if (sum(stats::complete.cases(m)) < 10L)
    return(.vr("unverifiable", NA_character_, "Fewer than 10 complete cases."))
  NULL
}
.vr_num <- function(x) suppressWarnings(as.numeric(x))

.vr_subject_id <- function(pos, data, vi, panel) {
  x <- data[[pos]]; nn <- x[!is.na(x)]
  uniq <- length(unique(nn)) == length(nn)
  st <- sprintf("unique=%s", uniq)
  single <- !is.null(panel) && !isTRUE(panel$ambiguous)
  if (single && identical(panel$cluster, names(data)[pos]))
    return(.vr("consistent", st, "The column is the entity cluster of a detected long panel."))
  if (!uniq) {
    if (!is.null(panel) && isTRUE(panel$ambiguous))
      return(.vr("unverifiable", st, "Balanced panel: entity/time axis ambiguous."))
    if (single)   # a panel was detected and this column is not its cluster: a falsification
      return(.vr("contradicted", st, "Values repeat (unique=FALSE) and the column is not the entity cluster of the detected long panel, so it cannot identify subjects."))
    return(.vr("unverifiable", st, "Values repeat; not a row identifier; may be the cluster of an undetected long panel."))
  }
  if (vi$user_type[pos] == "ID" || .id_like(x))
    return(.vr("consistent", st, "Every observed value is distinct and the column has an identifier signature."))
  .vr("unverifiable", st, "Values are distinct but the column has no identifier signature.")
}

.vr_pair <- function(role, pos, data) {
  a <- data[[pos[1]]]; b <- data[[pos[2]]]
  if (!is.numeric(a) || !is.numeric(b))
    return(.vr("unverifiable", "type=non-numeric", "Pair rules apply to numeric columns only."))
  thin <- .vr_thin(a, b); if (!is.null(thin)) return(thin)
  if (stats::var(a, na.rm = TRUE) == 0 || stats::var(b, na.rm = TRUE) == 0)
    return(.vr("contradicted", "unique=1", "A column is constant (variance 0), so it cannot be a member of a measurement pair."))
  p <- .pair_admitted(a, b)
  st <- sprintf("bits=%.1f", p$bits)
  if (!p$admitted) {
    if (isTRUE(p$bits > 0) && !isTRUE(p$agree))
      return(.vr("contradicted", st, sprintf("The columns differ more than they agree (var(a-b) >= var((a+b)/2); bits=%.1f).", p$bits)))
    return(.vr("contradicted", st, sprintf("The two columns share no information beyond the independence noise floor (bits=%.1f).", p$bits)))
  }
  if (p$route == "dcor-fallback") {   # T = dep * sqrt(m): the 2/sqrt(m) floor of .pair_admitted, in sigma units
    tt <- .dep(a, b) * sqrt(sum(is.finite(a) & is.finite(b)))
    st <- sprintf("T=%.2f", tt)
    msg <- sprintf("The two columns are dependent and their difference is smaller than their mean (T=%.2f).", tt)
  } else msg <- sprintf("The two columns share information and their difference is smaller than their mean (bits=%.1f).", p$bits)
  note <- ""
  if (role != "coupled_pair") {
    shift <- .mdl_shift_bits(a, b)
    if ((role == "paired_pairs") != (shift > 0))
      note <- sprintf(" Shift bits suggest %s; the label is n-dependent.",
                      if (role == "paired_pairs") "agreement_pairs" else "paired_pairs")
  }
  .vr("consistent", st, paste0(msg, note))
}

.vr_survival <- function(role, pos, data, claim) {
  tp <- claim[["time_variable"]]; ep <- claim[["event_variable"]]
  if (is.null(tp) || is.null(ep))
    return(.vr("unverifiable", NA_character_, "Survival roles are verified as a pair; claim both."))
  if (!is.numeric(data[[tp]])) return(.vr("unverifiable", "type=non-numeric", "Time is not numeric."))
  thin <- .vr_thin(data[[tp]], data[[ep]]); if (!is.null(thin)) return(thin)
  tv <- .vr_num(data[[tp]]); ev01 <- .normalize_event(data[[ep]])
  if (is.null(ev01)) return(.vr("contradicted", "type=not-event", "The event column is not a two-level or {0,1,2} coding."))
  if (any(tv < 0, na.rm = TRUE)) return(.vr("contradicted", "type=negative", "Time has negative values."))
  nz <- sum(tv == 0, na.rm = TRUE)
  zero <- if (nz > 0) sprintf(" %d zero follow-up times.", nz) else ""
  sa <- .survival_assoc(tv, ev01, min_n = 9L)   # min_n = 9: same >= 10 complete cases as .vr_thin
  if (!is.finite(sa$bits)) return(.vr("unverifiable", NA_character_, paste0("Too few complete cases for the association test.", zero)))
  if (!sa$admitted)
    return(.vr("unverifiable", sprintf("bits=%.1f", sa$bits), paste0(sprintf(
      "No time\u2013event association detected (bits=%.1f, T=%.2f); compatible with survival under constant-hazard censoring.", sa$bits, sa$t), zero)))
  dec <- .censoring_decisive(tv, ev01)
  if (dec) .vr("consistent", "decisive=TRUE", paste0("Event is associated with time and the censoring geometry is that of a survival pair.", zero))
  else .vr("unverifiable", "decisive=FALSE", paste0("Survival vs cause/outcome reading is not decidable from censoring geometry.", zero))
}

.vr_repeated <- function(pos, data, panel) {
  cn <- names(data)[pos]
  if (!is.null(panel) && isTRUE(panel$ambiguous))
    return(.vr("unverifiable", NA_character_, "Balanced panel: entity/time axis ambiguous."))
  if (!is.null(panel)) {
    ok <- all(cn %in% panel$measures)
    return(.vr(if (ok) "consistent" else "contradicted", "type=long-panel",
               if (ok) "All claimed columns are measures of the detected long panel."
               else sprintf("The claimed columns are not the panel's measures (a long panel with %d measures was detected).", length(panel$measures))))
  }
  if (length(pos) < 3L) return(.vr("unverifiable", NA_character_, "Fewer than three columns; a common factor is not testable."))
  if (!all(vapply(pos, function(i) is.numeric(data[[i]]), logical(1L))))
    return(.vr("unverifiable", "type=non-numeric", "The set is not numeric."))
  thin <- .vr_thin(data[, pos, drop = FALSE]); if (!is.null(thin)) return(thin)
  M <- as.matrix(data[, pos, drop = FALSE])
  f <- .mdl_factor_bits(M); st <- sprintf("bits=%.1f", f$bits)
  if (f$bits <= 0) return(.vr("contradicted", st, sprintf("The columns share no common factor (bits=%.1f).", f$bits)))
  if (isTRUE(f$same_quantity)) .vr("consistent", st, "The columns share one factor on a common scale.")
  else .vr("unverifiable", st, "Shared factor but unequal scales.")
}

.vr_scale <- function(pos, data) {
  ty <- vapply(pos, function(i) .scale_item_type(data[[i]]), character(1L))
  if (length(pos) < 3L || anyNA(ty) || length(unique(ty)) != 1L)
    return(.vr("unverifiable", NA_character_, "Scale items need at least three columns of one integer-ordinal or binary type."))
  if (ty[1] == "ordinal") {
    thin <- .vr_thin(data[, pos, drop = FALSE]); if (!is.null(thin)) return(thin)
    f <- .mdl_factor_bits(as.matrix(data[, pos, drop = FALSE])); st <- sprintf("bits=%.1f", f$bits)
    return(if (f$bits > 0) .vr("consistent", st, "The items share one latent factor.")
           else .vr("contradicted", st, sprintf("The items share no common factor (bits=%.1f).", f$bits)))
  }
  m <- as.matrix(data[, pos, drop = FALSE]); m <- m[stats::complete.cases(m), , drop = FALSE]
  if (nrow(m) > 0 && max(rowSums(m)) <= 1)
    return(.vr("contradicted", "one_hot=TRUE", "The binary columns form a one-hot partition of one categorical variable, not scale items."))
  .vr("unverifiable", sprintf("items=%d", length(pos)),
      "Binary batteries have no derived test in this version; the detector's ported scorer is not treated as evidence.")
}

.vr_time_index <- function(pos, vi) {
  t <- vi$user_type[pos]; st <- sprintf("type=%s", t)
  if (t == "Datetime") .vr("consistent", st, "The column holds calendar values.")
  else .vr("contradicted", st, sprintf("The column is not a calendar column (type=%s).", t))
}

.vr_design <- function(role, pos, data, vi) {
  idf <- "identifiability"
  if (role == "covariate") {   # a set: one row, per-column types pasted
    if (any(colSums(!is.na(data[, pos, drop = FALSE])) == 0L))
      return(.vr("unverifiable", NA_character_, "No observed values in a claimed column."))
    return(.vr("unverifiable", paste0("type=", paste(vi$user_type[pos], collapse = ",")),
               "The data cannot confirm a column as a covariate.", idf))
  }
  x <- data[[pos]]; t <- vi$user_type[pos]; nn <- x[!is.na(x)]
  if (length(nn) == 0L) return(.vr("unverifiable", NA_character_, "No observed values in the claimed column."))
  switch(role,
    group_var = {
      if (t == "Datetime" || t == "ID")
        return(.vr("contradicted", sprintf("type=%s", t), sprintf("Every value is distinct or the column is a calendar column (type=%s), so it cannot be a grouping.", t)))
      gb <- if (is.character(x) || is.factor(x)) .group_balance(x) else NULL
      if (is.null(gb)) {
        why <- if (is.numeric(x)) "The grouping is numeric-coded; the blind detector excludes numeric groups by design (Lemma 2)."
               else "A level has fewer than three observations; balance is not computed."
        .vr("unverifiable", "B=NA", why, idf)
      } else .vr("unverifiable", sprintf("B=%.3f", gb$B), "The data cannot say whether a column is the grouping; the entropy balance is reported.", idf)
    },
    outcome_continuous = {
      if (!is.numeric(x)) return(.vr("contradicted", sprintf("type=%s", t), sprintf("The column is not numeric (type=%s).", t)))
      .vr("unverifiable", sprintf("type=%s", t), "The data cannot say which numeric column is the outcome.", idf)
    },
    outcome_binary = {
      u <- length(unique(nn))
      if (u != 2L) return(.vr("contradicted", sprintf("unique=%d", u), sprintf("The column has %d observed values, not two (unique=%d).", u, u)))
      .vr("unverifiable", "unique=2", "The data cannot say which two-valued column is the outcome.", idf)
    })
}

# `pos` are column positions; `data` is the relabeled frame (names col_k); `claim` maps
# role -> first claimed position (used by the survival pair).
.verify_one <- function(role, pos, data, var_info, blind, panel, claim = list()) {
  switch(role,
    subject_id = .vr_subject_id(pos, data, var_info, panel),
    paired_pairs = , agreement_pairs = .vr_pair(role, pos, data),
    coupled_pair = {
      a <- data[[pos[1]]]; b <- data[[pos[2]]]
      ia <- !is.null(.normalize_event(a)) && length(unique(a[!is.na(a)])) <= 3L
      ib <- !is.null(.normalize_event(b)) && length(unique(b[!is.na(b)])) <= 3L
      if (xor(ia, ib)) {
        cont <- if (ia) b else a
        if (!is.numeric(cont)) return(.vr("unverifiable", "type=non-numeric", "The continuous member is not numeric."))
        thin <- .vr_thin(a, b); if (!is.null(thin)) return(thin)
        ev <- .normalize_event(if (ia) a else b); tv <- .vr_num(cont)
        sa <- .survival_assoc(tv, ev, min_n = 9L)   # same >= 10 complete cases as .vr_thin
        if (!is.finite(sa$bits)) return(.vr("unverifiable", NA_character_, "Too few complete cases for the association test."))
        if (sa$admitted) {
          if (!sa$monotone) .vr("consistent", sprintf("T=%.2f", sa$t), sprintf("The binary and continuous members are dependent (T=%.2f).", sa$t))
          else .vr("consistent", sprintf("bits=%.1f", sa$bits), sprintf("The binary and continuous members are dependent (bits=%.1f).", sa$bits))
        } else .vr("contradicted", sprintf("bits=%.1f", sa$bits), sprintf("The binary and continuous members are independent (bits=%.1f, T=%.2f).", sa$bits, sa$t))
      } else .vr("unverifiable", NA_character_, "A coupled pair needs one binary and one continuous column.")
    },
    time_variable = , event_variable = .vr_survival(role, pos, data, claim),
    repeated_measures = .vr_repeated(pos, data, panel),
    scale_items = .vr_scale(pos, data),
    time_index = .vr_time_index(pos, var_info),
    .vr_design(role, pos, data, var_info))
}

# Blind detector's columns for `role` (as original names), or "silent".
.blind_cols <- function(blind, role, nm) {
  r <- blind$roles[[role]]
  if (is.null(r) || !isTRUE(r$found) || !length(r$columns)) return("silent")
  paste(nm[as.integer(sub("col_", "", r$columns))], collapse = ", ")
}

#' Verify a declared analysis plan against the data signature
#'
#' Checks a user-declared assignment of variable roles against what the data
#' themselves can support. The claim is mapped to column positions once and the
#' data are then examined name-blind (the same information-theoretic signature
#' used by [detect_roles()]). Each claimed role gets one verdict:
#' \code{contradicted} (a falsification backed by a value statistic),
#' \code{consistent} (the data support it) or \code{unverifiable}
#' (\code{sub_reason} \code{identifiability} or \code{coverage}). Structure the
#' blind detector finds that the claim does not mention is listed separately.
#'
#' @details
#' \code{contradicted} is issued only for these falsifications, each backed by a value
#' statistic: \code{subject_id} whose values repeat while a detected long panel has a
#' different entity cluster; a \code{paired_pairs}/\code{agreement_pairs} pair that
#' shares no information beyond the independence noise floor, whose members differ
#' more than they agree, or with a constant member; a \code{coupled_pair} whose members
#' are independent; an event column that is not a two-level or 0/1/2 coding, or
#' negative follow-up times; \code{repeated_measures} that are not the measures of a
#' detected long panel, or a wide set with no common factor; ordinal
#' \code{scale_items} with no common factor, or binary items forming a one-hot
#' partition; a \code{time_index} that is not a calendar column; a \code{group_var}
#' that is all-distinct or calendar; a non-numeric \code{outcome_continuous}; an
#' \code{outcome_binary} without exactly two observed values. Absence of a time-event
#' association is not a falsification (constant-hazard censoring makes the event
#' independent of the observed time) and is reported as
#' \code{unverifiable(coverage)}. A \code{consistent} pair means the two columns share
#' information and their difference is smaller than their mean. Each claimed pair is
#' judged on its own. A survival claim (\code{time_variable} plus
#' \code{event_variable}) is one claim shown as two rows and is counted once in
#' \code{counts}.
#'
#' @param data A \code{data.frame} with unique, non-empty column names.
#' @param claim A named list. Names are role keys: scalar roles
#'   (\code{subject_id}, \code{group_var}, \code{outcome_continuous},
#'   \code{outcome_binary}, \code{time_variable}, \code{event_variable},
#'   \code{time_index}) take one column name; set roles
#'   (\code{repeated_measures}, \code{scale_items}, \code{covariate}) take one or
#'   more; pair roles (\code{coupled_pair}, \code{paired_pairs},
#'   \code{agreement_pairs}) take a length-2 character vector or a list of them.
#'
#' @return An object of class \code{"role_verification"}: a list with
#'   \code{table} (role, columns, verdict, sub_reason, statistic, reason, blind),
#'   \code{unclaimed} (role, columns, statistic), \code{counts}
#'   (claims, testable, contradicted, consistent, unverifiable), \code{n_obs}
#'   and \code{n_var}.
#'
#' @examples
#' set.seed(1)
#' rct <- data.frame(id = 1:60, arm = rep(c("a", "b"), 30),
#'                   pre = rnorm(60), post = rnorm(60))
#' verify_roles(rct, list(subject_id = "id", group_var = "arm",
#'                        paired_pairs = c("pre", "post")))
#' @export
verify_roles <- function(data, claim) {
  .validate_claim(data, claim)
  nm <- names(data)
  rows <- .claim_rows(claim)
  pos_of <- function(cols) match(cols, nm)  # the only name -> position mapping

  dblind <- stats::setNames(data, paste0("col_", seq_along(data)))
  blind <- detect_roles(dblind)

  # everything below works on positions of the relabeled frame
  vi <- .build_var_info(dblind)
  panel <- .detect_long_panel(dblind, vi)
  pos_list <- lapply(rows$cols, pos_of)
  first_pos <- function(r) { k <- match(r, rows$role); if (is.na(k)) NULL else pos_list[[k]][1] }
  cl <- list(time_variable = first_pos("time_variable"), event_variable = first_pos("event_variable"))
  res <- lapply(seq_len(nrow(rows)), function(i) {
    role <- rows$role[i]; cols <- rows$cols[[i]]
    v <- .verify_one(role, pos_list[[i]], dblind, vi, blind, panel, cl)
    others <- unique(rows$role[-i][vapply(pos_list[-i], function(p) any(p %in% pos_list[[i]]), logical(1L))])
    others <- setdiff(others, role)
    if (length(others)) v$reason <- paste0(v$reason, " Also claimed as ", paste(others, collapse = ", "), ".")
    .verdict_row(role, cols, v$verdict, v$sub_reason, v$statistic, v$reason,
                 .blind_cols(blind, role, nm))
  })
  table <- if (length(res)) do.call(rbind, res) else
    data.frame(role = character(), columns = character(), verdict = character(),
               sub_reason = character(), statistic = character(), reason = character(),
               blind = character(), stringsAsFactors = FALSE)
  rownames(table) <- NULL

  # blind-found structure (covariate excluded) not claimed under the same role: a pair is
  # compared against every claimed pair, other roles against the claimed column set
  bpos <- function(r) as.integer(sub("col_", "", blind$roles[[r]]$columns))
  ur <- character(0); uc <- character(0); us <- character(0)
  for (r in setdiff(names(blind$roles), "covariate")) {
    x <- blind$roles[[r]]
    if (!isTRUE(x$found) || length(x$columns) == 0L) next
    bp <- bpos(r)
    units <- if (identical(unname(.ROLE_ARITY[r]), "pair")) split(bp, ceiling(seq_along(bp) / 2)) else list(bp)
    claimed_r <- lapply(pos_list[rows$role == r], sort)
    for (u in units) {
      if (any(vapply(claimed_r, identical, logical(1L), sort(u)))) next
      ur <- c(ur, r); uc <- c(uc, paste(nm[u], collapse = ", "))
      us <- c(us, if (is.na(x$max_score)) sprintf("bits=%.1f", x$score) else sprintf("pct=%.0f", x$pct))
    }
  }
  unclaimed <- data.frame(role = ur, columns = uc, statistic = us, stringsAsFactors = FALSE)

  # a survival claim (time + event) is ONE claim shown as two rows: count it once
  ct <- table[!(table$role == "event_variable" & "time_variable" %in% table$role), , drop = FALSE]
  n <- function(v) sum(ct$verdict == v)
  counts <- c(claims = nrow(ct), testable = n("contradicted") + n("consistent"),
              contradicted = n("contradicted"), consistent = n("consistent"),
              unverifiable = n("unverifiable"))
  structure(list(table = table, unclaimed = unclaimed, counts = counts,
                 n_obs = nrow(data), n_var = ncol(data)), class = "role_verification")
}

#' @param x,object A \code{role_verification} object.
#' @param ... Ignored.
#' @rdname verify_roles
#' @export
print.role_verification <- function(x, ...) {
  k <- x$counts
  cat(sprintf("<role_verification> %d claims: %d testable, %d contradicted, %d consistent, %d unverifiable\n",
              k[["claims"]], k[["testable"]], k[["contradicted"]], k[["consistent"]], k[["unverifiable"]]))
  print(x$table[, setdiff(names(x$table), "blind"), drop = FALSE], row.names = FALSE)
  cat("Structure not in the claim:\n")
  if (nrow(x$unclaimed) == 0L) cat("  none\n") else
    for (i in seq_len(nrow(x$unclaimed)))
      cat(sprintf("  %-18s %s (%s)\n", x$unclaimed$role[i], x$unclaimed$columns[i], x$unclaimed$statistic[i]))
  invisible(x)
}

#' @rdname verify_roles
#' @export
summary.role_verification <- function(object, ...) object$table
