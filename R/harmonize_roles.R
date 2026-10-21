.hz_stop <- function(...) stop("harmonize_roles(): ", ..., call. = FALSE)

# Detect roles name-blind on one site; pos = integer column positions per role.
.hz_detect <- function(site) {
  roles <- detect_roles(stats::setNames(site, paste0("col_", seq_along(site))))$roles
  pos <- lapply(roles, function(r) if (isTRUE(r$found)) as.integer(sub("col_", "", r$columns)) else integer(0))
  list(roles = roles, vi = .build_var_info(site), pos = pos)
}

# covariate is the residual role (leftovers); aligning it only adds noise, so it is excluded.
.hz_roles <- function(dets) setdiff(unique(unlist(lapply(dets, function(d) names(d$pos)), use.names = FALSE)), "covariate")

.hz_map <- function(dets, site_names, sites) {
  roles <- .hz_roles(dets)
  cells <- lapply(seq_along(dets), function(s) vapply(roles, function(r) {
    p <- dets[[s]]$pos[[r]]
    if (length(p)) paste(sort(names(sites[[s]])[p], method = "radix"), collapse = ", ") else NA_character_
  }, character(1L)))
  names(cells) <- site_names
  nfound <- rowSums(!is.na(as.data.frame(cells, check.names = FALSE)))
  out <- data.frame(role = roles,
                    status = ifelse(nfound == length(dets), "aligned", ifelse(nfound == 0L, "absent", "partial")),
                    stringsAsFactors = FALSE)
  for (s in site_names) out[[s]] <- cells[[s]]
  rownames(out) <- NULL   # vapply names would otherwise leak into row names
  out
}

.hz_jaccard <- function(sets) {
  if (length(sets) < 2L) return(NA_real_)
  mean(utils::combn(length(sets), 2L, function(i) {
    u <- union(sets[[i[1]]], sets[[i[2]]]); if (!length(u)) 1 else length(intersect(sets[[i[1]]], sets[[i[2]]])) / length(u)
  }))
}
.hz_overlap <- function(rng) {   # rng: list of c(min, max)
  mean(utils::combn(length(rng), 2L, function(i) {
    a <- rng[[i[1]]]; b <- rng[[i[2]]]
    un <- max(a[2], b[2]) - min(a[1], b[1])
    if (un <= 0) 1 else max(0, min(a[2], b[2]) - max(a[1], b[1])) / un
  }))
}

.hz_checks <- function(dets, sites, site_names = names(sites)) {
  rows <- list()
  add <- function(role, check, vals, found, flag = NA_character_) {
    r <- data.frame(role = role, check = check, stringsAsFactors = FALSE)
    for (k in seq_along(site_names)) r[[site_names[k]]] <- if (found[k]) vals[[k]] else NA_character_
    r$flag <- flag
    rows[[length(rows) + 1L]] <<- r
  }
  for (role in .hz_roles(dets)) {
    found <- vapply(dets, function(d) length(d$pos[[role]]) > 0L, logical(1L))
    if (sum(found) < 2L) next
    # first(s) assumes a scalar role (one column); only scalar roles reach the level/range checks
    first <- function(s) sites[[s]][[dets[[s]]$pos[[role]][1]]]   # scalar roles: the single column
    each <- function(f) lapply(seq_along(dets), function(s) if (found[s]) f(s) else NA_character_)
    differ <- function(v) length(unique(unlist(v[found]))) > 1L
    ty <- each(function(s) paste(sort(dets[[s]]$vi$user_type[dets[[s]]$pos[[role]]], method = "radix"), collapse = ","))
    add(role, "type", ty, found, if (differ(ty)) "type_mismatch" else NA_character_)
    lv <- function(s) { x <- first(s); paste(sort(unique(as.character(x[!is.na(x)])), method = "radix"), collapse = ", ") }
    if (role %in% c("event_variable", "outcome_binary")) {
      add(role, "levels", each(lv), found)
      # 2-level codings (TRUE/FALSE, {1,2}, {0,1}) normalize to {0,1}; a 3-level {0,1,2} keeps its
      # three-level scheme (0 = censor, >0 = event is a different coding), so it differs from {0,1}.
      cd <- each(function(s) { x <- first(s); e <- .normalize_event(x)
        if (is.null(e)) "none" else if (length(unique(x[!is.na(x)])) == 3L) "0,1,2" else "0,1" })
      add(role, "coding", cd, found, if (differ(cd)) "coding_mismatch" else NA_character_)
    } else if (role == "group_var") {
      add(role, "levels", each(lv), found)
      sets <- lapply(which(found), function(s) { x <- first(s); unique(as.character(x[!is.na(x)])) })
      add(role, "level_jaccard", rep(list(sprintf("%.2f", .hz_jaccard(sets))), length(dets)), found)
    } else if (role %in% c("outcome_continuous", "time_variable")) {
      num <- function(s) suppressWarnings(as.numeric(first(s)))
      rg <- each(function(s) { x <- num(s); x <- x[is.finite(x)]
        if (length(x)) sprintf("%s/%s/%s", format(min(x), digits = 4), format(stats::median(x), digits = 4), format(max(x), digits = 4)) else NA_character_ })
      add(role, "range", rg, found)
      rr <- lapply(which(found), function(s) { x <- num(s); x <- x[is.finite(x)]; if (length(x)) range(x) else c(NA_real_, NA_real_) })
      if (!anyNA(unlist(rr))) add(role, "range_overlap", rep(list(sprintf("%.2f", .hz_overlap(rr))), length(dets)), found)
    } else if (role %in% c("scale_items", "repeated_measures")) {
      sz <- each(function(s) as.character(length(dets[[s]]$pos[[role]])))
      add(role, "size", sz, found, if (differ(sz)) "size_mismatch" else NA_character_)
    }
  }
  if (!length(rows)) {
    out <- data.frame(role = character(0), check = character(0), stringsAsFactors = FALSE)
    for (s in site_names) out[[s]] <- character(0)
    out$flag <- character(0)
    return(out)
  }
  do.call(rbind, rows)
}

.hz_site_only <- function(map) {
  sc <- setdiff(names(map), c("role", "status"))
  rows <- lapply(seq_len(nrow(map)), function(i) {
    v <- unlist(map[i, sc, drop = TRUE]); ok <- !is.na(v)
    if (sum(ok) == 1L) data.frame(site = sc[ok], role = map$role[i], columns = unname(v[ok]), stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(data.frame(site = character(0), role = character(0), columns = character(0), stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

#' Align detected variable roles across sites
#'
#' Runs [detect_roles()] name-blind on each site and reports, per role, which
#' column each site holds, plus descriptive consistency checks. Alignment is by
#' role only; column names appear only in the returned cells. The residual `covariate`
#' role is excluded. `aligned` means the detector found the role at every site; it does not establish that the columns measure the same thing; read `checks` (levels, Jaccard, ranges) before pooling. Cell order is alphabetical and carries no pairing direction.
#'
#' @param sites A named list (unique, non-empty names) of at least two data.frames.
#' @return An object of class `"role_harmonization"`: `map`, `checks`,
#'   `site_only`, `n_sites`, `sites`. Flags are `type_mismatch`,
#'   `coding_mismatch`, `size_mismatch`; no numeric threshold decides a flag.
#' @param x,object A `role_harmonization` object.
#' @param ... Ignored.
#' @export
harmonize_roles <- function(sites) {
  if (!is.list(sites) || is.data.frame(sites) || length(sites) < 2L) .hz_stop("`sites` must be a list of at least two data.frames.")
  nm <- names(sites)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm))) .hz_stop("`sites` must have non-empty names.")
  if (anyDuplicated(nm)) .hz_stop("site names must be unique (", paste(unique(nm[duplicated(nm)]), collapse = ", "), ").")
  if (!all(vapply(sites, is.data.frame, logical(1L)))) .hz_stop("every element of `sites` must be a data.frame.")
  if (any(vapply(sites, ncol, integer(1L)) < 1L)) .hz_stop("every site needs at least one column.")
  dets <- lapply(sites, .hz_detect)
  map <- .hz_map(dets, nm, sites)
  structure(list(map = map, checks = .hz_checks(dets, sites, nm), site_only = .hz_site_only(map),
                 n_sites = length(sites), sites = nm), class = "role_harmonization")
}

#' @rdname harmonize_roles
#' @export
print.role_harmonization <- function(x, ...) {
  st <- table(factor(x$map$status, c("aligned", "partial", "absent")))
  cat(sprintf("<role_harmonization> %d sites: %d aligned, %d partial roles\n\n",
              x$n_sites, st[["aligned"]], st[["partial"]]))
  print(x$map[x$map$status != "absent", , drop = FALSE], row.names = FALSE)
  fl <- x$checks[!is.na(x$checks$flag) | x$checks$check %in% c("level_jaccard", "range_overlap"), , drop = FALSE]
  if (nrow(fl)) { cat("\nFlagged and overlap checks\n"); print(fl, row.names = FALSE) }
  if (nrow(x$site_only)) { cat("\nSite-only structure\n"); print(x$site_only, row.names = FALSE) }
  cat("\nDescriptive checks (levels, Jaccard, ranges) are in $checks; flags mark structural mismatches only.\n")
  invisible(x)
}

#' @rdname harmonize_roles
#' @export
summary.role_harmonization <- function(object, ...) object$map
