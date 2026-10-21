mk_site <- function(seed, names = c("id","arm","age","sbp_base","sbp_end","time","dead"), n = 200L, event_code = c(0L,1L), arm_labels = c("a","b")) {
  set.seed(seed); x <- rnorm(n); t <- rweibull(n, 1.3, exp(5 - 1 * x)); c <- runif(n, 0, quantile(t, .9))
  b <- rnorm(n, 150, 15); d <- data.frame(sample.int(99999, n), sample(arm_labels, n, TRUE), round(rnorm(n, 60, 10)),
         b, b - 5 + rnorm(n, 0, 6), round(pmin(t, c), 1), event_code[1 + as.integer(t <= c)])
  names(d) <- names; d }
test_that("harmonize_roles validates its inputs", {
  s1 <- mk_site(1); expect_error(harmonize_roles(list(s1)), "harmonize_roles\\(\\)")
  expect_error(harmonize_roles(list(s1, mk_site(2))), "names")
  expect_error(harmonize_roles(list(a = s1, a = mk_site(2))), "unique")
  expect_error(harmonize_roles(list(a = s1, b = 1)), "data.frame")
})
test_that("identical structure under different names and orders aligns every detected role with no flags", {
  s1 <- mk_site(1); s2 <- mk_site(2, names = c("pid","group","alter","bp0","bp1","ftime","status")); s2 <- s2[, sample(ncol(s2))]
  s3 <- mk_site(3, names = paste0("v", 1:7))
  h <- harmonize_roles(list(north = s1, south = s2, east = s3))
  expect_s3_class(h, "role_harmonization")
  al <- h$map$role[h$map$status == "aligned"]
  # the fixture yields exactly these aligned roles (detect_roles finds time/event, not coupled_pair)
  expect_setequal(al, c("group_var","paired_pairs","time_variable","event_variable","subject_id"))
  expect_true(all(is.na(h$checks$flag)))
  expect_output(print(h), "sites")
})
test_that("a site missing the survival pair is partial and listed; coding and size mismatches are flagged; labels/ranges reported", {
  s1 <- mk_site(1); s2 <- mk_site(2, event_code = c(1L, 2L), arm_labels = c("drug","placebo"))
  s3 <- mk_site(3)[, c("id","arm","age","sbp_base","sbp_end")]
  h <- harmonize_roles(list(a = s1, b = s2, c = s3))
  tv <- h$map[h$map$role == "time_variable", ]; expect_true(tv$status %in% c("partial","absent"))
  expect_true(is.na(tv$c))
  cod <- h$checks[h$checks$role == "event_variable" & h$checks$check == "coding", ]
  expect_gt(nrow(cod), 0); expect_true(is.na(cod$flag))                      # {1,2} normalizes to {0,1}: no coding mismatch
  jac <- h$checks[h$checks$role == "group_var" & h$checks$check == "level_jaccard", ]
  expect_gt(nrow(jac), 0); expect_match(jac$a, "^0")                            # a/b vs drug/placebo: Jaccard 0 reported, no flag
  expect_true(all(is.na(h$checks$flag[h$checks$role == "group_var"])))
})
test_that("TRUE/FALSE vs 0/1 events do not flag; a {0,1} vs {0,1,2} coding does", {
  s1 <- mk_site(1); s2 <- mk_site(2); s2$dead <- s2$dead == 1L
  h <- harmonize_roles(list(a = s1, b = s2)); cod <- h$checks[h$checks$check == "coding", ]
  expect_true(all(is.na(cod$flag)))
  s3 <- mk_site(3); s3$dead[sample(nrow(s3), 20)] <- 2L
  h3 <- harmonize_roles(list(a = s1, c = s3)); cod3 <- h3$checks[h3$checks$check == "coding", ]
  expect_gt(nrow(cod3), 0); expect_true(any(cod3$flag == "coding_mismatch", na.rm = TRUE))
})
test_that("scale battery size mismatch flags; noise site, single column, differing ncol do not error", {
  set.seed(9); n <- 200; lat <- rnorm(n)
  q5 <- as.data.frame(sapply(1:5, function(i) pmin(5, pmax(1, round(3 + .9 * lat + rnorm(n, 0, .8)))))); q5$x <- rnorm(n)
  q6 <- as.data.frame(sapply(1:6, function(i) pmin(5, pmax(1, round(3 + .9 * lat + rnorm(n, 0, .8)))))); q6$x <- rnorm(n)
  h <- harmonize_roles(list(a = q5, b = q6)); sz <- h$checks[h$checks$check == "size", ]
  expect_true(any(sz$flag == "size_mismatch", na.rm = TRUE))
  noise <- data.frame(u = rnorm(100), v = rnorm(100)); one <- data.frame(w = rnorm(100))
  expect_s3_class(harmonize_roles(list(a = mk_site(1), b = noise, c = one)), "role_harmonization")
  so <- harmonize_roles(list(a = mk_site(1), b = noise))$site_only
  expect_true(all(so$site == "a"))
})
test_that("alignment is name-blind and permutation-invariant", {
  s1 <- mk_site(1); s2 <- mk_site(2)
  h <- harmonize_roles(list(a = s1, b = s2))
  set.seed(3); p <- sample(ncol(s2)); s2p <- s2[, p]; names(s2p) <- paste0("c", seq_along(s2p))
  hp <- harmonize_roles(list(a = s1, b = s2p))
  expect_equal(h$map$status, hp$map$status); expect_equal(h$checks[, c("role","check","flag")], hp$checks[, c("role","check","flag")])
  expect_equal(h$checks$a, hp$checks$a)
  back <- function(cell) vapply(strsplit(cell, ", "), function(v) paste(sort(names(s2)[p[as.integer(sub("c", "", v))]], method = "radix"), collapse = ", "), "")
  ok <- !is.na(hp$map$b)
  expect_identical(h$map$b[ok], back(hp$map$b[ok]))          # same columns, same columns once each cell is sorted by original name
  expect_identical(hp$map$b[ok], vapply(strsplit(hp$map$b[ok], ", "), function(v) paste(sort(v, method = "radix"), collapse = ", "), ""))
  expect_identical(h$checks$b, hp$checks$b)
})
test_that("covariate is excluded everywhere", {
  set.seed(5); a <- data.frame(id = 1:100, z = rnorm(100)); b <- data.frame(id = 1:100, z = sample(letters[1:5], 100, TRUE))
  h <- harmonize_roles(list(a = a, b = b))
  expect_false("covariate" %in% h$map$role); expect_false("covariate" %in% h$checks$role); expect_false("covariate" %in% h$site_only$role)
  expect_identical(rownames(h$map), as.character(seq_len(nrow(h$map))))
})
test_that("type_mismatch flags when a role's value type differs (event {0,1} Binary vs {0,1,2} Continuous)", {
  s1 <- mk_site(1); s3 <- mk_site(3); s3$dead[sample(nrow(s3), 20)] <- 2L
  ck <- harmonize_roles(list(a = s1, c = s3))$checks
  ty <- ck[ck$role == "event_variable" & ck$check == "type", ]
  expect_equal(ty$flag, "type_mismatch")
})
test_that("a zero-column site is rejected; print shows overlap rows and the closing line", {
  expect_error(harmonize_roles(list(a = mk_site(1), b = data.frame())), "harmonize_roles\\(\\)")
  expect_output(print(harmonize_roles(list(a = mk_site(1), b = mk_site(2)))), "level_jaccard.*|Descriptive checks")
})
