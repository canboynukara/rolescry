# Behaviour introduced in 0.3.0. Constructions are synthetic and construction-true.

make_surv <- function(n, seed, time_fun) {
  set.seed(seed)
  t <- time_fun(n)
  ev <- rbinom(n, 1, plogis(2 - 4 * (rank(t) / n)))   # early times -> events
  data.frame(t = t, ev = ev)
}

test_that("an all-distinct integer follow-up time is a time candidate, not a subject id", {
  d <- make_surv(60, 11, function(n) sample(30:3000, n))   # whole days, all distinct
  d$x <- rnorm(60)
  expect_identical(.build_var_info(d)$user_type[1], "ID")  # value-types as an identifier
  r <- detect_roles(d)$roles
  expect_identical(r$time_variable$columns, "t")
  expect_identical(r$event_variable$columns, "ev")
  expect_false("t" %in% r$subject_id$columns)
})

test_that("an identifier independent of the event stays a subject id", {
  d <- make_surv(200, 12, function(n) round(rweibull(n, 1.2, 400), 1))
  set.seed(13); d$id <- sample.int(50 * 200, 200)
  r <- detect_roles(d)$roles
  expect_identical(r$subject_id$columns, "id")
  expect_identical(r$time_variable$columns, "t")
})

test_that("survival pairs are compared on a common set of complete rows", {
  d <- make_surv(200, 14, function(n) round(rweibull(n, 1.2, 400), 1))
  set.seed(15)
  d$g <- rbinom(200, 1, 0.5)
  d$m <- rexp(200, 1 / 50)
  d$m[sample(200, 60)] <- NA                               # distractor column with missingness
  r <- detect_roles(d)$roles
  expect_identical(r$time_variable$columns, "t")
  expect_identical(r$event_variable$columns, "ev")
})

test_that("covariate is the residual role and never overlaps another role", {
  set.seed(16); n <- 200; pre <- rnorm(n, 10, 2)
  d <- data.frame(arm = rep(c("a", "b"), each = n / 2), pre = pre,
                  post = pre + rnorm(n, 1, 1), z1 = rnorm(n), z2 = rnorm(n))
  d$resp <- rbinom(n, 1, plogis(1.2 * d$z1 + 1.2 * d$z2))
  r <- detect_roles(d)$roles
  others <- unlist(lapply(setdiff(names(r), "covariate"), function(k) unlist(r[[k]]$columns)))
  expect_length(intersect(r$covariate$columns, others), 0)
  expect_setequal(r$covariate$columns, setdiff(names(d), others))
})

test_that("tied continuous values give finite dependence and no warnings", {
  set.seed(17)
  x <- round(rnorm(150), 0); y <- x + round(rnorm(150), 0)  # heavy ties
  expect_true(is.finite(.dep(x + 0.0, y + 0.0)))
  d <- data.frame(a = x, b = y, c = round(runif(150) * 3), w = rbinom(150, 1, 0.5))
  expect_no_warning(detect_roles(d))
})

test_that("a logical TRUE/FALSE event column is a survival candidate", {
  d <- make_surv(200, 21, function(n) round(rweibull(n, 1.2, 400), 1))
  d$ev <- d$ev == 1L                                       # logical coding
  r <- detect_roles(d)$roles
  expect_identical(r$time_variable$columns, "t")
  expect_identical(r$event_variable$columns, "ev")
})

test_that("an RCT-like table yields coupled_pair, never a survival pair", {
  set.seed(22); n <- 300
  arm <- rep(0:1, each = n / 2); b <- rnorm(n, 150, 15)
  e <- b - 10 * arm + rnorm(n, 0, 8)
  d <- data.frame(arm = arm, sbp_base = b, sbp_end = e,
                  y = rbinom(n, 1, plogis(0.05 * (e - 140))))
  r <- detect_roles(d)$roles
  expect_false(isTRUE(r$time_variable$found))
  expect_false(isTRUE(r$event_variable$found))
  expect_true(isTRUE(r$coupled_pair$found))
  expect_length(unlist(r$coupled_pair$columns), 2L)
  expect_length(r$coupled_pair$interpretations, 2L)
})

test_that("censoring geometry decisive -> survival roles; interval-scale mimic -> not", {
  d <- make_surv(300, 23, function(n) round(rweibull(n, 1.3, 400), 1))
  expect_true(.censoring_decisive(d$t, d$ev))
  set.seed(24); sbp <- rnorm(300, 150, 15)                 # interval-scale, far from zero
  arm <- as.integer(sbp < stats::median(sbp)) ; arm[sample(300, 60)] <- 1 - arm[sample(300, 60)]
  expect_false(.censoring_decisive(sbp, arm))
})

test_that("sink veto: a collider is not swallowed by repeated_measures", {
  set.seed(25); n <- 300
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- d$x1 + d$x2 + rnorm(n, 0, 0.5)
  r <- detect_roles(d)$roles
  expect_false("y" %in% unlist(r$repeated_measures$columns))
})

test_that("genuine wide repeated measures survive the sink veto", {
  set.seed(26); n <- 150; s <- rnorm(n, 50, 10)
  d <- data.frame(id = sample.int(9999, n),
                  v1 = s + rnorm(n, 0, 3), v2 = s + 2 + rnorm(n, 0, 3),
                  v3 = s + 4 + rnorm(n, 0, 3), v4 = s + 6 + rnorm(n, 0, 3))
  r <- detect_roles(d)$roles
  expect_true(all(paste0("v", 1:4) %in% unlist(r$repeated_measures$columns)))
})

test_that("Datetime columns are typed and reported as time_index", {
  d <- data.frame(dt = as.Date("2024-01-01") + 0:199,
                  price = cumsum(rnorm(200)) + 500, vol = rlnorm(200, 11, 0.5))
  vi <- .build_var_info(d)
  expect_identical(vi$user_type[vi$column == "dt"], "Datetime")
  r <- detect_roles(d)$roles
  expect_identical(r$time_index$columns, "dt")
  expect_false("dt" %in% r$covariate$columns)
  # character-coded dates too
  d2 <- data.frame(day = format(as.Date("2024-01-01") + 0:49, "%Y-%m-%d"), x = rnorm(50))
  expect_identical(detect_roles(d2)$roles$time_index$columns, "day")
})

test_that("panel entity axis is picked by serial dependence; measures intact", {
  set.seed(41); firms <- 40; mo <- 12
  d <- do.call(rbind, lapply(1:firms, function(f) {
    p <- numeric(mo); p[1] <- rnorm(1, 100, 5)
    for (t in 2:mo) p[t] <- p[t - 1] + rnorm(1, 0, 1)      # AR within-firm price path
    data.frame(firm = 10000 + f, month = 1:mo, price = p, vol = rlnorm(mo, 12, 0.8)) }))
  r <- detect_roles(d)$roles
  expect_identical(r$subject_id$columns, "firm")
  expect_true("price" %in% r$repeated_measures$columns)
})

test_that("a value-symmetric balanced panel is reported as axis-ambiguous or resolved, never mislabelled", {
  set.seed(44); firms <- 40; mo <- 12
  d <- do.call(rbind, lapply(1:firms, function(f)
    data.frame(firm = 10000 + f, month = 1:mo, ret = rnorm(mo, 0, 0.05), vol = rlnorm(mo, 12, 0.8))))
  r <- detect_roles(d)$roles
  sid <- r$subject_id$columns
  # contract: either the entity is resolved to firm, or subject_id is silent with
  # the axis readings reported on repeated_measures -- month is never the entity
  expect_false(identical(sid, "month"))
  if (length(sid) == 0) expect_true(length(r$repeated_measures$panel_axes) >= 2)
})

test_that("a 4-item rating battery is reported as scale_items", {
  set.seed(42); n <- 300; q <- rnorm(n)
  mk <- function() pmin(9, pmax(1, round(5 + 1.5 * q + rnorm(n, 0, 0.8))))
  d <- data.frame(r1 = mk(), r2 = mk(), r3 = mk(), r4 = mk(), e = rlnorm(n, 16, 1))
  expect_setequal(detect_roles(d)$roles$scale_items$columns, paste0("r", 1:4))
  # 3 mutually independent integer columns must NOT become a scale
  set.seed(43)
  d2 <- data.frame(a = sample(1:5, 300, TRUE), b = sample(1:5, 300, TRUE),
                   c = sample(1:5, 300, TRUE), x = rnorm(300))
  expect_false(isTRUE(detect_roles(d2)$roles$scale_items$found))
})

test_that("finance name-bonus dictionary is available and capped channel only", {
  nb <- rolescry_default_name_bonus("finance")
  expect_true("default" %in% nb$event_variable)
  expect_setequal(names(rolescry_default_name_bonus("all")),
                  union(names(rolescry_default_name_bonus()), names(nb)))
})

test_that("two disjoint decisive pairs without a dominance margin stay a coupled_pair", {
  a <- make_surv(300, 51, function(n) round(rweibull(n, 1.3, 400), 1))
  b <- make_surv(300, 52, function(n) round(rweibull(n, 1.3, 400), 1))
  d <- data.frame(t1 = a$t, e1 = a$ev, t2 = b$t, e2 = b$ev)   # two independent survival processes
  r <- detect_roles(d)$roles
  expect_false(isTRUE(r$time_variable$found))
  expect_true(isTRUE(r$coupled_pair$found))
  expect_match(r$coupled_pair$reason, "competing pair")
})

test_that("a same-time rival event does not veto the survival reading", {
  d <- make_surv(300, 51, function(n) round(rweibull(n, 1.3, 400), 1))
  set.seed(52)
  d$hep <- as.integer(d$ev == 1 | rbinom(300, 1, 0.2) == 1)   # correlate of the event, not a reading
  r <- detect_roles(d)$roles
  expect_identical(r$time_variable$columns, "t")
  expect_identical(r$event_variable$columns, "ev")
})

test_that("a dominant survival pair survives a weakly-associated decisive rival", {
  d <- make_surv(300, 53, function(n) round(rweibull(n, 1.3, 400), 1))
  set.seed(54)
  d$x <- rnorm(300, 10, 2)
  r <- detect_roles(d)$roles
  expect_identical(r$time_variable$columns, "t")
  expect_identical(r$event_variable$columns, "ev")
})

test_that("a real logistic signal is detected as the binary outcome", {
  set.seed(55); n <- 300; x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  d <- data.frame(x1, x2, x3, y = rbinom(n, 1, plogis(1.5 * x1)))
  expect_identical(detect_roles(d)$roles$outcome_binary$columns, "y")
})

test_that("group_var is the most balanced categorical; level-count points no longer decide", {
  set.seed(41L); n <- 240L
  d <- data.frame(arm = sample(c("A","B","C"), n, TRUE),                       # balanced 3 arms
                  region = sample(c("n","s"), n, TRUE, prob = c(0.8, 0.2)),    # 2 levels, unbalanced
                  y = rnorm(n))
  r <- detect_roles(d)$roles$group_var
  expect_identical(r$columns, "arm")
  expect_equal(r$components[[1]]$name, "normalized entropy")
})
test_that("a categorical with a level seen fewer than 3 times is not a group candidate", {
  d <- data.frame(g = c(rep("a", 50), rep("b", 50), "c", "c"), y = rnorm(102))
  expect_false(isTRUE(detect_roles(d)$roles$group_var$found))
})
test_that("an all-NA or constant categorical does not break group detection", {
  d <- data.frame(g = rep(c("x","y"), 30), z = rep(NA_character_, 60), w = rep("k", 60), y = rnorm(60),
                  stringsAsFactors = FALSE)
  expect_identical(detect_roles(d)$roles$group_var$columns, "g")
})

test_that("equal B: fewer levels wins; value-identical categoricals are silenced", {
  d <- data.frame(g3 = rep(c("a","b","c"), 40), g2 = rep(c("x","y"), 60), y = rnorm(120),
                  stringsAsFactors = FALSE)
  expect_identical(detect_roles(d)$roles$group_var$columns, "g2")
  d2 <- data.frame(g = rep(c("x","y"), 30), h = rep(c("p","q"), 30), y = rnorm(60),
                   stringsAsFactors = FALSE)
  expect_false(isTRUE(detect_roles(d2)$roles$group_var$found))
})

test_that("a pre/post pair with a systematic shift is paired_pairs and not agreement_pairs", {
  set.seed(51L); n <- 200L; pre <- rnorm(n, 30, 5); post <- pre + 4 + rnorm(n, 0, 2)
  r <- detect_roles(data.frame(pre, post, age = rnorm(n, 50, 8)))$roles
  expect_setequal(r$paired_pairs$columns, c("pre", "post"))
  expect_false(isTRUE(r$agreement_pairs$found))
})
test_that("two readings without bias are agreement_pairs and not paired_pairs", {
  set.seed(52L); n <- 200L; tr <- rnorm(n, 30, 5); a <- tr + rnorm(n, 0, 1.5); b <- tr + rnorm(n, 0, 1.5)
  r <- detect_roles(data.frame(a, b, age = rnorm(n, 50, 8)))$roles
  expect_setequal(r$agreement_pairs$columns, c("a", "b"))
  expect_false(isTRUE(r$paired_pairs$found))
})
test_that(".mdl_shift_bits is finite-safe on an exact duplicate column", {
  a <- rnorm(50); expect_identical(rolescry:::.mdl_shift_bits(a, a), -Inf)
  r <- detect_roles(data.frame(a, b = a, c = rnorm(50)))$roles
  expect_true(isTRUE(r$agreement_pairs$found))
  expect_false(isTRUE(r$paired_pairs$found))
})
test_that("an exact constant offset is paired_pairs and agreement is silent", {
  set.seed(53L); a <- rnorm(60, 30, 5)
  r <- detect_roles(data.frame(a, b = a + 4, c = rnorm(60)))$roles
  expect_true(isTRUE(r$paired_pairs$found))
  expect_false(isTRUE(r$agreement_pairs$found))
})

test_that("a 6-item Likert battery is scale_items via the factor code and alpha points are gone", {
  set.seed(81L); n <- 200L; lat <- rnorm(n)
  items <- as.data.frame(sapply(1:6, function(i) pmin(5, pmax(1, round(3 + 0.9 * lat + rnorm(n, 0, .8))))))
  names(items) <- paste0("q", 1:6); items$age <- rnorm(n, 40, 10)
  r <- detect_roles(items)$roles$scale_items
  expect_setequal(r$columns, paste0("q", 1:6))
  expect_equal(r$components[[1]]$name, "one-factor MDL block")
})
test_that("a one-hot partition of eight is not a scale_items battery", {
  set.seed(82L); n <- 200L
  oh <- as.data.frame(model.matrix(~ 0 + f, data.frame(f = factor(sample(1:8, n, TRUE))))); names(oh) <- paste0("h", 1:8)
  oh$x <- rnorm(n)
  expect_false(isTRUE(detect_roles(oh)$roles$scale_items$found))
})
