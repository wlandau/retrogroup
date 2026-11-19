test_that("group_analysis_probability() has the expected structure", {
  data <- probability_example()
  x <- group_analysis_probability(data$patient, data$summary)
  expect_type(x, "list")
  expect_equal(names(x), c("age", "node4", "sex"))
  for (variable in names(x)) {
    p <- x[[variable]]
    expect_true(is.matrix(p))
    expect_type(p, "double")
    expect_equal(nrow(p), nrow(data$patient))
    expect_identical(
      colnames(p),
      data$summary$level[data$summary$variable == variable]
    )
    expect_false(anyNA(p))
    expect_true(all(p > 0 & p < 1))
    expect_equal(unname(rowSums(p)), rep(1, nrow(p)))
  }
})

test_that("group_analysis_probability() column order follows summary", {
  data <- probability_small()
  summary <- data$summary[c(4, 2, 5, 1, 3), ]
  x <- group_analysis_probability(data$patient, summary)
  expect_equal(names(x), c("sex", "age"))
  expect_identical(colnames(x$sex), c("female", "other", "male"))
  expect_identical(colnames(x$age), c("old", "young"))
  y <- group_analysis_probability(data$patient, data$summary)
  expect_equal(x$sex, y$sex[, c("female", "other", "male")])
  expect_equal(x$age, y$age[, c("old", "young")])
})

test_that("group_analysis_probability() rows follow the rows of patient", {
  data <- probability_small()
  order <- c(3, 6, 1, 5, 2, 4)
  x <- group_analysis_probability(data$patient[order, ], data$summary)
  y <- group_analysis_probability(data$patient, data$summary)
  expect_equal(x$age, y$age[order, ])
  expect_equal(x$sex, y$sex[order, ])
})

test_that("group_analysis_probability() matches a direct loop over the formula", {
  data <- probability_small()
  patient <- data$patient
  summary <- data$summary
  z <- patient$treatment
  t <- patient$time
  d <- patient$event
  m <- nrow(patient)
  r <- vapply(seq_len(m), \(i) sum(z == 0 & t >= t[i]), numeric(1))
  Lambda_0 <- vapply(
    seq_len(m),
    \(i) sum(ifelse(z == 0 & t <= t[i], d / pmax(r, 1), 0)),
    numeric(1)
  )
  x <- group_analysis_probability(patient, summary)
  for (variable in unique(summary$variable)) {
    rows <- summary$variable == variable
    n <- summary$patients[rows]
    beta <- summary$log_hazard_ratio[rows]
    expected <- matrix(0, m, sum(rows))
    for (i in seq_len(m)) {
      w <- log(n) + d[i] * beta * z[i] - Lambda_0[i] * exp(beta * z[i])
      expected[i, ] <- exp(w) / sum(exp(w))
    }
    colnames(expected) <- summary$level[rows]
    expect_equal(x[[variable]], expected)
  }
})

test_that("Nelson-Aalen estimate agrees with survival::survfit()", {
  set.seed(0)
  n <- 200L
  z <- rep(0:1, length.out = n)
  t <- round(stats::rexp(n), 1)
  d <- stats::rbinom(n, 1, 0.7)
  out <- retrogroup:::probability_nelson_aalen(z = z, t = t, d = d)
  fit <- survival::survfit(
    survival::Surv(t[z == 0], d[z == 0]) ~ 1,
    ctype = 1
  )
  step <- stats::stepfun(fit$time, c(0, fit$cumhaz))
  expect_equal(out, step(t))
})

test_that("group_analysis_probability() validates patient", {
  data <- probability_small()
  expect_error(
    group_analysis_probability(as.list(data$patient), data$summary),
    "data frame",
    class = "retrogroup_error"
  )
  expect_error(
    group_analysis_probability(
      data$patient[, c("time", "event")],
      data$summary
    ),
    "missing column",
    class = "retrogroup_error"
  )
  patient <- data$patient
  patient$treatment[1] <- 2L
  expect_error(
    group_analysis_probability(patient, data$summary),
    "only 0 and 1",
    class = "retrogroup_error"
  )
  patient <- data$patient
  patient$treatment <- 1L
  expect_error(
    group_analysis_probability(patient, data$summary),
    "patient is invalid",
    class = "retrogroup_error"
  )
  patient <- data$patient
  patient$time[1] <- -1
  expect_error(
    group_analysis_probability(patient, data$summary),
    "group_data_patient",
    class = "retrogroup_error"
  )
})

test_that("group_analysis_probability() validates summary", {
  data <- probability_small()
  expect_error(
    group_analysis_probability(data$patient, as.list(data$summary)),
    "data frame",
    class = "retrogroup_error"
  )
  expect_error(
    group_analysis_probability(data$patient, data$summary[, -4]),
    "missing column",
    class = "retrogroup_error"
  )
  summary <- data$summary
  summary$patients[1] <- 0L
  expect_error(
    group_analysis_probability(data$patient, summary),
    "summary is invalid",
    class = "retrogroup_error"
  )
})

test_that("group_analysis_probability() checks patient and summary together", {
  data <- probability_small()
  summary <- data$summary
  summary$patients <- c(3L, 4L, 2L, 2L, 3L)
  expect_error(
    group_analysis_probability(data$patient, summary),
    "must equal the number of rows of patient",
    class = "retrogroup_error"
  )
  expect_error(
    group_analysis_probability(data$patient, summary),
    "nrow\\(patient\\) = 6. Totals by subgroup variable: age = 7, sex = 7",
    class = "retrogroup_error"
  )
})

test_that("group_analysis_probability() rejects numerically zero probabilities", {
  data <- probability_small()
  summary <- data$summary
  summary$log_hazard_ratio[1] <- 50
  expect_error(
    group_analysis_probability(data$patient, summary),
    "sqrt\\(.Machine\\$double.eps\\)",
    class = "retrogroup_error"
  )
})
