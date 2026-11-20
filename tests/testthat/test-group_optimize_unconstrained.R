test_that("group_optimize_unconstrained() has the expected structure", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_unconstrained(data$patient, data$summary)
  sample$extra <- seq_len(nrow(sample))
  x <- suppressWarnings(group_optimize_unconstrained(
    sample,
    data$summary,
    iterations = 20L
  ))
  expect_s3_class(x, "tbl_df")
  expect_identical(colnames(x), colnames(sample))
  expect_equal(nrow(x), nrow(sample))
  other <- c(colnames(data$patient), "extra")
  expect_identical(x[, other], sample[, other])
  for (variable in unique(data$summary$variable)) {
    levels <- data$summary$level[data$summary$variable == variable]
    expect_s3_class(x[[variable]], "factor")
    expect_identical(levels(x[[variable]]), levels)
    expect_false(anyNA(x[[variable]]))
    expect_true(all(x[[variable]] %in% levels))
  }
})

test_that("group_optimize_unconstrained() reduces the loss", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_unconstrained(data$patient, data$summary)
  x <- suppressWarnings(group_optimize_unconstrained(
    sample,
    data$summary,
    iterations = 50L
  ))
  for (variable in unique(data$summary$variable)) {
    before <- optimize_loss_reference(sample, data$summary, variable)
    after <- optimize_loss_reference(x, data$summary, variable)
    expect_lt(after, before)
  }
})

test_that("group_optimize_unconstrained() stops early", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_unconstrained(data$patient, data$summary)
  expect_identical(
    suppressWarnings(
      group_optimize_unconstrained(sample, data$summary, iterations = 0L)
    ),
    sample
  )
  expect_identical(
    group_optimize_unconstrained(sample, data$summary, tolerance = Inf),
    sample
  )
})

test_that("group_optimize_unconstrained() is reproducible with set.seed()", {
  data <- probability_example()
  summary <- data$summary[data$summary$variable == "age", ]
  set.seed(0)
  sample <- group_sample_unconstrained(data$patient, summary)
  set.seed(1)
  x <- suppressWarnings(group_optimize_unconstrained(
    sample,
    summary,
    iterations = 20L
  ))
  set.seed(1)
  y <- suppressWarnings(group_optimize_unconstrained(
    sample,
    summary,
    iterations = 20L
  ))
  expect_identical(x, y)
})

test_that("group_optimize_unconstrained() handles factor levels in patient", {
  data <- probability_example()
  summary <- data$summary[data$summary$variable == "age", ]
  set.seed(0)
  sample <- group_sample_unconstrained(data$patient, summary)
  factor <- sample
  factor$age <- factor(factor$age)
  set.seed(1)
  x <- suppressWarnings(group_optimize_unconstrained(
    sample,
    summary,
    iterations = 5L
  ))
  set.seed(1)
  y <- suppressWarnings(group_optimize_unconstrained(
    factor,
    summary,
    iterations = 5L
  ))
  expect_identical(x, y)
})

test_that("optimize_loss_term() is on the hazard ratio scale", {
  expect_equal(retrogroup:::optimize_loss_term(log(1.5), 1.2), 0.09)
  expect_equal(retrogroup:::optimize_loss_term(0, 1), 0)
  expect_identical(retrogroup:::optimize_loss_term(NA_real_, 1.2), Inf)
})

test_that("optimize_loss_terms() matches the reference loss", {
  data <- probability_example()
  summary <- data$summary
  set.seed(0)
  sample <- group_sample_unconstrained(data$patient, summary)
  for (variable in unique(summary$variable)) {
    rows <- summary$variable == variable
    terms <- retrogroup:::optimize_loss_terms(
      x = match(as.character(sample[[variable]]), summary$level[rows]),
      z = sample$treatment,
      t = sample$time,
      d = sample$event,
      hazard_ratio = summary$hazard_ratio[rows]
    )
    expect_equal(
      sqrt(sum(terms)),
      optimize_loss_reference(sample, summary, variable)
    )
  }
})

test_that("optimize_unconstrained_variable() never increases the loss", {
  data <- probability_small()
  patient <- data$patient
  set.seed(0)
  index <- sample(1:2, nrow(patient), replace = TRUE)
  hazard_ratio <- exp(c(-0.5, 0.3))
  args <- list(z = patient$treatment, t = patient$time, d = patient$event)
  loss <- function(index) {
    sqrt(sum(do.call(
      retrogroup:::optimize_loss_terms,
      c(list(x = index), args, list(hazard_ratio = hazard_ratio))
    )))
  }
  before <- loss(index)
  for (iterations in c(1L, 5L, 20L)) {
    result <- do.call(
      retrogroup:::optimize_unconstrained_variable,
      c(
        list(x = index),
        args,
        list(hazard_ratio = hazard_ratio),
        list(iterations = iterations, tolerance = 0)
      )
    )
    after <- loss(result$x)
    expect_identical(result$loss, after)
    expect_false(isTRUE(after > before))
  }
})

test_that("group_optimize_unconstrained() warns if tolerance is not reached", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_unconstrained(data$patient, data$summary)
  expect_warning(
    group_optimize_unconstrained(sample, data$summary, iterations = 2L),
    "did not fall below tolerance = 0.001 within iterations = 2 for 3",
    class = "retrogroup_warning"
  )
  warning <- tryCatch(
    group_optimize_unconstrained(sample, data$summary, iterations = 2L),
    retrogroup_warning = \(warning) warning
  )
  message <- conditionMessage(warning)
  expect_match(message, "\"age\" = ", fixed = TRUE)
  expect_match(message, "\"node4\" = ", fixed = TRUE)
  expect_match(message, "\"sex\" = ", fixed = TRUE)
  expect_match(message, "Increase iterations", fixed = TRUE)
  summary <- data$summary[data$summary$variable == "age", ]
  expect_warning(
    group_optimize_unconstrained(sample, summary, iterations = 0L),
    "for 1 subgroup variable",
    class = "retrogroup_warning"
  )
})

test_that("group_optimize_unconstrained() does not warn if tolerance is reached", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_unconstrained(data$patient, data$summary)
  expect_no_warning(
    group_optimize_unconstrained(sample, data$summary, tolerance = Inf)
  )
  expect_no_warning(
    group_optimize_unconstrained(
      sample,
      data$summary,
      iterations = 0L,
      tolerance = Inf
    )
  )
})

test_that("warn_optimize_tolerance() only lists unconverged variables", {
  loss <- c(age = 0.5, sex = 1e-8, node4 = Inf)
  warning <- tryCatch(
    retrogroup:::warn_optimize_tolerance(loss, 10L, 1e-6),
    retrogroup_warning = \(warning) warning
  )
  message <- conditionMessage(warning)
  expect_match(message, "for 2 subgroup variable", fixed = TRUE)
  expect_match(message, "\"age\" = ", fixed = TRUE)
  expect_match(message, "\"node4\" = ", fixed = TRUE)
  expect_false(grepl("\"sex\"", message, fixed = TRUE))
  expect_match(message, "linear hazard ratio scale", fixed = TRUE)
  expect_no_warning(
    retrogroup:::warn_optimize_tolerance(c(age = 1e-8), 10L, 1e-6)
  )
})

test_that("both optimizers default to tolerance = 1e-3", {
  expect_identical(formals(group_optimize_unconstrained)$tolerance, 1e-3)
  expect_identical(formals(group_optimize_constrained)$tolerance, 1e-3)
})

test_that("a loss below tolerance bounds every hazard ratio error", {
  set.seed(0)
  for (rep in seq_len(100L)) {
    observed <- exp(rnorm(4L, sd = 0.5))
    assumed <- exp(rnorm(4L, sd = 0.5))
    terms <- vapply(
      seq_along(observed),
      \(k) retrogroup:::optimize_loss_term(log(observed[k]), assumed[k]),
      numeric(1L)
    )
    loss <- sqrt(sum(terms))
    expect_true(all(abs(observed - assumed) <= loss))
  }
})

test_that("group_optimize_unconstrained() validates its arguments", {
  data <- probability_example()
  summary <- data$summary[data$summary$variable == "age", ]
  set.seed(0)
  sample <- group_sample_unconstrained(data$patient, summary)
  expect_error(
    group_optimize_unconstrained(data$patient, summary),
    "Missing subgroup column",
    class = "retrogroup_error"
  )
  bad <- sample
  bad$age <- as.character(bad$age)
  bad$age[1] <- "nonexistent"
  expect_error(
    group_optimize_unconstrained(bad, summary),
    "only levels",
    class = "retrogroup_error"
  )
  bad$age[1] <- NA_character_
  expect_error(
    group_optimize_unconstrained(bad, summary),
    "only levels",
    class = "retrogroup_error"
  )
  for (iterations in list(-1, 1.5, NA, c(1, 2), "1", Inf)) {
    expect_error(
      group_optimize_unconstrained(sample, summary, iterations = iterations),
      "iterations",
      class = "retrogroup_error"
    )
  }
  for (tolerance in list(-1, NA_real_, c(1, 2), "1")) {
    expect_error(
      group_optimize_unconstrained(sample, summary, tolerance = tolerance),
      "tolerance",
      class = "retrogroup_error"
    )
  }
  bad <- summary
  bad$patients[1] <- 0L
  expect_error(
    group_optimize_unconstrained(sample, bad),
    "summary is invalid",
    class = "retrogroup_error"
  )
  expect_error(
    group_optimize_unconstrained(as.list(sample), summary),
    "data frame",
    class = "retrogroup_error"
  )
})
