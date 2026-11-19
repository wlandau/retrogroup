test_that("group_optimize_constrained() has the expected structure", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_constrained(data$patient, data$summary)
  sample$extra <- seq_len(nrow(sample))
  x <- suppressWarnings(group_optimize_constrained(
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

test_that("group_optimize_constrained() preserves the published counts", {
  data <- probability_example()
  for (seed in seq_len(20L)) {
    set.seed(seed)
    sample <- group_sample_constrained(data$patient, data$summary)
    x <- suppressWarnings(group_optimize_constrained(
      sample,
      data$summary,
      iterations = 100L
    ))
    for (variable in unique(data$summary$variable)) {
      rows <- data$summary$variable == variable
      expect_equal(
        as.integer(table(x[[variable]])),
        data$summary$patients[rows],
        info = paste("seed", seed, "variable", variable)
      )
    }
  }
})

test_that("group_optimize_constrained() reduces the loss", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_constrained(data$patient, data$summary)
  x <- suppressWarnings(group_optimize_constrained(
    sample,
    data$summary,
    iterations = 200L
  ))
  for (variable in unique(data$summary$variable)) {
    before <- optimize_loss_reference(sample, data$summary, variable)
    after <- optimize_loss_reference(x, data$summary, variable)
    expect_lt(after, before)
  }
})

test_that("group_optimize_constrained() reduces the loss without SEs", {
  data <- probability_example()
  summary <- data$summary
  summary$standard_error <- NA_real_
  set.seed(0)
  sample <- group_sample_constrained(data$patient, summary)
  x <- suppressWarnings(group_optimize_constrained(
    sample,
    summary,
    iterations = 200L
  ))
  for (variable in unique(summary$variable)) {
    before <- optimize_loss_reference(sample, summary, variable)
    after <- optimize_loss_reference(x, summary, variable)
    expect_lt(after, before)
  }
})

test_that("group_optimize_constrained() stops early", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_constrained(data$patient, data$summary)
  expect_identical(
    suppressWarnings(
      group_optimize_constrained(sample, data$summary, iterations = 0L)
    ),
    sample
  )
  expect_identical(
    group_optimize_constrained(sample, data$summary, tolerance = Inf),
    sample
  )
})

test_that("group_optimize_constrained() is reproducible with set.seed()", {
  data <- probability_example()
  summary <- data$summary[data$summary$variable == "age", ]
  set.seed(0)
  sample <- group_sample_constrained(data$patient, summary)
  set.seed(1)
  x <- suppressWarnings(group_optimize_constrained(
    sample,
    summary,
    iterations = 20L
  ))
  set.seed(1)
  y <- suppressWarnings(group_optimize_constrained(
    sample,
    summary,
    iterations = 20L
  ))
  expect_identical(x, y)
})

test_that("group_optimize_constrained() handles factor levels in patient", {
  data <- probability_example()
  summary <- data$summary[data$summary$variable == "age", ]
  set.seed(0)
  sample <- group_sample_constrained(data$patient, summary)
  character <- sample
  character$age <- as.character(character$age)
  set.seed(1)
  x <- suppressWarnings(group_optimize_constrained(
    sample,
    summary,
    iterations = 5L
  ))
  set.seed(1)
  y <- suppressWarnings(group_optimize_constrained(
    character,
    summary,
    iterations = 5L
  ))
  expect_identical(x, y)
})

test_that("optimize_constrained_variable() never increases the loss", {
  data <- probability_small()
  patient <- data$patient
  set.seed(0)
  index <- sample(rep(1:3, times = c(1L, 2L, 3L)))
  beta <- c(-0.5, 0.3, 0)
  sigma <- c(NA_real_, NA_real_, NA_real_)
  args <- list(z = patient$treatment, t = patient$time, d = patient$event)
  loss <- function(index) {
    sum(do.call(
      retrogroup:::optimize_loss_terms,
      c(list(x = index), args, list(beta = beta, sigma = sigma))
    ))
  }
  before <- loss(index)
  for (iterations in c(1L, 5L, 20L)) {
    result <- do.call(
      retrogroup:::optimize_constrained_variable,
      c(
        list(x = index),
        args,
        list(beta = beta, sigma = sigma),
        list(iterations = iterations, tolerance = 0)
      )
    )
    after <- loss(result$x)
    expect_identical(result$loss, after)
    expect_false(isTRUE(after > before))
    expect_identical(tabulate(result$x, 3L), tabulate(index, 3L))
  }
})

test_that("optimize_constrained_variable() only swaps patients", {
  data <- probability_example()
  patient <- data$patient
  summary <- data$summary[data$summary$variable == "sex", ]
  set.seed(0)
  x <- sample(rep(1:2, times = summary$patients))
  y <- x
  for (iteration in seq_len(10L)) {
    result <- retrogroup:::optimize_constrained_variable(
      x = y,
      z = patient$treatment,
      t = patient$time,
      d = patient$event,
      beta = summary$log_hazard_ratio,
      sigma = summary$standard_error,
      iterations = 1L,
      tolerance = 0
    )
    expect_true(sum(result$x != y) %in% c(0L, 2L))
    y <- result$x
  }
  expect_identical(tabulate(y, 2L), tabulate(x, 2L))
})

test_that("optimize_constrained_variable() handles a single level", {
  data <- probability_small()
  patient <- data$patient
  x <- rep(1L, nrow(patient))
  result <- retrogroup:::optimize_constrained_variable(
    x = x,
    z = patient$treatment,
    t = patient$time,
    d = patient$event,
    beta = 0,
    sigma = NA_real_,
    iterations = 10L,
    tolerance = 0
  )
  expect_identical(result$x, x)
  expect_true(is.finite(result$loss))
})

test_that("group_optimize_constrained() warns if tolerance is not reached", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_constrained(data$patient, data$summary)
  expect_warning(
    group_optimize_constrained(sample, data$summary, iterations = 2L),
    "did not fall below tolerance = 1e-06 within iterations = 2 for 3",
    class = "retrogroup_warning"
  )
  warning <- tryCatch(
    group_optimize_constrained(sample, data$summary, iterations = 2L),
    retrogroup_warning = \(warning) warning
  )
  message <- conditionMessage(warning)
  expect_match(message, "\"age\" = ", fixed = TRUE)
  expect_match(message, "\"node4\" = ", fixed = TRUE)
  expect_match(message, "\"sex\" = ", fixed = TRUE)
  expect_match(message, "Increase iterations", fixed = TRUE)
  summary <- data$summary[data$summary$variable == "age", ]
  expect_warning(
    group_optimize_constrained(sample, summary, iterations = 0L),
    "for 1 subgroup variable",
    class = "retrogroup_warning"
  )
})

test_that("group_optimize_constrained() does not warn if tolerance is reached", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_constrained(data$patient, data$summary)
  expect_no_warning(
    group_optimize_constrained(sample, data$summary, tolerance = Inf)
  )
  expect_no_warning(
    group_optimize_constrained(
      sample,
      data$summary,
      iterations = 0L,
      tolerance = Inf
    )
  )
})

test_that("group_optimize_constrained() validates its arguments", {
  data <- probability_example()
  summary <- data$summary[data$summary$variable == "age", ]
  set.seed(0)
  sample <- group_sample_constrained(data$patient, summary)
  expect_error(
    group_optimize_constrained(data$patient, summary),
    "Missing subgroup column",
    class = "retrogroup_error"
  )
  bad <- sample
  bad$age <- as.character(bad$age)
  bad$age[1] <- "nonexistent"
  expect_error(
    group_optimize_constrained(bad, summary),
    "only levels",
    class = "retrogroup_error"
  )
  bad$age[1] <- NA_character_
  expect_error(
    group_optimize_constrained(bad, summary),
    "only levels",
    class = "retrogroup_error"
  )
  bad <- sample
  levels <- summary$level
  bad$age[bad$age == levels[1L]][1L] <- levels[2L]
  expect_error(
    group_optimize_constrained(bad, summary),
    "Mismatched subgroup variable\\(s\\): \"age\"",
    class = "retrogroup_error"
  )
  for (iterations in list(-1, 1.5, NA, c(1, 2), "1", Inf)) {
    expect_error(
      group_optimize_constrained(sample, summary, iterations = iterations),
      "iterations",
      class = "retrogroup_error"
    )
  }
  for (tolerance in list(-1, NA_real_, c(1, 2), "1")) {
    expect_error(
      group_optimize_constrained(sample, summary, tolerance = tolerance),
      "tolerance",
      class = "retrogroup_error"
    )
  }
  bad <- summary
  bad$patients[1] <- 0L
  expect_error(
    group_optimize_constrained(sample, bad),
    "summary is invalid",
    class = "retrogroup_error"
  )
  expect_error(
    group_optimize_constrained(as.list(sample), summary),
    "data frame",
    class = "retrogroup_error"
  )
})

test_that("assert_optimize_counts() lists only mismatched variables", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_constrained(data$patient, data$summary)
  expect_null(retrogroup:::assert_optimize_counts(sample, data$summary))
  bad <- sample
  levels <- levels(bad$sex)
  bad$sex[bad$sex == levels[1L]][1L] <- levels[2L]
  error <- tryCatch(
    retrogroup:::assert_optimize_counts(bad, data$summary),
    retrogroup_error = \(error) error
  )
  message <- conditionMessage(error)
  expect_match(message, "\"sex\"", fixed = TRUE)
  expect_false(grepl("\"age\"", message, fixed = TRUE))
  expect_false(grepl("\"node4\"", message, fixed = TRUE))
})
