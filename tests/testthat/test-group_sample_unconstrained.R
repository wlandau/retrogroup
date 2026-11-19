test_that("group_sample_unconstrained() has the expected structure", {
  data <- probability_example()
  set.seed(0)
  x <- group_sample_unconstrained(data$patient, data$summary)
  probability <- group_analysis_probability(data$patient, data$summary)
  expect_s3_class(x, "tbl_df")
  expect_equal(nrow(x), nrow(data$patient))
  expect_identical(
    colnames(x),
    c(colnames(data$patient), names(probability))
  )
  expect_identical(
    x[, colnames(data$patient)],
    tibble::as_tibble(data$patient)
  )
  for (variable in names(probability)) {
    levels <- data$summary$level[data$summary$variable == variable]
    expect_s3_class(x[[variable]], "factor")
    expect_identical(levels(x[[variable]]), levels)
    expect_false(anyNA(x[[variable]]))
    expect_true(all(x[[variable]] %in% levels))
  }
})

test_that("group_sample_unconstrained() copies labels from summary", {
  data <- probability_small()
  summary <- data$summary[c(4, 2, 5, 1, 3), ]
  set.seed(0)
  x <- group_sample_unconstrained(data$patient, summary)
  expect_identical(colnames(x), c("treatment", "time", "event", "sex", "age"))
  expect_true(all(x$sex %in% c("male", "female", "other")))
  expect_true(all(x$age %in% c("young", "old")))
})

test_that("group_sample_unconstrained() is reproducible with set.seed()", {
  data <- probability_example()
  set.seed(1)
  x <- group_sample_unconstrained(data$patient, data$summary)
  set.seed(1)
  y <- group_sample_unconstrained(data$patient, data$summary)
  set.seed(2)
  z <- group_sample_unconstrained(data$patient, data$summary)
  expect_identical(x, y)
  expect_false(identical(x, z))
})

test_that("sampled frequencies agree with the posterior probabilities", {
  data <- probability_small()
  probability <- group_analysis_probability(data$patient, data$summary)
  reps <- 4000L
  patient <- data$patient[rep(seq_len(nrow(data$patient)), times = reps), ]
  set.seed(0)
  for (variable in names(probability)) {
    p <- probability[[variable]]
    index <- retrogroup:::sample_unconstrained_variable(
      p[rep(seq_len(nrow(p)), times = reps), , drop = FALSE]
    )
    patient_id <- rep(seq_len(nrow(p)), times = reps)
    frequency <- table(
      factor(patient_id, levels = seq_len(nrow(p))),
      factor(index, levels = seq_len(ncol(p)))
    ) /
      reps
    expect_equal(
      unname(matrix(frequency, nrow(p))),
      unname(p),
      tolerance = 0.05
    )
  }
})

test_that("sample_unconstrained_variable() handles extreme probabilities", {
  set.seed(0)
  epsilon <- 1e-7
  p <- matrix(c(1 - 2 * epsilon, epsilon, epsilon), 1000L, 3L, byrow = TRUE)
  index <- retrogroup:::sample_unconstrained_variable(p)
  expect_true(all(index %in% 1:3))
  expect_gt(mean(index == 1L), 0.99)
  p <- p[, 3:1]
  index <- retrogroup:::sample_unconstrained_variable(p)
  expect_true(all(index %in% 1:3))
  expect_gt(mean(index == 3L), 0.99)
})

test_that("group_sample_unconstrained() validates its arguments", {
  data <- probability_small()
  expect_error(
    group_sample_unconstrained(as.list(data$patient), data$summary),
    "data frame",
    class = "retrogroup_error"
  )
  summary <- data$summary
  summary$patients[1] <- 0L
  expect_error(
    group_sample_unconstrained(data$patient, summary),
    "summary is invalid",
    class = "retrogroup_error"
  )
  summary <- data$summary
  summary$variable[summary$variable == "age"] <- "time"
  expect_error(
    group_sample_unconstrained(data$patient, summary),
    "clashing names: \"time\"",
    class = "retrogroup_error"
  )
})
