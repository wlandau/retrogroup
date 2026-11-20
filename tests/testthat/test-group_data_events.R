events_example <- function() {
  patient <- group_data_patient(
    group_example_patient(),
    arm = "arm",
    treatment = "Lev+5FU",
    time = "time",
    event = "event"
  )
  data <- group_example_summary()
  summary <- group_data_summary(
    data,
    variable = "variable",
    level = "level",
    patients = "patients",
    events = NULL,
    hazard_ratio = "hazard_ratio"
  )
  list(patient = patient, data = data, summary = summary)
}

test_that("events_approximate() matches the formula", {
  hazard_ratio <- c(0.8, 1.2, 0.5)
  bound <- c(1.1, 1.6, 0.3)
  proportion <- 0.4
  c <- qnorm(0.975)
  expected <- c^2 /
    (proportion * (1 - proportion) * (log(bound) - log(hazard_ratio))^2)
  expect_equal(
    retrogroup:::events_approximate(hazard_ratio, bound, proportion, 0.95),
    expected
  )
  expect_equal(
    retrogroup:::events_approximate(1, exp(c / 10), 0.5, 0.95),
    400
  )
})

test_that("events_approximate() gives the same result for either bound", {
  hazard_ratio <- c(0.8, 1.2)
  width <- c(0.3, 0.1)
  lower <- exp(log(hazard_ratio) - width)
  upper <- exp(log(hazard_ratio) + width)
  expect_equal(
    retrogroup:::events_approximate(hazard_ratio, lower, 0.5, 0.95),
    retrogroup:::events_approximate(hazard_ratio, upper, 0.5, 0.95)
  )
})

test_that("events_approximate() scales with the squared normal quantile", {
  ratio <- retrogroup:::events_approximate(0.8, 1.1, 0.5, 0.9) /
    retrogroup:::events_approximate(0.8, 1.1, 0.5, 0.95)
  expect_equal(ratio, (qnorm(0.95) / qnorm(0.975))^2)
})

test_that("events_approximate() is close to the true counts of the example", {
  data <- events_example()
  estimate <- retrogroup:::events_approximate(
    data$data$hazard_ratio,
    data$data$upper,
    mean(data$patient$treatment),
    0.95
  )
  expect_true(all(abs(estimate - data$data$events) / data$data$events < 0.1))
})

test_that("events_adjust() keeps counts that already add up", {
  expect_identical(
    retrogroup:::events_adjust(c(10.2, 4.6, 2.9), n = c(20L, 9L, 9L), 18L),
    c(10L, 5L, 3L)
  )
})

test_that("events_adjust() subtracts one per level in order", {
  expect_identical(
    retrogroup:::events_adjust(c(10, 5, 3), n = c(20L, 9L, 9L), 16L),
    c(9L, 4L, 3L)
  )
  expect_identical(
    retrogroup:::events_adjust(c(10, 5, 3), n = c(20L, 9L, 9L), 14L),
    c(8L, 4L, 2L)
  )
})

test_that("events_adjust() adds one per level in order", {
  expect_identical(
    retrogroup:::events_adjust(c(10, 5, 3), n = c(20L, 9L, 9L), 20L),
    c(11L, 6L, 3L)
  )
  expect_identical(
    retrogroup:::events_adjust(c(10, 5, 3), n = c(20L, 9L, 9L), 22L),
    c(12L, 6L, 4L)
  )
})

test_that("events_adjust() makes several passes and skips full levels", {
  expect_identical(
    retrogroup:::events_adjust(c(1, 1, 1), n = c(10L, 2L, 10L), 12L),
    c(5L, 2L, 5L)
  )
  expect_identical(
    retrogroup:::events_adjust(c(9, 2, 9), n = c(10L, 2L, 10L), 4L),
    c(1L, 1L, 2L)
  )
})

test_that("events_adjust() clamps counts above the number of patients", {
  expect_identical(
    retrogroup:::events_adjust(c(50, 3), n = c(10L, 10L), 13L),
    c(10L, 3L)
  )
  expect_identical(
    retrogroup:::events_adjust(c(50, 50), n = c(10L, 10L), 20L),
    c(10L, 10L)
  )
  expect_identical(
    retrogroup:::events_adjust(c(50, 3), n = c(10L, 10L), 15L),
    c(10L, 5L)
  )
})

test_that("events_adjust() raises counts below 1 to 1", {
  expect_identical(
    retrogroup:::events_adjust(c(0.2, 1e-8, 5), n = c(10L, 10L, 10L), 7L),
    c(1L, 1L, 5L)
  )
  expect_identical(
    retrogroup:::events_adjust(c(0.2, 1e-8, 5), n = c(10L, 10L, 10L), 3L),
    c(1L, 1L, 1L)
  )
})

test_that("events_adjust() always satisfies the constraints", {
  set.seed(0)
  for (rep in seq_len(200L)) {
    K <- sample(2:5, 1L)
    n <- sample(1:30, K, replace = TRUE)
    total <- sample(seq(K, sum(n)), 1L)
    estimate <- exp(rnorm(K, mean = log(total / K), sd = 1.5))
    events <- retrogroup:::events_adjust(estimate, n, total)
    expect_type(events, "integer")
    expect_equal(sum(events), total)
    expect_true(all(events >= 1L & events <= n))
  }
})

test_that("events_adjust() rejects impossible totals", {
  expect_error(
    retrogroup:::events_adjust(c(1, 1, 1), n = c(5L, 5L, 5L), 2L, "age"),
    "\"age\" has 3 levels and 15 patients, but there are 2 events",
    class = "retrogroup_error"
  )
  expect_error(
    retrogroup:::events_adjust(c(1, 1), n = c(2L, 2L), 5L),
    "at least the number of levels",
    class = "retrogroup_error"
  )
})

test_that("group_data_events() returns the structure of group_data_summary()", {
  data <- events_example()
  out <- group_data_events(data$patient, data$summary, data$data$upper)
  expect_s3_class(out, "tbl_df")
  expect_identical(
    colnames(out),
    c("variable", "level", "patients", "events", "hazard_ratio")
  )
  expect_type(out$events, "integer")
  expect_identical(out[, -4L], data$summary[, -4L])
  expect_identical(
    group_data_summary(
      out,
      variable = "variable",
      level = "level",
      patients = "patients",
      events = "events",
      hazard_ratio = "hazard_ratio"
    ),
    out
  )
  expect_no_error(retrogroup:::assert_probability_summary(out))
  expect_no_error(retrogroup:::assert_probability_joint(data$patient, out))
  set.seed(0)
  sample <- group_sample_constrained(data$patient, out)
  expect_equal(
    as.integer(tapply(sample$event, sample$age, sum)),
    out$events[out$variable == "age"]
  )
})

test_that("group_data_events() reconciles each subgroup variable", {
  data <- events_example()
  out <- group_data_events(data$patient, data$summary, data$data$lower)
  for (variable in unique(out$variable)) {
    rows <- out$variable == variable
    expect_equal(sum(out$events[rows]), sum(data$patient$event))
  }
  expect_true(all(out$events >= 1L & out$events <= out$patients))
  expect_true(all(abs(out$events - data$data$events) <= 10L))
})

test_that("group_data_events() is deterministic and uses either bound", {
  data <- events_example()
  upper <- group_data_events(data$patient, data$summary, data$data$upper)
  expect_identical(
    group_data_events(data$patient, data$summary, data$data$upper),
    upper
  )
  even <- seq_len(nrow(data$data)) %% 2L == 0L
  mixed <- ifelse(even, data$data$lower, data$data$upper)
  expect_type(
    group_data_events(data$patient, data$summary, mixed)$events,
    "integer"
  )
})

test_that("group_data_events() honors confidence", {
  data <- events_example()
  bound <- data$data$upper
  expect_false(identical(
    group_data_events(data$patient, data$summary, bound, confidence = 0.95),
    group_data_events(data$patient, data$summary, bound, confidence = 0.9)
  ))
})

test_that("group_data_events() ignores the events in summary", {
  data <- events_example()
  expected <- group_data_events(data$patient, data$summary, data$data$upper)
  valid <- group_data_summary(
    data$data,
    variable = "variable",
    level = "level",
    patients = "patients",
    events = "events",
    hazard_ratio = "hazard_ratio"
  )
  expect_identical(
    group_data_events(data$patient, valid, data$data$upper),
    expected
  )
  invalid <- data$summary
  invalid$events <- c(0L, 1000L, 3L, -1L, 2.5, NA)
  expect_identical(
    group_data_events(data$patient, invalid, data$data$upper),
    expected
  )
  missing <- data$summary[, setdiff(colnames(data$summary), "events")]
  expect_identical(
    group_data_events(data$patient, missing, data$data$upper),
    expected
  )
})

test_that("group_data_events() validates bound", {
  data <- events_example()
  bound <- data$data$upper
  for (bad in list(bound[-1L], as.character(bound), NULL)) {
    expect_error(
      group_data_events(data$patient, data$summary, bad),
      "one element per row of summary",
      class = "retrogroup_error"
    )
  }
  for (value in c(NA, Inf)) {
    bad <- bound
    bad[2L] <- value
    expect_error(
      group_data_events(data$patient, data$summary, bad),
      "finite",
      class = "retrogroup_error"
    )
  }
  for (value in c(0, -0.5)) {
    bad <- bound
    bad[2L] <- value
    expect_error(
      group_data_events(data$patient, data$summary, bad),
      "strictly positive",
      class = "retrogroup_error"
    )
  }
  bad <- bound
  bad[c(2L, 5L)] <- data$summary$hazard_ratio[c(2L, 5L)]
  expect_error(
    group_data_events(data$patient, data$summary, bad),
    "age = \"60 or over\", sex = \"female\"",
    class = "retrogroup_error"
  )
})

test_that("group_data_events() validates confidence", {
  data <- events_example()
  for (confidence in list(0, 1, -0.5, 1.5, NA_real_, c(0.9, 0.95), "0.95")) {
    expect_error(
      group_data_events(
        data$patient,
        data$summary,
        data$data$upper,
        confidence = confidence
      ),
      "confidence",
      class = "retrogroup_error"
    )
  }
})

test_that("group_data_events() validates patient and summary", {
  data <- events_example()
  expect_error(
    group_data_events(group_example_patient(), data$summary, data$data$upper),
    class = "retrogroup_error"
  )
  expect_error(
    group_data_events(data$patient[-1L, ], data$summary, data$data$upper),
    "must equal the number of rows of patient",
    class = "retrogroup_error"
  )
  bad <- data$summary
  bad$hazard_ratio[1L] <- -1
  expect_error(
    group_data_events(data$patient, bad, data$data$upper),
    "summary is invalid",
    class = "retrogroup_error"
  )
  expect_error(
    group_data_events(data$patient, as.list(data$summary), data$data$upper),
    "data frame",
    class = "retrogroup_error"
  )
})

test_that("downstream functions reject a summary without events", {
  data <- events_example()
  expect_error(
    group_sample_constrained(data$patient, data$summary),
    "summary is invalid",
    class = "retrogroup_error"
  )
  expect_error(
    group_analysis_probability(data$patient, data$summary),
    "summary is invalid",
    class = "retrogroup_error"
  )
})
