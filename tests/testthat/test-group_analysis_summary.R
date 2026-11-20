analysis_colon <- function() {
  data <- survival::colon
  data <- data[data$etype == 2 & data$rx %in% c("Obs", "Lev+5FU"), ]
  tibble::tibble(
    treatment = as.integer(data$rx == "Lev+5FU"),
    time = as.numeric(data$time),
    event = as.integer(data$status),
    age = factor(
      ifelse(data$age < 60, "under 60", "60 or over"),
      levels = c("under 60", "60 or over")
    ),
    sex = ifelse(data$sex == 1, "male", "female")
  )
}

test_that("group_analysis_summary() has the expected format", {
  patient <- analysis_colon()
  x <- group_analysis_summary(patient)
  expect_s3_class(x, "tbl_df")
  expect_identical(
    colnames(x),
    c("variable", "level", "patients", "events", "hazard_ratio")
  )
  expect_type(x$variable, "character")
  expect_type(x$level, "character")
  expect_type(x$patients, "integer")
  expect_type(x$events, "integer")
  expect_type(x$hazard_ratio, "double")
  expect_identical(
    group_data_summary(
      x,
      variable = "variable",
      level = "level",
      patients = "patients",
      events = "events",
      hazard_ratio = "hazard_ratio"
    ),
    x
  )
})

test_that("group_analysis_summary() orders variables and levels", {
  patient <- analysis_colon()
  x <- group_analysis_summary(patient)
  expect_identical(x$variable, c("age", "age", "sex", "sex"))
  expect_identical(x$level, c("under 60", "60 or over", "female", "male"))
  x <- group_analysis_summary(patient, c("sex", "age"))
  expect_identical(x$variable, c("sex", "sex", "age", "age"))
  x <- group_analysis_summary(patient, "sex")
  expect_identical(x$variable, c("sex", "sex"))
  patient$age <- factor(patient$age, levels = c("unused", levels(patient$age)))
  x <- group_analysis_summary(patient, "age")
  expect_identical(x$level, c("under 60", "60 or over"))
})

test_that("group_analysis_summary() agrees with survival::coxph()", {
  patient <- analysis_colon()
  x <- group_analysis_summary(patient)
  for (row in seq_len(nrow(x))) {
    subset <- patient[patient[[x$variable[row]]] == x$level[row], ]
    fit <- survival::coxph(
      survival::Surv(time, event) ~ treatment,
      data = subset
    )
    expect_identical(x$patients[row], nrow(subset))
    expect_identical(x$events[row], sum(subset$event))
    expect_equal(x$hazard_ratio[row], exp(unname(stats::coef(fit))))
  }
})

test_that("group_analysis_summary() counts constrained samples exactly", {
  data <- probability_example()
  set.seed(0)
  sample <- group_sample_constrained(data$patient, data$summary)
  x <- group_analysis_summary(sample, c("age", "node4", "sex"))
  expect_identical(x$patients[order(x$variable, x$level)], {
    a <- data$summary
    a$patients[order(a$variable, a$level)]
  })
  for (row in seq_len(nrow(x))) {
    members <- sample[[x$variable[row]]] == x$level[row]
    expect_identical(x$events[row], sum(sample$event[members]))
  }
  totals <- tapply(x$events, x$variable, sum)
  expect_true(all(totals == sum(sample$event)))
})

test_that("group_analysis_summary() gives NA when the model cannot be fit", {
  patient <- tibble::tibble(
    treatment = c(0L, 1L, 0L, 1L, 1L, 1L),
    time = c(1, 2, 3, 4, 5, 6),
    event = c(1L, 1L, 1L, 0L, 1L, 1L),
    group = c("a", "a", "a", "a", "b", "b")
  )
  x <- group_analysis_summary(patient)
  expect_identical(x$level, c("a", "b"))
  expect_identical(x$patients, c(4L, 2L))
  expect_identical(x$events, c(3L, 2L))
  expect_true(is.finite(x$hazard_ratio[1L]))
  expect_true(is.na(x$hazard_ratio[2L]))
})

test_that("group_analysis_summary() counts zero events in a level", {
  patient <- tibble::tibble(
    treatment = c(0L, 1L, 0L, 1L, 0L, 1L),
    time = c(1, 2, 3, 4, 5, 6),
    event = c(1L, 1L, 1L, 1L, 0L, 0L),
    group = factor(c("a", "a", "a", "a", "b", "b"), levels = c("a", "b"))
  )
  x <- group_analysis_summary(patient)
  expect_identical(x$level, c("a", "b"))
  expect_identical(x$patients, c(4L, 2L))
  expect_identical(x$events, c(4L, 0L))
  expect_true(is.finite(x$hazard_ratio[1L]))
  expect_true(is.na(x$hazard_ratio[2L]))
})

test_that("analysis_fit() agrees with survival::coxph()", {
  data <- probability_example()$patient
  fit <- survival::coxph(survival::Surv(time, event) ~ treatment, data)
  expect_equal(
    retrogroup:::analysis_fit(data$treatment, data$time, data$event),
    unname(stats::coef(fit))
  )
})

test_that("analysis_fit() is NA when the model cannot be fit", {
  expect_identical(
    retrogroup:::analysis_fit(integer(0L), numeric(0L), integer(0L)),
    NA_real_
  )
  expect_identical(
    retrogroup:::analysis_fit(c(0L, 0L, 1L), c(1, 2, 3), c(1L, 1L, 0L)),
    NA_real_
  )
  expect_identical(
    retrogroup:::analysis_fit(c(1L, 1L), c(1, 2), c(1L, 1L)),
    NA_real_
  )
})

test_that("analysis_fit() is NA when survival::coxph.fit() errors", {
  expect_identical(
    retrogroup:::analysis_fit(c(0L, 1L), c("a", "b"), c(1L, 1L)),
    NA_real_
  )
})

test_that("analysis_fit() is NA when the estimates are not finite", {
  expect_identical(
    retrogroup:::analysis_fit(c(0, 1, NaN), c(1, 2, 3), c(1L, 1L, 1L)),
    NA_real_
  )
})

test_that("group_analysis_summary() validates its arguments", {
  patient <- analysis_colon()
  expect_error(
    group_analysis_summary(as.list(patient), "age"),
    "data frame",
    class = "retrogroup_error"
  )
  expect_error(
    group_analysis_summary(patient[, c("treatment", "time", "event")]),
    "variables must be",
    class = "retrogroup_error"
  )
  for (variables in list(character(0L), NA_character_, "", 1L)) {
    expect_error(
      group_analysis_summary(patient, variables),
      "variables must be",
      class = "retrogroup_error"
    )
  }
  expect_error(
    group_analysis_summary(patient, c("age", "age")),
    "unique",
    class = "retrogroup_error"
  )
  expect_error(
    group_analysis_summary(patient, c("age", "time")),
    "must not include",
    class = "retrogroup_error"
  )
  expect_error(
    group_analysis_summary(patient, "nonexistent"),
    "not in patient",
    class = "retrogroup_error"
  )
  bad <- patient
  bad$age[1] <- NA
  expect_error(
    group_analysis_summary(bad, "age"),
    "missing values",
    class = "retrogroup_error"
  )
  bad <- patient
  bad$number <- seq_len(nrow(bad))
  expect_error(
    group_analysis_summary(bad, "number"),
    "character or factor",
    class = "retrogroup_error"
  )
})
