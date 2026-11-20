summary_data <- function() {
  data.frame(
    variable = c("age", "age", "sex", "sex"),
    level = c("young", "old", "female", "male"),
    patients = c(3L, 7L, 4L, 6L),
    events = c(2L, 4L, 3L, 3L),
    hazard_ratio = c(0.8, 1.1, 1, 1.3)
  )
}

summary_call <- function(data = summary_data(), ...) {
  args <- list(
    variable = "variable",
    level = "level",
    patients = "patients",
    events = "events",
    hazard_ratio = "hazard_ratio"
  )
  args[names(list(...))] <- list(...)
  do.call(group_data_summary, c(list(data), args))
}

test_that("group_data_summary() reproduces the example data", {
  data <- group_example_summary()
  x <- group_data_summary(
    data,
    variable = "variable",
    level = "level",
    patients = "patients",
    events = "events",
    hazard_ratio = "hazard_ratio"
  )
  expect_s3_class(x, "tbl_df")
  expect_identical(x, data)
})

test_that("group_data_summary() has the expected structure", {
  x <- summary_call()
  expect_equal(
    colnames(x),
    c("variable", "level", "patients", "events", "hazard_ratio")
  )
  expect_type(x$variable, "character")
  expect_type(x$level, "character")
  expect_type(x$patients, "integer")
  expect_type(x$events, "integer")
  expect_type(x$hazard_ratio, "double")
  expect_false(anyNA(x))
  expect_equal(x$level, summary_data()$level)
})

test_that("group_data_summary() accepts flexible input types", {
  data <- summary_data()
  data$patients <- as.numeric(data$patients)
  data$level <- factor(data$level)
  data$events <- as.numeric(data$events)
  data$hazard_ratio <- as.integer(c(1, 2, 1, 2))
  data$extra <- NA
  x <- summary_call(tibble::as_tibble(data))
  expect_type(x$patients, "integer")
  expect_equal(x$patients, c(3L, 7L, 4L, 6L))
  expect_type(x$level, "character")
  expect_equal(x$level, c("young", "old", "female", "male"))
  expect_type(x$events, "integer")
  expect_type(x$hazard_ratio, "double")
})

test_that("group_data_summary() accepts custom column names", {
  data <- summary_data()
  colnames(data) <- c("v", "l", "n", "e", "h")
  x <- group_data_summary(data, "v", "l", "n", "e", "h")
  expect_equal(
    colnames(x),
    c("variable", "level", "patients", "events", "hazard_ratio")
  )
  expect_equal(x$patients, c(3L, 7L, 4L, 6L))
  expect_equal(x$events, c(2L, 4L, 3L, 3L))
  expect_equal(x$hazard_ratio, c(0.8, 1.1, 1, 1.3))
})

test_that("group_data_summary() validates data", {
  expect_error(
    summary_call(as.list(summary_data())),
    "data frame",
    class = "retrogroup_error"
  )
  expect_error(
    summary_call(summary_data()[1, ]),
    "at least two rows",
    class = "retrogroup_error"
  )
})

test_that("group_data_summary() validates column name arguments", {
  expect_error(
    summary_call(variable = 1),
    "variable must be a single",
    class = "retrogroup_error"
  )
  expect_error(
    summary_call(hazard_ratio = NA_character_),
    "hazard_ratio must be a single",
    class = "retrogroup_error"
  )
  expect_error(
    summary_call(level = "variable"),
    "different columns",
    class = "retrogroup_error"
  )
  expect_error(
    summary_call(patients = "nope"),
    "not in data",
    class = "retrogroup_error"
  )
  expect_error(
    summary_call(events = NULL),
    "events must be a single",
    class = "retrogroup_error"
  )
})

test_that("group_data_summary() validates variable", {
  data <- summary_data()
  data$variable <- factor(data$variable)
  expect_error(summary_call(data), "character", class = "retrogroup_error")
  data <- summary_data()
  data$variable[1] <- NA
  expect_error(summary_call(data), "missing", class = "retrogroup_error")
  data$variable[1] <- ""
  expect_error(summary_call(data), "empty", class = "retrogroup_error")
})

test_that("group_data_summary() validates level", {
  data <- summary_data()
  data$level <- 1:4
  expect_error(
    summary_call(data),
    "character or factor",
    class = "retrogroup_error"
  )
  data <- summary_data()
  data$level[1] <- NA
  expect_error(summary_call(data), "missing", class = "retrogroup_error")
  data$level[1] <- ""
  expect_error(summary_call(data), "empty", class = "retrogroup_error")
  data <- summary_data()
  data$variable[2] <- "other"
  data$patients[2] <- 3L
  data$events[2] <- 2L
  expect_error(
    summary_call(data),
    "more than one level.*\"other\"",
    class = "retrogroup_error"
  )
  data <- summary_data()
  data$level[2] <- "young"
  expect_error(
    summary_call(data),
    "unique within.*\"age\"",
    class = "retrogroup_error"
  )
})

test_that("group_data_summary() validates patients", {
  data <- summary_data()
  data$patients <- as.character(data$patients)
  expect_error(summary_call(data), "numeric", class = "retrogroup_error")
  data <- summary_data()
  data$patients[1] <- NA
  expect_error(summary_call(data), "finite", class = "retrogroup_error")
  data$patients <- c(3.5, 6.5, 4, 6)
  expect_error(
    summary_call(data),
    "whole numbers",
    class = "retrogroup_error"
  )
  data$patients <- c(0L, 10L, 4L, 6L)
  data$events <- c(0L, 6L, 3L, 3L)
  expect_error(
    summary_call(data),
    "strictly positive",
    class = "retrogroup_error"
  )
  data$patients <- c(1, 2^31, 1, 2^31)
  expect_error(summary_call(data), "at most", class = "retrogroup_error")
})

test_that("group_data_summary() requires equal patient totals", {
  data <- summary_data()
  data$patients[4] <- 7L
  error <- tryCatch(summary_call(data), error = identity)
  expect_s3_class(error, "retrogroup_error")
  message <- conditionMessage(error)
  expect_match(message, "same for every subgroup variable")
  expect_match(message, "age = 10, sex = 11")
  expect_match(message, "disclose it transparently")
})

test_that("group_data_summary() validates events", {
  data <- summary_data()
  data$events <- as.character(data$events)
  expect_error(summary_call(data), "numeric", class = "retrogroup_error")
  data <- summary_data()
  data$events[1] <- NA
  expect_error(summary_call(data), "finite", class = "retrogroup_error")
  data$events[1] <- Inf
  expect_error(summary_call(data), "finite", class = "retrogroup_error")
  data$events <- c(2.5, 3.5, 3, 3)
  expect_error(
    summary_call(data),
    "whole numbers",
    class = "retrogroup_error"
  )
})

test_that("group_data_summary() requires 1 <= events <= patients", {
  data <- summary_data()
  data$events <- c(0L, 6L, 3L, 3L)
  error <- tryCatch(summary_call(data), error = identity)
  expect_s3_class(error, "retrogroup_error")
  message <- conditionMessage(error)
  expect_match(message, "at least 1 and at most the number of patients")
  expect_match(message, "age = \"young\" (0 events, 3 patients)", fixed = TRUE)
  expect_match(message, "no events cannot support a hazard ratio")
  data <- summary_data()
  data$events <- c(4L, 2L, 5L, 1L)
  error <- tryCatch(summary_call(data), error = identity)
  expect_s3_class(error, "retrogroup_error")
  message <- conditionMessage(error)
  expect_match(message, "age = \"young\" (4 events, 3 patients)", fixed = TRUE)
  expect_match(message, "sex = \"female\" (5 events, 4 patients)", fixed = TRUE)
  data <- summary_data()
  data$events <- c(3L, 7L, 4L, 6L)
  expect_identical(summary_call(data)$events, data$events)
})

test_that("group_data_summary() requires equal event totals", {
  data <- summary_data()
  data$events[4] <- 4L
  error <- tryCatch(summary_call(data), error = identity)
  expect_s3_class(error, "retrogroup_error")
  message <- conditionMessage(error)
  expect_match(message, "total number of events must be the same")
  expect_match(message, "age = 6, sex = 7")
  expect_match(message, "event counts define the reconstructed subgroups")
})

test_that("group_data_summary() validates hazard_ratio", {
  data <- summary_data()
  data$hazard_ratio <- as.character(data$hazard_ratio)
  expect_error(summary_call(data), "numeric", class = "retrogroup_error")
  data <- summary_data()
  data$hazard_ratio[1] <- NA
  expect_error(summary_call(data), "finite", class = "retrogroup_error")
  data$hazard_ratio[1] <- Inf
  expect_error(summary_call(data), "finite", class = "retrogroup_error")
  for (value in c(0, -0.2)) {
    data$hazard_ratio[1] <- value
    expect_error(
      summary_call(data),
      "strictly positive",
      class = "retrogroup_error"
    )
  }
  error <- tryCatch(summary_call(data), error = identity)
  expect_match(conditionMessage(error), "log hazard ratios")
})
