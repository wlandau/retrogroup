summary_data <- function() {
  data.frame(
    variable = c("age", "age", "sex", "sex"),
    level = c("young", "old", "female", "male"),
    patients = c(3L, 7L, 4L, 6L),
    log_hazard_ratio = c(-0.2, 0.1, 0, 0.3),
    standard_error = c(0.5, 0.4, 0.45, 0.35)
  )
}

summary_call <- function(data = summary_data(), ...) {
  args <- list(
    variable = "variable",
    level = "level",
    patients = "patients",
    log_hazard_ratio = "log_hazard_ratio",
    standard_error = "standard_error"
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
    log_hazard_ratio = "log_hazard_ratio",
    standard_error = "standard_error"
  )
  expect_s3_class(x, "tbl_df")
  expect_identical(x, data)
})

test_that("group_data_summary() has the expected structure", {
  x <- summary_call()
  expect_equal(
    colnames(x),
    c("variable", "level", "patients", "log_hazard_ratio", "standard_error")
  )
  expect_type(x$variable, "character")
  expect_type(x$level, "character")
  expect_type(x$patients, "integer")
  expect_type(x$log_hazard_ratio, "double")
  expect_type(x$standard_error, "double")
  expect_false(anyNA(x))
  expect_equal(x$level, summary_data()$level)
})

test_that("group_data_summary() allows missing standard errors", {
  x <- summary_call(standard_error = NULL)
  expect_equal(
    colnames(x),
    c("variable", "level", "patients", "log_hazard_ratio", "standard_error")
  )
  expect_type(x$standard_error, "double")
  expect_true(all(is.na(x$standard_error)))
  expect_equal(x$log_hazard_ratio, summary_data()$log_hazard_ratio)
})

test_that("group_data_summary() accepts flexible input types", {
  data <- summary_data()
  data$patients <- as.numeric(data$patients)
  data$level <- factor(data$level)
  data$log_hazard_ratio <- as.integer(c(0, 1, 0, 1))
  data$extra <- NA
  x <- summary_call(tibble::as_tibble(data))
  expect_type(x$patients, "integer")
  expect_equal(x$patients, c(3L, 7L, 4L, 6L))
  expect_type(x$level, "character")
  expect_equal(x$level, c("young", "old", "female", "male"))
  expect_type(x$log_hazard_ratio, "double")
})

test_that("group_data_summary() accepts custom column names", {
  data <- summary_data()
  colnames(data) <- c("v", "l", "n", "b", "s")
  x <- group_data_summary(data, "v", "l", "n", "b", "s")
  expect_equal(
    colnames(x),
    c("variable", "level", "patients", "log_hazard_ratio", "standard_error")
  )
  expect_equal(x$patients, c(3L, 7L, 4L, 6L))
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
    summary_call(standard_error = NA_character_),
    "standard_error must be a single",
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

test_that("group_data_summary() validates log_hazard_ratio", {
  data <- summary_data()
  data$log_hazard_ratio <- as.character(data$log_hazard_ratio)
  expect_error(summary_call(data), "numeric", class = "retrogroup_error")
  data <- summary_data()
  data$log_hazard_ratio[1] <- NA
  expect_error(summary_call(data), "finite", class = "retrogroup_error")
  data$log_hazard_ratio[1] <- -Inf
  expect_error(summary_call(data), "finite", class = "retrogroup_error")
})

test_that("group_data_summary() validates standard_error", {
  data <- summary_data()
  data$standard_error <- as.character(data$standard_error)
  expect_error(summary_call(data), "numeric", class = "retrogroup_error")
  data <- summary_data()
  data$standard_error[1] <- NA
  expect_error(summary_call(data), "finite", class = "retrogroup_error")
  data$standard_error[1] <- Inf
  expect_error(summary_call(data), "finite", class = "retrogroup_error")
  data$standard_error[1] <- 0
  expect_error(
    summary_call(data),
    "strictly positive",
    class = "retrogroup_error"
  )
})
