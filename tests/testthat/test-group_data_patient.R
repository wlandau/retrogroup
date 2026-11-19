patient_data <- function() {
  data.frame(
    arm = c("a", "b", "a", "b"),
    time = c(1, 2.5, 3, 0),
    event = c(1L, 0L, 1L, 0L)
  )
}

patient_call <- function(data = patient_data(), ...) {
  args <- list(arm = "arm", treatment = "a", time = "time", event = "event")
  args[names(list(...))] <- list(...)
  do.call(group_data_patient, c(list(data), args))
}

test_that("group_data_patient() has the expected structure", {
  data <- group_example_patient()
  x <- group_data_patient(
    data,
    arm = "arm",
    treatment = "Lev+5FU",
    time = "time",
    event = "event"
  )
  expect_s3_class(x, "tbl_df")
  expect_equal(colnames(x), c("treatment", "time", "event"))
  expect_type(x$treatment, "integer")
  expect_type(x$time, "double")
  expect_type(x$event, "integer")
  expect_equal(nrow(x), 619L)
  expect_false(anyNA(x))
})

test_that("group_data_patient() has the expected contents", {
  data <- group_example_patient()
  x <- group_data_patient(
    data,
    arm = "arm",
    treatment = "Lev+5FU",
    time = "time",
    event = "event"
  )
  expect_equal(x$treatment, as.integer(data$arm == "Lev+5FU"))
  expect_equal(x$time, data$time)
  expect_equal(x$event, data$event)
  expect_equal(sum(x$treatment), 304L)
  expect_equal(sum(x$event), 291L)
})

test_that("group_data_patient() treatment argument sets the coding", {
  data <- group_example_patient()
  x <- group_data_patient(
    data,
    arm = "arm",
    treatment = "Obs",
    time = "time",
    event = "event"
  )
  expect_equal(x$treatment, as.integer(data$arm == "Obs"))
  expect_equal(sum(x$treatment), 315L)
})

test_that("group_data_patient() accepts flexible input types", {
  data <- patient_data()
  data$arm <- factor(data$arm, levels = c("a", "b", "unused"))
  data$time <- as.integer(c(1, 2, 3, 0))
  data$event <- as.logical(data$event)
  data$extra <- NA
  x <- patient_call(data)
  expect_equal(x$treatment, c(1L, 0L, 1L, 0L))
  expect_type(x$time, "double")
  expect_equal(x$time, c(1, 2, 3, 0))
  expect_equal(x$event, c(1L, 0L, 1L, 0L))
})

test_that("group_data_patient() accepts numeric arms and custom names", {
  data <- patient_data()
  colnames(data) <- c("z", "t", "d")
  data$z <- c(1, 0, 1, 0)
  x <- group_data_patient(data, "z", 0, "t", "d")
  expect_equal(colnames(x), c("treatment", "time", "event"))
  expect_equal(x$treatment, c(0L, 1L, 0L, 1L))
  expect_equal(x$time, c(1, 2.5, 3, 0))
})

test_that("group_data_patient() validates data", {
  expect_error(
    patient_call(as.list(patient_data())),
    "data frame",
    class = "retrogroup_error"
  )
  expect_error(
    patient_call(patient_data()[0, ]),
    "at least one row",
    class = "retrogroup_error"
  )
})

test_that("group_data_patient() validates column name arguments", {
  expect_error(
    patient_call(arm = NA_character_),
    "arm must be a single",
    class = "retrogroup_error"
  )
  expect_error(
    patient_call(arm = ""),
    "arm must be a single",
    class = "retrogroup_error"
  )
  expect_error(
    patient_call(time = 1),
    "time must be a single",
    class = "retrogroup_error"
  )
  expect_error(
    patient_call(event = c("a", "b")),
    "event must be a single",
    class = "retrogroup_error"
  )
  expect_error(
    patient_call(event = "time"),
    "different columns",
    class = "retrogroup_error"
  )
  expect_error(
    patient_call(arm = "nope"),
    "not in data",
    class = "retrogroup_error"
  )
})

test_that("group_data_patient() validates arm", {
  data <- patient_data()
  data$arm <- as.list(data$arm)
  expect_error(patient_call(data), "atomic", class = "retrogroup_error")
  data <- patient_data()
  data$arm[1] <- NA
  expect_error(patient_call(data), "missing", class = "retrogroup_error")
  data$arm[1] <- "c"
  expect_error(patient_call(data), "exactly two", class = "retrogroup_error")
  data$arm <- "a"
  expect_error(patient_call(data), "exactly two", class = "retrogroup_error")
})

test_that("group_data_patient() validates treatment", {
  expect_error(
    patient_call(treatment = c("a", "b")),
    "single",
    class = "retrogroup_error"
  )
  expect_error(
    patient_call(treatment = NA),
    "single",
    class = "retrogroup_error"
  )
  expect_error(
    patient_call(treatment = "c"),
    "one of the values",
    class = "retrogroup_error"
  )
})

test_that("group_data_patient() validates time", {
  data <- patient_data()
  data$time <- as.character(data$time)
  expect_error(
    patient_call(data),
    "must be numeric",
    class = "retrogroup_error"
  )
  data <- patient_data()
  data$time[1] <- NA
  expect_error(patient_call(data), "finite", class = "retrogroup_error")
  data$time[1] <- Inf
  expect_error(patient_call(data), "finite", class = "retrogroup_error")
  data$time[1] <- -1
  expect_error(patient_call(data), "non-negative", class = "retrogroup_error")
})

test_that("group_data_patient() validates event", {
  data <- patient_data()
  data$event[1] <- NA
  expect_error(patient_call(data), "missing", class = "retrogroup_error")
  data$event[1] <- 2L
  expect_error(patient_call(data), "0/1", class = "retrogroup_error")
  data$event <- as.character(patient_data()$event)
  expect_error(patient_call(data), "0/1", class = "retrogroup_error")
  data$event <- 0L
  expect_error(
    patient_call(data),
    "at least one event",
    class = "retrogroup_error"
  )
})
