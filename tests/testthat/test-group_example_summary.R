test_that("group_example_summary() has the expected structure", {
  x <- group_example_summary()
  expect_s3_class(x, "tbl_df")
  expect_equal(nrow(x), 6L)
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
})

test_that("group_example_summary() has the expected contents", {
  x <- group_example_summary()
  expect_equal(x$variable, c("age", "age", "node4", "node4", "sex", "sex"))
  expect_equal(
    x$level,
    c("under 60", "60 or over", "no", "yes", "female", "male")
  )
  n <- nrow(group_example_patient())
  for (variable in unique(x$variable)) {
    expect_equal(sum(x$patients[x$variable == variable]), n)
  }
  expect_true(all(x$patients > 0L))
  expect_true(all(x$standard_error > 0))
})

test_that("group_example_summary() agrees with direct Cox fits", {
  x <- group_example_summary()
  data <- survival::colon
  data <- data[data$etype == 2 & data$rx %in% c("Obs", "Lev+5FU"), ]
  data$treatment <- as.integer(data$rx == "Lev+5FU")
  data$age <- ifelse(data$age < 60, "under 60", "60 or over")
  data$node4 <- ifelse(data$node4 == 1, "yes", "no")
  data$sex <- ifelse(data$sex == 1, "male", "female")
  for (row in seq_len(nrow(x))) {
    variable <- x$variable[row]
    level <- x$level[row]
    subset <- data[data[[variable]] == level, ]
    fit <- survival::coxph(
      survival::Surv(time, status) ~ treatment,
      data = subset
    )
    expect_equal(x$patients[row], nrow(subset))
    expect_equal(x$log_hazard_ratio[row], unname(stats::coef(fit)))
    expect_equal(
      x$standard_error[row],
      unname(sqrt(diag(stats::vcov(fit))))
    )
  }
})
