probability_example <- function() {
  patient <- group_data_patient(
    group_example_patient(),
    arm = "arm",
    treatment = "Lev+5FU",
    time = "time",
    event = "event"
  )
  summary <- group_data_summary(
    group_example_summary(),
    variable = "variable",
    level = "level",
    patients = "patients",
    events = "events",
    hazard_ratio = "hazard_ratio"
  )
  list(patient = patient, summary = summary)
}

probability_small <- function() {
  patient <- tibble::tibble(
    treatment = c(0L, 1L, 0L, 1L, 0L, 1L),
    time = c(1, 2, 2, 3, 4, 0.5),
    event = c(1L, 1L, 0L, 0L, 1L, 1L)
  )
  summary <- tibble::tibble(
    variable = c("age", "age", "sex", "sex", "sex"),
    level = c("young", "old", "male", "female", "other"),
    patients = c(2L, 4L, 1L, 2L, 3L),
    events = c(1L, 3L, 1L, 1L, 2L),
    hazard_ratio = exp(c(-0.5, 0.3, 0.2, -0.1, 0))
  )
  list(patient = patient, summary = summary)
}
