#' @title Patient-level survival data
#' @export
#' @family data
#' @description Convert a patient-level right-censored survival dataset
#'   from a randomized 2-arm trial into the standardized format
#'   used by `retrogroup`.
#' @return A tibble with one row per patient and the following columns:
#'   * `treatment`: integer, `1` for the treatment arm and `0` for the
#'     control arm.
#'   * `time`: numeric, time since randomization of the event or censoring.
#'   * `event`: integer, `1` if the event occurred and `0` if censored.
#' @param data Data frame with one row per patient. Must have the columns
#'   named by `arm`, `time`, and `event`. Other columns are ignored.
#' @param arm Character string, name of the column in `data` with the
#'   treatment arm. Must have no missing values and exactly two
#'   unique values.
#' @param treatment The value of `data[[arm]]` that denotes the
#'   treatment arm. The other value denotes the control arm.
#' @param time Character string, name of the column in `data` with the
#'   time since randomization of the event or censoring. The column must be
#'   numeric, finite, non-negative, and free of missing values.
#' @param event Character string, name of the column in `data` with the
#'   event indicator: `1` or `TRUE` for an event, `0` or `FALSE` for
#'   right censoring. Must have no missing values and at least one event.
#' @examples
#'   data <- group_example_patient()
#'   group_data_patient(
#'     data,
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
group_data_patient <- function(data, arm, treatment, time, event) {
  assert_data_patient(data, arm, treatment, time, event)
  tibble::tibble(
    treatment = as.integer(data[[arm]] == treatment),
    time = as.numeric(data[[time]]),
    event = as.integer(data[[event]])
  )
}

#' @title Validate the arguments of `group_data_patient()`
#' @noRd
#' @return `NULL` (invisibly). Throws an error if any argument is invalid.
#' @inheritParams group_data_patient
#' @examples
#'   data <- group_example_patient()
#'   assert_data_patient(data, "arm", "Lev+5FU", "time", "event")
#'   try(assert_data_patient(data, "arm", "placebo", "time", "event"))
assert_data_patient <- function(data, arm, treatment, time, event) {
  # data
  if (!is.data.frame(data)) {
    throw_error("data must be a data frame.")
  }
  if (nrow(data) < 1L) {
    throw_error("data must have at least one row.")
  }
  # column name arguments
  assert_columns(data, list(arm = arm, time = time, event = event))
  # arm
  values <- data[[arm]]
  if (!(is.atomic(values) || is.factor(values))) {
    throw_error(
      paste0("column \"", arm, "\" (arm) must be an atomic vector or factor.")
    )
  }
  if (anyNA(values)) {
    throw_error(
      paste0("column \"", arm, "\" (arm) must have no missing values.")
    )
  }
  if (length(unique(values)) != 2L) {
    throw_error(
      paste0("column \"", arm, "\" (arm) must have exactly two unique values.")
    )
  }
  # treatment
  if (length(treatment) != 1L || is.na(treatment)) {
    throw_error("treatment must be a single non-missing value.")
  }
  if (!treatment %in% values) {
    throw_error(
      paste0(
        "treatment must be one of the values in column \"",
        arm,
        "\" (arm)."
      )
    )
  }
  # time
  values <- data[[time]]
  if (!is.numeric(values)) {
    throw_error(paste0("column \"", time, "\" (time) must be numeric."))
  }
  if (anyNA(values) || !all(is.finite(values))) {
    throw_error(
      paste0(
        "column \"",
        time,
        "\" (time) must be finite with no missing values."
      )
    )
  }
  if (any(values < 0)) {
    throw_error(paste0("column \"", time, "\" (time) must be non-negative."))
  }
  # event
  values <- data[[event]]
  if (anyNA(values)) {
    throw_error(
      paste0("column \"", event, "\" (event) must have no missing values.")
    )
  }
  if (!(is.logical(values) || (is.numeric(values) && all(values %in% 0:1)))) {
    throw_error(
      paste0(
        "column \"",
        event,
        "\" (event) must contain only 0/1 or TRUE/FALSE."
      )
    )
  }
  if (sum(values) < 1L) {
    throw_error(
      paste0("column \"", event, "\" (event) must have at least one event.")
    )
  }
  invisible()
}
