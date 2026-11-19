#' @title Summary subgroup statistics
#' @export
#' @family data
#' @description Convert published subgroup statistics into the
#'   standardized format used by `retrogroup`.
#' @details Each row of `data` is one level of a subgroup variable.
#'   A subgroup variable must have at least two levels, the levels must be
#'   unique within the subgroup variable, and every level must have a
#'   strictly positive number of patients.
#'
#'   The patient counts must have the same total for every subgroup
#'   variable, because every subgroup variable partitions the same
#'   patients. The counts are the sizes of the levels of the reconstructed
#'   subgroup variables, so they need not copy the publication verbatim.
#'   Published counts are the natural starting point, but they sometimes
#'   do not add up to the total number of patients in the study, for
#'   example because of missing subgroup data or typographical errors. In
#'   that case, choosing counts that add up to the total number of patients
#'   is a decision about the reconstruction, not an attempt to
#'   recover the true counts from the clinical study. Please disclose that
#'   decision transparently wherever the reconstructed data or downstream
#'   results are reported.
#' @return A tibble with one row per subgroup level, in the same order as
#'   `data`, and the following columns:
#'   * `variable`: character, name of the subgroup variable.
#'   * `level`: character, level of the subgroup variable.
#'   * `patients`: integer, number of patients in the subgroup level.
#'   * `log_hazard_ratio`: numeric, published log hazard ratio of treatment
#'     versus control in the subgroup level.
#'   * `standard_error`: numeric, published standard error of the log
#'     hazard ratio, or `NA` if `standard_error` is `NULL`.
#' @param data Data frame with one row per subgroup level. Must have the
#'   columns named by the other arguments. Other columns are ignored.
#' @param variable Character string, name of the column in `data` with the
#'   name of the subgroup variable. The column must be character with
#'   non-empty, non-missing values.
#' @param level Character string, name of the column in `data` with the
#'   level of the subgroup variable. The column must be character or
#'   factor with non-empty, non-missing values.
#' @param patients Character string, name of the column in `data` with the
#'   number of patients to assign to each level of the reconstructed
#'   subgroup variable, usually taken from the publication. See the
#'   details section for when these counts may differ from the publication.
#'   The column must contain
#'   strictly positive whole numbers with no missing values.
#' @param log_hazard_ratio Character string, name of the column in `data`
#'   with the published log hazard ratio of treatment versus control in
#'   each subgroup level. The column must be numeric and finite with no
#'   missing values.
#' @param standard_error Character string or `NULL`. If a character string,
#'   name of the column in `data` with the published standard error of the
#'   log hazard ratio in each subgroup level. The column must be numeric,
#'   finite, and strictly positive with no missing values. If `NULL`,
#'   the standard errors are unavailable, and the `standard_error` column
#'   of the output is `NA`.
#' @examples
#'   data <- group_example_summary()
#'   group_data_summary(
#'     data,
#'     variable = "variable",
#'     level = "level",
#'     patients = "patients",
#'     log_hazard_ratio = "log_hazard_ratio",
#'     standard_error = "standard_error"
#'   )
#'   group_data_summary(
#'     data,
#'     variable = "variable",
#'     level = "level",
#'     patients = "patients",
#'     log_hazard_ratio = "log_hazard_ratio"
#'   )
group_data_summary <- function(
  data,
  variable,
  level,
  patients,
  log_hazard_ratio,
  standard_error = NULL
) {
  assert_data_summary(
    data,
    variable,
    level,
    patients,
    log_hazard_ratio,
    standard_error
  )
  tibble::tibble(
    variable = as.character(data[[variable]]),
    level = as.character(data[[level]]),
    patients = as.integer(data[[patients]]),
    log_hazard_ratio = as.numeric(data[[log_hazard_ratio]]),
    standard_error = if (is.null(standard_error)) {
      NA_real_
    } else {
      as.numeric(data[[standard_error]])
    }
  )
}

#' @title Validate the arguments of `group_data_summary()`
#' @noRd
#' @return `NULL` (invisibly). Throws an error if any argument is invalid.
#' @inheritParams group_data_summary
#' @examples
#'   data <- group_example_summary()
#'   assert_data_summary(
#'     data,
#'     variable = "variable",
#'     level = "level",
#'     patients = "patients",
#'     log_hazard_ratio = "log_hazard_ratio",
#'     standard_error = "standard_error"
#'   )
#'   try(
#'     assert_data_summary(
#'       data,
#'       variable = "variable",
#'       level = "level",
#'       patients = "missing_column",
#'       log_hazard_ratio = "log_hazard_ratio",
#'       standard_error = NULL
#'     )
#'   )
assert_data_summary <- function(
  data,
  variable,
  level,
  patients,
  log_hazard_ratio,
  standard_error
) {
  # data
  if (!is.data.frame(data)) {
    throw_error("data must be a data frame.")
  }
  if (nrow(data) < 2L) {
    throw_error("data must have at least two rows.")
  }
  # column name arguments
  if (!is.null(standard_error)) {
    assert_string(standard_error, "standard_error")
  }
  assert_columns(
    data,
    list(
      variable = variable,
      level = level,
      patients = patients,
      log_hazard_ratio = log_hazard_ratio,
      standard_error = standard_error
    )
  )
  # variable
  variables <- data[[variable]]
  label <- paste0("column \"", variable, "\" (variable)")
  if (!is.character(variables)) {
    throw_error(paste(label, "must be character."))
  }
  if (anyNA(variables) || !all(nzchar(variables))) {
    throw_error(paste(label, "must have no missing or empty values."))
  }
  # level
  levels <- data[[level]]
  label <- paste0("column \"", level, "\" (level)")
  if (!(is.character(levels) || is.factor(levels))) {
    throw_error(paste(label, "must be character or factor."))
  }
  levels <- as.character(levels)
  if (anyNA(levels) || !all(nzchar(levels))) {
    throw_error(paste(label, "must have no missing or empty values."))
  }
  counts <- table(variables)
  single <- names(counts)[counts < 2L]
  if (length(single)) {
    throw_error(
      c(
        "Each subgroup variable must have more than one level.",
        x = paste0(
          "Subgroup variable(s) with only one level: ",
          paste0("\"", single, "\"", collapse = ", "),
          "."
        )
      )
    )
  }
  duplicated <- unique(variables[duplicated(data.frame(variables, levels))])
  if (length(duplicated)) {
    throw_error(
      c(
        "Levels must be unique within each subgroup variable.",
        x = paste0(
          "Subgroup variable(s) with duplicated levels: ",
          paste0("\"", duplicated, "\"", collapse = ", "),
          "."
        )
      )
    )
  }
  # patients
  values <- data[[patients]]
  label <- paste0("column \"", patients, "\" (patients)")
  if (!is.numeric(values)) {
    throw_error(paste(label, "must be numeric."))
  }
  if (anyNA(values) || !all(is.finite(values))) {
    throw_error(paste(label, "must be finite with no missing values."))
  }
  if (any(values != round(values))) {
    throw_error(paste(label, "must contain whole numbers."))
  }
  if (any(values <= 0)) {
    throw_error(paste(label, "must be strictly positive."))
  }
  if (any(values > .Machine$integer.max)) {
    throw_error(
      paste(label, "must be at most", .Machine$integer.max, "patients.")
    )
  }
  totals <- tapply(values, variables, sum)[unique(variables)]
  if (length(unique(totals)) > 1L) {
    throw_error(
      c(
        paste(
          "The total number of patients must be the same",
          "for every subgroup variable."
        ),
        x = paste0(
          "Totals by subgroup variable: ",
          paste(names(totals), "=", totals, collapse = ", "),
          "."
        ),
        i = message_count_adjustment()
      )
    )
  }
  # log_hazard_ratio
  values <- data[[log_hazard_ratio]]
  label <- paste0("column \"", log_hazard_ratio, "\" (log_hazard_ratio)")
  if (!is.numeric(values)) {
    throw_error(paste(label, "must be numeric."))
  }
  if (anyNA(values) || !all(is.finite(values))) {
    throw_error(
      c(
        paste(label, "must be finite with no missing values."),
        i = paste(
          "Non-finite log hazard ratios usually come from a typo",
          "or from taking the log of a hazard ratio of 0."
        )
      )
    )
  }
  # standard_error
  if (!is.null(standard_error)) {
    values <- data[[standard_error]]
    label <- paste0("column \"", standard_error, "\" (standard_error)")
    if (!is.numeric(values)) {
      throw_error(paste(label, "must be numeric."))
    }
    if (anyNA(values) || !all(is.finite(values))) {
      throw_error(paste(label, "must be finite with no missing values."))
    }
    if (any(values <= 0)) {
      throw_error(paste(label, "must be strictly positive."))
    }
  }
  invisible()
}
