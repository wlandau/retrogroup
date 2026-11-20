#' @title Summary subgroup statistics
#' @export
#' @family data
#' @description Convert assumed subgroup statistics into the
#'   standardized format used by `retrogroup`.
#' @details Each row of `data` is one level of a subgroup variable, with
#'   the number of patients \eqn{n_k}, the number of events \eqn{e_k}, and
#'   the hazard ratio \eqn{\exp(\beta_k)} of treatment versus control
#'   defined in `vignette("framing", package = "retrogroup")`.
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
#'
#'   Likewise, the event counts must have the same total for every
#'   subgroup variable, and each level must satisfy
#'   \eqn{1 \le e_k \le n_k}. See the "Event counts" section of
#'   `vignette("framing", package = "retrogroup")` for how to choose event
#'   counts, and [group_data_events()] for approximating them from
#'   confidence intervals when the publication does not report them.
#'   For that purpose, `events = NULL` skips the checks on event counts
#'   and sets every `events` value to `NA`. Such a summary is only valid
#'   input to [group_data_events()]. Every other function rejects it.
#' @return A tibble with one row per subgroup level, in the same order as
#'   `data`, and the following columns:
#'   * `variable`: character, name of the subgroup variable.
#'   * `level`: character, level of the subgroup variable.
#'   * `patients`: integer, number of patients \eqn{n_k} in the subgroup
#'     level.
#'   * `events`: integer, number of events \eqn{e_k} in the subgroup level,
#'     or `NA` if `events` is `NULL`.
#'   * `hazard_ratio`: numeric, assumed hazard ratio
#'     \eqn{\exp(\beta_k)} of treatment versus control in the subgroup
#'     level.
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
#' @param events Character string, name of the column in `data` with the
#'   number of events to assign to each level of the reconstructed
#'   subgroup variable. The column must contain whole numbers with no
#'   missing values, each at least 1 and at most the number of patients
#'   in the level. The event counts must have the same total for every
#'   subgroup variable. Alternatively, `NULL` to set every `events` value
#'   to `NA` and approximate the event counts later with
#'   [group_data_events()].
#' @param hazard_ratio Character string, name of the column in `data`
#'   with the assumed hazard ratio of treatment versus control in each
#'   subgroup level, on the natural scale (not the log scale). The column
#'   must be numeric, finite, and strictly positive with no missing values.
#' @examples
#'   group_data_summary(
#'     group_example_summary(),
#'     variable = "variable",
#'     level = "level",
#'     patients = "patients",
#'     events = "events",
#'     hazard_ratio = "hazard_ratio"
#'   )
group_data_summary <- function(
  data,
  variable,
  level,
  patients,
  events,
  hazard_ratio
) {
  assert_data_summary(
    data,
    variable,
    level,
    patients,
    events,
    hazard_ratio
  )
  tibble::tibble(
    variable = as.character(data[[variable]]),
    level = as.character(data[[level]]),
    patients = as.integer(data[[patients]]),
    events = if (is.null(events)) {
      rep(NA_integer_, nrow(data))
    } else {
      as.integer(data[[events]])
    },
    hazard_ratio = as.numeric(data[[hazard_ratio]])
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
#'     events = "events",
#'     hazard_ratio = "hazard_ratio"
#'   )
#'   try(
#'     assert_data_summary(
#'       data,
#'       variable = "variable",
#'       level = "level",
#'       patients = "missing_column",
#'       events = "events",
#'       hazard_ratio = "hazard_ratio"
#'     )
#'   )
assert_data_summary <- function(
  data,
  variable,
  level,
  patients,
  events,
  hazard_ratio
) {
  # data
  if (!is.data.frame(data)) {
    throw_error("data must be a data frame.")
  }
  if (nrow(data) < 2L) {
    throw_error("data must have at least two rows.")
  }
  # column name arguments
  columns <- list(
    variable = variable,
    level = level,
    patients = patients,
    events = events,
    hazard_ratio = hazard_ratio
  )
  assert_columns(data, Filter(Negate(is.null), columns))
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
  n <- values <- data[[patients]]
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
  # events
  if (!is.null(events)) {
    assert_summary_events(data[[events]], events, variables, levels, n)
  }
  # hazard_ratio
  values <- data[[hazard_ratio]]
  label <- paste0("column \"", hazard_ratio, "\" (hazard_ratio)")
  if (!is.numeric(values)) {
    throw_error(paste(label, "must be numeric."))
  }
  if (anyNA(values) || !all(is.finite(values))) {
    throw_error(paste(label, "must be finite with no missing values."))
  }
  if (any(values <= 0)) {
    throw_error(
      c(
        paste(label, "must be strictly positive."),
        i = paste(
          "Hazard ratios must be on the natural scale.",
          "A value of 0 or below usually comes from a typo",
          "or from supplying log hazard ratios."
        )
      )
    )
  }
  invisible()
}

#' @title Validate the event counts of `group_data_summary()`
#' @noRd
#' @return `NULL` (invisibly). Throws an error if the event counts are
#'   invalid.
#' @param values The column of event counts.
#' @param events Character string, name of the column of event counts.
#' @param variables Character vector, subgroup variable of each row.
#' @param levels Character vector, level of each row.
#' @param n Numeric vector, number of patients of each row.
#' @examples
#'   assert_summary_events(
#'     c(2L, 1L, 1L, 2L),
#'     "events",
#'     c("age", "age", "sex", "sex"),
#'     c("young", "old", "female", "male"),
#'     c(3L, 2L, 2L, 3L)
#'   )
#'   try(
#'     assert_summary_events(
#'       c(0L, 3L, 1L, 2L),
#'       "events",
#'       c("age", "age", "sex", "sex"),
#'       c("young", "old", "female", "male"),
#'       c(3L, 2L, 2L, 3L)
#'     )
#'   )
assert_summary_events <- function(values, events, variables, levels, n) {
  label <- paste0("column \"", events, "\" (events)")
  if (!is.numeric(values)) {
    throw_error(paste(label, "must be numeric."))
  }
  if (anyNA(values) || !all(is.finite(values))) {
    throw_error(paste(label, "must be finite with no missing values."))
  }
  if (any(values != round(values))) {
    throw_error(paste(label, "must contain whole numbers."))
  }
  invalid <- values < 1 | values > n
  if (any(invalid)) {
    throw_error(
      c(
        paste(
          label,
          "must be at least 1 and at most the number of patients",
          "in each subgroup level."
        ),
        x = paste0(
          "Invalid subgroup level(s): ",
          paste0(
            variables[invalid],
            " = \"",
            levels[invalid],
            "\" (",
            values[invalid],
            " events, ",
            n[invalid],
            " patients)",
            collapse = ", "
          ),
          "."
        ),
        i = "A subgroup level with no events cannot support a hazard ratio."
      )
    )
  }
  totals <- tapply(values, variables, sum)[unique(variables)]
  if (length(unique(totals)) > 1L) {
    throw_error(
      c(
        paste(
          "The total number of events must be the same",
          "for every subgroup variable."
        ),
        x = paste0(
          "Totals by subgroup variable: ",
          paste(names(totals), "=", totals, collapse = ", "),
          "."
        ),
        i = message_count_adjustment("event")
      )
    )
  }
  invisible()
}
