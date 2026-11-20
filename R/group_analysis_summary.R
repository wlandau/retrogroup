#' @title Subgroup-level statistics
#' @export
#' @family analysis
#' @description Compute the number of patients, the number of events, and
#'   the hazard ratio of treatment versus control in each level of each
#'   subgroup variable of a patient-level dataset.
#' @details For each level of each subgroup variable,
#'   `group_analysis_summary()` counts the patients and events in that
#'   level. These are the quantities \eqn{n_k} and \eqn{e_k} of
#'   `vignette("framing", package = "retrogroup")`, computed from the
#'   subgroup assignments in `patient` instead of assumed. Unlike the
#'   assumed \eqn{e_k}, the number of events of a level can be 0.
#'
#'   `group_analysis_summary()` also fits a Cox proportional hazards model
#'   by partial likelihood to the patients in each level, with treatment as
#'   the only covariate and the Efron method for ties, as in
#'   [survival::coxph()]. The hazard ratio is the exponential of the
#'   estimated coefficient of treatment versus control, the estimate of
#'   \eqn{\exp(\beta_k)} of `vignette("framing", package = "retrogroup")`.
#'   The log hazard ratio is the observed statistic \eqn{\hat\beta_k(x)} of
#'   the "Loss function" section of
#'   `vignette("optimization", package = "retrogroup")`.
#'
#'   If the model cannot be fit for some subgroup level, for example
#'   because one of its arms has no patients or no events, the
#'   `hazard_ratio`, `log_hazard_ratio`, and `standard_error` of that level
#'   are `NA`.
#'
#'   The levels of a factor subgroup column are its factor levels, in
#'   order, excluding levels with no patients. The levels of a character
#'   subgroup column are its sorted unique values.
#' @return A tibble with one row per subgroup level, ordered by subgroup
#'   variable as in `variables` and then by level, and the following
#'   columns:
#'   * `variable`: character, name of the subgroup variable.
#'   * `level`: character, level of the subgroup variable.
#'   * `patients`: integer, number of patients in the subgroup level.
#'   * `events`: integer, number of events in the subgroup level.
#'   * `hazard_ratio`: numeric, estimated hazard ratio of treatment
#'     versus control in the subgroup level.
#'   * `log_hazard_ratio`: numeric, estimated log hazard ratio of treatment
#'     versus control in the subgroup level. Legacy column, to be removed
#'     in a future version.
#'   * `standard_error`: numeric, standard error of the log hazard ratio.
#'     Legacy column, to be removed in a future version.
#' @param patient Patient-level survival data with subgroup columns,
#'   usually the output of [group_sample_unconstrained()],
#'   [group_sample_constrained()], [group_optimize_unconstrained()], or
#'   [group_optimize_constrained()]. Must have the columns of the output of
#'   [group_data_patient()] and the subgroup columns named in `variables`.
#' @param variables Character vector of unique names of subgroup columns
#'   in `patient`. Each subgroup column must be character or factor with
#'   no missing values. Defaults to every column of `patient` other than
#'   `treatment`, `time`, and `event`.
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   summary <- group_data_summary(
#'     group_example_summary(),
#'     variable = "variable",
#'     level = "level",
#'     patients = "patients",
#'     log_hazard_ratio = "log_hazard_ratio",
#'     standard_error = "standard_error"
#'   )
#'   set.seed(0)
#'   sample <- group_sample_constrained(patient, summary)
#'   group_analysis_summary(sample)
#'   summary
group_analysis_summary <- function(
  patient,
  variables = setdiff(colnames(patient), c("treatment", "time", "event"))
) {
  assert_probability_patient(patient)
  assert_analysis_variables(patient, variables)
  rows <- lapply(variables, \(variable) {
    values <- patient[[variable]]
    levels <- if (is.factor(values)) {
      intersect(levels(values), as.character(values))
    } else {
      sort(unique(values))
    }
    values <- as.character(values)
    fits <- vapply(
      levels,
      \(level) {
        members <- values == level
        c(
          sum(members),
          sum(patient$event[members]),
          analysis_fit(
            patient$treatment[members],
            patient$time[members],
            patient$event[members]
          )
        )
      },
      numeric(4L)
    )
    tibble::tibble(
      variable = variable,
      level = levels,
      patients = as.integer(fits[1L, ]),
      events = as.integer(fits[2L, ]),
      hazard_ratio = exp(unname(fits[3L, ])),
      log_hazard_ratio = unname(fits[3L, ]),
      standard_error = unname(fits[4L, ])
    )
  })
  do.call(rbind, rows)
}

#' @title Fit the Cox model of one subgroup level
#' @noRd
#' @description Fit a Cox proportional hazards model by partial likelihood
#'   with treatment as the only covariate and the Efron method for ties,
#'   as in [survival::coxph()].
#' @return Numeric vector of length 2: the estimated log hazard ratio
#'   and its standard error. Both elements are `NA` if the model cannot be
#'   fit, for example because the level is empty or one of its arms has
#'   no events.
#' @param z Integer vector, treatment indicator of each patient in the
#'   subgroup level.
#' @param t Numeric vector, time of each patient in the subgroup level.
#' @param d Integer vector, event indicator of each patient in the
#'   subgroup level.
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   analysis_fit(patient$treatment, patient$time, patient$event)
#'   analysis_fit(z = c(0L, 1L), t = c(1, 2), d = c(0L, 0L))
analysis_fit <- function(z, t, d) {
  if (!any(d[z == 1L] == 1L) || !any(d[z == 0L] == 1L)) {
    return(c(NA_real_, NA_real_))
  }
  fit <- tryCatch(
    suppressWarnings(
      survival::coxph.fit(
        x = matrix(as.numeric(z), ncol = 1L),
        y = survival::Surv(t, d),
        strata = NULL,
        offset = NULL,
        init = NULL,
        control = survival::coxph.control(),
        weights = NULL,
        method = "efron",
        rownames = NULL,
        resid = FALSE
      )
    ),
    error = \(error) NULL
  )
  if (is.null(fit)) {
    return(c(NA_real_, NA_real_))
  }
  out <- c(unname(fit$coefficients[1L]), sqrt(fit$var[1L, 1L]))
  if (!all(is.finite(out))) {
    return(c(NA_real_, NA_real_))
  }
  out
}

#' @title Validate the `variables` argument of `group_analysis_summary()`
#' @noRd
#' @return `NULL` (invisibly). Throws an error if `variables` is invalid.
#' @inheritParams group_analysis_summary
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   patient$sex <- rep(c("female", "male"), length.out = nrow(patient))
#'   assert_analysis_variables(patient, "sex")
#'   try(assert_analysis_variables(patient, "age"))
assert_analysis_variables <- function(patient, variables) {
  if (
    !is.character(variables) ||
      length(variables) < 1L ||
      anyNA(variables) ||
      !all(nzchar(variables))
  ) {
    throw_error(
      "variables must be a character vector of non-missing non-empty names."
    )
  }
  if (anyDuplicated(variables)) {
    throw_error("variables must be unique.")
  }
  reserved <- intersect(variables, c("treatment", "time", "event"))
  if (length(reserved)) {
    throw_error(
      paste0(
        "variables must not include \"treatment\", \"time\", or \"event\": ",
        paste0("\"", reserved, "\"", collapse = ", "),
        "."
      )
    )
  }
  missing <- setdiff(variables, colnames(patient))
  if (length(missing)) {
    throw_error(
      paste0(
        "column(s) not in patient: ",
        paste0("\"", missing, "\"", collapse = ", "),
        "."
      )
    )
  }
  for (variable in variables) {
    values <- patient[[variable]]
    label <- paste0("column \"", variable, "\" of patient")
    if (!(is.character(values) || is.factor(values))) {
      throw_error(paste(label, "must be character or factor."))
    }
    if (anyNA(values)) {
      throw_error(paste(label, "must have no missing values."))
    }
  }
  invisible()
}
