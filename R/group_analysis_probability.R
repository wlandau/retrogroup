#' @title Posterior subgroup probabilities
#' @export
#' @family analysis
#' @description Estimate the posterior probability that each patient
#'   belongs to each level of each subgroup variable, given the
#'   patient-level survival data and the assumed subgroup statistics.
#' @details `group_analysis_probability()` computes the estimated posterior
#'   probabilities \eqn{p_{ik}} defined in the "Computation" section of
#'   `vignette("model", package = "retrogroup")`, separately for each
#'   subgroup variable. See that section for the formulas, including the
#'   error when any probability is numerically indistinguishable from zero.
#'
#'   The total number of patients of every subgroup variable in
#'   `summary` must equal the number of rows of `patient`, and the total
#'   number of events of every subgroup variable must equal the number of
#'   events in `patient`. See the details of [group_data_summary()] for how
#'   to choose the counts. The model uses the patient counts and the log of
#'   the hazard ratios in `summary`, not the event counts.
#' @return A named list with one element per subgroup variable, in the
#'   order the subgroup variables first appear in `summary$variable`.
#'   Each element is a numeric matrix of the estimated posterior
#'   probabilities \eqn{p_{ik}} of that subgroup variable:
#'   * Rows correspond to patients, in the same order as the rows of
#'     `patient`.
#'   * Columns correspond to the levels of the subgroup variable. The
#'     column names are exactly the levels, in the same order as in
#'     `summary`.
#'   * Each row sums to 1.
#' @param patient Patient-level survival data, the output of
#'   [group_data_patient()].
#' @param summary Summary subgroup statistics, the output of
#'   [group_data_summary()].
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
#'     events = "events",
#'     hazard_ratio = "hazard_ratio"
#'   )
#'   probability <- group_analysis_probability(patient, summary)
#'   names(probability)
#'   head(probability$age)
group_analysis_probability <- function(patient, summary) {
  assert_probability_patient(patient)
  assert_probability_summary(summary)
  assert_probability_joint(patient, summary)
  z <- patient$treatment
  t <- patient$time
  d <- patient$event
  Lambda_0 <- probability_nelson_aalen(z = z, t = t, d = d)
  variables <- unique(summary$variable)
  out <- lapply(variables, \(variable) {
    rows <- summary$variable == variable
    probability_variable(
      variable = variable,
      level = summary$level[rows],
      n = summary$patients[rows],
      beta = log(summary$hazard_ratio[rows]),
      z = z,
      d = d,
      Lambda_0 = Lambda_0
    )
  })
  names(out) <- variables
  out
}

#' @title Posterior probabilities of one subgroup variable
#' @noRd
#' @return Numeric matrix of \eqn{p_{ik}} with one row per patient and one
#'   column per level, with column names `level`.
#' @param variable Character string, name of the subgroup variable.
#' @param level Character vector, levels of the subgroup variable.
#' @param n Integer vector, assumed number of patients in each level.
#' @param beta Numeric vector, log of the assumed hazard ratio in each
#'   level.
#' @param z Integer vector, treatment indicator of each patient.
#' @param d Integer vector, event indicator of each patient.
#' @param Lambda_0 Numeric vector, estimated baseline cumulative hazard
#'   \eqn{\hat{\Lambda}_0(t_i)} of each patient.
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   summary <- group_example_summary()
#'   sex <- summary[summary$variable == "sex", ]
#'   z <- patient$treatment
#'   t <- patient$time
#'   d <- patient$event
#'   p <- probability_variable(
#'     variable = "sex",
#'     level = sex$level,
#'     n = sex$patients,
#'     beta = log(sex$hazard_ratio),
#'     z = z,
#'     d = d,
#'     Lambda_0 = probability_nelson_aalen(z, t, d)
#'   )
#'   head(p)
probability_variable <- function(
  variable,
  level,
  n,
  beta,
  z,
  d,
  Lambda_0
) {
  w <- matrix(log(n), nrow = length(z), ncol = length(n), byrow = TRUE) +
    outer(d * z, beta) -
    Lambda_0 * exp(outer(z, beta))
  w <- w - apply(w, 1L, max)
  p <- exp(w) / rowSums(exp(w))
  colnames(p) <- level
  assert_probability_positive(p, variable)
  p
}

#' @title Nelson-Aalen baseline cumulative hazard from the control arm
#' @noRd
#' @return Numeric vector with \eqn{\hat{\Lambda}_0(t_i)} for each patient.
#' @param z Integer vector, treatment indicator of each patient.
#' @param t Numeric vector, time of each patient.
#' @param d Integer vector, event indicator of each patient.
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   Lambda_0 <- probability_nelson_aalen(
#'     z = patient$treatment,
#'     t = patient$time,
#'     d = patient$event
#'   )
#'   head(Lambda_0)
probability_nelson_aalen <- function(z, t, d) {
  control <- z == 0L
  t_c <- t[control]
  d_c <- d[control]
  order <- order(t_c)
  t_c <- t_c[order]
  d_c <- d_c[order]
  r <- length(t_c) - findInterval(t_c, t_c, left.open = TRUE)
  cumulative <- c(0, cumsum(d_c / r))
  cumulative[findInterval(t, t_c) + 1L]
}

#' @title Validate the `patient` argument of `group_analysis_probability()`
#' @noRd
#' @return `NULL` (invisibly). Throws an error if `patient` is invalid.
#' @inheritParams group_analysis_probability
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   assert_probability_patient(patient)
#'   try(assert_probability_patient(group_example_patient()))
assert_probability_patient <- function(patient) {
  hint <- "patient must be the output of group_data_patient()."
  if (!is.data.frame(patient)) {
    throw_error(c("patient must be a data frame.", i = hint))
  }
  missing <- setdiff(c("treatment", "time", "event"), colnames(patient))
  if (length(missing)) {
    throw_error(
      c(
        paste0(
          "patient is missing column(s): ",
          paste0("\"", missing, "\"", collapse = ", "),
          "."
        ),
        i = hint
      )
    )
  }
  values <- patient$treatment
  if (!is.numeric(values) || anyNA(values) || !all(values %in% 0:1)) {
    throw_error(
      c(
        "column \"treatment\" of patient must contain only 0 and 1.",
        i = hint
      )
    )
  }
  tryCatch(
    assert_data_patient(
      patient,
      arm = "treatment",
      treatment = 1L,
      time = "time",
      event = "event"
    ),
    retrogroup_error = \(error) {
      throw_error(c("patient is invalid.", i = hint), parent = error)
    }
  )
  invisible()
}

#' @title Validate the `summary` argument of `group_analysis_probability()`
#' @noRd
#' @return `NULL` (invisibly). Throws an error if `summary` is invalid.
#' @inheritParams group_analysis_probability
#' @param events Logical, whether to check the event counts. `FALSE`
#'   does not require an `events` column and ignores it, for
#'   [group_data_events()].
#' @examples
#'   assert_probability_summary(group_example_summary())
#'   try(assert_probability_summary(group_example_patient()))
assert_probability_summary <- function(summary, events = TRUE) {
  hint <- "summary must be the output of group_data_summary()."
  if (!is.data.frame(summary)) {
    throw_error(c("summary must be a data frame.", i = hint))
  }
  columns <- c("variable", "level", "patients", "events", "hazard_ratio")
  if (!events) {
    columns <- setdiff(columns, "events")
  }
  missing <- setdiff(columns, colnames(summary))
  if (length(missing)) {
    throw_error(
      c(
        paste0(
          "summary is missing column(s): ",
          paste0("\"", missing, "\"", collapse = ", "),
          "."
        ),
        i = hint
      )
    )
  }
  tryCatch(
    assert_data_summary(
      summary,
      variable = "variable",
      level = "level",
      patients = "patients",
      events = if (events) "events",
      hazard_ratio = "hazard_ratio"
    ),
    retrogroup_error = \(error) {
      throw_error(c("summary is invalid.", i = hint), parent = error)
    }
  )
  invisible()
}

#' @title Validate `patient` and `summary` together
#' @noRd
#' @description Checks that require both `patient` and `summary`.
#'   Assumes `assert_probability_patient()` and
#'   `assert_probability_summary()` already passed.
#' @return `NULL` (invisibly). Throws an error if `patient` and `summary`
#'   are inconsistent with each other.
#' @inheritParams group_analysis_probability
#' @param events Logical, whether to check the total number of events.
#'   `FALSE` skips that check, for [group_data_events()].
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   summary <- group_example_summary()
#'   assert_probability_joint(patient, summary)
#'   try(assert_probability_joint(patient[1:10, ], summary))
assert_probability_joint <- function(patient, summary, events = TRUE) {
  variables <- unique(summary$variable)
  totals <- tapply(summary$patients, summary$variable, sum)[variables]
  if (any(totals != nrow(patient))) {
    throw_error(
      c(
        paste(
          "The total number of patients of every subgroup variable",
          "in summary must equal the number of rows of patient."
        ),
        x = paste0(
          "nrow(patient) = ",
          nrow(patient),
          ". Totals by subgroup variable: ",
          paste(names(totals), "=", totals, collapse = ", "),
          "."
        ),
        i = message_count_adjustment()
      )
    )
  }
  if (!events) {
    return(invisible())
  }
  events <- tapply(summary$events, summary$variable, sum)[variables]
  if (any(events != sum(patient$event))) {
    throw_error(
      c(
        paste(
          "The total number of events of every subgroup variable",
          "in summary must equal the number of events in patient."
        ),
        x = paste0(
          "sum(patient$event) = ",
          sum(patient$event),
          ". Totals by subgroup variable: ",
          paste(names(events), "=", events, collapse = ", "),
          "."
        ),
        i = message_count_adjustment("event")
      )
    )
  }
  invisible()
}

#' @title Validate that posterior probabilities are numerically positive
#' @noRd
#' @return `NULL` (invisibly). Throws an error if any probability is less
#'   than `sqrt(.Machine$double.eps)`.
#' @param p Numeric matrix of posterior probabilities of one subgroup
#'   variable, with one column per level.
#' @param variable Character string, name of the subgroup variable.
#' @examples
#'   p <- matrix(c(0.7, 0.3, 0.4, 0.6), nrow = 2, byrow = TRUE)
#'   colnames(p) <- c("female", "male")
#'   assert_probability_positive(p, "sex")
#'   p[1, ] <- c(1, 0)
#'   try(assert_probability_positive(p, "sex"))
assert_probability_positive <- function(p, variable) {
  threshold <- sqrt(.Machine$double.eps)
  small <- !(p >= threshold)
  if (any(small)) {
    throw_error(
      c(
        paste0(
          "Posterior probabilities of subgroup variable \"",
          variable,
          "\" must be at least sqrt(.Machine$double.eps) = ",
          format(threshold),
          "."
        ),
        x = paste0(
          sum(rowSums(small) > 0L),
          " patient(s) have a numerically zero probability of level(s): ",
          paste0("\"", colnames(p)[colSums(small) > 0L], "\"", collapse = ", "),
          "."
        ),
        i = paste(
          "This usually comes from an extreme assumed hazard ratio or a",
          "very small assumed number of patients in a level."
        )
      )
    )
  }
  invisible()
}
