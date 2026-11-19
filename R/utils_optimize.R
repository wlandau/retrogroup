#' @title Loss term of one subgroup level
#' @noRd
#' @return Numeric scalar, the contribution of the level to the loss
#'   \eqn{\Delta(x)}: \eqn{(\hat\beta_k(x) - \beta_k)^2}, plus
#'   \eqn{(\hat\sigma_k(x) - \sigma_k)^2} if `sigma` is available.
#'   `Inf` if the model could not be fit.
#' @param fit Numeric vector of length 2 from `analysis_fit()`.
#' @param beta Numeric scalar, published log hazard ratio of the level.
#' @param sigma Numeric scalar, published standard error of the level,
#'   or `NA` if unavailable.
#' @examples
#'   optimize_loss_term(c(-0.3, 0.15), beta = -0.25, sigma = 0.1)
#'   optimize_loss_term(c(-0.3, 0.15), beta = -0.25, sigma = NA)
#'   optimize_loss_term(c(NA, NA), beta = -0.25, sigma = 0.1)
optimize_loss_term <- function(fit, beta, sigma) {
  if (anyNA(fit)) {
    return(Inf)
  }
  out <- (fit[1L] - beta)^2
  if (!is.na(sigma)) {
    out <- out + (fit[2L] - sigma)^2
  }
  out
}

#' @title Loss terms of all levels of one subgroup variable
#' @noRd
#' @return Numeric vector with one loss term per level. The loss
#'   \eqn{\Delta(x)} is the sum.
#' @param x Integer vector with one element per patient: the index of
#'   the current subgroup level \eqn{x_i} of each patient.
#' @param z Integer vector, treatment indicator of each patient.
#' @param t Numeric vector, time of each patient.
#' @param d Integer vector, event indicator of each patient.
#' @param beta Numeric vector, published log hazard ratio of each level.
#' @param sigma Numeric vector, published standard error of each level,
#'   with `NA` where unavailable.
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
#'   x <- sample.int(2L, nrow(patient), replace = TRUE)
#'   optimize_loss_terms(
#'     x = x,
#'     z = patient$treatment,
#'     t = patient$time,
#'     d = patient$event,
#'     beta = sex$log_hazard_ratio,
#'     sigma = sex$standard_error
#'   )
optimize_loss_terms <- function(x, z, t, d, beta, sigma) {
  vapply(
    seq_along(beta),
    \(k) {
      members <- x == k
      fit <- analysis_fit(z[members], t[members], d[members])
      optimize_loss_term(fit, beta[k], sigma[k])
    },
    numeric(1L)
  )
}

#' @title Validate the arguments of the optimization functions
#' @noRd
#' @return `NULL` (invisibly). Throws an error if any argument is invalid.
#' @param patient Patient-level data with subgroup columns, the output of
#'   [group_sample_unconstrained()] or [group_sample_constrained()].
#' @param summary Summary subgroup statistics, the output of
#'   [group_data_summary()].
#' @param iterations Maximum number of iterations per subgroup variable.
#' @param tolerance Loss tolerance.
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   summary <- group_example_summary()
#'   sampled <- group_sample_unconstrained(patient, summary)
#'   assert_optimize_arguments(sampled, summary, iterations = 10, tolerance = 1e-6)
#'   try(assert_optimize_arguments(sampled, summary, iterations = -1, tolerance = 1e-6))
assert_optimize_arguments <- function(
  patient,
  summary,
  iterations,
  tolerance
) {
  assert_probability_patient(patient)
  assert_probability_summary(summary)
  assert_probability_joint(patient, summary)
  variables <- unique(summary$variable)
  missing <- setdiff(variables, colnames(patient))
  if (length(missing)) {
    throw_error(
      c(
        "patient must have a column for every subgroup variable in summary.",
        x = paste0(
          "Missing subgroup column(s): ",
          paste0("\"", missing, "\"", collapse = ", "),
          "."
        ),
        i = paste(
          "patient is usually the output of",
          "group_sample_unconstrained() or group_sample_constrained()."
        )
      )
    )
  }
  for (variable in variables) {
    values <- patient[[variable]]
    levels <- summary$level[summary$variable == variable]
    if (anyNA(values) || !all(as.character(values) %in% levels)) {
      throw_error(
        paste0(
          "column \"",
          variable,
          "\" of patient must contain only levels of subgroup variable \"",
          variable,
          "\" in summary, with no missing values."
        )
      )
    }
  }
  if (
    !is.numeric(iterations) ||
      length(iterations) != 1L ||
      is.na(iterations) ||
      !is.finite(iterations) ||
      iterations < 0 ||
      iterations != round(iterations)
  ) {
    throw_error("iterations must be a single non-negative whole number.")
  }
  if (
    !is.numeric(tolerance) ||
      length(tolerance) != 1L ||
      is.na(tolerance) ||
      tolerance < 0
  ) {
    throw_error("tolerance must be a single non-negative number.")
  }
  invisible()
}

#' @title Warn if the optimization did not reach the tolerance
#' @noRd
#' @return `NULL` (invisibly). Throws a warning if any element of `loss`
#'   is not below `tolerance`.
#' @param loss Named numeric vector with the final loss of each subgroup
#'   variable.
#' @param iterations Maximum number of iterations per subgroup variable.
#' @param tolerance Loss tolerance.
#' @examples
#'   warn_optimize_tolerance(
#'     loss = c(age = 1e-8, sex = 0.02),
#'     iterations = 100,
#'     tolerance = 1e-6
#'   )
warn_optimize_tolerance <- function(loss, iterations, tolerance) {
  unconverged <- loss[!(loss < tolerance)]
  if (length(unconverged)) {
    throw_warning(
      c(
        paste0(
          "The loss did not fall below tolerance = ",
          format(tolerance),
          " within iterations = ",
          format(iterations),
          " for ",
          length(unconverged),
          " subgroup variable(s)."
        ),
        x = paste0(
          "Final loss by subgroup variable: ",
          paste0(
            "\"",
            names(unconverged),
            "\" = ",
            format(unconverged, digits = 3L),
            collapse = ", "
          ),
          "."
        ),
        i = paste(
          "Increase iterations, rerun from freshly sampled subgroup levels,",
          "or increase tolerance."
        )
      )
    )
  }
  invisible()
}
