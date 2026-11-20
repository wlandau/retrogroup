#' @title Unconstrained subgroup optimization
#' @export
#' @family optimization
#' @description Refine sampled subgroup levels so that the subgroup
#'   statistics observed in the data match the assumed subgroup
#'   statistics in `summary`. Patient counts are not preserved.
#' @details `group_optimize_unconstrained()` implements the
#'   "Optimization for unconstrained sampling" section of
#'   `vignette("optimization", package = "retrogroup")`, separately for each
#'   subgroup variable. The loss \eqn{\Delta(x)} is the Euclidean
#'   distance between the observed and assumed hazard ratios and event
#'   fractions of the levels of the subgroup variable, as defined in the
#'   "Loss function" section of that vignette. The observed event fraction
#'   of a level is its observed number of events divided by its assumed
#'   number of patients in `summary$patients`, and the assumed event
#'   fraction is `summary$events / summary$patients`. The denominator only
#'   puts the event counts on a scale comparable to the hazard ratios, so
#'   the event fraction terms are a soft constraint that steers the event
#'   counts toward `summary$events` without forcing them to match. The event
#'   fraction terms do not depend on the observed number of patients in each
#'   level. The observed hazard ratios come from a Cox proportional hazards
#'   model of each subgroup level with treatment as the only
#'   covariate, as in [survival::coxph()]. If the model cannot be fit for
#'   some subgroup level, for example because the level is empty or one of
#'   its arms has no events, the loss is infinite.
#'
#'   Each iteration takes a flip step with probability 1/2 and a swap step
#'   otherwise, as defined in the "Moves" section of the vignette. If no
#'   swap is possible, the iteration takes a flip step instead. A flip step
#'   draws a patient \eqn{i} uniformly at random and moves that patient to
#'   the subgroup level \eqn{k} that minimizes \eqn{\Delta(x^{(i \to k)})},
#'   or leaves the patient in place if the current level ties for the
#'   minimum. A swap step draws a pair of patients \eqn{(i, j)} uniformly
#'   at random from the pairs with \eqn{x_i \ne x_j} and exchanges their
#'   subgroup levels if the exchange reduces the loss. Neither step ever
#'   increases the loss. The search
#'   for each subgroup variable stops when the loss falls below `tolerance`
#'   or after `iterations` iterations,
#'   whichever comes first. If the loss of any subgroup variable is still
#'   at least `tolerance` after `iterations` iterations,
#'   `group_optimize_unconstrained()` throws a warning of class
#'   `"retrogroup_warning"` that lists those subgroup variables and their
#'   final losses.
#'
#'   The search is greedy, so it may settle in a local minimum. The number
#'   of patients and the number of events in each subgroup level are not
#'   guaranteed to equal the assumed counts in `summary`. Use
#'   [group_optimize_constrained()] on the output of
#'   [group_sample_constrained()] to preserve the assumed number of
#'   patients.
#'
#'   The timing of events within each subgroup level is not reproduced
#'   from the original data. We make no claims about the prognostic
#'   information of optimized subgroups, but
#'   `vignette("prognostic", package = "retrogroup")` examines it
#'   empirically. The statistician who uses the optimized data in an
#'   analysis, for example for covariate adjustment, is responsible for
#'   the consequences.
#'
#'   In addition to the requirements of [group_analysis_probability()], `patient`
#'   must have one column per subgroup variable in `summary`, named after
#'   the subgroup variable, with no missing values and only the levels of
#'   that subgroup variable in `summary`. The output of
#'   [group_sample_unconstrained()] satisfies these requirements. These
#'   columns are the starting subgroup levels, and all other columns of
#'   `patient` are returned unchanged.
#'
#'   Optimization uses R's random number generator. Call [set.seed()]
#'   beforehand for reproducible results.
#' @return A tibble with the same rows and columns as `patient`, in the
#'   same order. Each subgroup column is a factor of the optimized
#'   levels whose levels are those of the subgroup variable in
#'   `summary$level`, in the same order. All other columns are unchanged.
#' @inheritParams group_analysis_probability
#' @param iterations Non-negative whole number, the maximum number of
#'   iterations for each subgroup variable. A flip step fits one Cox model
#'   per subgroup level, and a swap step fits two Cox models.
#' @param tolerance Non-negative number on the linear (natural) scale of
#'   the hazard ratio, not the log scale, and on the scale of the event
#'   fractions. The search for each subgroup variable stops when the loss
#'   \eqn{\Delta(x)} falls below `tolerance`. Because \eqn{\Delta(x)} is
#'   the Euclidean distance between the observed and assumed hazard ratios
#'   and event fractions themselves, a loss below `tolerance` guarantees
#'   that the observed hazard ratio and event fraction of every subgroup
#'   level are within `tolerance` of their assumed values.
#'   The default of `1e-3` is much finer than the two decimal places
#'   hazard ratios are usually published with. However, the event fraction
#'   of level \eqn{k} moves in steps of \eqn{1 / n_k}, and the hazard ratio
#'   and event fraction terms may pull in different directions, so a small
#'   `tolerance` may be unreachable, in which case the search runs for all
#'   `iterations`. Larger values such as `1e-2` give faster but coarser
#'   fits.
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
#'   summary <- summary[summary$variable == "age", ]
#'   set.seed(0)
#'   sample <- group_sample_unconstrained(patient, summary)
#'   optimized <- group_optimize_unconstrained(
#'     sample,
#'     summary,
#'     tolerance = 1e-2,
#'     iterations = 100L
#'   )
#'   summary
#'   group_analysis_summary(sample, "age")
#'   group_analysis_summary(optimized, "age")
group_optimize_unconstrained <- function(
  patient,
  summary,
  iterations = 10000L,
  tolerance = 1e-3
) {
  assert_optimize_arguments(patient, summary, iterations, tolerance)
  out <- tibble::as_tibble(patient)
  variables <- unique(summary$variable)
  loss <- numeric(length(variables))
  names(loss) <- variables
  for (variable in variables) {
    rows <- summary$variable == variable
    levels <- summary$level[rows]
    result <- optimize_unconstrained_variable(
      x = match(as.character(patient[[variable]]), levels),
      z = patient$treatment,
      t = patient$time,
      d = patient$event,
      hazard_ratio = summary$hazard_ratio[rows],
      patients = summary$patients[rows],
      events = summary$events[rows],
      iterations = iterations,
      tolerance = tolerance
    )
    out[[variable]] <- factor(levels[result$x], levels = levels)
    loss[variable] <- result$loss
  }
  warn_optimize_tolerance(loss, iterations, tolerance)
  out
}

#' @title Unconstrained optimization of one subgroup variable
#' @noRd
#' @return A list with elements:
#'   * `x`: integer vector with one element per patient, the index of the
#'     optimized subgroup level of each patient.
#'   * `loss`: numeric scalar, the final loss \eqn{\Delta(x)}.
#' @param x Integer vector with one element per patient: the index of
#'   the starting subgroup level \eqn{x_i} of each patient.
#' @param z Integer vector, treatment indicator of each patient.
#' @param t Numeric vector, time of each patient.
#' @param d Integer vector, event indicator of each patient.
#' @param hazard_ratio Numeric vector, assumed hazard ratio of each level.
#' @param patients Integer vector, assumed number of patients of each level.
#' @param events Integer vector, assumed number of events of each level.
#' @param iterations Maximum number of iterations.
#' @param tolerance Loss tolerance on the scale of the hazard ratios and
#'   event fractions.
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
#'   out <- optimize_unconstrained_variable(
#'     x = x,
#'     z = patient$treatment,
#'     t = patient$time,
#'     d = patient$event,
#'     hazard_ratio = sex$hazard_ratio,
#'     patients = sex$patients,
#'     events = sex$events,
#'     iterations = 20,
#'     tolerance = 1e-3
#'   )
#'   table(out$x)
#'   out$loss
optimize_unconstrained_variable <- function(
  x,
  z,
  t,
  d,
  hazard_ratio,
  patients,
  events,
  iterations,
  tolerance
) {
  problem <- optimize_problem(z, t, d, hazard_ratio, patients, events)
  # Loss terms of x, one per subgroup level, and the loss Delta(x).
  terms <- optimize_loss_terms(x, problem)
  loss <- sqrt(sum(terms))
  iteration <- 0L
  while (iteration < iterations && !(loss < tolerance)) {
    iteration <- iteration + 1L
    # Take a flip step with probability 1/2 and a swap step otherwise.
    # If no swap is possible, take a flip step instead.
    swap <- stats::runif(1L) < 0.5 && optimize_swap_possible(x)
    move <- if (swap) {
      optimize_swap(x, terms, loss, problem)
    } else {
      optimize_flip(x, terms, loss, problem)
    }
    x[move$index] <- move$level
    terms <- move$terms
    loss <- move$loss
  }
  list(x = x, loss = unname(loss))
}
