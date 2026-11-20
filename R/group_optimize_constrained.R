#' @title Constrained subgroup optimization
#' @export
#' @family optimization
#' @description Refine sampled subgroup levels so that the subgroup
#'   statistics observed in the data match the assumed subgroup
#'   statistics in `summary`, while preserving the assumed number of
#'   patients in each subgroup level.
#' @details `group_optimize_constrained()` implements the
#'   "Optimization for constrained sampling" section of
#'   `vignette("optimization", package = "retrogroup")`, separately for each
#'   subgroup variable. The loss \eqn{\Delta(x)} is the same as in
#'   [group_optimize_unconstrained()]: the Euclidean distance between the
#'   observed and assumed hazard ratios and event fractions. The observed hazard ratios come
#'   from a Cox proportional hazards model of each
#'   subgroup level with treatment as the only covariate, as in
#'   [survival::coxph()]. If the model cannot be fit for some subgroup
#'   level, for example because one of its arms has no events, the loss is
#'   infinite.
#'
#'   Each iteration takes a swap step, as defined in the "Moves" section of
#'   the vignette. A swap step draws a pair of patients \eqn{(i, j)}
#'   uniformly at random from the pairs with \eqn{x_i \ne x_j} and
#'   exchanges their subgroup levels if the exchange reduces the loss.
#'   Otherwise, the subgroup levels stay the same, so the loss never
#'   increases. If no such pair exists,
#'   no swap is possible and the subgroup levels are returned unchanged.
#'   The search for each subgroup variable stops when the loss falls below
#'   `tolerance` or after `iterations` iterations,
#'   whichever comes first. If the loss of any subgroup variable is still
#'   at least `tolerance` after `iterations` iterations,
#'   `group_optimize_constrained()` throws a warning of class
#'   `"retrogroup_warning"` that lists those subgroup variables and their
#'   final losses.
#'
#'   The search is greedy, so it may settle in a local minimum. Swaps
#'   preserve the number of patients in each subgroup level, so the
#'   optimized counts equal `summary$patients`.
#'   Use [group_optimize_unconstrained()] on the output of
#'   [group_sample_unconstrained()] to optimize without the count
#'   constraint.
#'
#'   The event fraction terms of the loss are only a soft constraint, so
#'   the number of events in each subgroup level is not guaranteed to
#'   equal `summary$events`. In addition, the timing of events within each
#'   subgroup level is not reproduced from the original data. We make no
#'   claims about the prognostic information of optimized subgroups, but
#'   `vignette("prognostic", package = "retrogroup")` examines it
#'   empirically. The statistician who uses the optimized data in an analysis, for example
#'   for covariate adjustment, is responsible for the consequences.
#'
#'   In addition to the requirements of [group_analysis_probability()],
#'   `patient` must have one column per subgroup variable in `summary`,
#'   named after the subgroup variable, with no missing values and only the
#'   levels of that subgroup variable in `summary`. The number of patients
#'   in each subgroup level of these columns must equal
#'   `summary$patients`. The output of
#'   [group_sample_constrained()] satisfies these requirements. These
#'   columns are the starting subgroup levels, and all other columns of
#'   `patient` are returned unchanged.
#'
#'   Optimization uses R's random number generator. Call [set.seed()]
#'   beforehand for reproducible results.
#' @return A tibble with the same rows and columns as `patient`, in the
#'   same order. Each subgroup column is a factor of the optimized
#'   levels whose levels are those of the subgroup variable in
#'   `summary$level`, in the same order. The number of patients in each
#'   level equals `summary$patients`. All other columns are unchanged.
#' @inheritParams group_optimize_unconstrained
#' @param iterations Non-negative whole number, the maximum number of
#'   iterations for each subgroup variable. Each iteration draws one pair
#'   of patients and fits two Cox models.
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
#'   sample <- group_sample_constrained(patient, summary)
#'   optimized <- group_optimize_constrained(
#'     sample,
#'     summary,
#'     tolerance = 1e-2,
#'     iterations = 500L
#'   )
#'   summary
#'   group_analysis_summary(sample, "age")
#'   group_analysis_summary(optimized, "age")
group_optimize_constrained <- function(
  patient,
  summary,
  iterations = 10000L,
  tolerance = 1e-3
) {
  assert_optimize_arguments(patient, summary, iterations, tolerance)
  assert_optimize_counts(patient, summary)
  out <- tibble::as_tibble(patient)
  variables <- unique(summary$variable)
  loss <- numeric(length(variables))
  names(loss) <- variables
  for (variable in variables) {
    rows <- summary$variable == variable
    levels <- summary$level[rows]
    result <- optimize_constrained_variable(
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

#' @title Constrained optimization of one subgroup variable
#' @noRd
#' @return A list with elements:
#'   * `x`: integer vector with one element per patient, the index of the
#'     optimized subgroup level of each patient. Each level has the same
#'     number of patients as in the starting `x`.
#'   * `loss`: numeric scalar, the final loss \eqn{\Delta(x)}.
#' @inheritParams optimize_unconstrained_variable
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
#'   x <- sample(rep(1:2, times = sex$patients))
#'   out <- optimize_constrained_variable(
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
optimize_constrained_variable <- function(
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
  # Swaps preserve the subgroup counts, so if no swap is possible now,
  # no swap is ever possible.
  if (!optimize_swap_possible(x)) {
    return(list(x = x, loss = unname(loss)))
  }
  iteration <- 0L
  while (iteration < iterations && !(loss < tolerance)) {
    iteration <- iteration + 1L
    # Every iteration takes a swap step.
    move <- optimize_swap(x, terms, loss, problem)
    x[move$index] <- move$level
    terms <- move$terms
    loss <- move$loss
  }
  list(x = x, loss = unname(loss))
}

#' @title Validate the subgroup counts of constrained optimization
#' @noRd
#' @return `NULL` (invisibly). Throws an error if the number of patients
#'   in any subgroup level of `patient` differs from the assumed number
#'   of patients in `summary`.
#' @param patient Patient-level data with subgroup columns, the output of
#'   [group_sample_constrained()]. Must already satisfy
#'   `assert_optimize_arguments()`.
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
#'   summary <- group_example_summary()
#'   assert_optimize_counts(group_sample_constrained(patient, summary), summary)
#'   sampled <- group_sample_unconstrained(patient, summary)
#'   try(assert_optimize_counts(sampled, summary))
assert_optimize_counts <- function(patient, summary) {
  variables <- unique(summary$variable)
  mismatched <- Filter(
    \(variable) {
      rows <- summary$variable == variable
      levels <- summary$level[rows]
      observed <- tabulate(
        match(as.character(patient[[variable]]), levels),
        nbins = length(levels)
      )
      any(observed != summary$patients[rows])
    },
    variables
  )
  if (length(mismatched)) {
    throw_error(
      c(
        paste(
          "The number of patients in each subgroup level of patient",
          "must equal the assumed number of patients in summary."
        ),
        x = paste0(
          "Mismatched subgroup variable(s): ",
          paste0("\"", mismatched, "\"", collapse = ", "),
          "."
        ),
        i = paste(
          "patient is usually the output of group_sample_constrained().",
          "Use group_optimize_unconstrained() for the output of",
          "group_sample_unconstrained()."
        )
      )
    )
  }
  invisible()
}
