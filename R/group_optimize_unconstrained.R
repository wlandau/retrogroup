#' @title Unconstrained subgroup optimization
#' @export
#' @family optimization
#' @description Refine sampled subgroup levels so that the subgroup
#'   statistics observed in the data match the published subgroup
#'   statistics in `summary`. Patient counts are not preserved.
#' @details `group_optimize_unconstrained()` implements the
#'   "Optimization for unconstrained sampling" section of
#'   `vignette("optimization", package = "retrogroup")`, separately for each
#'   subgroup variable. The loss \eqn{\Delta(x)} is the Euclidean
#'   distance between the observed and assumed hazard ratios of the levels
#'   of the subgroup variable, as defined in the "Loss function" section of
#'   that vignette. The observed hazard ratios come from a Cox proportional
#'   hazards model of each subgroup level with treatment as the only
#'   covariate, as in [survival::coxph()]. If the model cannot be fit for
#'   some subgroup level, for example because the level is empty or one of
#'   its arms has no events, the loss is infinite.
#'
#'   Each iteration draws a patient \eqn{i} uniformly at random and moves
#'   that patient to the subgroup level \eqn{k} that minimizes
#'   \eqn{\Delta(x^{(i \to k)})}. The patient stays put if the current
#'   level ties for the minimum, so the loss never increases. The search
#'   for each subgroup variable stops when the loss falls below `tolerance`
#'   or after `iterations` iterations, whichever comes first. If the loss
#'   of any subgroup variable is still at least `tolerance` after
#'   `iterations` iterations, `group_optimize_unconstrained()` throws a
#'   warning of class `"retrogroup_warning"` that lists those subgroup
#'   variables and their final losses.
#'
#'   The search is greedy, so it may settle in a local minimum. The number
#'   of patients in each subgroup level is not guaranteed to equal the
#'   published number of patients in `summary`. Use
#'   [group_optimize_constrained()] on the output of
#'   [group_sample_constrained()] to preserve the published counts.
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
#'   iterations for each subgroup variable. Each iteration draws one
#'   patient and fits one Cox model per subgroup level.
#' @param tolerance Non-negative number. The search for each subgroup
#'   variable stops when the loss \eqn{\Delta(x)} falls below `tolerance`.
#'   The loss is on the hazard ratio scale, so the default of `1e-3` stops
#'   when the observed and assumed hazard ratios agree to roughly `1e-3`.
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
  iterations = 1000L,
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
#' @param iterations Maximum number of iterations.
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
#'   sex <- summary[summary$variable == "sex", ]
#'   x <- sample.int(2L, nrow(patient), replace = TRUE)
#'   out <- optimize_unconstrained_variable(
#'     x = x,
#'     z = patient$treatment,
#'     t = patient$time,
#'     d = patient$event,
#'     hazard_ratio = sex$hazard_ratio,
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
  iterations,
  tolerance
) {
  n <- length(x)
  K <- length(hazard_ratio)
  # Loss term of subgroup level k if its patients are the ones in members.
  loss_term <- function(members, k) {
    fit <- analysis_fit(z[members], t[members], d[members])
    optimize_loss_term(fit, hazard_ratio[k])
  }
  # Loss terms of x, one per subgroup level, and the loss Delta(x).
  terms <- optimize_loss_terms(x, z, t, d, hazard_ratio)
  loss <- sqrt(sum(terms))
  iteration <- 0L
  while (iteration < iterations && !(loss < tolerance)) {
    iteration <- iteration + 1L
    # Step 1: draw a patient i uniformly at random from 1, ..., n.
    i <- sample.int(n, 1L)
    x_i <- x[i]
    # Step 2: compute the counterfactual loss Delta(x^(i -> k)) for every
    # subgroup level k = 1, ..., K. A flip from level x_i to level k changes
    # the patients in those two levels only, so every other loss term
    # carries over from x. Column k of counterfactual_terms holds the loss
    # terms of x^(i -> k), and column x_i holds those of x itself.
    # Each iteration needs K model fits: one for level x_i without
    # patient i, and one for each other level k with patient i.
    members <- x == x_i
    members[i] <- FALSE
    term_without_i <- loss_term(members, x_i)
    counterfactual_terms <- vapply(
      seq_len(K),
      \(k) {
        if (k == x_i) {
          return(terms)
        }
        members <- x == k
        members[i] <- TRUE
        flipped <- terms
        flipped[x_i] <- term_without_i
        flipped[k] <- loss_term(members, k)
        flipped
      },
      numeric(K)
    )
    counterfactual_loss <- sqrt(colSums(counterfactual_terms))
    # Step 3: update x <- x^(i -> k_star), where k_star minimizes the
    # counterfactual loss. The candidates include x^(i -> x_i) = x, so the
    # loss never increases. If the current level x_i ties for the minimum,
    # the patient stays put.
    k_star <- which.min(counterfactual_loss)
    if (counterfactual_loss[k_star] < loss) {
      x[i] <- k_star
      terms <- counterfactual_terms[, k_star]
      loss <- counterfactual_loss[k_star]
    }
  }
  list(x = x, loss = unname(loss))
}
