#' @title Constrained subgroup optimization
#' @export
#' @family optimization
#' @description Refine sampled subgroup levels so that the subgroup
#'   statistics observed in the data match the published subgroup
#'   statistics in `summary`, while preserving the published number of
#'   patients in each subgroup level.
#' @details `group_optimize_constrained()` implements the
#'   "Optimization for constrained sampling" section of
#'   `vignette("optimization", package = "retrogroup")`, separately for each
#'   subgroup variable. The loss \eqn{\Delta(x)} is the same as in
#'   [group_optimize_unconstrained()]: the squared Euclidean distance
#'   between the observed and published log hazard ratios, plus the squared
#'   distance between the observed and published standard errors for the
#'   subgroup levels whose `standard_error` in `summary` is not `NA`. The
#'   observed statistics come from a Cox proportional hazards model of each
#'   subgroup level with treatment as the only covariate, as in
#'   [survival::coxph()]. If the model cannot be fit for some subgroup
#'   level, for example because one of its arms has no events, the loss is
#'   infinite.
#'
#'   Each iteration draws a pair of patients \eqn{(i, j)} uniformly at
#'   random from the pairs with \eqn{x_i \ne x_j} and swaps their subgroup
#'   levels if the swap reduces the loss. Otherwise, the subgroup levels
#'   stay the same, so the loss never increases. The search for each
#'   subgroup variable stops when the loss falls below `tolerance` or after
#'   `iterations` iterations, whichever comes first. If the loss of any
#'   subgroup variable is still at least `tolerance` after `iterations`
#'   iterations, `group_optimize_constrained()` throws a warning of class
#'   `"retrogroup_warning"` that lists those subgroup variables and their
#'   final losses.
#'
#'   The search is greedy, so it may settle in a local minimum. Swaps
#'   preserve the number of patients in each subgroup level, so the
#'   optimized counts equal the published number of patients in `summary`.
#'   Use [group_optimize_unconstrained()] on the output of
#'   [group_sample_unconstrained()] to optimize without the count
#'   constraint.
#'
#'   In addition to the requirements of [group_analysis_probability()],
#'   `patient` must have one column per subgroup variable in `summary`,
#'   named after the subgroup variable, with no missing values and only the
#'   levels of that subgroup variable in `summary`. The number of patients
#'   in each subgroup level of these columns must equal the published
#'   number of patients in `summary`. The output of
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
#'   level equals the published number of patients in `summary`. All other
#'   columns are unchanged.
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
#'     log_hazard_ratio = "log_hazard_ratio",
#'     standard_error = "standard_error"
#'   )
#'   summary <- summary[summary$variable == "age", ]
#'   set.seed(0)
#'   sample <- group_sample_constrained(patient, summary)
#'   optimized <- group_optimize_constrained(
#'     sample,
#'     summary,
#'     tolerance = 1e-4,
#'     iterations = 500L
#'   )
#'   summary
#'   group_analysis_summary(sample, "age")
#'   group_analysis_summary(optimized, "age")
group_optimize_constrained <- function(
  patient,
  summary,
  iterations = 1000L,
  tolerance = 1e-6
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
      beta = summary$log_hazard_ratio[rows],
      sigma = summary$standard_error[rows],
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
#'     beta = sex$log_hazard_ratio,
#'     sigma = sex$standard_error,
#'     iterations = 20,
#'     tolerance = 1e-6
#'   )
#'   table(out$x)
#'   out$loss
optimize_constrained_variable <- function(
  x,
  z,
  t,
  d,
  beta,
  sigma,
  iterations,
  tolerance
) {
  K <- length(beta)
  # Loss term of subgroup level k if its patients are the ones in members.
  loss_term <- function(members, k) {
    fit <- analysis_fit(z[members], t[members], d[members])
    optimize_loss_term(fit, beta[k], sigma[k])
  }
  # Loss terms of x, one per subgroup level, and the loss Delta(x).
  terms <- optimize_loss_terms(x, z, t, d, beta, sigma)
  loss <- sum(terms)
  # Swaps preserve the number of patients n_k in each subgroup level k,
  # so the number of pairs (i, j) with x_i = k and x_j = l is n_k * n_l
  # in every iteration. Each unordered pair of distinct levels k < l is
  # one row of level_pairs, with weight n_k * n_l.
  n <- tabulate(x, nbins = K)
  level_pairs <- which(upper.tri(diag(K)), arr.ind = TRUE)
  weights <- n[level_pairs[, 1L]] * n[level_pairs[, 2L]]
  # Without any pair of patients in different levels, no swap is possible.
  if (!any(weights > 0)) {
    return(list(x = x, loss = unname(loss)))
  }
  iteration <- 0L
  while (iteration < iterations && !(loss < tolerance)) {
    iteration <- iteration + 1L
    # Step 1: draw a pair of patients (i, j) uniformly at random from the
    # pairs with x_i != x_j. First draw the pair of levels (x_i, x_j) with
    # probability proportional to the number of such pairs n_k * n_l, then
    # draw i and j uniformly from the patients in those levels.
    pair <- level_pairs[sample.int(nrow(level_pairs), 1L, prob = weights), ]
    x_i <- pair[[1L]]
    x_j <- pair[[2L]]
    patients_i <- which(x == x_i)
    patients_j <- which(x == x_j)
    i <- patients_i[sample.int(length(patients_i), 1L)]
    j <- patients_j[sample.int(length(patients_j), 1L)]
    # Step 2: compute the loss Delta(x^(i <-> j)) of the swap. A swap changes
    # the patients in levels x_i and x_j only, so every other loss term
    # carries over from x. Each iteration needs two model fits: one for
    # level x_i with j in place of i, and one for level x_j with i in place
    # of j.
    members_i <- x == x_i
    members_i[c(i, j)] <- c(FALSE, TRUE)
    members_j <- x == x_j
    members_j[c(i, j)] <- c(TRUE, FALSE)
    swapped_terms <- terms
    swapped_terms[x_i] <- loss_term(members_i, x_i)
    swapped_terms[x_j] <- loss_term(members_j, x_j)
    swapped_loss <- sum(swapped_terms)
    # Step 3: update x <- x^(i <-> j) if the swap reduces the loss, and
    # otherwise leave x unchanged, so the loss never increases.
    if (swapped_loss < loss) {
      x[c(i, j)] <- c(x_j, x_i)
      terms <- swapped_terms
      loss <- swapped_loss
    }
  }
  list(x = x, loss = unname(loss))
}

#' @title Validate the subgroup counts of constrained optimization
#' @noRd
#' @return `NULL` (invisibly). Throws an error if the number of patients
#'   in any subgroup level of `patient` differs from the published number
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
          "must equal the published number of patients in summary."
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
