#' @title Constrained subgroup optimization
#' @export
#' @family optimization
#' @description Refine sampled subgroup levels so that the subgroup
#'   statistics observed in the data match the published subgroup
#'   statistics in `summary`, while preserving the assumed number of
#'   patients and number of events in each subgroup level.
#' @details `group_optimize_constrained()` implements the
#'   "Optimization for constrained sampling" section of
#'   `vignette("optimization", package = "retrogroup")`, separately for each
#'   subgroup variable. The loss \eqn{\Delta(x)} is the same as in
#'   [group_optimize_unconstrained()]: the Euclidean distance between the
#'   observed and assumed hazard ratios. The observed hazard ratios come
#'   from a Cox proportional hazards model of each
#'   subgroup level with treatment as the only covariate, as in
#'   [survival::coxph()]. If the model cannot be fit for some subgroup
#'   level, for example because one of its arms has no events, the loss is
#'   infinite.
#'
#'   Each iteration draws a pair of patients \eqn{(i, j)} uniformly at
#'   random from the pairs with \eqn{x_i \ne x_j} and \eqn{d_i = d_j},
#'   where \eqn{d_i} is the event indicator of patient \eqn{i}, and swaps
#'   their subgroup levels if the swap reduces the loss. Otherwise, the
#'   subgroup levels stay the same, so the loss never increases. If no
#'   such pair exists, no swap is possible and the subgroup levels are
#'   returned unchanged. The search for each
#'   subgroup variable stops when the loss falls below `tolerance` (on the
#'   linear hazard ratio scale) or after `iterations` iterations,
#'   whichever comes first. If the loss of any subgroup variable is still
#'   at least `tolerance` after `iterations` iterations,
#'   `group_optimize_constrained()` throws a warning of class
#'   `"retrogroup_warning"` that lists those subgroup variables and their
#'   final losses.
#'
#'   The search is greedy, so it may settle in a local minimum. Swaps
#'   exchange two patients with the same event status, so they preserve
#'   both the number of patients and the number of events in each
#'   subgroup level. The optimized counts equal `summary$patients` and
#'   `summary$events`.
#'   Use [group_optimize_unconstrained()] on the output of
#'   [group_sample_unconstrained()] to optimize without the count
#'   constraint.
#'
#'   The loss only targets the hazard ratios of treatment versus control.
#'   The preserved event counts carry some general prognostic information,
#'   but the timing of events within each subgroup level is not reproduced
#'   from the original data, so some prognostic information is lost. The
#'   statistician who uses the optimized data in an analysis, for example
#'   for covariate adjustment, is responsible for the consequences.
#'
#'   In addition to the requirements of [group_analysis_probability()],
#'   `patient` must have one column per subgroup variable in `summary`,
#'   named after the subgroup variable, with no missing values and only the
#'   levels of that subgroup variable in `summary`. The number of patients
#'   and the number of events in each subgroup level of these columns must
#'   equal `summary$patients` and `summary$events`. The output of
#'   [group_sample_constrained()] satisfies these requirements. These
#'   columns are the starting subgroup levels, and all other columns of
#'   `patient` are returned unchanged.
#'
#'   Optimization uses R's random number generator. Call [set.seed()]
#'   beforehand for reproducible results.
#' @return A tibble with the same rows and columns as `patient`, in the
#'   same order. Each subgroup column is a factor of the optimized
#'   levels whose levels are those of the subgroup variable in
#'   `summary$level`, in the same order. The number of patients and the
#'   number of events in each level equal `summary$patients` and
#'   `summary$events`. All other columns are unchanged.
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
  iterations = 1000L,
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
#'     number of patients and the same number of patients with `d == 1`
#'     as in the starting `x`.
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
  iterations,
  tolerance
) {
  K <- length(hazard_ratio)
  # Loss term of subgroup level k if its patients are the ones in members.
  loss_term <- function(members, k) {
    fit <- analysis_fit(z[members], t[members], d[members])
    optimize_loss_term(fit, hazard_ratio[k])
  }
  # Loss terms of x, one per subgroup level, and the loss Delta(x).
  terms <- optimize_loss_terms(x, z, t, d, hazard_ratio)
  loss <- sqrt(sum(terms))
  # Swaps only exchange patients with the same event status, so they
  # preserve the number n_ks of patients with event status s in each
  # subgroup level k. The number of pairs (i, j) with x_i = k, x_j = l, and
  # d_i = d_j = s is therefore n_ks * n_ls in every iteration. Each
  # combination of an unordered pair of distinct levels k < l and an event
  # status s is one row of candidates, with weight n_ks * n_ls.
  counts <- cbind(
    tabulate(x[d == 0L], nbins = K),
    tabulate(x[d == 1L], nbins = K)
  )
  level_pairs <- which(upper.tri(diag(K)), arr.ind = TRUE)
  candidates <- cbind(
    k = rep(level_pairs[, 1L], times = 2L),
    l = rep(level_pairs[, 2L], times = 2L),
    s = rep(c(0L, 1L), each = nrow(level_pairs))
  )
  status <- candidates[, "s"] + 1L
  weights <- counts[cbind(candidates[, "k"], status)] *
    counts[cbind(candidates[, "l"], status)]
  # Without any pair of patients in different levels with the same event
  # status, no swap is possible.
  if (!any(weights > 0)) {
    return(list(x = x, loss = unname(loss)))
  }
  iteration <- 0L
  while (iteration < iterations && !(loss < tolerance)) {
    iteration <- iteration + 1L
    # Step 1: draw a pair of patients (i, j) uniformly at random from the
    # pairs with x_i != x_j and d_i = d_j. First draw the levels (x_i, x_j)
    # and the event status s with probability proportional to the number
    # of such pairs n_ks * n_ls, then draw i and j uniformly from the
    # patients with event status s in those levels.
    candidate <- candidates[sample.int(nrow(candidates), 1L, prob = weights), ]
    x_i <- candidate[["k"]]
    x_j <- candidate[["l"]]
    s <- candidate[["s"]]
    patients_i <- which(x == x_i & d == s)
    patients_j <- which(x == x_j & d == s)
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
    swapped_loss <- sqrt(sum(swapped_terms))
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
#'   or the number of events in any subgroup level of `patient` differs
#'   from the assumed counts in `summary`.
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
  mismatch <- vapply(
    variables,
    \(variable) {
      rows <- summary$variable == variable
      levels <- summary$level[rows]
      index <- match(as.character(patient[[variable]]), levels)
      patients <- tabulate(index, nbins = length(levels))
      events <- tabulate(index[patient$event == 1L], nbins = length(levels))
      counts <- c(
        patients = any(patients != summary$patients[rows]),
        events = any(events != summary$events[rows])
      )
      paste(names(counts)[counts], collapse = ", ")
    },
    character(1L)
  )
  mismatch <- mismatch[nzchar(mismatch)]
  if (length(mismatch)) {
    throw_error(
      c(
        paste(
          "The number of patients and the number of events in each",
          "subgroup level of patient must equal the assumed counts",
          "in summary."
        ),
        x = paste0(
          "Mismatched subgroup variable(s): ",
          paste0("\"", names(mismatch), "\" (", mismatch, ")", collapse = ", "),
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
