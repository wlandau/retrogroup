#' @title Loss term of one subgroup level
#' @noRd
#' @return Numeric scalar, the contribution of the level to the squared
#'   loss \eqn{\Delta(x)^2}:
#'   \eqn{(\exp(\hat{\beta}_k(x)) - \exp(\beta_k))^2 + (\hat{f}_k(x) - f_k)^2},
#'   where \eqn{\hat{f}_k(x)} is the observed number of events divided by
#'   the assumed number of patients \eqn{n_k}, and \eqn{f_k = e_k / n_k}.
#'   `Inf` if the model could not be fit.
#' @param fit Numeric scalar from `analysis_fit()`, the observed log hazard
#'   ratio \eqn{\hat{\beta}_k(x)} of the level, or `NA`.
#' @param observed Integer scalar, observed number of events in the level.
#' @param hazard_ratio Numeric scalar, assumed hazard ratio
#'   \eqn{\exp(\beta_k)} of the level.
#' @param patients Integer scalar, assumed number of patients \eqn{n_k}
#'   of the level.
#' @param events Integer scalar, assumed number of events \eqn{e_k} of
#'   the level.
#' @examples
#'   optimize_loss_term(-0.3, 40L, hazard_ratio = 0.8, patients = 100L, events = 45L)
#'   optimize_loss_term(NA_real_, 40L, hazard_ratio = 0.8, patients = 100L, events = 45L)
optimize_loss_term <- function(fit, observed, hazard_ratio, patients, events) {
  if (is.na(fit)) {
    return(Inf)
  }
  # Both event fractions divide by the assumed number of patients n_k,
  # which only puts the event counts on the scale of the hazard ratios.
  (exp(fit) - hazard_ratio)^2 + ((observed - events) / patients)^2
}

#' @title Loss term of one subgroup level with the given patients
#' @noRd
#' @return Numeric scalar, the loss term of subgroup level `k` if its
#'   patients are those in `members`, from `optimize_loss_term()`.
#'   Needs one model fit.
#' @param members Logical vector with one element per patient, `TRUE` for
#'   the patients in subgroup level `k`.
#' @param k Integer scalar, index of the subgroup level.
#' @param problem List of the fixed inputs of the search, from
#'   `optimize_problem()`.
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
#'   problem <- optimize_problem(
#'     z = patient$treatment,
#'     t = patient$time,
#'     d = patient$event,
#'     hazard_ratio = sex$hazard_ratio,
#'     patients = sex$patients,
#'     events = sex$events
#'   )
#'   x <- sample.int(2L, nrow(patient), replace = TRUE)
#'   optimize_level_term(x == 1L, 1L, problem)
optimize_level_term <- function(members, k, problem) {
  d <- problem$d[members]
  fit <- analysis_fit(problem$z[members], problem$t[members], d)
  optimize_loss_term(
    fit = fit,
    observed = sum(d),
    hazard_ratio = problem$hazard_ratio[k],
    patients = problem$patients[k],
    events = problem$events[k]
  )
}

#' @title Fixed inputs of the search for one subgroup variable
#' @noRd
#' @return A list with the arguments as elements. The list refers to the
#'   argument vectors without copying them.
#' @param z Integer vector, treatment indicator of each patient.
#' @param t Numeric vector, time of each patient.
#' @param d Integer vector, event indicator of each patient.
#' @param hazard_ratio Numeric vector, assumed hazard ratio of each level.
#' @param patients Integer vector, assumed number of patients of each level.
#' @param events Integer vector, assumed number of events of each level.
#' @examples
#'   optimize_problem(
#'     z = c(0L, 1L, 0L, 1L),
#'     t = c(1, 2, 3, 4),
#'     d = c(1L, 1L, 0L, 1L),
#'     hazard_ratio = c(0.8, 1.2),
#'     patients = c(2L, 2L),
#'     events = c(1L, 2L)
#'   )
optimize_problem <- function(z, t, d, hazard_ratio, patients, events) {
  list(
    z = z,
    t = t,
    d = d,
    hazard_ratio = hazard_ratio,
    patients = patients,
    events = events
  )
}

#' @title Loss terms of all levels of one subgroup variable
#' @noRd
#' @return Numeric vector with one loss term per level. The loss
#'   \eqn{\Delta(x)} is the square root of the sum.
#' @param x Integer vector with one element per patient: the index of
#'   the current subgroup level \eqn{x_i} of each patient.
#' @inheritParams optimize_level_term
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
#'   problem <- optimize_problem(
#'     z = patient$treatment,
#'     t = patient$time,
#'     d = patient$event,
#'     hazard_ratio = sex$hazard_ratio,
#'     patients = sex$patients,
#'     events = sex$events
#'   )
#'   x <- sample.int(2L, nrow(patient), replace = TRUE)
#'   optimize_loss_terms(x, problem)
optimize_loss_terms <- function(x, problem) {
  vapply(
    seq_along(problem$hazard_ratio),
    \(k) optimize_level_term(x == k, k, problem),
    numeric(1L)
  )
}

#' @title Flip step
#' @noRd
#' @description Implements the flip step in the "Moves" section of
#'   `vignette("optimization", package = "retrogroup")`. The step reads
#'   the current state of the search and returns the move without
#'   modifying `x`.
#' @return A move: a list with elements
#'   * `index`: integer scalar, the patient \eqn{i} drawn in step 1.
#'   * `level`: integer scalar, the level \eqn{k^*} of patient \eqn{i}
#'     after step 3. Equals `x[index]` if the patient stays put.
#'   * `terms`: numeric vector, the loss terms of the assignments after
#'     step 3.
#'   * `loss`: numeric scalar, the loss \eqn{\Delta} after step 3, never
#'     greater than `loss`.
#'
#'   Apply the move with `x[move$index] <- move$level`.
#' @param x Integer vector with one element per patient: the index of
#'   the current subgroup level \eqn{x_i} of each patient.
#' @param terms Numeric vector, the loss terms of `x` from
#'   `optimize_loss_terms()`.
#' @param loss Numeric scalar, the loss \eqn{\Delta(x)}, equal to
#'   `sqrt(sum(terms))`.
#' @inheritParams optimize_level_term
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
#'   problem <- optimize_problem(
#'     z = patient$treatment,
#'     t = patient$time,
#'     d = patient$event,
#'     hazard_ratio = sex$hazard_ratio,
#'     patients = sex$patients,
#'     events = sex$events
#'   )
#'   x <- sample.int(2L, nrow(patient), replace = TRUE)
#'   terms <- optimize_loss_terms(x, problem)
#'   move <- optimize_flip(x, terms, sqrt(sum(terms)), problem)
#'   x[move$index] <- move$level
optimize_flip <- function(x, terms, loss, problem) {
  K <- length(terms)
  # Step 1: draw a patient i uniformly at random from 1, ..., n.
  i <- sample.int(length(x), 1L)
  x_i <- x[i]
  # Step 2: compute the counterfactual loss Delta(x^(i -> k)) for every
  # subgroup level k = 1, ..., K. A flip from level x_i to level k changes
  # the patients in those two levels only, so every other loss term
  # carries over from x. Column k of counterfactual_terms holds the loss
  # terms of x^(i -> k), and column x_i holds those of x itself.
  # The step needs K model fits: one for level x_i without patient i,
  # and one for each other level k with patient i.
  members <- x == x_i
  members[i] <- FALSE
  term_without_i <- optimize_level_term(members, x_i, problem)
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
      flipped[k] <- optimize_level_term(members, k, problem)
      flipped
    },
    numeric(K)
  )
  # vapply() returns a vector instead of a matrix when K = 1.
  dim(counterfactual_terms) <- c(K, K)
  counterfactual_loss <- sqrt(colSums(counterfactual_terms))
  # Step 3: update x <- x^(i -> k_star), where k_star minimizes the
  # counterfactual loss. The candidates include x^(i -> x_i) = x, so the
  # loss never increases. If the current level x_i ties for the minimum,
  # the patient stays put.
  k_star <- which.min(counterfactual_loss)
  if (!(counterfactual_loss[k_star] < loss)) {
    k_star <- x_i
  }
  list(
    index = i,
    level = k_star,
    terms = counterfactual_terms[, k_star],
    loss = unname(counterfactual_loss[k_star])
  )
}

#' @title Check whether a swap is possible
#' @noRd
#' @return `TRUE` if some pair of patients \eqn{(i, j)} has
#'   \eqn{x_i \ne x_j}, `FALSE` otherwise.
#' @inheritParams optimize_flip
#' @examples
#'   optimize_swap_possible(c(1L, 2L, 1L))
#'   optimize_swap_possible(c(2L, 2L, 2L))
optimize_swap_possible <- function(x) {
  any(x != x[1L])
}

#' @title Swap step
#' @noRd
#' @description Implements the swap step in the "Moves" section of
#'   `vignette("optimization", package = "retrogroup")`. The step reads
#'   the current state of the search and returns the move without
#'   modifying `x`. Requires `optimize_swap_possible(x)`.
#' @return A move: a list with elements
#'   * `index`: integer vector of length 2, the patients \eqn{(i, j)}
#'     drawn in step 1.
#'   * `level`: integer vector of length 2, the levels of patients
#'     \eqn{i} and \eqn{j} after step 3: \eqn{(x_j, x_i)} if the swap is
#'     accepted and \eqn{(x_i, x_j)} otherwise.
#'   * `terms`: numeric vector, the loss terms of the assignments after
#'     step 3.
#'   * `loss`: numeric scalar, the loss \eqn{\Delta} after step 3, never
#'     greater than `loss`.
#'
#'   Apply the move with `x[move$index] <- move$level`.
#' @inheritParams optimize_flip
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
#'   problem <- optimize_problem(
#'     z = patient$treatment,
#'     t = patient$time,
#'     d = patient$event,
#'     hazard_ratio = sex$hazard_ratio,
#'     patients = sex$patients,
#'     events = sex$events
#'   )
#'   x <- sample.int(2L, nrow(patient), replace = TRUE)
#'   terms <- optimize_loss_terms(x, problem)
#'   move <- optimize_swap(x, terms, sqrt(sum(terms)), problem)
#'   x[move$index] <- move$level
optimize_swap <- function(x, terms, loss, problem) {
  K <- length(terms)
  # Step 1: draw a pair of patients (i, j) uniformly at random from the
  # pairs with x_i != x_j. The number of pairs with x_i = k and x_j = l is
  # n_k * n_l, where n_k is the current number of patients in level k.
  # First draw the pair of levels (x_i, x_j) with probability proportional
  # to n_k * n_l, then draw i and j uniformly from the patients in those
  # levels. Each unordered pair of distinct levels k < l is one row of
  # level_pairs.
  n <- tabulate(x, nbins = K)
  level_pairs <- which(upper.tri(diag(K)), arr.ind = TRUE)
  weights <- n[level_pairs[, 1L]] * n[level_pairs[, 2L]]
  pair <- level_pairs[sample.int(nrow(level_pairs), 1L, prob = weights), ]
  x_i <- pair[[1L]]
  x_j <- pair[[2L]]
  patients_i <- which(x == x_i)
  patients_j <- which(x == x_j)
  i <- patients_i[sample.int(length(patients_i), 1L)]
  j <- patients_j[sample.int(length(patients_j), 1L)]
  # Step 2: compute the loss Delta(x^(i <-> j)) of the swap. A swap changes
  # the patients in levels x_i and x_j only, so every other loss term
  # carries over from x. The step needs two model fits: one for level x_i
  # with j in place of i, and one for level x_j with i in place of j.
  members_i <- x == x_i
  members_i[c(i, j)] <- c(FALSE, TRUE)
  members_j <- x == x_j
  members_j[c(i, j)] <- c(TRUE, FALSE)
  swapped_terms <- terms
  swapped_terms[x_i] <- optimize_level_term(members_i, x_i, problem)
  swapped_terms[x_j] <- optimize_level_term(members_j, x_j, problem)
  swapped_loss <- sqrt(sum(swapped_terms))
  # Step 3: update x <- x^(i <-> j) if the swap reduces the loss, and
  # otherwise leave x unchanged, so the loss never increases.
  if (swapped_loss < loss) {
    list(
      index = c(i, j),
      level = c(x_j, x_i),
      terms = swapped_terms,
      loss = swapped_loss
    )
  } else {
    list(index = c(i, j), level = c(x_i, x_j), terms = terms, loss = loss)
  }
}


#' @title Validate the arguments of the optimization functions
#' @noRd
#' @return `NULL` (invisibly). Throws an error if any argument is invalid.
#' @param patient Patient-level data with subgroup columns, the output of
#'   [group_sample_unconstrained()] or [group_sample_constrained()].
#' @param summary Summary subgroup statistics, the output of
#'   [group_data_summary()].
#' @param iterations Maximum number of iterations per subgroup variable.
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
#'   sampled <- group_sample_unconstrained(patient, summary)
#'   assert_optimize_arguments(sampled, summary, iterations = 10, tolerance = 1e-3)
#'   try(assert_optimize_arguments(sampled, summary, iterations = -1, tolerance = 1e-3))
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
#' @param tolerance Loss tolerance on the scale of the hazard ratios and
#'   event fractions.
#' @examples
#'   warn_optimize_tolerance(
#'     loss = c(age = 1e-8, sex = 0.02),
#'     iterations = 100,
#'     tolerance = 1e-3
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
          "tolerance applies to the hazard ratios on the linear scale",
          "and to the event fractions.",
          "Increase iterations, rerun from freshly sampled subgroup levels,",
          "or increase tolerance."
        )
      )
    )
  }
  invisible()
}
