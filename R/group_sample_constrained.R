#' @title Constrained subgroup sampling
#' @export
#' @family sampling
#' @description Sample a level of each subgroup variable for each patient
#'   from the posterior probabilities estimated by [group_analysis_probability()],
#'   such that the number of patients in each subgroup level equals the
#'   assumed number of patients in `summary`.
#' @details `group_sample_constrained()` implements the
#'   "Constrained sampling" section of
#'   `vignette("sampling", package = "retrogroup")`. For each subgroup
#'   variable, it iterates through a random permutation of patients and
#'   draws the level of each patient from a categorical distribution with
#'   probabilities proportional to the \eqn{p_{ik}} from
#'   [group_analysis_probability()], restricted to the subgroup levels that
#'   have not yet reached their assumed counts. All the requirements of
#'   [group_analysis_probability()] on `patient` and `summary` apply. In addition,
#'   no subgroup variable in `summary` may have the same name as a column
#'   of `patient`.
#'
#'   The sampled number of patients in each subgroup level always equals
#'   `summary$patients`. However, the result is only an approximate draw from independent categorical
#'   distributions conditional on those counts. Use
#'   [group_sample_unconstrained()] for independent draws without the
#'   count constraint.
#'
#'   Sampling does not target the assumed event counts in
#'   `summary$events`. The model behind [group_analysis_probability()] has
#'   no prognostic effect, so the sampled subgroups carry no prognostic
#'   information and should not be used for covariate adjustment.
#'   The loss functions of [group_optimize_constrained()]
#'   and [group_optimize_unconstrained()] steer the event counts toward
#'   `summary$events`, but only as a soft constraint.
#'   The statistician who uses the sampled data in an analysis is
#'   responsible for the consequences.
#'
#'   Sampling uses R's random number generator. Call [set.seed()]
#'   beforehand for reproducible results.
#' @return A tibble with one row per patient, in the same order as the
#'   rows of `patient`. The columns are all the columns of `patient`,
#'   followed by one column per subgroup variable, in the order the
#'   subgroup variables first appear in `summary$variable`. Each subgroup
#'   column is a factor of the sampled levels whose levels are those of the
#'   subgroup variable in `summary$level`, in the same order.
#' @inheritParams group_analysis_probability
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
#'   set.seed(0)
#'   sample <- group_sample_constrained(patient, summary)
#'   sample
#'   table(sample$age)
#'   summary[summary$variable == "age", c("level", "patients")]
group_sample_constrained <- function(patient, summary) {
  probability <- group_analysis_probability(patient, summary)
  variables <- names(probability)
  assert_sample_names(patient, variables)
  labels <- lapply(variables, \(variable) {
    rows <- summary$variable == variable
    levels <- summary$level[rows]
    index <- sample_constrained_variable(
      p = probability[[variable]],
      n = summary$patients[rows]
    )
    factor(levels[index], levels = levels)
  })
  names(labels) <- variables
  tibble::as_tibble(c(as.list(patient), labels))
}

#' @title Constrained sampling of one subgroup variable
#' @noRd
#' @return Integer vector with one element per row of `p`: the column index
#'   of the sampled level of each patient. Level `k` appears exactly
#'   `n[k]` times.
#' @param p Numeric matrix of posterior probabilities of one subgroup
#'   variable, with one row per patient and one column per level.
#'   Each row sums to 1, and every element is strictly positive.
#' @param n Integer vector with one element per column of `p`: the
#'   assumed number of patients in each level. Must sum to `nrow(p)`.
#' @examples
#'   p <- matrix(c(0.9, 0.1, 0.5, 0.5, 0.2, 0.8), ncol = 2, byrow = TRUE)
#'   colnames(p) <- c("female", "male")
#'   p
#'   sample_constrained_variable(p, n = c(1L, 2L))
sample_constrained_variable <- function(p, n) {
  # Since sum(n) = nrow(p), the remaining slots always equal the number of
  # unassigned patients, so some level always has room.
  remaining <- as.integer(n)
  index <- integer(nrow(p))
  u <- stats::runif(nrow(p))
  for (i in sample.int(nrow(p))) {
    open <- which(remaining > 0L)
    # Choosing among open levels only guards against rounding error.
    cumulative <- cumsum(p[i, open])
    threshold <- u[i] * cumulative[length(open)]
    k <- open[1L + sum(threshold > cumulative[-length(open)])]
    index[i] <- k
    remaining[k] <- remaining[k] - 1L
  }
  index
}
