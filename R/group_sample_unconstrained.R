#' @title Unconstrained subgroup sampling
#' @export
#' @family sampling
#' @description Sample a level of each subgroup variable for each patient,
#'   independently across patients, from the posterior probabilities
#'   estimated by [group_analysis_probability()].
#' @details `group_sample_unconstrained()` implements the
#'   "Unconstrained sampling" section of
#'   `vignette("sampling", package = "retrogroup")`. For each subgroup
#'   variable and each patient \eqn{i}, it draws
#'   \eqn{x_i \sim \text{Categorical}(p_{i1}, \ldots, p_{iK})}
#'   independently, where the \eqn{p_{ik}} come from [group_analysis_probability()].
#'   All the requirements of [group_analysis_probability()] on `patient` and
#'   `summary` apply. In addition, no subgroup variable in `summary`
#'   may have the same name as a column of `patient`.
#'
#'   The sampled number of patients in each subgroup level is not
#'   guaranteed to equal the assumed number of patients in `summary`.
#'   Use [group_sample_constrained()] to preserve the assumed counts.
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
#'   sample <- group_sample_unconstrained(patient, summary)
#'   sample
#'   table(sample$age)
group_sample_unconstrained <- function(patient, summary) {
  probability <- group_analysis_probability(patient, summary)
  variables <- names(probability)
  assert_sample_names(patient, variables)
  labels <- lapply(variables, \(variable) {
    levels <- summary$level[summary$variable == variable]
    factor(
      levels[sample_unconstrained_variable(probability[[variable]])],
      levels = levels
    )
  })
  names(labels) <- variables
  tibble::as_tibble(c(as.list(patient), labels))
}

#' @title Unconstrained sampling of one subgroup variable
#' @noRd
#' @return Integer vector with one element per row of `p`: the column index
#'   of the sampled level of each patient.
#' @param p Numeric matrix of posterior probabilities of one subgroup
#'   variable, with one row per patient and one column per level.
#'   Each row sums to 1.
#' @examples
#'   p <- matrix(c(0.9, 0.1, 0.5, 0.5, 0.2, 0.8), ncol = 2, byrow = TRUE)
#'   colnames(p) <- c("female", "male")
#'   p
#'   sample_unconstrained_variable(p)
sample_unconstrained_variable <- function(p) {
  levels <- ncol(p)
  cumulative <- matrix(
    t(apply(p, 1L, cumsum)),
    nrow = nrow(p),
    ncol = levels
  )
  u <- stats::runif(nrow(p))
  # Omitting the last cumulative column guards against rounding error.
  1L + as.integer(rowSums(u > cumulative[, -levels, drop = FALSE]))
}

#' @title Validate that subgroup variable names do not clash with `patient`
#' @noRd
#' @return `NULL` (invisibly). Throws an error if any subgroup variable has
#'   the same name as a column of `patient`.
#' @inheritParams group_analysis_probability
#' @param variables Character vector, names of the subgroup variables.
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   assert_sample_names(patient, c("age", "sex"))
#'   try(assert_sample_names(patient, "time"))
assert_sample_names <- function(patient, variables) {
  clash <- intersect(variables, colnames(patient))
  if (length(clash)) {
    throw_error(
      c(
        "Subgroup variables must not have the same names as columns of patient.",
        x = paste0(
          "Subgroup variable(s) with clashing names: ",
          paste0("\"", clash, "\"", collapse = ", "),
          "."
        )
      )
    )
  }
  invisible()
}
