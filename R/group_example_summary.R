#' @title Example summary subgroup statistics
#' @export
#' @family example
#' @description Return example summary subgroup statistics in the format
#'   expected by [group_data_summary()].
#' @details The statistics are computed from the same patients as
#'   [group_example_patient()]: the death endpoint (`etype == 2`) of the
#'   adjuvant colon cancer trial of levamisole plus fluorouracil
#'   (Moertel et al., 1990), as distributed in [survival::colon],
#'   restricted to the observation (`"Obs"`, control) and levamisole plus
#'   fluorouracil (`"Lev+5FU"`, treatment) arms.
#'
#'   There are three subgroup variables:
#'   * `"age"`: age at randomization, with levels `"under 60"`
#'     (younger than 60 years) and `"60 or over"` (60 years or older).
#'   * `"node4"`: more than 4 positive lymph nodes, with levels
#'     `"no"` and `"yes"`.
#'   * `"sex"`: with levels `"female"` and `"male"`.
#'
#'   Within each subgroup level, [group_analysis_summary()] counts the
#'   patients and events and fits a Cox proportional hazards model of
#'   survival on treatment (`"Lev+5FU"` versus `"Obs"`) to get the hazard
#'   ratio. The columns `lower` and `upper` are the bounds of the 95%
#'   Wald confidence interval of the hazard ratio from the same model, as
#'   a publication would report them. [group_data_summary()] ignores
#'   these columns, but [group_data_events()] can use them to approximate
#'   the event counts.
#' @return A tibble with one row per subgroup level and the following columns:
#'   * `variable`: character, name of the subgroup variable.
#'   * `level`: character, level of the subgroup variable.
#'   * `patients`: integer, number of patients in the subgroup level.
#'   * `events`: integer, number of events in the subgroup level.
#'   * `hazard_ratio`: numeric, estimated hazard ratio of treatment
#'     versus control in the subgroup level.
#'   * `lower`: numeric, lower bound of the 95% confidence interval of the
#'     hazard ratio.
#'   * `upper`: numeric, upper bound of the 95% confidence interval of the
#'     hazard ratio.
#' @references Moertel CG, Fleming TR, Macdonald JS, et al. (1990).
#'   Levamisole and fluorouracil for adjuvant therapy of resected colon
#'   carcinoma. New England Journal of Medicine, 322(6), 352-358.
#' @examples
#'   group_example_summary()
group_example_summary <- function() {
  data <- example_colon()
  patient <- tibble::tibble(
    treatment = as.integer(data$rx == "Lev+5FU"),
    time = as.numeric(data$time),
    event = as.integer(data$status),
    age = factor(
      ifelse(data$age < 60, "under 60", "60 or over"),
      levels = c("under 60", "60 or over")
    ),
    node4 = factor(
      ifelse(data$node4 == 1, "yes", "no"),
      levels = c("no", "yes")
    ),
    sex = factor(
      ifelse(data$sex == 1, "male", "female"),
      levels = c("female", "male")
    )
  )
  out <- group_analysis_summary(patient, c("age", "node4", "sex"))
  bounds <- vapply(
    seq_len(nrow(out)),
    \(row) {
      members <- as.character(patient[[out$variable[row]]]) == out$level[row]
      example_interval(patient[members, ])
    },
    numeric(2L)
  )
  out$lower <- bounds[1L, ]
  out$upper <- bounds[2L, ]
  out
}

#' @title 95% confidence interval of the hazard ratio of one subgroup level
#' @noRd
#' @return Numeric vector with the lower and upper bounds of the 95% Wald
#'   confidence interval of the hazard ratio of treatment versus control,
#'   from a Cox proportional hazards model with treatment as the only
#'   covariate and the Efron method for ties.
#' @param patient Data frame with columns `treatment`, `time`, and `event`
#'   for the patients of one subgroup level.
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   example_interval(patient)
example_interval <- function(patient) {
  fit <- survival::coxph(
    survival::Surv(time, event) ~ treatment,
    data = patient,
    ties = "efron"
  )
  unname(exp(stats::confint(fit, level = 0.95))[1L, ])
}
