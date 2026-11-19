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
#'   Within each subgroup level, [group_analysis_summary()] fits a Cox
#'   proportional hazards model of survival on treatment (`"Lev+5FU"`
#'   versus `"Obs"`) to get the log hazard ratio and its standard error.
#'   The standard error is on the log hazard ratio scale,
#'   not the hazard ratio scale.
#' @return A tibble with one row per subgroup level and the following columns:
#'   * `variable`: character, name of the subgroup variable.
#'   * `level`: character, level of the subgroup variable.
#'   * `patients`: integer, number of patients in the subgroup level.
#'   * `log_hazard_ratio`: numeric, estimated log hazard ratio of treatment
#'     versus control in the subgroup level.
#'   * `standard_error`: numeric, standard error of the log hazard ratio.
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
  group_analysis_summary(patient, c("age", "node4", "sex"))
}
