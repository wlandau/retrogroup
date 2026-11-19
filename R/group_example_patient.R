#' @title Example patient-level survival data
#' @export
#' @family example
#' @description Return example patient-level data in the format expected
#'   by [group_data_patient()].
#' @details The data come from the adjuvant colon cancer trial of
#'   levamisole plus fluorouracil (Moertel et al., 1990), as distributed in
#'   [survival::colon]. The subset keeps the death endpoint (`etype == 2`)
#'   and two arms: observation (`"Obs"`, control) and levamisole plus
#'   fluorouracil (`"Lev+5FU"`, treatment). The levamisole-only arm is
#'   excluded. See [survival::colon] for the full dataset and its
#'   documentation.
#' @return A tibble with one row per patient and the following columns:
#'   * `arm`: character, study arm, either `"Obs"` or `"Lev+5FU"`.
#'   * `time`: numeric, days from randomization to death or censoring.
#'   * `event`: integer, `1` if the patient died and `0` if censored.
#' @references Moertel CG, Fleming TR, Macdonald JS, et al. (1990).
#'   Levamisole and fluorouracil for adjuvant therapy of resected colon
#'   carcinoma. New England Journal of Medicine, 322(6), 352-358.
#' @examples
#' group_example_patient()
group_example_patient <- function() {
  data <- example_colon()
  tibble::tibble(
    arm = as.character(data$rx),
    time = as.numeric(data$time),
    event = as.integer(data$status)
  )
}
