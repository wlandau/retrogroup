#' @title Example colon cancer trial subset
#' @noRd
#' @description Subset of [survival::colon] shared by the example datasets:
#'   the death endpoint (`etype == 2`) in the observation (`"Obs"`) and
#'   levamisole plus fluorouracil (`"Lev+5FU"`) arms.
#' @return A data frame with one row per patient and all the columns
#'   of [survival::colon].
#' @examples
#'   head(example_colon())
example_colon <- function() {
  data <- survival::colon
  data[data$etype == 2 & data$rx %in% c("Obs", "Lev+5FU"), ]
}
