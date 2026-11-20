#' @title Approximate subgroup event counts
#' @export
#' @family data
#' @description Approximate the number of events in each subgroup level
#'   from the confidence intervals of the assumed hazard ratios, for
#'   publications that do not report subgroup event counts.
#' @details `group_data_events()` implements
#'   `vignette("events", package = "retrogroup")`, separately for each
#'   subgroup variable. Assuming the confidence interval of each hazard
#'   ratio is a Wald interval that is symmetric on the log scale, and using
#'   the variance approximation of Schoenfeld (1983) as in Parmar et al.
#'   (1998), the approximate
#'   number of events in subgroup level \eqn{k} is
#'   \deqn{\hat{e}_k = \frac{1}{\bar{z} (1 - \bar{z})}
#'   \frac{c^2}{(\log(b_k) - \beta_k)^2}}
#'   where \eqn{\exp(\beta_k)} is `summary$hazard_ratio`, \eqn{b_k} is
#'   `bound`, \eqn{c} is the \eqn{1 - \alpha / 2} quantile of the standard
#'   normal distribution for `confidence` \eqn{= 1 - \alpha}, and
#'   \eqn{\bar{z}} is the proportion of patients in `patient` randomized to
#'   treatment.
#'
#'   The approximate counts do not necessarily add up to the total number
#'   of events in `patient`, and they may fall below 1 or exceed the number
#'   of patients in the level. So `group_data_events()` rounds each
#'   \eqn{\hat{e}_k} to the nearest whole number and clamps it to between 1
#'   and \eqn{n_k}. Then it makes passes over the levels in the order of
#'   `summary`, adding 1 to (or subtracting 1 from) each level that has
#'   room, until the counts of each subgroup variable add up to the total
#'   number of events in `patient`. The result satisfies all the
#'   requirements on event counts of [group_data_summary()]. It is an
#'   error if the total number of events is less than the number of levels
#'   of a subgroup variable, because each level needs at least one event.
#'
#'   This is only one of many possible ways to choose event counts. The
#'   choice is a decision about the reconstruction, and it is your
#'   responsibility. Please disclose the event counts used and how they
#'   were chosen wherever the reconstructed data or downstream results are
#'   reported.
#'
#'   `summary` must satisfy the requirements of
#'   [group_analysis_probability()], except that its event counts are
#'   ignored. Usually `summary` is the output of [group_data_summary()]
#'   with `events = NULL`.
#' @return A tibble with the same rows and columns as the output of
#'   [group_data_summary()], in the same order as `summary`. The `events`
#'   column has the approximate event counts. All other columns are the
#'   same as in `summary`.
#' @references
#'   Parmar MKB, Torri V, Stewart L (1998). Extracting summary statistics
#'   to perform meta-analyses of the published literature for survival
#'   endpoints. Statistics in Medicine, 17(24), 2815-2834.
#'   \doi{10.1002/(SICI)1097-0258(19981230)17:24<2815::AID-SIM110>3.0.CO;2-8}
#'
#'   Schoenfeld DA (1983). Sample-size formula for the
#'   proportional-hazards regression model. Biometrics, 39(2), 499-503.
#'   \doi{10.2307/2531021}
#' @inheritParams group_analysis_probability
#' @param summary Summary subgroup statistics, the output of
#'   [group_data_summary()], usually with `events = NULL`. The event counts
#'   in `summary`, if any, are ignored.
#' @param bound Numeric vector with one element per row of `summary`: the
#'   lower or upper bound \eqn{b_k} of the two-sided confidence interval of
#'   the hazard ratio of each subgroup level, on the natural scale (not
#'   the log scale). The lower and upper bounds give the same result, so
#'   each element may be either one. Every element must be finite,
#'   strictly positive, and different from the hazard ratio of its level.
#' @param confidence Number strictly between 0 and 1, the confidence level
#'   \eqn{1 - \alpha} of the confidence intervals.
#' @examples
#'   patient <- group_data_patient(
#'     group_example_patient(),
#'     arm = "arm",
#'     treatment = "Lev+5FU",
#'     time = "time",
#'     event = "event"
#'   )
#'   data <- group_example_summary()
#'   summary <- group_data_summary(
#'     data,
#'     variable = "variable",
#'     level = "level",
#'     patients = "patients",
#'     events = NULL,
#'     hazard_ratio = "hazard_ratio"
#'   )
#'   summary
#'   group_data_events(patient, summary, bound = data$upper)
#'   data[, c("variable", "level", "events")]
group_data_events <- function(patient, summary, bound, confidence = 0.95) {
  assert_probability_patient(patient)
  assert_probability_summary(summary, events = FALSE)
  assert_probability_joint(patient, summary, events = FALSE)
  assert_data_events(summary, bound, confidence)
  estimate <- events_approximate(
    hazard_ratio = summary$hazard_ratio,
    bound = bound,
    proportion = mean(patient$treatment),
    confidence = confidence
  )
  events <- integer(nrow(summary))
  total <- sum(patient$event)
  for (variable in unique(summary$variable)) {
    rows <- summary$variable == variable
    events[rows] <- events_adjust(
      estimate = estimate[rows],
      n = summary$patients[rows],
      total = total,
      variable = variable
    )
  }
  tibble::tibble(
    variable = as.character(summary$variable),
    level = as.character(summary$level),
    patients = as.integer(summary$patients),
    events = events,
    hazard_ratio = as.numeric(summary$hazard_ratio)
  )
}

#' @title Approximate event counts before adjustment
#' @noRd
#' @return Numeric vector with one element per subgroup level: the
#'   approximate number of events \eqn{\hat{e}_k} from the Schoenfeld
#'   variance approximation, before rounding and adjustment.
#' @param hazard_ratio Numeric vector, assumed hazard ratio of each level.
#' @param bound Numeric vector, lower or upper confidence bound of the
#'   hazard ratio of each level.
#' @param proportion Number, the proportion \eqn{\bar{z}} of patients
#'   randomized to treatment.
#' @param confidence Number, the confidence level of the intervals.
#' @examples
#'   events_approximate(
#'     hazard_ratio = c(0.8, 1.2),
#'     bound = c(1.1, 1.6),
#'     proportion = 0.5,
#'     confidence = 0.95
#'   )
events_approximate <- function(hazard_ratio, bound, proportion, confidence) {
  c <- stats::qnorm(1 - (1 - confidence) / 2)
  c^2 / (proportion * (1 - proportion) * (log(bound) - log(hazard_ratio))^2)
}

#' @title Reconcile approximate event counts
#' @noRd
#' @return Integer vector with one element per level of one subgroup
#'   variable: the event counts, with sum `total` and each element between
#'   1 and the corresponding element of `n`.
#' @details Rounds `estimate` to whole numbers, clamps each element to
#'   between 1 and `n`, then makes passes over the levels in order. Each
#'   pass adds 1 to (or subtracts 1 from) each level with room, in order,
#'   until the sum equals `total`.
#' @param estimate Numeric vector, approximate event count of each level.
#' @param n Integer vector, number of patients of each level.
#' @param total Whole number, the total number of events.
#' @param variable Character string, name of the subgroup variable for
#'   error messages.
#' @examples
#'   events_adjust(c(10.4, 3.6, 0.2), n = c(20L, 4L, 6L), total = 12L)
#'   events_adjust(c(10.4, 3.6, 0.2), n = c(20L, 4L, 6L), total = 20L)
events_adjust <- function(estimate, n, total, variable = "variable") {
  n <- as.integer(n)
  if (total < length(n) || total > sum(n)) {
    throw_error(
      c(
        paste(
          "The total number of events must be at least the number of",
          "levels and at most the number of patients of each subgroup",
          "variable, so that every level can have at least one event and",
          "at most as many events as patients."
        ),
        x = paste0(
          "Subgroup variable \"",
          variable,
          "\" has ",
          length(n),
          " levels and ",
          sum(n),
          " patients, but there are ",
          total,
          " events."
        )
      )
    )
  }
  events <- as.integer(pmin(pmax(round(estimate), 1), n))
  difference <- total - sum(events)
  while (difference != 0L) {
    open <- if (difference > 0L) {
      which(events < n)
    } else {
      which(events > 1L)
    }
    change <- open[seq_len(min(abs(difference), length(open)))]
    events[change] <- events[change] + as.integer(sign(difference))
    difference <- total - sum(events)
  }
  events
}

#' @title Validate the arguments of `group_data_events()`
#' @noRd
#' @return `NULL` (invisibly). Throws an error if `bound` or `confidence`
#'   is invalid.
#' @inheritParams group_data_events
#' @examples
#'   data <- group_example_summary()
#'   assert_data_events(data, data$upper, 0.95)
#'   try(assert_data_events(data, data$hazard_ratio, 0.95))
assert_data_events <- function(summary, bound, confidence) {
  if (!is.numeric(bound) || length(bound) != nrow(summary)) {
    throw_error(
      c(
        "bound must be a numeric vector with one element per row of summary.",
        x = paste0(
          "length(bound) = ",
          length(bound),
          " and nrow(summary) = ",
          nrow(summary),
          "."
        )
      )
    )
  }
  if (anyNA(bound) || !all(is.finite(bound))) {
    throw_error("bound must be finite with no missing values.")
  }
  if (any(bound <= 0)) {
    throw_error(
      c(
        "bound must be strictly positive.",
        i = paste(
          "Confidence bounds must be on the natural scale of the hazard",
          "ratio, not the log scale."
        )
      )
    )
  }
  same <- bound == summary$hazard_ratio
  if (any(same)) {
    throw_error(
      c(
        paste(
          "Each element of bound must differ from the hazard ratio of",
          "its subgroup level. Otherwise, the approximate number of",
          "events is infinite."
        ),
        x = paste0(
          "Subgroup level(s) with bound equal to the hazard ratio: ",
          paste0(
            summary$variable[same],
            " = \"",
            summary$level[same],
            "\"",
            collapse = ", "
          ),
          "."
        )
      )
    )
  }
  if (
    !is.numeric(confidence) ||
      length(confidence) != 1L ||
      is.na(confidence) ||
      confidence <= 0 ||
      confidence >= 1
  ) {
    throw_error("confidence must be a single number strictly between 0 and 1.")
  }
  invisible()
}
