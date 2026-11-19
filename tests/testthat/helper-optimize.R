optimize_loss_reference <- function(patient, summary, variable) {
  rows <- summary$variable == variable
  levels <- summary$level[rows]
  beta <- summary$log_hazard_ratio[rows]
  sigma <- summary$standard_error[rows]
  terms <- vapply(
    seq_along(levels),
    \(k) {
      data <- patient[patient[[variable]] == levels[k], ]
      fit <- survival::coxph(survival::Surv(time, event) ~ treatment, data)
      out <- (unname(stats::coef(fit)) - beta[k])^2
      if (!is.na(sigma[k])) {
        out <- out + (sqrt(stats::vcov(fit)[1L, 1L]) - sigma[k])^2
      }
      out
    },
    numeric(1L)
  )
  sum(terms)
}
