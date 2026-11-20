optimize_loss_reference <- function(patient, summary, variable) {
  rows <- summary$variable == variable
  levels <- summary$level[rows]
  hazard_ratio <- summary$hazard_ratio[rows]
  terms <- vapply(
    seq_along(levels),
    \(k) {
      data <- patient[patient[[variable]] == levels[k], ]
      fit <- survival::coxph(survival::Surv(time, event) ~ treatment, data)
      (exp(unname(stats::coef(fit))) - hazard_ratio[k])^2
    },
    numeric(1L)
  )
  sqrt(sum(terms))
}
