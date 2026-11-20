optimize_test_problem <- function(variable = "node4") {
  data <- probability_example()
  rows <- data$summary$variable == variable
  summary <- data$summary[rows, ]
  problem <- retrogroup:::optimize_problem(
    z = data$patient$treatment,
    t = data$patient$time,
    d = data$patient$event,
    hazard_ratio = summary$hazard_ratio,
    patients = summary$patients,
    events = summary$events
  )
  list(data = data, summary = summary, problem = problem)
}

test_that("optimize_level_term() agrees with optimize_loss_terms()", {
  setup <- optimize_test_problem("age")
  set.seed(0)
  K <- nrow(setup$summary)
  x <- sample.int(K, nrow(setup$data$patient), replace = TRUE)
  terms <- retrogroup:::optimize_loss_terms(x, setup$problem)
  for (k in seq_len(K)) {
    expect_identical(
      retrogroup:::optimize_level_term(x == k, k, setup$problem),
      terms[k]
    )
  }
})

test_that("optimize_swap_possible() detects pairs in different levels", {
  expect_true(retrogroup:::optimize_swap_possible(c(1L, 2L, 1L)))
  expect_true(retrogroup:::optimize_swap_possible(c(3L, 3L, 1L)))
  expect_false(retrogroup:::optimize_swap_possible(c(2L, 2L, 2L)))
  expect_false(retrogroup:::optimize_swap_possible(1L))
})

test_that("optimize_flip() takes the argmin over all K levels", {
  setup <- optimize_test_problem("age")
  problem <- setup$problem
  K <- nrow(setup$summary)
  set.seed(0)
  x <- sample.int(K, nrow(setup$data$patient), replace = TRUE)
  for (rep in seq_len(10L)) {
    terms <- retrogroup:::optimize_loss_terms(x, problem)
    loss <- sqrt(sum(terms))
    move <- retrogroup:::optimize_flip(x, terms, loss, problem)
    expect_length(move$index, 1L)
    expect_length(move$level, 1L)
    expect_length(move$terms, K)
    # Brute force: the loss of x^(i -> k) for every level k.
    candidates <- vapply(
      seq_len(K),
      \(k) {
        y <- x
        y[move$index] <- k
        sqrt(sum(retrogroup:::optimize_loss_terms(y, problem)))
      },
      numeric(1L)
    )
    if (min(candidates) < loss) {
      expect_identical(move$level, which.min(candidates))
    } else {
      expect_identical(move$level, x[move$index])
    }
    expect_equal(move$loss, min(candidates))
    expect_lte(move$loss, loss)
    x[move$index] <- move$level
    expect_equal(move$terms, retrogroup:::optimize_loss_terms(x, problem))
    expect_equal(move$loss, sqrt(sum(move$terms)))
  }
})

test_that("optimize_flip() leaves the patient in place on a tie", {
  data <- probability_small()
  patient <- data$patient
  # With one level, the only candidate is x itself.
  problem <- retrogroup:::optimize_problem(
    z = patient$treatment,
    t = patient$time,
    d = patient$event,
    hazard_ratio = 1,
    patients = nrow(patient),
    events = sum(patient$event)
  )
  x <- rep(1L, nrow(patient))
  terms <- retrogroup:::optimize_loss_terms(x, problem)
  loss <- sqrt(sum(terms))
  set.seed(0)
  move <- retrogroup:::optimize_flip(x, terms, loss, problem)
  expect_identical(move$level, 1L)
  expect_identical(move$terms, terms)
  expect_identical(move$loss, loss)
})

test_that("optimize_swap() accepts exactly the swaps that reduce the loss", {
  setup <- optimize_test_problem("node4")
  problem <- setup$problem
  K <- nrow(setup$summary)
  set.seed(0)
  x <- sample.int(K, nrow(setup$data$patient), replace = TRUE)
  accepted <- 0L
  rejected <- 0L
  for (rep in seq_len(30L)) {
    terms <- retrogroup:::optimize_loss_terms(x, problem)
    loss <- sqrt(sum(terms))
    move <- retrogroup:::optimize_swap(x, terms, loss, problem)
    i <- move$index[1L]
    j <- move$index[2L]
    expect_length(move$index, 2L)
    expect_false(x[i] == x[j])
    y <- x
    y[c(i, j)] <- x[c(j, i)]
    swapped_terms <- retrogroup:::optimize_loss_terms(y, problem)
    swapped_loss <- sqrt(sum(swapped_terms))
    if (swapped_loss < loss) {
      accepted <- accepted + 1L
      expect_identical(move$level, x[c(j, i)])
      expect_equal(move$terms, swapped_terms)
      expect_equal(move$loss, swapped_loss)
    } else {
      rejected <- rejected + 1L
      expect_identical(move$level, x[c(i, j)])
      expect_identical(move$terms, terms)
      expect_identical(move$loss, loss)
    }
    counts <- tabulate(x, K)
    x[move$index] <- move$level
    expect_identical(tabulate(x, K), counts)
  }
  expect_gt(accepted, 0L)
  expect_gt(rejected, 0L)
})
