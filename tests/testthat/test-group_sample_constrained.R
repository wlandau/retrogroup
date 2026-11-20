test_that("group_sample_constrained() has the expected structure", {
  data <- probability_example()
  set.seed(0)
  x <- group_sample_constrained(data$patient, data$summary)
  probability <- group_analysis_probability(data$patient, data$summary)
  expect_s3_class(x, "tbl_df")
  expect_equal(nrow(x), nrow(data$patient))
  expect_identical(
    colnames(x),
    c(colnames(data$patient), names(probability))
  )
  expect_identical(
    x[, colnames(data$patient)],
    tibble::as_tibble(data$patient)
  )
  for (variable in names(probability)) {
    levels <- data$summary$level[data$summary$variable == variable]
    expect_s3_class(x[[variable]], "factor")
    expect_identical(levels(x[[variable]]), levels)
    expect_false(anyNA(x[[variable]]))
    expect_true(all(x[[variable]] %in% levels))
  }
})

test_that("group_sample_constrained() copies labels from summary", {
  data <- probability_small()
  summary <- data$summary[c(4, 2, 5, 1, 3), ]
  set.seed(0)
  x <- group_sample_constrained(data$patient, summary)
  expect_identical(colnames(x), c("treatment", "time", "event", "sex", "age"))
  expect_true(all(x$sex %in% c("male", "female", "other")))
  expect_true(all(x$age %in% c("young", "old")))
})

test_that("group_sample_constrained() is reproducible with set.seed()", {
  data <- probability_example()
  set.seed(1)
  x <- group_sample_constrained(data$patient, data$summary)
  set.seed(1)
  y <- group_sample_constrained(data$patient, data$summary)
  set.seed(2)
  z <- group_sample_constrained(data$patient, data$summary)
  expect_identical(x, y)
  expect_false(identical(x, z))
})

test_that("group_sample_constrained() preserves subgroup counts", {
  set.seed(0)
  for (data in list(probability_example(), probability_small())) {
    for (rep in seq_len(25L)) {
      x <- group_sample_constrained(data$patient, data$summary)
      for (variable in unique(data$summary$variable)) {
        rows <- data$summary$variable == variable
        levels <- data$summary$level[rows]
        counts <- table(factor(x[[variable]], levels = levels))
        expect_identical(
          as.integer(counts),
          data$summary$patients[rows]
        )
      }
    }
  }
})

test_that("group_sample_constrained() preserves subgroup event counts", {
  set.seed(0)
  for (data in list(probability_example(), probability_small())) {
    for (rep in seq_len(200L)) {
      x <- group_sample_constrained(data$patient, data$summary)
      for (variable in unique(data$summary$variable)) {
        rows <- data$summary$variable == variable
        levels <- data$summary$level[rows]
        values <- factor(x[[variable]], levels = levels)
        expect_identical(
          as.integer(table(values)),
          data$summary$patients[rows]
        )
        expect_identical(
          as.integer(tapply(x$event, values, sum, default = 0L)),
          data$summary$events[rows]
        )
      }
    }
  }
})

test_that("sample_constrained_variable() preserves counts", {
  set.seed(0)
  p <- matrix(stats::runif(300L), 100L, 3L)
  p <- p / rowSums(p)
  n <- c(17L, 50L, 33L)
  e <- c(5L, 30L, 10L)
  for (rep in seq_len(50L)) {
    d <- sample(rep(c(1L, 0L), c(sum(e), sum(n) - sum(e))))
    index <- retrogroup:::sample_constrained_variable(p, n, e, d)
    expect_identical(tabulate(index, nbins = 3L), n)
    expect_identical(tabulate(index[d == 1L], nbins = 3L), e)
  }
})

test_that("sample_constrained_variable() never gets stuck", {
  # Random configurations, including e_k = n_k and e_k = 0 or 1, with
  # probabilities that push patients toward levels with few slots.
  set.seed(0)
  for (rep in seq_len(500L)) {
    K <- sample(2:5, 1L)
    n <- sample.int(8L, K, replace = TRUE)
    e <- vapply(n, \(n_k) sample(0:n_k, 1L), integer(1L))
    e[1L] <- n[1L]
    if (K > 2L) e[2L] <- min(1L, n[2L])
    N <- sum(n)
    d <- sample(rep(c(1L, 0L), c(sum(e), N - sum(e))))
    epsilon <- 1e-7
    favorite <- ifelse(d == 1L, which.min(replace(e, e == 0L, Inf)), K)
    p <- matrix(epsilon, N, K)
    p[cbind(seq_len(N), favorite)] <- 1
    p <- p / rowSums(p)
    index <- retrogroup:::sample_constrained_variable(p, n, e, d)
    expect_false(anyNA(index))
    expect_true(all(index %in% seq_len(K)))
    expect_identical(tabulate(index, nbins = K), n)
    expect_identical(tabulate(index[d == 1L], nbins = K), e)
  }
})

test_that("sample_constrained_variable() handles extreme probabilities", {
  set.seed(0)
  epsilon <- 1e-7
  p <- matrix(c(1 - 2 * epsilon, epsilon, epsilon), 1000L, 3L, byrow = TRUE)
  n <- c(10L, 490L, 500L)
  e <- c(10L, 1L, 289L)
  for (rep in seq_len(10L)) {
    d <- sample(rep(c(1L, 0L), c(sum(e), sum(n) - sum(e))))
    index <- retrogroup:::sample_constrained_variable(p, n, e, d)
    expect_identical(tabulate(index, nbins = 3L), n)
    expect_identical(tabulate(index[d == 1L], nbins = 3L), e)
  }
})

test_that("sample_constrained_variable() follows the probabilities", {
  # Patients with high probability of level 1 should usually get level 1.
  set.seed(0)
  p <- rbind(
    matrix(c(0.98, 0.01, 0.01), 50L, 3L, byrow = TRUE),
    matrix(c(0.01, 0.98, 0.01), 50L, 3L, byrow = TRUE),
    matrix(c(0.01, 0.01, 0.98), 50L, 3L, byrow = TRUE)
  )
  n <- c(50L, 50L, 50L)
  truth <- rep(1:3, each = 50L)
  d <- rep(rep(c(1L, 0L), c(20L, 30L)), 3L)
  e <- c(20L, 20L, 20L)
  agreement <- replicate(
    20L,
    mean(retrogroup:::sample_constrained_variable(p, n, e, d) == truth)
  )
  expect_gt(mean(agreement), 0.9)
})

test_that("group_sample_constrained() validates its arguments", {
  data <- probability_small()
  expect_error(
    group_sample_constrained(as.list(data$patient), data$summary),
    "data frame",
    class = "retrogroup_error"
  )
  summary <- data$summary
  summary$patients[1] <- 0L
  expect_error(
    group_sample_constrained(data$patient, summary),
    "summary is invalid",
    class = "retrogroup_error"
  )
  summary <- data$summary
  summary$variable[summary$variable == "age"] <- "time"
  expect_error(
    group_sample_constrained(data$patient, summary),
    "clashing names: \"time\"",
    class = "retrogroup_error"
  )
})
