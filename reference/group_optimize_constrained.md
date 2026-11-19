# Constrained subgroup optimization

Refine sampled subgroup levels so that the subgroup statistics observed
in the data match the published subgroup statistics in `summary`, while
preserving the published number of patients in each subgroup level.

## Usage

``` r
group_optimize_constrained(
  patient,
  summary,
  iterations = 1000L,
  tolerance = 1e-06
)
```

## Arguments

- patient:

  Patient-level survival data, the output of
  [`group_data_patient()`](https://wlandau.github.io/retrogroup/reference/group_data_patient.md).

- summary:

  Summary subgroup statistics, the output of
  [`group_data_summary()`](https://wlandau.github.io/retrogroup/reference/group_data_summary.md).

- iterations:

  Non-negative whole number, the maximum number of iterations for each
  subgroup variable. Each iteration draws one pair of patients and fits
  two Cox models.

- tolerance:

  Non-negative number. The search for each subgroup variable stops when
  the loss \\\Delta(x)\\ falls below `tolerance`. The loss is a sum of
  squares, so the default of `1e-6` stops when the observed and
  published statistics agree to roughly `1e-3`.

## Value

A tibble with the same rows and columns as `patient`, in the same order.
Each subgroup column is a factor of the optimized levels whose levels
are those of the subgroup variable in `summary$level`, in the same
order. The number of patients in each level equals the published number
of patients in `summary`. All other columns are unchanged.

## Details

`group_optimize_constrained()` implements the "Optimization for
constrained sampling" section of
[`vignette("optimization", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/optimization.md),
separately for each subgroup variable. The loss \\\Delta(x)\\ is the
same as in
[`group_optimize_unconstrained()`](https://wlandau.github.io/retrogroup/reference/group_optimize_unconstrained.md):
the squared Euclidean distance between the observed and published log
hazard ratios, plus the squared distance between the observed and
published standard errors for the subgroup levels whose `standard_error`
in `summary` is not `NA`. The observed statistics come from a Cox
proportional hazards model of each subgroup level with treatment as the
only covariate, as in
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html). If
the model cannot be fit for some subgroup level, for example because one
of its arms has no events, the loss is infinite.

Each iteration draws a pair of patients \\(i, j)\\ uniformly at random
from the pairs with \\x_i \ne x_j\\ and swaps their subgroup levels if
the swap reduces the loss. Otherwise, the subgroup levels stay the same,
so the loss never increases. The search for each subgroup variable stops
when the loss falls below `tolerance` or after `iterations` iterations,
whichever comes first. If the loss of any subgroup variable is still at
least `tolerance` after `iterations` iterations,
`group_optimize_constrained()` throws a warning of class
`"retrogroup_warning"` that lists those subgroup variables and their
final losses.

The search is greedy, so it may settle in a local minimum. Swaps
preserve the number of patients in each subgroup level, so the optimized
counts equal the published number of patients in `summary`. Use
[`group_optimize_unconstrained()`](https://wlandau.github.io/retrogroup/reference/group_optimize_unconstrained.md)
on the output of
[`group_sample_unconstrained()`](https://wlandau.github.io/retrogroup/reference/group_sample_unconstrained.md)
to optimize without the count constraint.

In addition to the requirements of
[`group_analysis_probability()`](https://wlandau.github.io/retrogroup/reference/group_analysis_probability.md),
`patient` must have one column per subgroup variable in `summary`, named
after the subgroup variable, with no missing values and only the levels
of that subgroup variable in `summary`. The number of patients in each
subgroup level of these columns must equal the published number of
patients in `summary`. The output of
[`group_sample_constrained()`](https://wlandau.github.io/retrogroup/reference/group_sample_constrained.md)
satisfies these requirements. These columns are the starting subgroup
levels, and all other columns of `patient` are returned unchanged.

Optimization uses R's random number generator. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand for
reproducible results.

## See also

Other optimization:
[`group_optimize_unconstrained()`](https://wlandau.github.io/retrogroup/reference/group_optimize_unconstrained.md)

## Examples

``` r
  patient <- group_data_patient(
    group_example_patient(),
    arm = "arm",
    treatment = "Lev+5FU",
    time = "time",
    event = "event"
  )
  summary <- group_data_summary(
    group_example_summary(),
    variable = "variable",
    level = "level",
    patients = "patients",
    log_hazard_ratio = "log_hazard_ratio",
    standard_error = "standard_error"
  )
  summary <- summary[summary$variable == "age", ]
  set.seed(0)
  sample <- group_sample_constrained(patient, summary)
  optimized <- group_optimize_constrained(
    sample,
    summary,
    tolerance = 1e-4,
    iterations = 500L
  )
  summary
#> # A tibble: 2 × 5
#>   variable level      patients log_hazard_ratio standard_error
#>   <chr>    <chr>         <int>            <dbl>          <dbl>
#> 1 age      under 60        280           -0.169          0.178
#> 2 age      60 or over      339           -0.546          0.160
  group_analysis_summary(sample, "age")
#> # A tibble: 2 × 5
#>   variable level      patients log_hazard_ratio standard_error
#>   <chr>    <chr>         <int>            <dbl>          <dbl>
#> 1 age      under 60        280           -0.198          0.174
#> 2 age      60 or over      339           -0.517          0.163
  group_analysis_summary(optimized, "age")
#> # A tibble: 2 × 5
#>   variable level      patients log_hazard_ratio standard_error
#>   <chr>    <chr>         <int>            <dbl>          <dbl>
#> 1 age      under 60        280           -0.169          0.174
#> 2 age      60 or over      339           -0.542          0.163
```
