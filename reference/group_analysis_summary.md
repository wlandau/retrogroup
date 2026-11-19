# Subgroup-level hazard ratios

Compute the number of patients, the log hazard ratio, and its standard
error in each level of each subgroup variable of a patient-level
dataset.

## Usage

``` r
group_analysis_summary(
  patient,
  variables = setdiff(colnames(patient), c("treatment", "time", "event"))
)
```

## Arguments

- patient:

  Patient-level survival data with subgroup columns, usually the output
  of
  [`group_sample_unconstrained()`](https://wlandau.github.io/retrogroup/reference/group_sample_unconstrained.md),
  [`group_sample_constrained()`](https://wlandau.github.io/retrogroup/reference/group_sample_constrained.md),
  [`group_optimize_unconstrained()`](https://wlandau.github.io/retrogroup/reference/group_optimize_unconstrained.md),
  or
  [`group_optimize_constrained()`](https://wlandau.github.io/retrogroup/reference/group_optimize_constrained.md).
  Must have the columns of the output of
  [`group_data_patient()`](https://wlandau.github.io/retrogroup/reference/group_data_patient.md)
  and the subgroup columns named in `variables`.

- variables:

  Character vector of unique names of subgroup columns in `patient`.
  Each subgroup column must be character or factor with no missing
  values. Defaults to every column of `patient` other than `treatment`,
  `time`, and `event`.

## Value

A tibble in the format of the output of
[`group_data_summary()`](https://wlandau.github.io/retrogroup/reference/group_data_summary.md),
with one row per subgroup level, ordered by subgroup variable as in
`variables` and then by level, and the following columns:

- `variable`: character, name of the subgroup variable.

- `level`: character, level of the subgroup variable.

- `patients`: integer, number of patients in the subgroup level.

- `log_hazard_ratio`: numeric, estimated log hazard ratio of treatment
  versus control in the subgroup level.

- `standard_error`: numeric, standard error of the log hazard ratio.

## Details

For each level of each subgroup variable, `group_analysis_summary()`
fits a Cox proportional hazards model by partial likelihood to the
patients in that level, with treatment as the only covariate and the
Efron method for ties, as in
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html). The
log hazard ratio is the estimated coefficient of treatment versus
control, and the standard error is on the log hazard ratio scale. These
are the observed statistics \\\hat\beta_k(x)\\ and \\\hat\sigma_k(x)\\
of the "Loss function" section of
[`vignette("optimization", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/optimization.md).

If the model cannot be fit for some subgroup level, for example because
one of its arms has no patients or no events, the `log_hazard_ratio` and
`standard_error` of that level are `NA`.

The levels of a factor subgroup column are its factor levels, in order,
excluding levels with no patients. The levels of a character subgroup
column are its sorted unique values.

## See also

Other analysis:
[`group_analysis_probability()`](https://wlandau.github.io/retrogroup/reference/group_analysis_probability.md)

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
  set.seed(0)
  sample <- group_sample_constrained(patient, summary)
  group_analysis_summary(sample)
#> # A tibble: 6 × 5
#>   variable level      patients log_hazard_ratio standard_error
#>   <chr>    <chr>         <int>            <dbl>          <dbl>
#> 1 age      under 60        280           -0.198          0.174
#> 2 age      60 or over      339           -0.517          0.163
#> 3 node4    no              453           -0.291          0.139
#> 4 node4    yes             166           -0.587          0.230
#> 5 sex      female          312           -0.188          0.160
#> 6 sex      male            307           -0.587          0.179
  summary
#> # A tibble: 6 × 5
#>   variable level      patients log_hazard_ratio standard_error
#>   <chr>    <chr>         <int>            <dbl>          <dbl>
#> 1 age      under 60        280           -0.169          0.178
#> 2 age      60 or over      339           -0.546          0.160
#> 3 node4    no              453           -0.417          0.153
#> 4 node4    yes             166           -0.312          0.190
#> 5 sex      female          312           -0.147          0.162
#> 6 sex      male            307           -0.656          0.179
```
