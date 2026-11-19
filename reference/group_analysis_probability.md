# Posterior subgroup probabilities

Estimate the posterior probability that each patient belongs to each
level of each subgroup variable, given the patient-level survival data
and the published subgroup statistics.

## Usage

``` r
group_analysis_probability(patient, summary)
```

## Arguments

- patient:

  Patient-level survival data, the output of
  [`group_data_patient()`](https://wlandau.github.io/retrogroup/reference/group_data_patient.md).

- summary:

  Summary subgroup statistics, the output of
  [`group_data_summary()`](https://wlandau.github.io/retrogroup/reference/group_data_summary.md).

## Value

A named list with one element per subgroup variable, in the order the
subgroup variables first appear in `summary$variable`. Each element is a
numeric matrix of the estimated posterior probabilities \\p\_{ik}\\ of
that subgroup variable:

- Rows correspond to patients, in the same order as the rows of
  `patient`.

- Columns correspond to the levels of the subgroup variable. The column
  names are exactly the levels, in the same order as in `summary`.

- Each row sums to 1.

## Details

`group_analysis_probability()` computes the estimated posterior
probabilities \\p\_{ik}\\ defined in the "Computation" section of
[`vignette("model", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/model.md),
separately for each subgroup variable. See that section for the
formulas, including the error when any probability is numerically
indistinguishable from zero.

The total number of patients of every subgroup variable in `summary`
must equal the number of rows of `patient`. See the details of
[`group_data_summary()`](https://wlandau.github.io/retrogroup/reference/group_data_summary.md)
for how to choose the counts. The `standard_error` column of `summary`
is not used.

## See also

Other analysis:
[`group_analysis_summary()`](https://wlandau.github.io/retrogroup/reference/group_analysis_summary.md)

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
  probability <- group_analysis_probability(patient, summary)
  names(probability)
#> [1] "age"   "node4" "sex"  
  head(probability$age)
#>       under 60 60 or over
#> [1,] 0.5082235  0.4917765
#> [2,] 0.3946053  0.6053947
#> [3,] 0.4523425  0.5476575
#> [4,] 0.5422896  0.4577104
#> [5,] 0.4523425  0.5476575
#> [6,] 0.5055396  0.4944604
```
