# Summary subgroup statistics

Convert published subgroup statistics into the standardized format used
by `retrogroup`.

## Usage

``` r
group_data_summary(
  data,
  variable,
  level,
  patients,
  log_hazard_ratio,
  standard_error = NULL
)
```

## Arguments

- data:

  Data frame with one row per subgroup level. Must have the columns
  named by the other arguments. Other columns are ignored.

- variable:

  Character string, name of the column in `data` with the name of the
  subgroup variable. The column must be character with non-empty,
  non-missing values.

- level:

  Character string, name of the column in `data` with the level of the
  subgroup variable. The column must be character or factor with
  non-empty, non-missing values.

- patients:

  Character string, name of the column in `data` with the number of
  patients to assign to each level of the reconstructed subgroup
  variable, usually taken from the publication. See the details section
  for when these counts may differ from the publication. The column must
  contain strictly positive whole numbers with no missing values.

- log_hazard_ratio:

  Character string, name of the column in `data` with the published log
  hazard ratio of treatment versus control in each subgroup level. The
  column must be numeric and finite with no missing values.

- standard_error:

  Character string or `NULL`. If a character string, name of the column
  in `data` with the published standard error of the log hazard ratio in
  each subgroup level. The column must be numeric, finite, and strictly
  positive with no missing values. If `NULL`, the standard errors are
  unavailable, and the `standard_error` column of the output is `NA`.

## Value

A tibble with one row per subgroup level, in the same order as `data`,
and the following columns:

- `variable`: character, name of the subgroup variable.

- `level`: character, level of the subgroup variable.

- `patients`: integer, number of patients in the subgroup level.

- `log_hazard_ratio`: numeric, published log hazard ratio of treatment
  versus control in the subgroup level.

- `standard_error`: numeric, published standard error of the log hazard
  ratio, or `NA` if `standard_error` is `NULL`.

## Details

Each row of `data` is one level of a subgroup variable. A subgroup
variable must have at least two levels, the levels must be unique within
the subgroup variable, and every level must have a strictly positive
number of patients.

The patient counts must have the same total for every subgroup variable,
because every subgroup variable partitions the same patients. The counts
are the sizes of the levels of the reconstructed subgroup variables, so
they need not copy the publication verbatim. Published counts are the
natural starting point, but they sometimes do not add up to the total
number of patients in the study, for example because of missing subgroup
data or typographical errors. In that case, choosing counts that add up
to the total number of patients is a decision about the reconstruction,
not an attempt to recover the true counts from the clinical study.
Please disclose that decision transparently wherever the reconstructed
data or downstream results are reported.

## See also

Other data:
[`group_data_patient()`](https://wlandau.github.io/retrogroup/reference/group_data_patient.md)

## Examples

``` r
  data <- group_example_summary()
  group_data_summary(
    data,
    variable = "variable",
    level = "level",
    patients = "patients",
    log_hazard_ratio = "log_hazard_ratio",
    standard_error = "standard_error"
  )
#> # A tibble: 6 × 5
#>   variable level      patients log_hazard_ratio standard_error
#>   <chr>    <chr>         <int>            <dbl>          <dbl>
#> 1 age      under 60        280           -0.169          0.178
#> 2 age      60 or over      339           -0.546          0.160
#> 3 node4    no              453           -0.417          0.153
#> 4 node4    yes             166           -0.312          0.190
#> 5 sex      female          312           -0.147          0.162
#> 6 sex      male            307           -0.656          0.179
  group_data_summary(
    data,
    variable = "variable",
    level = "level",
    patients = "patients",
    log_hazard_ratio = "log_hazard_ratio"
  )
#> # A tibble: 6 × 5
#>   variable level      patients log_hazard_ratio standard_error
#>   <chr>    <chr>         <int>            <dbl>          <dbl>
#> 1 age      under 60        280           -0.169             NA
#> 2 age      60 or over      339           -0.546             NA
#> 3 node4    no              453           -0.417             NA
#> 4 node4    yes             166           -0.312             NA
#> 5 sex      female          312           -0.147             NA
#> 6 sex      male            307           -0.656             NA
```
