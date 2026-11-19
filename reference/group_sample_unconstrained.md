# Unconstrained subgroup sampling

Sample a level of each subgroup variable for each patient, independently
across patients, from the posterior probabilities estimated by
[`group_analysis_probability()`](https://wlandau.github.io/retrogroup/reference/group_analysis_probability.md).

## Usage

``` r
group_sample_unconstrained(patient, summary)
```

## Arguments

- patient:

  Patient-level survival data, the output of
  [`group_data_patient()`](https://wlandau.github.io/retrogroup/reference/group_data_patient.md).

- summary:

  Summary subgroup statistics, the output of
  [`group_data_summary()`](https://wlandau.github.io/retrogroup/reference/group_data_summary.md).

## Value

A tibble with one row per patient, in the same order as the rows of
`patient`. The columns are all the columns of `patient`, followed by one
column per subgroup variable, in the order the subgroup variables first
appear in `summary$variable`. Each subgroup column is a factor of the
sampled levels whose levels are those of the subgroup variable in
`summary$level`, in the same order.

## Details

`group_sample_unconstrained()` implements the "Unconstrained sampling"
section of
[`vignette("sampling", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/sampling.md).
For each subgroup variable and each patient \\i\\, it draws \\x_i \sim
\text{Categorical}(p\_{i1}, \ldots, p\_{iK})\\ independently, where the
\\p\_{ik}\\ come from
[`group_analysis_probability()`](https://wlandau.github.io/retrogroup/reference/group_analysis_probability.md).
All the requirements of
[`group_analysis_probability()`](https://wlandau.github.io/retrogroup/reference/group_analysis_probability.md)
on `patient` and `summary` apply. In addition, no subgroup variable in
`summary` may have the same name as a column of `patient`.

The sampled number of patients in each subgroup level is not guaranteed
to equal the published number of patients in `summary`. Use
[`group_sample_constrained()`](https://wlandau.github.io/retrogroup/reference/group_sample_constrained.md)
to preserve the published counts.

Sampling uses R's random number generator. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand for
reproducible results.

## See also

Other sampling:
[`group_sample_constrained()`](https://wlandau.github.io/retrogroup/reference/group_sample_constrained.md)

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
  sample <- group_sample_unconstrained(patient, summary)
  sample
#> # A tibble: 619 × 6
#>    treatment  time event age        node4 sex   
#>        <int> <dbl> <int> <fct>      <fct> <fct> 
#>  1         1  1521     1 60 or over no    female
#>  2         1  3087     0 under 60   no    female
#>  3         0   963     1 under 60   yes   male  
#>  4         1   293     1 60 or over no    female
#>  5         0   659     1 60 or over yes   female
#>  6         1  1767     1 under 60   no    female
#>  7         0  3192     0 60 or over no    female
#>  8         1  3308     0 60 or over no    male  
#>  9         1  3309     0 60 or over no    male  
#> 10         0  2085     1 60 or over no    male  
#> # ℹ 609 more rows
  table(sample$age)
#> 
#>   under 60 60 or over 
#>        276        343 
```
