# Patient-level survival data

Convert a patient-level right-censored survival dataset from a
randomized 2-arm trial into the standardized format used by
`retrogroup`.

## Usage

``` r
group_data_patient(data, arm, treatment, time, event)
```

## Arguments

- data:

  Data frame with one row per patient. Must have the columns named by
  `arm`, `time`, and `event`. Other columns are ignored.

- arm:

  Character string, name of the column in `data` with the treatment arm.
  Must have no missing values and exactly two unique values.

- treatment:

  The value of `data[[arm]]` that denotes the treatment arm. The other
  value denotes the control arm.

- time:

  Character string, name of the column in `data` with the time since
  randomization of the event or censoring. The column must be numeric,
  finite, non-negative, and free of missing values.

- event:

  Character string, name of the column in `data` with the event
  indicator: `1` or `TRUE` for an event, `0` or `FALSE` for right
  censoring. Must have no missing values and at least one event.

## Value

A tibble with one row per patient and the following columns:

- `treatment`: integer, `1` for the treatment arm and `0` for the
  control arm.

- `time`: numeric, time since randomization of the event or censoring.

- `event`: integer, `1` if the event occurred and `0` if censored.

## See also

Other data:
[`group_data_summary()`](https://wlandau.github.io/retrogroup/reference/group_data_summary.md)

## Examples

``` r
  data <- group_example_patient()
  group_data_patient(
    data,
    arm = "arm",
    treatment = "Lev+5FU",
    time = "time",
    event = "event"
  )
#> # A tibble: 619 × 3
#>    treatment  time event
#>        <int> <dbl> <int>
#>  1         1  1521     1
#>  2         1  3087     0
#>  3         0   963     1
#>  4         1   293     1
#>  5         0   659     1
#>  6         1  1767     1
#>  7         0  3192     0
#>  8         1  3308     0
#>  9         1  3309     0
#> 10         0  2085     1
#> # ℹ 609 more rows
```
