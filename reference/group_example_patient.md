# Example patient-level survival data

Return example patient-level data in the format expected by
[`group_data_patient()`](https://wlandau.github.io/retrogroup/reference/group_data_patient.md).

## Usage

``` r
group_example_patient()
```

## Value

A tibble with one row per patient and the following columns:

- `arm`: character, study arm, either `"Obs"` or `"Lev+5FU"`.

- `time`: numeric, days from randomization to death or censoring.

- `event`: integer, `1` if the patient died and `0` if censored.

## Details

The data come from the adjuvant colon cancer trial of levamisole plus
fluorouracil (Moertel et al., 1990), as distributed in
[survival::colon](https://rdrr.io/pkg/survival/man/colon.html). The
subset keeps the death endpoint (`etype == 2`) and two arms: observation
(`"Obs"`, control) and levamisole plus fluorouracil (`"Lev+5FU"`,
treatment). The levamisole-only arm is excluded. See
[survival::colon](https://rdrr.io/pkg/survival/man/colon.html) for the
full dataset and its documentation.

## References

Moertel CG, Fleming TR, Macdonald JS, et al. (1990). Levamisole and
fluorouracil for adjuvant therapy of resected colon carcinoma. New
England Journal of Medicine, 322(6), 352-358.

## See also

Other example:
[`group_example_summary()`](https://wlandau.github.io/retrogroup/reference/group_example_summary.md)

## Examples

``` r
group_example_patient()
#> # A tibble: 619 × 3
#>    arm      time event
#>    <chr>   <dbl> <int>
#>  1 Lev+5FU  1521     1
#>  2 Lev+5FU  3087     0
#>  3 Obs       963     1
#>  4 Lev+5FU   293     1
#>  5 Obs       659     1
#>  6 Lev+5FU  1767     1
#>  7 Obs      3192     0
#>  8 Lev+5FU  3308     0
#>  9 Lev+5FU  3309     0
#> 10 Obs      2085     1
#> # ℹ 609 more rows
```
