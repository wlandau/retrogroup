# The motivating problem

``` r

library(retrogroup)
```

This vignette frames the problem `retrogroup` tries to solve.
Throughout, it introduces concepts and notation that the other vignettes
use.

## The goal

`retrogroup` takes as input:

1.  A patient-level survival dataset with columns for treatment, time,
    and event.
2.  A published table giving the log hazard ratio, number of patients,
    and number of events for each level of a subgroup variable.

The goal is to reconstruct the subgroup variable with patient-level
assignments that are consistent with the published subgroup statistics.
Reconstructed subgroup variables become part of datasets suitable for
subgroup analyses in clinical trial design and health technology
assessment.

The following sections describe the structure of the survival dataset
and the published subgroup statistics that `retrogroup` uses as input.

## Patient-level data

We begin with a right-censored patient-level survival dataset from a
randomized controlled parallel 2-arm clinical trial. For each patient
\\i = 1, \ldots, n\\, the dataset has one row per patient and columns
with the following values:

- \\z_i\\: 1 if randomized to the treatment arm, 0 if randomized to the
  control arm.
- \\t_i\\: time since randomization of the event or censoring.
- \\d_i\\: 1 if the event occurred, 0 if censored.

Our running example is the adjuvant colon cancer trial of Moertel et al.
(1990), as distributed in
[`survival::colon`](https://rdrr.io/pkg/survival/man/colon.html). We
keep the death endpoint and compare levamisole plus fluorouracil
(treatment) against observation (control), excluding the levamisole-only
arm.
[`group_example_patient()`](https://wlandau.github.io/retrogroup/reference/group_example_patient.md)
returns these patients, and
[`group_data_patient()`](https://wlandau.github.io/retrogroup/reference/group_data_patient.md)
converts them into the standardized format of `retrogroup`. Below,
variables `treatment`, `time` (days), and `event` correspond to \\z_i\\,
\\t_i\\, and \\d_i\\, respectively.

``` r

patient <- group_data_patient(
  group_example_patient(),
  arm = "arm",
  treatment = "Lev+5FU",
  time = "time",
  event = "event"
)
patient
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

In practice, a dataset may come from a clinical study directly, or it
may be reconstructed from published Kaplan-Meier curves and summary
statistics, as in Figure 1A of Gandhi et al. (2018). Tools like
`retroglyph` (Landau 2026), KM-GPT (Zhao et al. 2025), or
`SurvdigitizeR` (Zhang et al. 2024) can assist in this reconstruction.

## Summary-level data

Clinical studies often look for subgroups of patients with enhanced
treatment benefit relative to the overall population. A subgroup
variable is categorical with two or more levels, not randomized, and not
time-varying: for example, age group, smoking status, or biomarker
status. Oncology publications often report inferential and descriptive
statistics for each subgroup level, as in Figure 1B of Gandhi et al.
(2018).

For each level of the subgroup variable, we assume we know the estimated
log hazard ratio of treatment versus control and the number of patients.
For the colon cancer trial,
[`group_example_summary()`](https://wlandau.github.io/retrogroup/reference/group_example_summary.md)
provides these statistics from the patients above for subgroup variables
age, more than 4 positive lymph nodes (`node4`), and sex.
[`group_data_summary()`](https://wlandau.github.io/retrogroup/reference/group_data_summary.md)
converts them into the standardized format of `retrogroup`.

``` r

summary <- group_data_summary(
  group_example_summary(),
  variable = "variable",
  level = "level",
  patients = "patients",
  log_hazard_ratio = "log_hazard_ratio",
  standard_error = "standard_error"
)
```

| Variable | Level      | Patients | Log hazard ratio | Standard error |
|:---------|:-----------|---------:|-----------------:|---------------:|
| age      | under 60   |      280 |           -0.169 |          0.178 |
| age      | 60 or over |      339 |           -0.546 |          0.160 |
| node4    | no         |      453 |           -0.417 |          0.153 |
| node4    | yes        |      166 |           -0.312 |          0.190 |
| sex      | female     |      312 |           -0.147 |          0.162 |
| sex      | male       |      307 |           -0.656 |          0.179 |

For a given subgroup variable with \\K\\ levels, define:

- \\k = 1, \ldots, K\\ indexes the levels of the subgroup variable.
- \\n_k\\: number of patients in subgroup level \\k\\.
- \\\beta_k\\: published log hazard ratio of treatment versus control in
  subgroup level \\k\\.
- \\\sigma_k\\: published standard error of \\\beta_k\\, if available.

For the age subgroup variable above, we have \\K = 2\\ with levels \\k =
1\\ for under 60 years and \\k = 2\\ for 60 years or over. And we have
\\n_1 = 280\\, \\n_2 = 339\\, \\\beta_1 = -0.169\\, \\\beta_2 =
-0.546\\, \\\sigma_1 = 0.178\\, and \\\sigma_2 = 0.160\\.

We assume \\n = \sum\_{k = 1}^K n_k\\: in other words, the sum of the
subgroup counts equals the total number of patients. Each \\n_k\\ should
be strictly positive: unrepresented subgroup levels must be excluded
from the data reconstruction.

Strictly speaking, \\n_k\\ is the number of patients assigned to level
\\k\\ of the reconstructed subgroup variable. Published subgroup counts
are the natural starting point, but they do not always add up to the
total number of patients in the study, for example because of missing
subgroup data or reporting errors. In that case, choosing values of
\\n_k\\ that satisfy \\n = \sum\_{k = 1}^K n_k\\ is a decision about the
reconstruction, not an attempt to recover the true counts from the
clinical study. Likewise, if the publication does not report the \\n_k\\
values, a reasonable approach is to supply roughly equal counts across
levels, or to try several configurations in a sensitivity analysis.
Either way, any downstream subgroup analyses should explicitly state the
subgroup counts used and how they were chosen.

The standard errors \\\sigma_k\\ of the log hazard ratios may or may not
be reported in the publication. If they are not reported, they can often
be obtained from published confidence intervals after transforming the
bounds to the log scale, or the number of events in the subgroup. If
available, they aid the optimization of subgroup assignments in
[`vignette("optimization", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/optimization.md).
However, they are not strictly required.

The methods in `retrogroup` focus on a single subgroup variable at a
time. To generate multiple subgroup variables, the package applies the
method separately for each variable.

## References

Gandhi, Leena, Delvys Rodríguez-Abreu, Shirish Gadgeel, et al. 2018.
“Pembrolizumab Plus Chemotherapy in Metastatic Non–Small-Cell Lung
Cancer.” *New England Journal of Medicine* 378 (22): 2078–92.
<https://doi.org/10.1056/NEJMoa1801005>.

Landau, William Michael. 2026. *retroglyph: A Trusted Mini-Agent to
Reconstruct Individual Patient Data from Kaplan-Meier Plots*.
<https://wlandau.github.io/retroglyph/>.

Moertel, Charles G., Thomas R. Fleming, John S. Macdonald, et al. 1990.
“Levamisole and Fluorouracil for Adjuvant Therapy of Resected Colon
Carcinoma.” *New England Journal of Medicine* 322 (6): 352–58.
<https://doi.org/10.1056/NEJM199002083220602>.

Zhang, Jasper Zhongyuan, Juan David Rios, Tilemanchos Pechlivanoglou, et
al. 2024. “SurvdigitizeR: An Algorithm for Automated Survival Curve
Digitization.” *BMC Medical Research Methodology* 24 (1): 147.
<https://doi.org/10.1186/s12874-024-02273-8>.

Zhao, Yao, Haoyue Sun, Yantian Ding, and Yanxun Xu. 2025. *KM-GPT: An
Automated Pipeline for Reconstructing Individual Patient Data from
Kaplan-Meier Plots*. <https://arxiv.org/abs/2509.18141>.
