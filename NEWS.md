# retrogroup 0.0.3

## Breaking changes

* `group_data_summary()` takes the number of events (`events`) and the hazard ratio on the natural scale (`hazard_ratio`) instead of `log_hazard_ratio` and `standard_error`.
* `group_example_summary()` and `group_analysis_summary()` report `events` and `hazard_ratio` instead of `log_hazard_ratio` and `standard_error`.
  `group_example_summary()` also reports the bounds `lower` and `upper` of the 95% confidence interval of each hazard ratio.
* The optimization loss is now the Euclidean distance between the observed and assumed hazard ratios on the linear scale, without standard errors.
  `tolerance` is on the linear scale of the hazard ratio.
* `group_sample_constrained()` and `group_optimize_constrained()` preserve the number of events in each subgroup level as well as the number of patients.

## New features

* Add `group_data_events()` to approximate subgroup event counts from the confidence intervals of the hazard ratios (Schoenfeld 1983; Parmar et al. 1998).
  `group_data_summary(events = NULL)` prepares its input.
* Add the appendix vignette "Approximating event counts using confidence intervals of assumed hazard ratios".

## Documentation

* Clarify that the event counts carry some general prognostic information, but the timing of events within subgroup levels is not reproduced, and the statistician who uses a reconstructed dataset is responsible for the consequences.

# retrogroup 0.0.2

* Prepare for open-source release: adopt the MIT license and move the source code and documentation website to <https://github.com/wlandau/retrogroup> and <https://wlandau.github.io/retrogroup/>.

# retrogroup 0.0.1

* First version.
