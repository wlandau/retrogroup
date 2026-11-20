# retrogroup 0.0.3

## Breaking changes

* `group_data_summary()` takes the number of events (`events`) and the hazard ratio on the natural scale (`hazard_ratio`) instead of `log_hazard_ratio` and `standard_error`.
* `group_example_summary()` and `group_analysis_summary()` report `events` and `hazard_ratio` instead of `log_hazard_ratio` and `standard_error`.
  `group_example_summary()` also reports the bounds `lower` and `upper` of the 95% confidence interval of each hazard ratio.
* The optimization loss is now the Euclidean distance between the observed and assumed hazard ratios on the linear scale, without standard errors, and between the observed and assumed event fractions.
  The event fraction of a subgroup level is its number of events divided by its assumed number of patients in `summary$patients`, so the event fraction terms steer the event counts toward `summary$events` as a soft constraint.
  `tolerance` is on the linear scale of the hazard ratio and on the scale of the event fractions.
* The default `iterations` of `group_optimize_constrained()` and `group_optimize_unconstrained()` is now 10000 instead of 1000.
  Iterations are cheap, and 1000 was often too few for the colon cancer example to converge.
* `group_sample_constrained()` and `group_optimize_constrained()` preserve the number of patients in each subgroup level.
  `group_optimize_constrained()` and `group_optimize_unconstrained()` steer the number of events in each subgroup level toward `summary$events`, but only as a soft constraint.
* Each iteration of `group_optimize_unconstrained()` now takes a flip step with probability 1/2 and a swap step otherwise, instead of always taking a flip step.
  With flips alone, the search often stalled in a local minimum above the tolerance, for example for `node4` in the colon cancer example.
  Results for a given seed differ from previous versions.

## New features

* Add `group_data_events()` to approximate subgroup event counts from the confidence intervals of the hazard ratios (Schoenfeld 1983; Parmar et al. 1998).
  `group_data_summary(events = NULL)` prepares its input.
* Add the appendix vignette "Approximating event counts using confidence intervals of assumed hazard ratios".

## Documentation

* Clarify that the event counts carry some general prognostic information, but the timing of events within subgroup levels is not reproduced, and the statistician who uses a reconstructed dataset is responsible for the consequences.
* Add the verification vignette "Prognostic benchmarks".
  It tests how well sampled and optimized subgroups recover the true prognostic hazard ratios of the colon cancer subgroups.

# retrogroup 0.0.2

* Prepare for open-source release: adopt the MIT license and move the source code and documentation website to <https://github.com/wlandau/retrogroup> and <https://wlandau.github.io/retrogroup/>.

# retrogroup 0.0.1

* First version.
