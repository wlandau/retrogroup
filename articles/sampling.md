# Sampling subgroups

``` r

library(retrogroup)
```

[`vignette("model", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/model.md)
derives \\p\_{ik}\\, an estimate of the posterior probability that
patient \\i\\ belongs to subgroup level \\k\\ (\\i = 1, \ldots, n\\, \\k
= 1, \ldots, K\\). This vignette explains a method to sample subgroup
levels \\x_1, \ldots, x_n\\ given fixed posterior probability estimates
\\p\_{ik}\\.

Over a large number of simulations, subgroup log hazard ratios computed
on the simulated datasets should match the published log hazard ratios
on average. However, each individual simulated dataset is subject to
random variation and may not exactly match the published statistics.
[`vignette("optimization", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/optimization.md)
defines methods to refine an individual dataset so the observed and
published subgroup log hazard ratios match.

## Unconstrained sampling

If we treat subgroup levels as unconstrained, each \\x_i\\ has its own
independent categorical distribution with probabilities \\p\_{i1},
\ldots, p\_{iK}\\. Independently for each patient \\i = 1, \ldots, n\\,
we draw:

\\ x_i \sim \text{Categorical}(p\_{i1}, \ldots, p\_{iK}) \\

Unconstrained sampling is simple, but it does not guarantee that the
resulting subgroup assignments match \\n_k\\, the published number of
patients in each subgroup level \\k\\ from
[`vignette("framing", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/framing.md).
For that, we need constrained sampling.

## Constrained sampling

Constrained sampling forces the sampled subgroup levels to match the
published counts \\n_k\\ for each subgroup level \\k\\. To preserve the
counts, we iterate through a random permutation of patients and assign
each one to a random subgroup level with fewer than \\n_k\\ patients.

Let \\\rho_1, \ldots, \rho_n\\ be a random permutation of patient
indexes \\1, \ldots, n\\. Every \\n_k\\ is strictly positive, so every
subgroup level has room for the first patient, and we draw
\\x\_{\rho_1}\\ from a categorical distribution:

\\ x\_{\rho_1} \sim \text{Categorical}(p\_{{\rho_1}1}, \ldots,
p\_{{\rho_1}K}) \\

At arbitrary patient \\\rho_j\\, patients \\\rho_1, \ldots,
\rho\_{j-1}\\ have already been assigned to subgroup levels, with
\\n^\*\_{jk} := \sum\_{h = 1}^{j - 1} I(x\_{\rho_h} = k)\\ patients in
subgroup level \\k\\ so far. We draw \\x\_{\rho_j}\\ from a categorical
distribution that excludes subgroup levels that have already reached
their published counts:

\\ x\_{\rho_j} \sim \text{Categorical}(q\_{{\rho_j}1}, \ldots,
q\_{{\rho_j}K}) \\

where:

\\ q\_{{\rho_j}k} := \begin{cases} \dfrac{p\_{{\rho_j}k}}{\sum\_{l =
1}^K I(n^\*\_{jl} \< n_l) p\_{{\rho_j}l}}, & \text{if } n^\*\_{jk} \<
n_k \\ 0, & \text{if } n^\*\_{jk} = n_k \end{cases} \\

We assume \\\sum\_{k = 1}^K n_k = n\\, so at least one subgroup level
has room at every step. Every \\p\_{ik}\\ is strictly positive, as shown
in
[`vignette("model", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/model.md),
so the denominator above is positive. This scheme guarantees the
published counts, but it is only an approximate draw from independent
categorical distributions conditional on those counts, because patients
late in the permutation can only join subgroup levels that still have
room.
