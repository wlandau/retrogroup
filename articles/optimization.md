# Optimizing subgroup levels

``` r

library(retrogroup)
```

[`vignette("sampling", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/sampling.md)
defines methods for sampling subgroup levels based on posterior
probabilities. Over a large number of simulations, subgroup log hazard
ratios computed on the simulated datasets should match the published log
hazard ratios on average. However, each individual simulated dataset is
subject to random variation and may not exactly match the published
statistics. This vignette defines methods to refine an individual
simulated dataset so that the observed subgroup log hazard ratios match
the published log hazard ratios.

We define two optimization methods: one for subgroup variables produced
with unconstrained sampling, and another for those produced with
constrained sampling. Both methods are the same greedy local search.
Starting from sampled subgroup levels \\x = (x_1, \ldots, x_n)\\, each
iteration proposes a small random change to \\x\\ and keeps the change
only if it reduces a loss function. The loss function is shared. The
methods differ only in the kind of change they propose: a flip for
unconstrained sampling and a swap for constrained sampling.

## Loss function

The loss function measures the distance between the published subgroup
statistics and the same statistics observed under the current subgroup
levels \\x\\.

To get the observed statistics for subgroup level \\k\\, we fit a Cox
proportional hazards model to the patients currently assigned to that
level, \\\\i : x_i = k\\\\, with treatment \\z_i\\ as the only
covariate. The model is fit by partial likelihood, which leaves the
baseline hazard unspecified (Cox 1972, 1975), as in
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html). Let
\\\hat\beta_k(x)\\ be the estimated log hazard ratio of treatment versus
control and \\\hat\sigma_k(x)\\ its standard error.

If the publication reports only log hazard ratios, the loss is the
squared Euclidean distance between observed and published log hazard
ratios:

\\ \Delta(x) = \sum\_{k = 1}^K \left( \hat\beta_k(x) - \beta_k \right)^2
\\

If the standard errors \\\sigma_k\\ of the published log hazard ratios
from
[`vignette("framing", package = "retrogroup")`](https://wlandau.github.io/retrogroup/articles/framing.md)
are also available, the loss includes them:

\\ \Delta(x) = \sum\_{k = 1}^K \left\[ \left( \hat\beta_k(x) - \beta_k
\right)^2 + \left( \hat\sigma_k(x) - \sigma_k \right)^2 \right\] \\

The standard error depends mostly on the number of events in each arm,
so matching it steers the event counts in each subgroup level toward
those behind the published statistics.

If the model cannot be fit for some subgroup level, for example because
the level is empty or one of its arms has no events, we set \\\Delta(x)
= \infty\\.

## Optimization for unconstrained sampling

Unconstrained sampling draws subgroup levels independently, and we
expect the observed number of patients in each subgroup level to deviate
from the published subgroup counts \\n_1, \ldots, n_K\\. The
optimization does not preserve the counts either: each change is a flip,
which reassigns one patient to a different subgroup level, independently
of the other patients.

Let \\x^{(i \to k)}\\ denote the assignments that agree with \\x\\
except that patient \\i\\ belongs to subgroup level \\k\\. Each
iteration:

1.  Draws a patient \\i\\ uniformly at random from \\1, \ldots, n\\,
    independently of previous iterations and with replacement.
2.  Computes the counterfactual loss \\\Delta(x^{(i \to k)})\\ for every
    subgroup level \\k = 1, \ldots, K\\.
3.  Updates \\x \leftarrow x^{(i \to k^\*)}\\, where \\k^\* =
    \operatorname{argmin}\_k \Delta(x^{(i \to k)})\\.

The candidates include \\x^{(i \to x_i)}\\, which is just \\x\\
unchanged, so the loss never increases. If the current level ties for
the minimum, the patient stays put.

A flip from level \\x_i\\ to level \\k\\ changes the patients in those
two levels only, so every other term of the loss carries over from the
previous iteration. Each iteration needs \\K\\ model fits (step 2).

## Optimization for constrained sampling

Constrained sampling ensures that the observed number of patients in
each subgroup matches the published subgroup counts \\n_1, \ldots,
n_K\\. The optimization method for constrained sampling adjusts subgroup
assignments to approximate the published log hazard ratios while
maintaining the published subgroup counts. Subgroup counts are invariant
because each reassignment is a swap between subgroup levels of two
different patients.

Let \\x^{(i \leftrightarrow j)}\\ denote the assignments that agree with
\\x\\ except that patients \\i\\ and \\j\\ have exchanged subgroup
levels. Each iteration:

1.  Draws a pair of patients \\(i, j)\\ uniformly at random from the
    pairs with \\x_i \ne x_j\\, independently of previous iterations and
    with replacement.
2.  Computes the loss \\\Delta(x^{(i \leftrightarrow j)})\\ of the swap.
3.  Updates \\x \leftarrow x^{(i \leftrightarrow j)}\\ if \\\Delta(x^{(i
    \leftrightarrow j)}) \< \Delta(x)\\, and otherwise leaves \\x\\
    unchanged.

As with flips, the loss never increases. A swap changes the patients in
levels \\x_i\\ and \\x_j\\ only, so an iteration needs two model fits.

## Stopping

Both methods stop when the loss falls below a tolerance or after a
maximum number of iterations, whichever comes first. The loss is
non-increasing and bounded below by zero, so it converges. The search is
greedy, however, so it may settle in a local minimum where no single
flip or swap improves the loss. Rerunning from freshly sampled subgroup
levels gives a different starting point.

## References

Cox, D. R. 1972. “Regression Models and Life-Tables.” *Journal of the
Royal Statistical Society: Series B (Methodological)* 34 (2): 187–220.
<https://doi.org/10.1111/j.2517-6161.1972.tb00899.x>.

Cox, D. R. 1975. “Partial Likelihood.” *Biometrika* 62 (2): 269–76.
<https://doi.org/10.1093/biomet/62.2.269>.
