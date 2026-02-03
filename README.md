
# bcmvr

**bcmvr** provides treatment effect inference for longitudinal outcomes
with skewed distributions using a Box–Cox multivariate regression
(BCMVR) model. Models are fitted separately by group, and prespecified
group comparisons are performed using:

- covariate-adjusted marginal medians (with Wald-type confidence
  intervals),
- median differences for ordered group pairs (CI and p-values),
- probability-based effect measures comparing full marginal
  distributions, reported as `P(Y_{g1,t} < Y_{g2,t})` with inference on
  the logit scale.

> **Direction matters.** In `pairs = list(c(g1, g2), ...)`, `g1` is the
> **test** group and `g2` is the **control** group.  
> Median differences are computed as `median(g1) - median(g2)`.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("kzkzmr/bcmvr")
```

## Quick start

This example uses `aidscd4` from **bcmixed**.

``` r
data(aidscd4, package = "bcmixed")
set.seed(1)

# sample ~140 subjects (IDs), keep all rows for those IDs
id_keep <- sample(unique(aidscd4$id), size = 140, replace = FALSE)
dat <- subset(aidscd4, id %in% id_keep)

# make treatment a factor so group labels are stable/readable
dat$treatment <- factor(dat$treatment)

# ordered pairs: g1=test, g2=control (here: compare 2/3/4 vs 1)
pairs <- list(
  c("2", "1"),
  c("3", "1"),
  c("4", "1")
)

# CD4 is typically "larger is better" -> smaller = FALSE
res <- bcmvr(
  dat,
  outcome = cd4,
  id      = id,
  group   = treatment,
  time    = weekc,
  cov     = c("age", "cd4.bl"),
  pairs   = pairs,
  smaller = FALSE,
  t_eval  = 1:4,
  robust  = TRUE
)

head(res$median)
head(res$meddiff)
head(res$prob)
```

## Notes

- `data` must be in **long format**, with one row per subject–time
  combination.
- `outcome` must be **positive**, as required by the Box–Cox
  transformation.
- `cov` must be **time-invariant within subject**.
- Covariates that are **constant within any group** are dropped with a
  warning to ensure comparability of covariate-adjusted effects across
  groups.
- Group comparisons are **directional**: in each pair `(g1, g2)`, `g1`
  is treated as the test group and `g2` as the control group.
- The probability-based effect measure is interpreted as
  `P(Y_{g1,t} < Y_{g2,t})` when `smaller = TRUE`, and
  `P(Y_{g1,t} > Y_{g2,t})` when `smaller = FALSE`.

## References

Maruo K et al. *Inference on treatment effects for longitudinal outcomes
with skewed distributions.* Submitting (available from author).
