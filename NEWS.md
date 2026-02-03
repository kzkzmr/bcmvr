# bcmvr 0.1.0

## Initial release

- Initial CRAN-style release of the **bcmvr** package.
- Implements treatment effect inference for longitudinal outcomes using
the Box–Cox multivariate regression (BCMVR) model.
- Provides covariate-adjusted inference based on:
  - group-specific marginal medians with Wald-type confidence intervals,
- median differences between prespecified group pairs,
- probability-based treatment effect measures
\( P(Y_{g1,t} < Y_{g2,t}) \) with logit-scale inference.
- Supports robust (sandwich) and model-based variance estimation.
- Includes a single high-level user-facing function `bcmvr()` that performs
model fitting, inference, and group comparisons in one step.
- Example dataset usage demonstrated via `aidscd4` from the **bcmixed** package.
