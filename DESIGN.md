# Bayesian MSE design

## Aim

Develop a hierarchical Bayesian-threshold estimator for multiple-systems estimation, inspired by Section 6 of Silverman (2020), but designed from scratch to use the machinery and ideas of MultipleSystemsEstimation wherever appropriate.

## Initial principles

- Poisson log-linear model for capture-pattern counts.
- Main effects always included.
- Interactions allowed up to a user-specified `maxorder`.
- Shrinkage priors on interaction terms.
- Threshold interactions using posterior evidence.
- Preserve hierarchy when retaining higher-order interactions.
- Refit after thresholding.
- Treat sparse/zero-overlap structures using the existence and hierarchy machinery already developed for MultipleSystemsEstimation where possible.
- Return posterior population estimates, credible intervals, retained model structure, and diagnostics.
