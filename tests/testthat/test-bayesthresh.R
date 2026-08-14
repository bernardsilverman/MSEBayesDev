test_that("bayesian threshold reproduces UK pairwise structure", {

  data("UKdat", package = "MultipleSystemsEstimation")
  data("Kosovo", package = "MultipleSystemsEstimation")

  res <- estimate_population_bayesthresh(
    UKdat,
    prior = "proper",
    prior_variance = 1,
    threshold = 2,
    burnin = 1000,
    mcmc = 5000,
    thin = 1,
    seed = 12345
  )

  expect_equal(
    res$model,
    "[12,13,25,35,36,45]"
  )

  expect_setequal(
    res$retained_interactions,
    c(
      "LA:NG",
      "LA:PF",
      "NG:GP",
      "PF:GP",
      "GO:GP",
      "PF:NCA"
    )
  )
})

test_that("bayesian threshold reproduces Kosovo pairwise structure", {

  res <- estimate_population_bayesthresh(
    Kosovo,
    prior = "proper",
    prior_variance = 1,
    threshold = 2,
    burnin = 1000,
    mcmc = 5000,
    thin = 1,
    seed = 12345
  )

  expect_equal(
    res$model,
    "[12,13,14,23,34]"
  )
})
