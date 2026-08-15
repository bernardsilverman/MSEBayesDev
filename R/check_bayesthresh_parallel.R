# Parallel checks for the simplified Bayesian-threshold implementation
#
# Source this only during development.  It does not form part of the
# package API.  The same MCMC settings and seed should be used for the
# existing and simplified functions.

compare_bayesthresh <- function(zdat, ..., tolerance = 1e-8) {
    old <- estimate_population_bayesthresh(zdat, ...)
    new <- estimate_population_bayesthresh_simple(zdat, ...)

    old_pairs <- sort(old$retained_interactions)
    new_pairs <- sort(new$retained_interactions)

    old_triples <- if (is.null(old$retained_triples))
        character(0) else sort(old$retained_triples)
    new_triples <- sort(new$retained_triples)

    list(
        popest_old = old$popest,
        popest_new = new$popest,
        popest_difference = new$popest - old$popest,
        quantiles_old = old$quantiles,
        quantiles_new = new$quantiles,
        quantiles_equal = isTRUE(all.equal(
            old$quantiles, new$quantiles, tolerance = tolerance
        )),
        pairs_old = old_pairs,
        pairs_new = new_pairs,
        pairs_equal = identical(old_pairs, new_pairs),
        triples_old = old_triples,
        triples_new = new_triples,
        triples_equal = identical(old_triples, new_triples),
        removed_effects_new = new$removed_effects
    )
}


# Suggested first checks:
#
# set.seed(1234)
# compare_bayesthresh(
#     Korea,
#     prior = "proper",
#     maxorder = 2,
#     mcmc = 10000,
#     burnin = 1000,
#     thin = 10,
#     seed = 1234
# )
#
# Then repeat on a 4+-list dataset for maxorder = 3, and on data with
# zero sufficient statistics using prior = "improper".
#
# IMPORTANT:
# MCMCpack's seed argument, rather than set.seed() alone, should be kept
# identical between the two calls when exact draw-by-draw reproducibility
# is being investigated.
