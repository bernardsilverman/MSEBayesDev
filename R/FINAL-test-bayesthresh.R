test_that("Bayesian threshold regression agrees on UKdat_5", {
    skip_if_not_installed("MCMCpack")

    fit <- estimate_population_bayesthresh(
        UKdat_5,
        prior = "improper",
        maxorder = 3,
        mcmc = 10000,
        burnin = 1000,
        thin = 10,
        seed = 1234
    )

    expect_equal(
        fit$popest,
        12163.46,
        tolerance = 0.01
    )

    expect_equal(
        fit$quantiles,
        c(
            `2.5%` = 10609.25,
            `10%` = 11083.69,
            `50%` = 12163.46,
            `90%` = 13266.99,
            `97.5%` = 13694.40
        ),
        tolerance = 0.02
    )

    expect_identical(
        sort(fit$retained_interactions),
        sort(c(
            "GO:GP",
            "LA:NG",
            "LA:PFNCA",
            "NG:GP",
            "PFNCA:GP"
        ))
    )

    expect_length(fit$retained_triples, 0)
    expect_identical(fit$removed_pairs, "LA:GP")
})


test_that("zero sufficient statistic three-way effect is removed", {
    z <- expand.grid(
        A = 0:1,
        B = 0:1,
        C = 0:1,
        D = 0:1
    )

    z <- z[rowSums(z) > 0, , drop = FALSE]
    z$count <- seq_len(nrow(z)) + 3L

    # A:B:C has zero sufficient statistic, while all pairwise
    # sufficient statistics remain positive.
    z$count[
        z$A == 1 &
        z$B == 1 &
        z$C == 1
    ] <- 0L

    pairs <- .bayesthresh_all_effects(4, 2)
    triples <- .bayesthresh_all_effects(4, 3)

    pair_stats <- vapply(
        pairs,
        function(e)
            .bayesthresh_sufficient_stat(z, e),
        numeric(1)
    )

    triple_stats <- vapply(
        triples,
        function(e)
            .bayesthresh_sufficient_stat(z, e),
        numeric(1)
    )

    names(triple_stats) <- .bayesthresh_pretty(
        triples,
        z
    )

    expect_true(all(pair_stats > 0))
    expect_equal(triple_stats["A:B:C"], 0)

    setup <- .bayesthresh_remove_zero_effects(
        z,
        c(pairs, triples)
    )

    removed <- .bayesthresh_pretty(
        setup$removed,
        z
    )

    expect_true("A:B:C" %in% removed)
    expect_false(any(
        c("A:B", "A:C", "B:C") %in% removed
    ))

    expect_false(any(
        setup$data$A == 1 &
        setup$data$B == 1 &
        setup$data$C == 1
    ))
})
