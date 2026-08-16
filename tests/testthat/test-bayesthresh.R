test_that("Bayesian threshold regression agrees on UKdat_5", {
    skip_if_not_installed("MCMCpack")

  data("UKdat_5", package = "MultipleSystemsEstimation")

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
    expect_identical(fit$minus_infinite_estimated_effects, "LA:GP")
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

    # A:B:C has zero sufficient statistic, while all its pairwise
    # sufficient statistics remain positive.
    z$count[
        z$A == 1 &
        z$B == 1 &
        z$C == 1
    ] <- 0L

    sufficient_stat <- function(vars) {
        active <- apply(
            as.matrix(z[, vars, drop = FALSE]),
            1,
            prod
        )
        sum(z$count * active)
    }

    expect_true(all(c(
        sufficient_stat(c("A", "B")),
        sufficient_stat(c("A", "C")),
        sufficient_stat(c("B", "C"))
    ) > 0))

    expect_equal(
        sufficient_stat(c("A", "B", "C")),
        0
    )

    pairs <- .bayesthresh_all_effects(4, 2)
    triples <- .bayesthresh_all_effects(4, 3)

    setup <- .bayesthresh_remove_zero_effects(
        z,
        c(pairs, triples)
    )

    expect_true("x1:x2:x3" %in% setup$removed)

    expect_false(any(
        c("x1:x2", "x1:x3", "x2:x3") %in% setup$removed
    ))

    expect_false(any(
        setup$data$A == 1 &
        setup$data$B == 1 &
        setup$data$C == 1
    ))
})
