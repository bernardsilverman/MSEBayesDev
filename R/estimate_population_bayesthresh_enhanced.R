#' Estimate population size by Bayesian thresholding
#' with order-general improper-prior handling
#'
#' Development version extending estimate_population_bayesthresh_simple().
#' It agrees with the verified parallel implementation for pairwise fits,
#' but when three-way effects are introduced under an improper prior it
#' additionally removes any zero-sufficient-statistic triple effects and
#' the corresponding cells from the likelihood.
#'
#' @param zdat Multiple-systems data in the usual
#'   MultipleSystemsEstimation format.
#' @param prior Either "proper" or "improper".
#' @param prior_variance Prior variance for interaction terms when a proper
#'   prior is used.
#' @param main_variance Prior variance for the intercept and main effects when
#'   a proper prior is used.
#' @param threshold Threshold applied to abs(posterior mean / posterior SD).
#' @param maxorder Maximum interaction order, either 2 or 3.
#' @param ... Additional arguments passed to MCMCpack::MCMCpoisson().
#'
#' @export
estimate_population_bayesthresh_enhanced <- function(
    zdat,
    prior = c("proper", "improper"),
    prior_variance = 1,
    main_variance = 1e4,
    threshold = 2,
    maxorder = 2,
    ...
) {
    prior <- match.arg(prior)

    if (!maxorder %in% c(2, 3))
        stop("maxorder must be either 2 or 3.")

    zfull <- MultipleSystemsEstimation::tidy_lists(
        zdat,
        includezerocounts = TRUE
    )

    nlists <- ncol(zfull) - 1

    if (maxorder == 3 && nlists < 3)
        stop("Three-way interactions require at least three lists.")

    ## Stage 1: all pairwise interactions.
    pair_candidates <- .bayesthresh2_all_effects(nlists, 2)

    fit1 <- .bayesthresh2_fit(
        zfull,
        pair_candidates,
        prior,
        prior_variance,
        main_variance,
        prune_improper = TRUE,
        ...
    )

    pair_threshold <- .bayesthresh2_threshold(
        fit1$fit,
        fit1$effects,
        2,
        threshold
    )

    pairs <- pair_threshold$retained

    ## The pairwise refit uses the once-pruned table, exactly as in the
    ## verified parallel implementation.
    zrefit <- if (prior == "improper") fit1$data else zfull

    fit2 <- .bayesthresh2_fit(
        zrefit,
        pairs,
        prior,
        prior_variance,
        main_variance,
        prune_improper = FALSE,
        ...
    )

    final_fit <- fit2
    final_data <- zrefit

    triple_candidates <- character(0)
    triple_threshold <- NULL
    triples <- character(0)
    removed_triples <- character(0)
    fit3 <- NULL

    ## Stage 3: admissible three-way effects.
    if (maxorder == 3) {
        triple_candidates <- .bayesthresh2_admissible_triples(
            pairs,
            nlists
        )

        if (length(triple_candidates)) {

            ## This is the enhancement.  Under an improper prior, inspect
            ## the candidate three-way effects for zero sufficient statistics
            ## and prune them (and their active zero-count cells) before fit 3.
            fit3 <- .bayesthresh2_fit(
                zrefit,
                c(pairs, triple_candidates),
                prior,
                prior_variance,
                main_variance,
                prune_improper = (prior == "improper"),
                ...
            )

            removed_triples <- fit3$removed[
                .bayesthresh2_order(fit3$removed) == 3
            ]

            triple_threshold <- .bayesthresh2_threshold(
                fit3$fit,
                fit3$effects,
                3,
                threshold
            )

            triples <- triple_threshold$retained
            final_data <- fit3$data

            final_fit <- .bayesthresh2_fit(
                final_data,
                c(pairs, triples),
                prior,
                prior_variance,
                main_variance,
                prune_improper = FALSE,
                ...
            )
        }
    }

    population <- .bayesthresh2_population(
        final_fit$fit,
        zfull
    )

    list(
        popest = unname(population$quantiles["50%"]),
        quantiles = population$quantiles,
        posterior = population$total_population,

        retained_interactions =
            .bayesthresh2_pretty(pairs, zfull),

        retained_triples =
            .bayesthresh2_pretty(triples, zfull),

        threshold_statistics =
            .bayesthresh2_pretty_named(
                pair_threshold$ratios,
                zfull
            ),

        triple_threshold_statistics =
            if (is.null(triple_threshold))
                numeric(0)
            else
                .bayesthresh2_pretty_named(
                    triple_threshold$ratios,
                    zfull
                ),

        removed_pairs =
            .bayesthresh2_pretty(
                fit1$removed[
                    .bayesthresh2_order(fit1$removed) == 2
                ],
                zfull
            ),

        removed_triples =
            .bayesthresh2_pretty(
                removed_triples,
                zfull
            ),

        removed_effects =
            .bayesthresh2_pretty(
                unique(c(fit1$removed, removed_triples)),
                zfull
            ),

        admissible_triples =
            .bayesthresh2_pretty(
                triple_candidates,
                zfull
            ),

        # Temporary intermediate objects for checking.
        first_fit = fit1$fit,
        second_fit = fit2$fit,
        third_fit =
            if (is.null(fit3)) NULL else fit3$fit,
        final_fit = final_fit$fit
    )
}
