#' Estimate population size by Bayesian thresholding
#' (parallel simplified implementation)
#'
#' This function is kept under a separate name while it is checked against
#' the existing implementation.  It preserves the existing sequence of fits
#' and, for an improper prior, prunes the likelihood table once at the first
#' fit and reuses that same pruned table in all later refits.
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
estimate_population_bayesthresh_simple <- function(
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

    # Preserve the present implementation during parallel checking.
    if (maxorder == 3 && nlists < 4)
        stop("Three-way interactions require at least four lists.")

    ## First fit: all pairwise interactions.
    pair_candidates <- .bayesthresh2_all_effects(
        nlists,
        2
    )

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

    ## For an improper prior, the existing implementation carries forward
    ## the table pruned at the first fit.  Do exactly the same here.
    zrefit <- if (prior == "improper") fit1$data else zfull

    ## Second fit: retained pairwise interactions only.
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

    triple_candidates <- character(0)
    triple_threshold <- NULL
    triples <- character(0)
    fit3 <- NULL

    ## Optional three-way stage.
    if (maxorder == 3) {
        triple_candidates <- .bayesthresh2_admissible_triples(
            pairs,
            nlists
        )

        if (length(triple_candidates)) {
            fit3 <- .bayesthresh2_fit(
                zrefit,
                c(pairs, triple_candidates),
                prior,
                prior_variance,
                main_variance,
                prune_improper = FALSE,
                ...
            )

            triple_threshold <- .bayesthresh2_threshold(
                fit3$fit,
                fit3$effects,
                3,
                threshold
            )

            triples <- triple_threshold$retained

            final_fit <- .bayesthresh2_fit(
                zrefit,
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

        removed_effects =
            .bayesthresh2_pretty(
                fit1$removed,
                zfull
            ),

        admissible_triples =
            .bayesthresh2_pretty(
                triple_candidates,
                zfull
            ),

        # Temporary intermediate objects for parallel checking.
        first_fit = fit1$fit,
        second_fit = fit2$fit,
        third_fit =
            if (is.null(fit3)) NULL else fit3$fit,
        final_fit = final_fit$fit
    )
}
