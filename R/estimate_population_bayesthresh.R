#' Bayesian-threshold multiple-systems estimation
#'
#' Fits the Bayesian threshold estimator of Silverman (2020), with optional
#' three-way interactions.  Pairwise interactions are thresholded first.
#' A three-way interaction is considered only when all three constituent
#' pairwise interactions survive the pairwise threshold.
#'
#' With an improper prior, any interaction having zero sufficient statistic
#' is assigned the limiting value minus infinity, removed from the fitted
#' model, and all likelihood cells containing that interaction are removed.
#' The same rule applies to pairwise and three-way interactions.
#'
#' @param zdat Multiple-systems data in the usual
#'   MultipleSystemsEstimation format.
#' @param prior Either `"proper"` or `"improper"`.
#' @param prior_variance Prior variance for interaction terms when a proper
#'   prior is used.
#' @param main_variance Prior variance for the intercept and main effects when
#'   a proper prior is used.
#' @param threshold Threshold applied to the absolute posterior mean divided
#'   by posterior standard deviation for interaction parameters.
#' @param maxorder Maximum interaction order, either 2 or 3.
#' @param ... Additional arguments passed to
#'   \code{MCMCpack::MCMCpoisson()}.
#'
#' @return A list containing the posterior median population estimate,
#'   posterior quantiles and draws, retained pairwise and three-way
#'   interactions, threshold statistics, and effects removed under an
#'   improper prior.
#'
#' @export
estimate_population_bayesthresh <- function(
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

    # Stage 1: fit all pairwise interactions.
    pair_candidates <- .bayesthresh_all_effects(
        nlists,
        2
    )

    fit1 <- .bayesthresh_fit(
        zfull,
        pair_candidates,
        prior,
        prior_variance,
        main_variance,
        prune_improper = TRUE,
        ...
    )

    pair_threshold <- .bayesthresh_threshold(
        fit1$fit,
        fit1$effects,
        2,
        threshold
    )

    pairs <- pair_threshold$retained

    # For an improper prior, all later fits start from the likelihood
    # table pruned at the first stage.
    zrefit <- if (prior == "improper") fit1$data else zfull

    # Stage 2: refit the retained pairwise model.
    fit2 <- .bayesthresh_fit(
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
    removed_triples <- character(0)

    # Stage 3: consider admissible three-way effects.
    if (maxorder == 3) {
        triple_candidates <- .bayesthresh_admissible_triples(
            pairs,
            nlists
        )

        if (length(triple_candidates)) {
            fit3 <- .bayesthresh_fit(
                zrefit,
                c(pairs, triple_candidates),
                prior,
                prior_variance,
                main_variance,
                prune_improper = (prior == "improper"),
                ...
            )

            removed_triples <- fit3$removed[
                .bayesthresh_order(fit3$removed) == 3
            ]

            triple_threshold <- .bayesthresh_threshold(
                fit3$fit,
                fit3$effects,
                3,
                threshold
            )

            triples <- triple_threshold$retained

            final_fit <- .bayesthresh_fit(
                fit3$data,
                c(pairs, triples),
                prior,
                prior_variance,
                main_variance,
                prune_improper = FALSE,
                ...
            )
        }
    }

    population <- .bayesthresh_population(
        final_fit$fit,
        zfull
    )

    removed_pairs <- fit1$removed[
        .bayesthresh_order(fit1$removed) == 2
    ]

    list(
        popest = unname(population$quantiles["50%"]),
        quantiles = population$quantiles,
        posterior = population$total_population,

        retained_interactions =
            .bayesthresh_pretty(pairs, zfull),

        retained_triples =
            .bayesthresh_pretty(triples, zfull),

        threshold_statistics =
            .bayesthresh_pretty_named(
                pair_threshold$ratios,
                zfull
            ),

        triple_threshold_statistics =
            if (is.null(triple_threshold))
                numeric(0)
            else
                .bayesthresh_pretty_named(
                    triple_threshold$ratios,
                    zfull
                ),

        removed_pairs =
            .bayesthresh_pretty(
                removed_pairs,
                zfull
            ),

        removed_triples =
            .bayesthresh_pretty(
                removed_triples,
                zfull
            ),

        removed_effects =
            .bayesthresh_pretty(
                unique(c(removed_pairs, removed_triples)),
                zfull
            ),

        admissible_triples =
            .bayesthresh_pretty(
                triple_candidates,
                zfull
            )
    )
}
