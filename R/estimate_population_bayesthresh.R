#' Bayesian-threshold multiple-systems estimation
#'
#' @description
#' Fits the Bayesian threshold estimator of Silverman (2020), based on a
#' Poisson log-linear model for the capture-pattern counts.
#' A Bayesian/thresholding approach is adopted, with priors placed on the parameters
#' in the log-linear model. Improper priors are used for the intercept and
#' main effects.
#' For the Bayesian parts of the procedure, inference is carried out by Markov chain
#' Monte Carlo, calling \code{MCMCpack::MCMCpoisson()}.
#'
#' The method begins by including all pairwise interactions. The interactions are
#' then thresholded, by discarding all those whose posterior mean to standard
#' deviation ratio has absolute value below \code{threshold}. The model
#' containing the retained interactions is then re-estimated, and the posterior
#' distribution of the total population is obtained by adding the observed
#' population to the posterior estimate of the unobserved cell.
#'
#' @details
#' By default a proper normal prior is used for the log-linear parameters
#' corresponding to pairwise and higher-order effects.
#' Setting \code{prior = "improper"} instead uses an improper flat prior. In
#' that case interactions having zero sufficient statistic have posterior
#' distribution concentrated at minus infinity. Such effects are retained
#' in the model but are accounted for separately before carrying out the
#' actual MCMC, as set out in Silverman (2020).
#'
#' If \code{maxorder = 3}, three-way interactions are also considered.
#' Pairwise interactions are thresholded first. A
#' three-way interaction is then eligible for consideration only if all three of
#' its constituent pairwise interactions have been retained. The model
#' containing the retained pairs and all eligible three-way interactions is
#' then fitted, the three-way interactions are thresholded in the same way,
#' and the resulting hierarchical model is refitted.
#'
#' In all cases, the full pairwise model is checked for identifiability and
#' for the Fienberg-Rinaldo existence criterion before thresholding begins. If
#' it fails, the procedure stops. With \code{maxorder = 3}, if the model
#' containing all eligible three-way interactions fails the Fienberg-Rinaldo
#' criterion, the third-order fit is not carried out and the completed order-2
#' analysis is returned, together with the eligible three-way interactions.
#'
#' @param zdat Multiple-systems data in the usual
#'   MultipleSystemsEstimation format.
#' @param prior Either \code{"proper"} (the default) or \code{"improper"}.
#' @param prior_variance Prior variance for interaction terms when a proper
#'   prior is used.
#' @param main_variance Prior variance for the intercept and main effects when
#'   a proper prior is used.
#' @param threshold Threshold applied to the absolute posterior mean to
#'   standard deviation ratio for interaction parameters.
#' @param maxorder Maximum interaction order, either 2 or 3.
#' @param return_posterior Logical. If \code{TRUE}, include the full posterior
#'   sample of the total population in the returned object. The default is
#'   \code{FALSE}.
#' @param ... Additional arguments passed to \code{MCMCpack::MCMCpoisson()}.
#'   Useful arguments include:
#'   \itemize{
#'     \item \code{burnin}: number of burn-in iterations; default \code{1000}.
#'     \item \code{mcmc}: number of MCMC iterations retained after burn-in;
#'       default \code{10000}.
#'     \item \code{thin}: thinning interval; default \code{1}.
#'     \item \code{tune}: Metropolis tuning parameter; default \code{1.1}.
#'     \item \code{seed}: random-number seed. The default \code{NA} causes
#'       \code{MCMCpack::MCMCpoisson()} to use the Mersenne Twister generator
#'       with fixed seed 12345, so repeated calls with otherwise identical
#'       arguments are reproducible. Supply another seed to obtain a different
#'       MCMC run.
#'   }
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{call}}{
#'       The matched function call used to obtain the result.
#'     }
#'     \item{\code{popest}}{
#'       Posterior median estimate of the total population size.
#'     }
#'     \item{\code{quantiles}}{
#'       Posterior quantiles of the total population size.
#'     }
#'     \item{\code{retained_interactions}}{
#'       Pairwise interaction effects retained by the thresholding procedure
#'       and estimated by MCMC.
#'     }
#'     \item{\code{threshold_statistics}}{
#'       Absolute posterior mean to posterior standard deviation ratios for
#'       the pairwise interaction effects considered in the initial
#'       thresholding step.
#'     }
#'     \item{\code{eligible_triples}}{
#'       If \code{maxorder = 3}, the three-way interactions eligible for
#'       consideration because all three constituent pairwise interactions
#'       have been retained.
#'     }
#'     \item{\code{retained_triples}}{
#'       If the third-order thresholding step is carried out, the eligible
#'       three-way interactions retained by that thresholding step.
#'     }
#'     \item{\code{triple_threshold_statistics}}{
#'       If the third-order thresholding step is carried out, the threshold
#'       statistics for the eligible three-way interactions.
#'     }
#'     \item{\code{minus_infinite_estimated_effects}}{
#'       With an improper prior, interaction effects whose posterior
#'       distribution is concentrated at minus infinity. This component is
#'       present only if such effects occur. These effects are retained in
#'       the fitted model but are reported separately from
#'       \code{retained_interactions}.
#'     }
#'     \item{\code{posterior}}{
#'       If \code{return_posterior = TRUE}, the full posterior sample of the
#'       total population size.
#'     }
#'   }
#'
#' @references
#' Silverman, B. W. (2020). Multiple-systems analysis for the quantification of
#' modern slavery: classical and Bayesian approaches. \emph{Journal of the
#' Royal Statistical Society: Series A (Statistics in Society)}, 183, 691--736.
#'
#' Fienberg, S. E. and Rinaldo, A. (2012). Maximum likelihood estimation in
#' log-linear models. \emph{The Annals of Statistics}, 40, 996--1023.
#'
#' Martin, A. D., Quinn, K. M. and Park, J. H. (2011). MCMCpack: Markov Chain
#' Monte Carlo in R. \emph{Journal of Statistical Software}, 42(9), 1--21.
#'
#' @examples
#' data(Western, package = "MultipleSystemsEstimation")
#'
#' fit <- estimate_population_bayesthresh(
#'     Western,
#'     burnin = 100,
#'     mcmc = 1000,
#'     seed = 1234
#' )
#' fit$popest
#' fit$retained_interactions
#'
#' @export
estimate_population_bayesthresh <- function(
    zdat,
    prior = "proper",
    prior_variance = 1,
    main_variance = 1e4,
    threshold = 2,
    maxorder = 2,
    return_posterior = FALSE,
    ...
) {
  call <- match.call()
  prior <- match.arg(prior, c("proper", "improper"))

    if (!maxorder %in% c(2, 3))
        stop("maxorder must be either 2 or 3.")

    zfull <- MultipleSystemsEstimation::tidy_lists(
        zdat,
        includezerocounts = TRUE
    )

    nlists <- ncol(zfull) - 1

    if (maxorder == 3 && nlists < 3)
        stop("Three-way interactions require at least three lists.")

    .bayesthresh_check_pair_start(zfull)

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

    triple_candidates <- character(0)
    triple_threshold <- NULL
    triples <- character(0)
    removed_triples <- character(0)

    # For maxorder = 3, try the three-way extension before doing the
    # retained-pairwise refit.  This avoids an unnecessary MCMC fit when
    # the three-way extension is successfully carried through.
    if (maxorder == 3) {
        triple_candidates <- .bayesthresh_eligible_triples(
            pairs,
            nlists
        )
    }

    triple_start_ok <- FALSE

    if (maxorder == 3 && length(triple_candidates)) {
        triple_start_ok <- .bayesthresh_check_triple_start(
            zfull,
            c(pairs, triple_candidates)
        )

        if (!triple_start_ok) {
            warning(
                paste(
                    "The maximal three-way extension fails the",
                    "Fienberg-Rinaldo existence criterion;",
                    "reverting to maxorder = 2."
                ),
                call. = FALSE
            )
        }
    }

    if (maxorder == 3 && length(triple_candidates) && triple_start_ok) {
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
    } else {
        # Final order-2 fit: used for maxorder = 2, when there are no
        # eligible triples, or when the proposed three-way extension
        # fails the Fienberg-Rinaldo check.
        final_fit <- .bayesthresh_fit(
            zrefit,
            pairs,
            prior,
            prior_variance,
            main_variance,
            prune_improper = FALSE,
            ...
        )
    }

    population <- .bayesthresh_population(
        final_fit$fit,
        zfull
    )

    removed_pairs <- fit1$removed[
        .bayesthresh_order(fit1$removed) == 2
    ]

    ans <- list(
      call = call,
        popest = unname(population$quantiles["50%"]),
        quantiles = population$quantiles,

        retained_interactions =
            .bayesthresh_pretty(pairs, zfull)
    )

    if (maxorder == 3) {
        ans$eligible_triples <-
            .bayesthresh_pretty(
                triple_candidates,
                zfull
            )
    }

    if (!is.null(triple_threshold)) {
        ans$retained_triples <-
            .bayesthresh_pretty(triples, zfull)
    }

    ans$threshold_statistics <-
        .bayesthresh_pretty_named(
            pair_threshold$ratios,
            zfull
        )

    if (!is.null(triple_threshold)) {
        ans$triple_threshold_statistics <-
            .bayesthresh_pretty_named(
                triple_threshold$ratios,
                zfull
            )
    }

    if (prior == "improper") {
        minus_infinite_effects <-
            unique(c(removed_pairs, removed_triples))

        if (length(minus_infinite_effects)) {
            ans$minus_infinite_estimated_effects <-
                .bayesthresh_pretty(
                    minus_infinite_effects,
                    zfull
                )
        }
    }

    if (return_posterior)
        ans$posterior <- population$total_population

    ans
}
