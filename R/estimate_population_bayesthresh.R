#' Estimate population size by Bayesian thresholding
#'
#' Implements the Bayesian threshold method of Silverman (2020),
#' initially for pairwise interactions only.
#'
#' @param zdat Multiple-systems data in the usual
#'   MultipleSystemsEstimation format.
#' @param prior Either "proper" or "improper".
#' @param prior_variance Prior variance for pairwise interaction terms
#'   when a proper prior is used.
#' @param main_variance Prior variance for the intercept and main effects
#'   when a proper prior is used.
#' @param threshold Threshold applied to
#'   abs(posterior mean / posterior SD).
#' @param ... Additional arguments passed to MCMCpoisson().
#'
#' @export
estimate_population_bayesthresh <- function(
    zdat,
    prior = c("proper", "improper"),
    prior_variance = 1,
    main_variance = 1e4,
    threshold = 2,
    ...
) {

  prior <- match.arg(prior)

  # Full completed contingency table.
  # Keep this unchanged because it is needed for the observed population total.
  zfull <- .bayesthresh_prepare_data(zdat)

  nlists <- ncol(zfull) - 1

  # Construct the full pairwise model.
  model <- .bayesthresh_full_pairwise_model(zfull)

  terms <- .bayesthresh_model_terms(
    model,
    nlists = nlists
  )

  design <- .bayesthresh_design_matrix(
    zfull,
    terms
  )

  # By default the likelihood uses the full table.
  zfit <- zfull
  removed_pairs <- character(0)

  # For the improper prior, remove zero-overlap interactions
  # and all cells containing those pairs.
  if (prior == "improper") {

    improper_setup <- .bayesthresh_improper_setup(
      zfull,
      terms,
      design
    )

    zfit <- improper_setup$data
    model <- improper_setup$model
    terms <- improper_setup$terms
    design <- improper_setup$design
    removed_pairs <- improper_setup$removed_pairs
  }

  # First-stage MCMC fit.
  first_fit <- .bayesthresh_first_fit(
    zfit,
    design,
    terms,
    prior = prior,
    prior_variance = prior_variance,
    main_variance = main_variance,
    ...
  )

  # Threshold the pairwise interactions.
  threshold_result <- .bayesthresh_threshold_interactions(
    first_fit,
    design,
    terms,
    threshold = threshold
  )

  # Construct the thresholded hierarchical model.
  reduced_model <- .bayesthresh_reduced_model(
    terms,
    threshold_result,
    nlists
  )

  reduced_terms <- .bayesthresh_model_terms(
    reduced_model,
    nlists
  )

  # Important: for the improper case this must use zfit,
  # not the original full table.
  reduced_design <- .bayesthresh_design_matrix(
    zfit,
    reduced_terms
  )

  # Second-stage MCMC fit.
  second_fit <- .bayesthresh_first_fit(
    zfit,
    reduced_design,
    reduced_terms,
    prior = prior,
    prior_variance = prior_variance,
    main_variance = main_variance,
    ...
  )

  # Population total must use the ORIGINAL complete table,
  # not zfit, because all observed cases remain part of the population.
  population <- .bayesthresh_population_summary(
    second_fit,
    zfull
  )

  # Verbose return object while the method is under development.
  return(
    list(
      popest = unname(population$quantiles["50%"]),
      quantiles = population$quantiles,
      model = reduced_model,
      retained_interactions = threshold_result$retained,
      threshold_statistics = threshold_result$ratios,
      removed_pairs = removed_pairs,
      posterior = population$total_population
    )
  )
}

.bayesthresh_prepare_data <- function(zdat) {
  MultipleSystemsEstimation::tidy_lists(
    zdat,
    includezerocounts = TRUE
  )
}

.bayesthresh_pair_counts <- function(zfull) {

  nlists <- ncol(zfull) - 1
  list_names <- colnames(zfull)[seq_len(nlists)]
  counts <- zfull[[nlists + 1]]

  pairs <- combn(seq_len(nlists), 2)

  pair_counts <- apply(
    pairs,
    2,
    function(ij) {
      sum(
        counts *
          zfull[[ij[1]]] *
          zfull[[ij[2]]]
      )
    }
  )

  names(pair_counts) <- apply(
    pairs,
    2,
    function(ij) paste(list_names[ij], collapse = ":")
  )

  pair_counts
}

.bayesthresh_full_pairwise_model <- function(zfull) {

  nlists <- ncol(zfull) - 1

  main_effects <- vapply(
    seq_len(nlists),
    function(i) {
      z <- integer(nlists)
      z[i] <- 1L
      encode_capture(z)
    },
    numeric(1)
  )

  pairs <- combn(seq_len(nlists), 2)

  pair_effects <- apply(
    pairs,
    2,
    function(ij) {
      z <- integer(nlists)
      z[ij] <- 1L
      encode_capture(z)
    }
  )

  convert_to_hierarchy(
    c(main_effects, pair_effects),
    nlists = nlists
  )
}

.bayesthresh_model_terms <- function(model, nlists) {

  effects <- convert_from_hierarchy(model)

  decoded <- lapply(
    effects,
    decode_capture,
    nlists = nlists
  )

  orders <- vapply(decoded, sum, numeric(1))

  list(
    main = effects[orders == 1],
    interactions = effects[orders >= 2]
  )
}

.bayesthresh_design_matrix <- function(zfull, terms) {

  nlists <- ncol(zfull) - 1
  x <- as.matrix(zfull[, seq_len(nlists), drop = FALSE])

  effects <- c(terms$main, terms$interactions)

  effect_patterns <- lapply(
    effects,
    decode_capture,
    nlists = nlists
  )

  modmat <- vapply(
    effect_patterns,
    function(pattern) {
      rows_needed <- which(pattern)
      apply(x[, rows_needed, drop = FALSE], 1, prod)
    },
    numeric(nrow(x))
  )

  modmat <- cbind("(Intercept)" = 1, modmat)

  term_names <- vapply(
    effect_patterns,
    function(pattern) {
      paste(colnames(x)[which(pattern)], collapse = ":")
    },
    character(1)
  )

  colnames(modmat) <- c(
    "(Intercept)",
    term_names
  )

  modmat
}

.bayesthresh_mcmc_data <- function(zfull, design) {

  x <- design[, -1, drop = FALSE]

  colnames(x) <- paste0("x", seq_len(ncol(x)))

  out <- as.data.frame(x)
  out$y <- zfull[[ncol(zfull)]]

  out
}

.bayesthresh_prior_precision <- function(
    terms,
    prior_variance = 1,
    main_variance = 1e4
) {

  nmain <- length(terms$main)
  ninteractions <- length(terms$interactions)

  diag(c(
    rep(1 / main_variance, 1 + nmain),
    rep(1 / prior_variance, ninteractions)
  ))
}

.bayesthresh_start_values <- function(zfull, design) {

  nobs <- sum(zfull[[ncol(zfull)]])

  x <- design[, -1, drop = FALSE]
  nmain <- ncol(zfull) - 1

  # Intercept and main-effect starting values,
  # following the strategy used in the original implementation.
  main_counts <- as.numeric(
    t(as.matrix(zfull[, seq_len(nmain), drop = FALSE])) %*%
      zfull[[ncol(zfull)]]
  )

  start <- c(
    log(nobs * 5),
    log(main_counts / (nobs * 5))
  )

  # Start all interaction effects at zero.
  ninteractions <- ncol(x) - nmain

  c(
    start,
    rep(0, ninteractions)
  )
}

.bayesthresh_first_fit <- function(
    zfull,
    design,
    terms,
    prior = "proper",
    prior_variance = 1,
    main_variance = 1e4,
    ...
) {

  mcmc_data <- .bayesthresh_mcmc_data(
    zfull,
    design
  )

  predictor_names <- paste0(
    "x",
    seq_len(ncol(design) - 1)
  )

  form <- reformulate(
    predictor_names,
    response = "y"
  )

  if (prior == "improper") {

    B0 <- 0

  } else {

    B0 <- .bayesthresh_prior_precision(
      terms,
      prior_variance = prior_variance,
      main_variance = main_variance
    )
  }

  beta_start <- .bayesthresh_start_values(
    zfull,
    design
  )

  fit <- withCallingHandlers(
    MCMCpack::MCMCpoisson(
      formula = form,
      data = mcmc_data,
      B0 = B0,
      beta.start = beta_start,
      ...
    ),
    warning = function(w) {
      if (grepl(
        "fitted rates numerically 0 occurred",
        conditionMessage(w),
        fixed = TRUE
      )) {
        invokeRestart("muffleWarning")
      }
    }
  )

  fit
}

.bayesthresh_threshold_interactions <- function(
    first_fit,
    design,
    terms,
    threshold = 2
) {

  stats <- summary(first_fit)$statistics

  nmain <- length(terms$main)

  interaction_rows <- (nmain + 2):nrow(stats)

  interaction_stats <- stats[
    interaction_rows,
    ,
    drop = FALSE
  ]

  ratios <- abs(
    interaction_stats[, "Mean"] /
      interaction_stats[, "SD"]
  )

  interaction_names <- colnames(design)[interaction_rows]

  names(ratios) <- interaction_names

  keep <- ratios >= threshold

  list(
    ratios = ratios,
    retained = names(ratios)[keep],
    dropped = names(ratios)[!keep],
    retained_effects = terms$interactions[keep],
    dropped_effects = terms$interactions[!keep]
  )
}

.bayesthresh_reduced_model <- function(
    terms,
    threshold_result,
    nlists
) {

  convert_to_hierarchy(
    c(
      terms$main,
      threshold_result$retained_effects
    ),
    nlists = nlists
  )
}

.bayesthresh_population_summary <- function(
    second_fit,
    zfull,
    probs = c(0.025, 0.1, 0.5, 0.9, 0.975)
) {

  nobserved <- sum(zfull[[ncol(zfull)]])

  mu_draws <- second_fit[, "(Intercept)"]

  dark_figure <- exp(mu_draws)
  total_population <- nobserved + dark_figure

  list(
    nobserved = nobserved,
    dark_figure = dark_figure,
    total_population = total_population,
    quantiles = quantile(
      total_population,
      probs = probs
    )
  )
}

.bayesthresh_remove_empty_overlaps <- function(zfull) {

  nlists <- ncol(zfull) - 1
  list_names <- colnames(zfull)[seq_len(nlists)]

  pair_counts <- .bayesthresh_pair_counts(zfull)
  empty_pairs <- which(pair_counts == 0)

  if (length(empty_pairs) == 0) {
    return(
      list(
        data = zfull,
        removed_pairs = character(0),
        removed_pair_indices = matrix(integer(0), nrow = 2)
      )
    )
  }

  pairs <- combn(seq_len(nlists), 2)
  empty_pair_indices <- pairs[, empty_pairs, drop = FALSE]

  keep <- rep(TRUE, nrow(zfull))

  for (k in seq_len(ncol(empty_pair_indices))) {
    ij <- empty_pair_indices[, k]

    contains_pair <-
      zfull[[ij[1]]] == 1 &
      zfull[[ij[2]]] == 1

    keep[contains_pair] <- FALSE
  }

  removed_pairs <- apply(
    empty_pair_indices,
    2,
    function(ij) paste(list_names[ij], collapse = ":")
  )

  list(
    data = zfull[keep, , drop = FALSE],
    removed_pairs = removed_pairs,
    removed_pair_indices = empty_pair_indices
  )
}

.bayesthresh_improper_setup <- function(
    zfull,
    terms,
    design
) {

  empty <- .bayesthresh_remove_empty_overlaps(zfull)

  nlists <- ncol(zfull) - 1
  nmain <- length(terms$main)

  interaction_names <- colnames(design)[
    (nmain + 2):ncol(design)
  ]

  keep_interactions <-
    !interaction_names %in% empty$removed_pairs

  retained_effects <-
    terms$interactions[keep_interactions]

  model <- convert_to_hierarchy(
    c(
      terms$main,
      retained_effects
    ),
    nlists = nlists
  )

  new_terms <- .bayesthresh_model_terms(
    model,
    nlists
  )

  new_design <- .bayesthresh_design_matrix(
    empty$data,
    new_terms
  )

  list(
    data = empty$data,
    model = model,
    terms = new_terms,
    design = new_design,
    removed_pairs = empty$removed_pairs
  )
}

