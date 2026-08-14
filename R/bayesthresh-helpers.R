# Internal helpers for Bayesian threshold estimation

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

.bayesthresh_admissible_triples <- function(
    retained_effects,
    nlists
) {

  if (nlists < 3) {
    return(numeric(0))
  }

  triples <- combn(seq_len(nlists), 3)

  admissible <- apply(
    triples,
    2,
    function(ijk) {

      required_pairs <- combn(ijk, 2)

      pair_codes <- apply(
        required_pairs,
        2,
        function(ij) {
          z <- integer(nlists)
          z[ij] <- 1L
          encode_capture(z)
        }
      )

      all(pair_codes %in% retained_effects)
    }
  )

  if (!any(admissible)) {
    return(numeric(0))
  }

  apply(
    triples[, admissible, drop = FALSE],
    2,
    function(ijk) {
      z <- integer(nlists)
      z[ijk] <- 1L
      encode_capture(z)
    }
  )
}

.bayesthresh_third_fit_model <- function(
    pair_terms,
    admissible_triples,
    nlists
) {

  convert_to_hierarchy(
    c(
      pair_terms$main,
      pair_terms$interactions,
      admissible_triples
    ),
    nlists = nlists
  )
}

.bayesthresh_threshold_triples <- function(
    third_fit,
    third_design,
    third_terms,
    threshold = 2,
    nlists
) {

  stats <- summary(third_fit)$statistics

  interaction_effects <- third_terms$interactions

  interaction_patterns <- lapply(
    interaction_effects,
    decode_capture,
    nlists = nlists
  )

  interaction_orders <- vapply(
    interaction_patterns,
    sum,
    numeric(1)
  )

  triple_effects <- interaction_effects[
    interaction_orders == 3
  ]

  triple_names <- vapply(
    triple_effects,
    function(effect) {
      pattern <- decode_capture(
        effect,
        nlists = nlists
      )

      paste(
        colnames(third_design)[
          1 + which(pattern)
        ],
        collapse = ":"
      )
    },
    character(1)
  )

  term_names <- colnames(third_design)[-1]

  triple_rows <- match(
    triple_names,
    term_names
  ) + 1

  triple_stats <- stats[
    triple_rows,
    ,
    drop = FALSE
  ]

  ratios <- abs(
    triple_stats[, "Mean"] /
      triple_stats[, "SD"]
  )

  names(ratios) <- triple_names

  keep <- ratios >= threshold

  list(
    ratios = ratios,
    retained = names(ratios)[keep],
    dropped = names(ratios)[!keep],
    retained_effects = triple_effects[keep],
    dropped_effects = triple_effects[!keep]
  )
}

.bayesthresh_final_model <- function(
    reduced_terms,
    triple_threshold_result,
    nlists
) {

  convert_to_hierarchy(
    c(
      reduced_terms$main,
      reduced_terms$interactions,
      triple_threshold_result$retained_effects
    ),
    nlists = nlists
  )
}
