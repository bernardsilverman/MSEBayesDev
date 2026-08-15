# Internal helpers for Bayesian-threshold multiple-systems estimation


.bayesthresh_all_effects <- function(nlists, effect_order) {
    if (nlists < effect_order)
        return(character(0))

    if (effect_order == 1)
        return(paste0("x", seq_len(nlists)))

    ans <- character(0)

    # Preserve the coefficient ordering used by the existing implementation.
    # For pairs with five lists:
    # 12, 13, 23, 14, 24, 34, 15, 25, 35, 45.
    for (last in effect_order:nlists) {
        earlier <- utils::combn(
            seq_len(last - 1),
            effect_order - 1
        )

        new <- apply(
            earlier,
            2,
            function(ii)
                paste0("x", c(ii, last), collapse = ":")
        )

        ans <- c(ans, new)
    }

    ans
}


.bayesthresh_indices <- function(effect) {
    as.integer(
        sub(
            "^x",
            "",
            strsplit(effect, ":x", fixed = TRUE)[[1]]
        )
    )
}


.bayesthresh_order <- function(effect) {
    lengths(strsplit(effect, ":", fixed = TRUE))
}


.bayesthresh_pretty <- function(effects, zfull) {
    if (!length(effects))
        return(character(0))

    list_names <- colnames(zfull)[seq_len(ncol(zfull) - 1)]

    vapply(
        effects,
        function(effect)
            paste(
                list_names[.bayesthresh_indices(effect)],
                collapse = ":"
            ),
        character(1),
        USE.NAMES = FALSE
    )
}


.bayesthresh_pretty_named <- function(x, zfull) {
    if (!length(x))
        return(x)

    names(x) <- .bayesthresh_pretty(names(x), zfull)
    x
}


.bayesthresh_sufficient_stat <- function(zfull, effect) {
    ii <- .bayesthresh_indices(effect)

    active <- apply(
        as.matrix(zfull[, ii, drop = FALSE]),
        1,
        prod
    )

    sum(zfull[[ncol(zfull)]] * active)
}


.bayesthresh_remove_zero_effects <- function(zfull, effects) {
    if (!length(effects)) {
        return(list(
            data = zfull,
            effects = effects,
            removed = character(0)
        ))
    }

    zero <- vapply(
        effects,
        function(effect)
            .bayesthresh_sufficient_stat(zfull, effect) == 0,
        logical(1)
    )

    removed <- effects[zero]

    if (!length(removed)) {
        return(list(
            data = zfull,
            effects = effects,
            removed = character(0)
        ))
    }

    keep <- rep(TRUE, nrow(zfull))

    for (effect in removed) {
        ii <- .bayesthresh_indices(effect)

        contains <- apply(
            as.matrix(zfull[, ii, drop = FALSE]),
            1,
            prod
        ) == 1

        keep[contains] <- FALSE
    }

    retained <- effects[!zero]

    # Preserve hierarchy: if a lower-order effect is removed, any
    # higher-order effect containing it must also be removed.
    if (length(retained)) {
        drop_higher <- vapply(
            retained,
            function(effect) {
                jj <- .bayesthresh_indices(effect)

                any(vapply(
                    removed,
                    function(rem) {
                        ii <- .bayesthresh_indices(rem)
                        all(ii %in% jj)
                    },
                    logical(1)
                ))
            },
            logical(1)
        )

        removed <- unique(c(removed, retained[drop_higher]))
        retained <- retained[!drop_higher]
    }

    list(
        data = zfull[keep, , drop = FALSE],
        effects = retained,
        removed = removed
    )
}


.bayesthresh_B0 <- function(
    nlists,
    neffects,
    prior_variance,
    main_variance
) {
    diag(c(
        rep(1 / main_variance, 1 + nlists),
        rep(1 / prior_variance, neffects)
    ))
}


.bayesthresh_start <- function(zfull, neffects) {
    nlists <- ncol(zfull) - 1
    counts <- zfull[[ncol(zfull)]]
    nobs <- sum(counts)

    main_counts <- as.numeric(
        t(as.matrix(
            zfull[, seq_len(nlists), drop = FALSE]
        )) %*% counts
    )

    c(
        log(nobs * 5),
        log(main_counts / (nobs * 5)),
        rep(0, neffects)
    )
}


.bayesthresh_design <- function(zfull, effects) {
    nlists <- ncol(zfull) - 1

    mains <- as.matrix(
        zfull[, seq_len(nlists), drop = FALSE]
    )

    if (!length(effects))
        return(mains)

    interactions <- vapply(
        effects,
        function(effect) {
            ii <- .bayesthresh_indices(effect)

            apply(
                as.matrix(zfull[, ii, drop = FALSE]),
                1,
                prod
            )
        },
        numeric(nrow(zfull))
    )

    cbind(mains, as.matrix(interactions))
}


.bayesthresh_fit <- function(
    zfull,
    effects,
    prior,
    prior_variance,
    main_variance,
    prune_improper = TRUE,
    ...
) {
    nlists <- ncol(zfull) - 1
    zfit <- zfull
    removed <- character(0)

    if (prior == "improper") {
        B0 <- 0

        if (prune_improper) {
            setup <- .bayesthresh_remove_zero_effects(
                zfit,
                effects
            )
            zfit <- setup$data
            effects <- setup$effects
            removed <- setup$removed
        }
    } else {
        B0 <- .bayesthresh_B0(
            nlists,
            length(effects),
            prior_variance,
            main_variance
        )
    }

    # Build the numerical design explicitly and give MCMCpoisson an
    # additive formula.  This preserves fixed-seed reproducibility with
    # the previous implementation.
    design <- .bayesthresh_design(zfit, effects)

    dat <- as.data.frame(design)
    names(dat) <- paste0("x", seq_len(ncol(design)))
    dat$y <- zfit[[ncol(zfit)]]

    form <- stats::reformulate(
        paste0("x", seq_len(ncol(design))),
        response = "y"
    )

    beta_start <- .bayesthresh_start(
        zfit,
        length(effects)
    )

    fit <- withCallingHandlers(
        MCMCpack::MCMCpoisson(
            formula = form,
            data = dat,
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

    colnames(fit) <- c(
        "(Intercept)",
        paste0("x", seq_len(nlists)),
        effects
    )

    list(
        fit = fit,
        effects = effects,
        data = zfit,
        design = design,
        removed = removed
    )
}


.bayesthresh_threshold <- function(
    fit,
    effects,
    effect_order,
    threshold
) {
    candidates <- effects[
        .bayesthresh_order(effects) == effect_order
    ]

    if (!length(candidates)) {
        return(list(
            ratios = numeric(0),
            retained = character(0),
            dropped = character(0)
        ))
    }

    stats <- summary(fit)$statistics
    rows <- match(candidates, rownames(stats))
    present <- !is.na(rows)

    candidates <- candidates[present]
    rows <- rows[present]

    ratios <- abs(
        stats[rows, "Mean"] /
            stats[rows, "SD"]
    )
    names(ratios) <- candidates

    keep <- ratios >= threshold

    list(
        ratios = ratios,
        retained = candidates[keep],
        dropped = candidates[!keep]
    )
}


.bayesthresh_admissible_triples <- function(pairs, nlists) {
    triples <- .bayesthresh_all_effects(
        nlists,
        3
    )

    if (!length(triples))
        return(character(0))

    triples[vapply(
        triples,
        function(triple) {
            ii <- .bayesthresh_indices(triple)

            required <- apply(
                utils::combn(ii, 2),
                2,
                function(jj)
                    paste0("x", jj, collapse = ":")
            )

            all(required %in% pairs)
        },
        logical(1)
    )]
}


.bayesthresh_population <- function(
    fit,
    zfull,
    probs = c(0.025, 0.1, 0.5, 0.9, 0.975)
) {
    nobserved <- sum(zfull[[ncol(zfull)]])
    dark_figure <- exp(fit[, "(Intercept)"])
    total_population <- nobserved + dark_figure

    list(
        nobserved = nobserved,
        dark_figure = dark_figure,
        total_population = total_population,
        quantiles = stats::quantile(
            total_population,
            probs = probs
        )
    )
}
