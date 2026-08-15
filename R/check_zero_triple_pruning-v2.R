# Direct test of zero-three-way improper-prior handling
#
# This deliberately tests the pruning logic itself rather than forcing
# MCMCpoisson to fit a highly symmetric artificial improper-prior model.

make_zero_triple_data <- function() {
    z <- expand.grid(
        A = 0:1,
        B = 0:1,
        C = 0:1,
        D = 0:1
    )

    z <- z[rowSums(z) > 0, , drop = FALSE]

    # Unequal positive counts avoid accidental symmetries.
    z$count <- seq_len(nrow(z)) + 3L

    # No observed case is simultaneously in A, B and C.
    z$count[z$A == 1 & z$B == 1 & z$C == 1] <- 0L

    z
}


zero_triple_full <- make_zero_triple_data()

# Internal helper notation uses x1, x2, ...
pair_effects <- .bayesthresh2_all_effects(4, 2)
triple_effects <- .bayesthresh2_all_effects(4, 3)

pair_stats <- vapply(
    pair_effects,
    function(e)
        .bayesthresh2_sufficient_stat(zero_triple_full, e),
    numeric(1)
)

triple_stats <- vapply(
    triple_effects,
    function(e)
        .bayesthresh2_sufficient_stat(zero_triple_full, e),
    numeric(1)
)

names(pair_stats) <- .bayesthresh2_pretty(
    pair_effects,
    zero_triple_full
)

names(triple_stats) <- .bayesthresh2_pretty(
    triple_effects,
    zero_triple_full
)

cat("Pair sufficient statistics:\n")
print(pair_stats)

cat("\nTriple sufficient statistics:\n")
print(triple_stats)

stopifnot(all(pair_stats > 0))
stopifnot(triple_stats["A:B:C"] == 0)
stopifnot(all(triple_stats[names(triple_stats) != "A:B:C"] > 0))

setup <- .bayesthresh2_remove_zero_effects(
    zero_triple_full,
    c(pair_effects, triple_effects)
)

removed <- .bayesthresh2_pretty(
    setup$removed,
    zero_triple_full
)

cat("\nRemoved effects:\n")
print(removed)

stopifnot("A:B:C" %in% removed)
stopifnot(!any(c(
    "A:B", "A:C", "B:C"
) %in% removed))

# Every row retained after pruning must fail to contain A:B:C.
stopifnot(
    !any(
        setup$data$A == 1 &
        setup$data$B == 1 &
        setup$data$C == 1
    )
)

cat("\nZero-triple improper-prior pruning test: PASS\n")
