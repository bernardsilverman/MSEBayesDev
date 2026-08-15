# Deliberate zero-three-way test for estimate_population_bayesthresh_enhanced()
#
# Four lists, with positive counts in every two-way overlap but no observed
# case belonging simultaneously to A, B and C.  With threshold = 0 all
# estimable pair interactions are retained, so A:B:C is necessarily an
# admissible three-way candidate.  Under an improper prior it should then be
# recognised as a zero-sufficient-statistic triple and removed.

make_zero_triple_data <- function() {
    z <- expand.grid(
        A = 0:1,
        B = 0:1,
        C = 0:1,
        D = 0:1
    )

    z <- z[rowSums(z) > 0, , drop = FALSE]
    z$count <- 10L

    # Make the A:B:C sufficient statistic zero.
    z$count[z$A == 1 & z$B == 1 & z$C == 1] <- 0L

    # Supply only observed capture histories; tidy_lists() will complete
    # the table internally.
    z[z$count > 0, , drop = FALSE]
}


zero_triple_dat <- make_zero_triple_data()

zero_triple_fit <- estimate_population_bayesthresh_enhanced(
    zero_triple_dat,
    prior = "improper",
    maxorder = 3,
    threshold = 0,
    mcmc = 10000,
    burnin = 1000,
    thin = 10,
    seed = 1234
)

zero_triple_fit$removed_pairs
zero_triple_fit$admissible_triples
zero_triple_fit$removed_triples
zero_triple_fit$retained_triples

# Expected key result:
#   "A:B:C" appears in removed_triples,
#   and does not appear in retained_triples.
