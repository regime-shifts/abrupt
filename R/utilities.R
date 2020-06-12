##' Seed the RNG state possibly with a user supplied seed
##'
##' @param seed numeric seed for RNG. See [base::set.seed()].
##'
##' @return A list containing `seed`, the supplied seed, `initial_state`,
##'   the initial state of the RNG before the seed was set, and
##'   `kind`, the type of RNG used, as returned by [base::RNGkind()].
##'
##' @importFrom stats runif
`seed_rng` <- function(seed = NULL) {
    ## initialise seed if not set in session
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        runif(1)
    }
    rnd_kind <- as.list(RNGkind()) # need kind to be reproducible

    ## grab the current state of the RNG system
    ## want to return this so it can be reset on exit from
    ## calling function
    initial_state <- get(".Random.seed", envir = .GlobalEnv)

    ## if user provided a seed, set the seed now
    if (!is.null(seed)) {
        set.seed(seed)
    }

    ## return the seed and other info as a list
    list(seed = seed, initial_state = initial_state, kind = rnd_kind)
}
