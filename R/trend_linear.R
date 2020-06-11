##' Linear trend response model
##'
##' @param t numeric; vector of time points. 
##' @param start_value,end_value numeric vectors of length 1; the start and end
##'   values for the linear trend.
##' @param ... other arguments. Ignored here.
##'
##' @importFrom tibble tibble
##' @importFrom stats approx
##'
`linear_trend` <- function(t, start_value = 0, end_value = 1, ...) {
    ## is t in order?
    if (is.unsorted(t)) {
        stop("'t' must be in increasing order.")
    }
    
    nt <- length(t) # length of series

    ## linear sequence from start to end of the length of t
    trend <- seq(start_value, end_value, length.out = nt)

    ## if t is irregular, interpolate truth to the irregular t points
    irregular <- length(unique(diff(t))) > 1L
    if (irregular) {
        ## use approx to interpolate
        trend <- approx(x = seq(t[1], t[nt], length.out = nt),
                        y = trend, xout = t)$y
    }

    ## arrange in a tibble
    tibble(time = t, trend = trend)
}

##' Simulate data from a linear trend model
##'
##' @param t numeric; vector of time points. 
##' @param start_value,end_value numeric vectors of length 1; the start and end
##'   values for the linear trend.
##' @param sampling_distribution function; a random number generating function,
##'   which takes as it's first argument the number of observations to sample.
##'   The second argument should be the expected value. The default, if nothing
##'   is supplied, is [stats::rnorm()].
##' @param seed numeric; a seed for the simulation.
##' @param ... additional arguments that will be passed to
##'   `sampling_distribution`.
##' 
##' @importFrom stats approx
##' @importFrom tibble add_column
`simulate_linear_trend` <- function(t, start_value = 0, end_value = 1,
                                    sampling_distribution = NULL, seed = NULL,
                                    ...) {
    ## initialise the RNG, possibly with the user-supplied seed
    rng_state <- seed_rng(seed = seed)
    ## arrange for RNG state to be reset upon exit from function
    on.exit(assign(".Random.seed", rng_state$initial_state, envir = .GlobalEnv))

    ## match the sampling_distribution to a function
    fun <- if (is.null(sampling_distribution)) {
        stats::rnorm # use rnorm() for the default
    } else {
        match.fun(sampling_distribution)
    }
    
    nt <- length(t) # length of series

    ## generate linear trend
    out <- linear_trend(t = t, start_value = start_value, end_value = end_value)

    ## generate noisy values from trend
    out <- add_column(out, y = fun(nt, out$trend))
    attr(out, "rng_state") <- rng_state
    out
}

`linear_trend_fun` <- function(t, start_value = 0, end_value = 1,
                               sampling_distribution = NULL, ...) {
    .NotYetImplemented()
}
