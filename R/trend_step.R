##' Step change trend response model
##'
##' @param t numeric; vector of time points. 
##' @param change_points numeric; vector of change points, within `t`
##' @param means numeric; vector of means for the regimes implied by
##'   `change_points`. Must be of length `length(change_points) + 1`.
##' @param ... other arguments. Ignored here.
##' 
##' @importFrom tibble tibble
##' @importFrom stats approx
##'
##' @export
##' 
##' @examples
##' \dontshow{
##' set.seed(1)
##' op <- options(digits = 3, cli.unicode = FALSE)
##' }
##' sims <- step_trend(1:100, change_points = c(25, 75), means = c(2, 8, 4))
##' sims
##'
##' library("ggplot2")
##' ggplot(sims, aes(x = time, y = trend)) +
##'   geom_step()
##' \dontshow{options(op)}
`step_trend` <- function(t, change_points, means, ...) {
    ## is t in order?
    if (is.unsorted(t)) {
        stop("'t' must be in increasing order.")
    }
    
    nt <- length(t) # length of series
    n_pts <- length(change_points)
    n_mns <- length(means)
    if (!identical(n_pts, n_mns - 1L)) {
        stop("Number of change points must be one fewer than number of means.")
    }

    ## set trend vector to last mean
    trend <- rep(means[n_mns], nt)
    ## loop over the change points, repeat mean[i] for t < change_point[i]
    start <- t[1]
    for (i in seq_along(change_points)) {
        ind <-  (t >= start) & (t < change_points[i])
        trend[ind] <- means[i]
        start <- change_points[i]
    }
    
    ## linear sequence from start to end of the length of t
    ## trend <- seq(start_value, end_value, length.out = nt)

    ## if t is irregular, interpolate truth to the irregular t points
    irregular <- length(unique(diff(t))) > 1L
    if (irregular) {
        ## use approx to interpolate
        trend <- approx(x = seq(t[1], t[nt], length.out = nt),
                        y = trend, xout = t)$y
    }

    ## arrange in a tibble
    out <- tibble(t = t, trend = trend)
    class(out) <- c("step_trend", "abrupt_driver", class(out))
    out
}

##' Simulate data from a linear trend model
##'
##' @param sampling_distribution function; a random number generating function,
##'   which takes as it's first argument the number of observations to sample.
##'   The second argument should be the expected value. The default, if nothing
##'   is supplied, is [stats::rnorm()].
##' @param seed numeric; a seed for the simulation.
##' @param ... additional arguments that will be passed to
##'   `sampling_distribution`.
##'
##' @inheritParams step_trend
##' 
##' @importFrom stats approx
##' @importFrom tibble add_column
##'
##' @export
##'
##' @examples
##' \dontshow{
##' set.seed(1)
##' op <- options(digits = 3, cli.unicode = FALSE)
##' }
##' sims <- simulate_step_trend(1:100, change_points = c(25, 75),
##'                             means = c(2, 8, 4))
##' sims
##'
##' library("ggplot2")
##' ggplot(sims, aes(x = time, y = y)) +
##'   geom_point() +
##'   geom_step(aes(y = trend))
##' \dontshow{options(op)}
`simulate_step_trend` <- function(t, change_points, means,
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
    out <- step_trend(t = t, change_points = change_points, means = means)

    ## generate noisy values from trend
    out <- add_column(out, y = fun(nt, out$trend))
    class(out) <- c("simulate_step_trend", "simulate_driver",
                    "step_trend", "abrupt_driver", class(out))
    attr(out, "rng_state") <- rng_state
    out
}
