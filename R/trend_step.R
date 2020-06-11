##' Step change trend response model
##'
##' @param t numeric; vector of time points. 
##' @param start_value,end_value numeric vectors of length 1; the start and end
##'   values for the linear trend.
##' @param ... other arguments. Ignored here.
##'
##' @importFrom tibble tibble
##' @importFrom stats approx
##'
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
    tibble(time = t, trend = trend)
}
