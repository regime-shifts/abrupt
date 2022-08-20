#' Generate a deterministic skeleton for a abrupt change driver
#' 
#' @param t numeric; either a length 1 vector indicating the number ot time
#'   points in the skeleton, or a vector of time points in the interval [0,1].
#' @param scenario
#' @param breaks
#' @param mean
#' @param ...
#'
#' @export
#'
#' @examples
#'
#' skeleton(t = 10, scenario = "constant", mean = 0.5)
#'
#' skeleton(t = 10, scenario = "stepped", breaks = c(0.3, 0.75),
#'          mean = c(2, 4, 6))
`skeleton` <- function(t = 100,
    scenario = c("constant", "linear", "segmented", "stepped", "unimodal",
                 "cyclic", "auto_regressive"),
    breaks = 0.5,
    mean = 0.5,
    ...
) {
    # handle scalar t
    if (length(t) == 1L) {
        if (t < 3L || isFALSE(t %% 1 == 0)) {
            stop("If scalar 't',  must be an integer > 2")
        }
        t <- seq(0, 1, length.out = t)
    }

    # handle scenario
    scenario <- match.arg(scenario)
    fun <-  paste0(scenario, "_skeleton")
    fun <- match.fun(fun)

    # call the actual skeleton function
    skeleton <- fun(t, breaks = breaks, mean = mean, ...)

    skeleton
}

#' @importFrom tibble tibble
#' @importFrom dplyr arrange
`constant_skeleton` <- function(t, mean, ...) {
    if (length(mean) > 1L) {
        stop("Single mean only")
    }
    skeleton <- tibble(t = t, value = rep_len(mean, length(t))) |>
        arrange(t)
    class(skeleton) <- append("skeleton", class(skeleton))
    skeleton
}

#' @importFrom tibble tibble
#' @importFrom dplyr arrange
`stepped_skeleton` <- function(t, breaks, mean, ...) {
    n_breaks <- length(breaks)
    n_means <- length(mean)
    if (isFALSE(n_breaks == (n_means - 1))) {
        stop("Number of breaks must be 1 fewer than the number of means")
    }
    cuts <- cut(t, c(0L, breaks, 1L), include.lowest = TRUE, right = FALSE)
    skeleton <- tibble(t = t, value = mean[cuts]) |>
        arrange(t)
    class(skeleton) <- append("skeleton", class(skeleton))
    skeleton
}
