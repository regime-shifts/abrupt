#' Estimate change points using a regression tree
#'
`abr_fit.rpart` <- function(data, ...) {
    # fail early if rpart and partykit are not avilable
    pkgs <- c("rpart", "partykit")
    load_required_packages(pkgs)
    fit <- rpart::rpart(y ~ t, data = data)
    fit <- partykit::as.party(fit)
    fit
}

#' @importFrom purrr map_lgl
#' @importFrom glue glue_collapse
`load_required_packages` <- function(packages, quietly = TRUE, ...) {
    loaded <- map_lgl(packages, load_required_package, quietly = quietly, ...)
    # report which failed
    if (any(!loaded)) {
        pkgs <- packages[!loaded]
        msg <- glue::glue_collapse(glue::glue("{pkgs}"), sep = ", ")
        msg <- glue::glue("Failed to load packages: {tmp}.

Install the missing packages and try again.")
        stop(msg, call. = FALSE)
    }
}

`load_required_package` <- function(package, quietly = TRUE, ...) {
    requireNamespace(package, quietly = quietly, ...)
}