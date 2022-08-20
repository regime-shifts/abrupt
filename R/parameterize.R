`parameterize` <- function(t = 100,
                           model = c("mean_plus_noise", "lottery", "ricker",
                               "trophic_triangle"),
                           scenario = c("constant", "linear", "segmented",
                               "stepped", "unimodal",
                               "cyclic", "auto_regressive"),
                           observation_dist = NULL,
                           parameters = NULL,
                           n_spp = 1,
                           ...) {

    model <- match.arg(model)
    scenario <- match.arg(scenario)

    par_fun <- paste("parameterize", model, sep = "_")
    par_fun <- match.fun(par_fun)

    pars <- par_fun(t = t, scenario = scenario, n_spp = n_spp, ...)
}

# vary noise, survival (u), density dependence (s), birth rate sd (sd)
`parameterize_lottery` <- function(t,
                                   scenario = "constant",
                                   n_spp = 2,
                                   breaks = NULL,
                                   survival = 0.95,
                                   density_dependence = 0.01,
                                   birth_sd = 0.8, ...) {
    birth <- 20 # this is parameter b, hard coded
    # birth needs to be in data frame of parameters

    # survival can be a list of length n_spp, each element is the mean for the
    # skeleton, each species follows the same skeleton/scenario

    # genearate data frame of parameters with length(t) rows

    # scenario affects survival *only*
    # survival gets passed to the mean argument of the skeleton fun

    # parameters is a data frame with columns:
    # t, spp, birth, survival, density_dependence, birth_sd

    # after generated parameters do parameter checks
}

# constant - mean
# linear - start and end values, linearlly interpolate
# segmented - specify internal breaks, mean = breaks = 2 (start and end) values
#           - linearlly interpolate between the start, breaks, end
# stepped - breaks and mean (length(mean) == length(breaks) + 1)
# unimdal - optima, tolerance assume a gaussian curve
# cyclic  - mean, amplitude, freq: output mean + (amplitude * sin(freq*t*2*pi)),
#         - min freq is 2
# auto_regressive mean rho t_1 sigma

# t <- 1:100
# t_1 <- 0
# mean <- 10
# rho <- 0.99
# sd <- 2
# x <- rnorm(length(t), mean = 0, sd = sd)
# x[1] <- t_1 - mean # (init)
# f <- stats::filter(x, filter = rho, method = "recursive", sides = 1)
# f + mean
# plot(f + mean)
