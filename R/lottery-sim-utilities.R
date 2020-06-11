library(dplyr)
library(ggplot2)
library(reshape2)
library(styler)


# All code for simulating the lottery model (refer to lottery_description for details)

############################################################
###### utility functions ###################################

# construct temporal variation in a parameter (one  change so far)
AbruptParam <- function(base_param, rel_change, rel_time, time) {
  result <- rep(base_param, time)
  ac <- floor(rel_time * time) # timing of the change
  result[ac:time] <- result[ac:time] * rel_change
  return(result)
}

# computes the density-dependent survival
SurvivalDD <- function(freq, u, s) {
  return(u / (u + (1 - u) * exp(-s * freq)))
}

# generates random recruitment for each grid cell
FecundityStoch <- function(size, param, st) {
  sg <- ifelse(param$rho == 0, param$sb, sqrt(st - (1 - param$rho)^2 * param$sb^2) / param$rho) # normalize noise
  temp <- exp(log(param$b) + param$rho * rnorm(1, 0, sg) + (1 - param$rho) * rnorm(size * size, 0, param$sb))
  return(array(temp, c(size, size)))
}

# Construct a dispersal kernel so that a proportion fdisp is distributed in the 8-cell neighbors
KernelConstruct <- function(size, fdisp) {
  kernel <- matrix(0, nrow = size, ncol = size)
  kernel[1, 1] <- 1 - fdisp
  kernel[c(size, 2), c(size, 1, 2)] <- fdisp / 8
  kernel[1, c(size, 2)] <- fdisp / 8
  return(kernel)
}

# Distributes the propagules given a kernel (using convolution product and fourier transform)
DispersalSeed <- function(kernel, BB, size) {
  fkernel <- fft(kernel)
  fBB <- fft(BB, inverse = TRUE) / (size * size)
  result <- round(Re(fft(fkernel * fBB)), 10) # numerical errors round to the 10th decimal
  return(result)
}