############################################################
##### main code to simulate the lottery model ##############
LotterySim <- function(nsp = nsp, time = time, size = size, param = param, driver = driver, seed = 1) {
  set.seed(seed)
  
  # contruct time-varying parameter
  UU <- AbruptParam(param$u, driver$ru, driver$rt, time)
  SS <- AbruptParam(param$s, driver$rs, driver$rt, time)
  ST <- AbruptParam(param$sT, driver$rsT, driver$rt, time)
  
  # initial conditions
  init <- rep(1 / nsp, nsp) # assumes same frequency
  
  # stores results
  res <- array(NA, c(nsp, size, size, time))
  for (x in 1:size) {
    for (y in 1:size) {
      res[, x, y, 1] <- init
    }
  }
  
  # generates recruitment for each grid cell and each time step (follows a lognormal distribution)
  BB <- array(NA, c(nsp, size, size, time))
  for (i in 1:nsp) {
    for (t in 1:time) {
      BB[i, , , t] <- FecundityStoch(size, param, ST[t])
    }
  }
  
  # dispersal kernel (only adjacent cells)
  kernel <- KernelConstruct(size, param$d)
  
  # per species
  occ_site <- array(NA, c(nsp, size, size)) # stores the proportion of each species in each cell
  fec_site <- occ_site # recruitment potential for each species in each cell
  
  # summed across species
  occ_tot <- array(0, c(size, size)) # stores the proportion of occuppied sites in each cell
  fec_tot <- occ_tot # stores the total recruitment in each cell
  
  for (t in 2:time) {
    # updating survival and fecundity for each species
    for (i in 1:nsp) {
      occ_site[i, , ] <- res[i, , , t - 1] * SurvivalDD(res[i, , , t - 1], UU[t], SS[t])
      temp <- res[i, , , t - 1] * BB[i, , , t] # seed produced (i.e., recruitment)
      fec_site[i, , ] <- DispersalSeed(kernel, temp, size) # distribute seeds
    }
    
    # updating total occupied site and total fecundity
    occ_tot <- colSums(occ_site, dims = 1) # sums across species for each cell
    fec_tot <- colSums(fec_site, dims = 1) # sums across species for each cell
    
    # updating frequency for each species
    for (i in 1:nsp) {
      res[i, , , t] <- occ_site[i, , ] + (1 - occ_tot) * fec_site[i, , ] / fec_tot
    }
  }
  
  df_param <- data.frame("nspecies" = nsp, "sim_time" = time, "sim_size" = size, "seed" = seed, "param" = param, "driver" = driver)
  result <- list("sim_param" = df_param, "sim_out" = res)
  return(result)
}