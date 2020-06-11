create_driver_function <- function(change_times,
                                   parm_values, 
                                   interpolation = c("linear","constant","spline")){
  interpolation <- match.arg(interpolation)
  if(length(change_times)+1 != length(parm_values)){
    stop("parm_values should be a vector 1 longer than change_times to account for the initial regime")
  }
  stopifnot(is.numeric(change_times))
  stopifnot(is.numeric(parm_values))
  
  interp_values = c(0, change_times)
  
  if(interpolation=="spline"){
    fun <- stats::splinefun(x= interp_values, 
                            y = parm_values,
                            method = "natural")
  }else{
    fun <- stats::approxfun(x = interp_values,
                            y= parm_values,
                            rule = 2,
                            method = interpolation)
  }
  
  fun
}

random_var_matrix <- function(n_sp,diag_min = 0, diag_max =1, off_diag_min = -1,
                              frac_connected = 0.25, 
                              off_diag_max = 1, check_eigen = TRUE){
  
  if(!check_eigen){
    n_iter = 1
    mat <- matrix(runif(n_sp^2, 
                        min = off_diag_min, 
                        max = off_diag_max) *
                    rbinom(n_sp^2,size = 1,prob = frac_connected),
                  nrow = n_sp)
    diag(mat) = runif(n_sp,diag_min,diag_max)
  } else{
    n_iter = 0
    while(check_eigen){
      n_iter = n_iter + 1
      mat <- matrix(runif(n_sp^2, 
                          min = off_diag_min, 
                          max = off_diag_max) *
                      rbinom(n_sp^2,size = 1,prob = frac_connected),
                    nrow = n_sp)
      diag(mat) = runif(n_sp,diag_min,diag_max)
      
      if(max(abs(eigen(mat)$value))<1) check_eigen = FALSE
    }
  }
  attr(mat, "n_iter") = n_iter
  return(mat)
}


sim_troph_triangle <- function(time_out = 100, 
                               n_steps = c(100), 
                               measurement_sigma = c(0,0,0),
                               harvest_rates = create_driver_function(change_times = 50,
                                                                        parm_values = c(0,0.25),
                                                                        interpolation = "constant"),
                               pred_overlap = create_driver_function(change_times = 50,
                                                                     parm_values = c(1,1),
                                                                     interpolation = "constant")
                               ) {
  
  #Describes the basic dynamic system 
  
  library(deSolve)
  library(dplyr)
  library(tidyr)
  
  model_parameters <- list(
    T_mat = 5,    #length of time it takes a juvenile predator to mature to an adult
    m = 0.025,    #mortality rate of adult and juvenile fish
    s = 0.1,     #The stocking rate (amount of new fish being added from outside) for adult predators
    
    harvest_rates = harvest_rates,
    pred_overlap = pred_overlap,

    a_PF_base = 0.1, 
    
    f = 0.5,      #amount of new offspring for each adult predator per unit time
    a_PJ  = 0.05, #Cannibalism rate of adult predators on juveniles
    a_FJ  = 0.1,  #attack rate of forage fish on juvenile predators
    
    r = 0.25,     #population growth rate of forage fish at low densities
    b = 0.005,    #density-dependence term for the forage fish
    a_PF_start = 0.1,   #attack rate of adult predators on forage fish when species fully overlap
    d = 1,       #Stocking rate for forage fish
    
    max_time = time_out
  )
  
  
  
  init_cond <- c(adult = 77, forage = 0.067, juv = 9.37)
  
  troph_tri_static <-  function(t,y, parms){
    #make the current model state variables also available to refer to by name
    adult = y["adult"]
    forage = y["forage"]
    juv = y["juv"]
    
    e <-  parms$harvest_rates(t)
    a_PF <-  parms$pred_overlap(t)*parms[["a_PF_base"]]
    
    #This next code calculates the derivatives at each point in time. 
    #the with(x,...) function here make the model parameters available by name
    #without having to type parms$e*adult + parms$s...
    d_adult <-  with(parms, juv/T_mat - m*adult - e*adult + s) 
    d_forage <-  with(parms, r*forage - b*forage^2 - a_PF*adult*forage + d)
    d_juv <-  with(parms, f*adult - juv/T_mat - m*juv - a_PJ*adult*juv - a_FJ*forage*juv)
    return(list(c(adult=d_adult,forage=d_forage, juv=d_juv)))
  }
  
  
  simulation  <-  ode(y=init_cond,
             times = seq(0,time_out,length.out = n_steps),
             func = troph_tri_static,
             parms = model_parameters)
  
  simulation <- as.data.frame(simulation)

  noisy_obs <- mutate(simulation,
                      adult = adult*rlnorm(n_steps,0, measurement_sigma[1]),
                      juv = juv*rlnorm(n_steps,0, measurement_sigma[2]),
                      forage = forage*rlnorm(n_steps,0, measurement_sigma[3]),
                      )
  noisy_obs <- gather(noisy_obs,key = species, value = abundance_obs, adult, juv,forage)
  
  
  simulation <- gather(simulation,key = species, value = abundance_true, adult, juv,forage)
  simulation <- left_join(simulation, noisy_obs)
  simulation <- mutate(simulation,
                       exploitation_rate = harvest_rates(time),
                       predator_prey_attack = model_parameters[["a_PF_base"]]*pred_overlap(time)
  )
  
  return(simulation)
  
}


sim_var <- function(n_steps = 50, 
                    n_species = 10,
                    regime_change_points = c(),
                    regime_means = "random",
                    regime_coefs = "random",
                    regime_vars  = 1
                    ) {
  
  n_regimes <-  length(regime_change_points)+1
  if(!(regime_means[[1]][[1]]=="random" |(is.list(regime_means)&length(regime_means) == n_regimes)|length(regime_means)==1|length(regime_means)==n_species)){
    stop("regime_means should be one of: 'random', a list of vectors equal to the number of regime change points plus 1, a single number, or a vector of length n_species")
  }
  
  if(!(regime_coefs[[1]][[1]]=="random" |(is.list(regime_coefs)&length(regime_coefs) == n_regimes)|length(regime_coefs)==1|length(regime_coefs)==n_species^2)){
    stop("regime_coefs should be one of: 'random', a list of matrices equal to the number of regime change points plus 1, a single number, or a matrix of length n_species")
  }
  
  if(!((is.list(regime_vars)&length(regime_vars) == n_regimes)|length(regime_vars)==1|length(regime_vars)==n_species)){
    stop("regime_vars should be one of:  a list of vectors equal to the number of regime change points plus 1, a single number, or a vector of length n_species")
  }
  
  regime_change_points <- c(0, sort(regime_change_points), n_steps+1)
  
  if(regime_means[1]=="random"){
    regime_means <-  list()
    for(i in 1:n_regimes){
      regime_means[[i]] <-  runif(n_species)
    }
  } else if(!is.list(regime_means)){
    new_means <- list()
    for(i in 1:n_regimes){
      new_means[[i]] <- regime_means
    }
    regime_means <- new_means
  }
  
  if(regime_coefs[1]=="random"){
    regime_coefs <-  list()
    for(i in 1:n_regimes){
      regime_coefs[[i]] <-  random_var_matrix(n_species)
    }
  } else if(!is.list(regime_coefs)){
    new_coefs <- list()
    for(i in 1:n_regimes){
      new_coefs[[i]] <- regime_coefs
    }
    regime_coefs <- new_coefs
  }
  
  if(!is.list(regime_vars)){
    new_vars <- list()
    for(i in 1:n_regimes){
      new_vars[[i]] <- regime_vars
    }
    regime_vars <- new_vars
  }
  
  
  
  
  
  out_data <- matrix(nrow = n_steps, ncol = n_species)
  out_data[1,] <- regime_means[[1]]
  regime <- 1
  
  for(i in 2:n_steps){
    if(i>regime_change_points[regime+1]) regime = regime +1
    out_data[i,] <- regime_means[[regime]] + regime_coefs[[regime]]%*%(out_data[i-1,]-regime_means[[regime]]) + rnorm(n_species)
  }
  
  out_data <- as.data.frame(out_data)
  out_data <- dplyr::mutate(out_data, 
                     time = 1:n_steps)
  out_data <- tidyr::gather(out_data,species,abundance,-time)
  
  out_list <- list(simulation = out_data,
                   n_regimes = n_regimes,
                   regime_change_points = regime_change_points,
                   regime_means = regime_means,
                   regime_coefs = regime_coefs)
  out_list
}
  