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
  
  