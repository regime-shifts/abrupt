

sim_troph_triangle <- function(time_out, 
                               n_steps, 
                               measurement_sigma = c(0,0,0),
                               harvest_start = 0, 
                               harvest_end = 0.25, 
                               harvest_shape = 0,
                               pred_overlap_start = 1,
                               pred_overlap_end   = 1,
                               pred_overlap_shape = 0
                               ) {
  
  #Describes the basic dynamic system 
  
  library(deSolve)
  library(dplyr)
  library(tidyr)
  
  model_parameters <- list(
    T_mat = 5,    #length of time it takes a juvenile predator to mature to an adult
    m = 0.025,    #mortality rate of adult and juvenile fish
    s = 0.05,     #The stocking rate (amount of new fish being added from outside) for adult predators
    
    e_start = harvest_start,
    e_end   = harvest_end,
    
    a_PF_start = pred_overlap_start,
    a_PF_end   = pred_overlap_end,
    
    
    f = 0.5,      #amount of new offspring for each adult predator per unit time
    a_PJ  = 0.05, #Cannibalism rate of adult predators on juveniles
    a_FJ  = 0.1,  #attack rate of forage fish on juvenile predators
    
    r = 0.25,     #population growth rate of forage fish at low densities
    b = 0.005,    #density-dependence term for the forage fish
    a_PF_start = 0.1,   #attack rate of adult predators on forage fish when species fully overlap
    d = 0.5,       #Stocking rate for forage fish
    
    max_time = time_out
  )
  
  
  
  init_cond <- c(adult = 77, forage = 0.067, juv = 9.37)
  
  troph_tri_static <-  function(t,y, parms){
    #make the current model state variables also available to refer to by name
    adult = y["adult"]
    forage = y["forage"]
    juv = y["juv"]
    
    e <-  parms[["e_start"]] + (parms[["e_end"]]-parms[["e_start"]])*t/parms[["max_time"]]
    a_PF <-  parms[["a_PF_start"]] + (parms[["a_PF_end"]]-parms[["a_PF_start"]])*t/parms[["max_time"]]
    
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
                       exploitation_rate = harvest_start + (harvest_end-harvest_start)*time/time_out,
                       predator_prey_attack = pred_overlap_start + (pred_overlap_end-pred_overlap_start)*time/time_out
  )
  
  return(simulation)
  
}
  
  