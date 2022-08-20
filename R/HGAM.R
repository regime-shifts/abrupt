##' Fits an HGAM model
##'
##' @param data Data frame used for fitting. Should be constructed by the
##'   abr_data function
##' @param type The type of regime detection model. One of "t", "st","s", "gt",
##'   "gst","gs". "t" indicates a temporal smoother, "s" indicates a spatial
##'   smoother, and "g" indicates a grouped (HGAM) smoother
##' @param bs_t/bs_s/bs_g basis types to be used for the temporal, spatial, and
##'   group-level smoothers
##' @param k_t/k_s/k_g number of basis functions to use for temporal, spatial,
##'   and group-level smoothers
##' @param family The family mgcv or brms will use to fit the model. 
##'
##'
##' @return A fitted object of class 'abrupt_model' and 'HGAM'. 
##'
##' @depend mgcv
##' @suggest brms
`abr_fit.hgam` <- function(data,type = "t", 
                          bs_t  = "tp", bs_s = "gp", bs_g = "re",
                          k_t = 20,k_s=20, k_g=NA,  
                          func = "gam", 
                          family = "gaussian",
                          calc_breaks = FALSE,
                          break_args = list(),
                          ... ) {
  
  #Currently, just loading the gam function via mgcv::gam does not work,
  #as it is unable to find other associated functions from the mgcv namespace
  #like ldTweedie
  library(mgcv)
  
  #is the data of the right type?
  if (!inherits(data,"abdata")) {
    stop("'data' must be constructed by the the 'abr_data' function.")
  }
  
  if(grepl("s", type,fixed = TRUE)){
    #Add code turning spatial column into either coordinates (for point-like 
    #data) or an adjacency matrix (for polygon-like data)
  }
  
  ## To do: add checks to ensure that basis types, families provided are valid
  
  if(type=="t"){
    model <- mgcv::gam(value~s(time, bs=bs_t), 
                       family  = family, 
                       method = "REML", 
                       data =data)
  }
  if(type=="gt"){
    model <- mgcv::gam(formula = value~t2(time, variable, 
                                          bs= c(bs_t,bs_g), 
                                          k=c(k_t,k_g)), 
                       family  = family, 
                       method = "REML", 
                       data =data) 
  }
  out <- list()
  out$data <- data
  out$type <- type
  out$model <- model
  out$break_args <- break_args
  class(out) <- c("abrupt_model","hgam", class(model))
  
  if(calc_breaks){
    out$breaks <- abr_breaks.hgam(out)
  }else{
    out$breaks <- NULL
  }
  
  out
}



##' Extracts breakpoints from an HGAM model
##'
##' @param abr_fitobj A previously fitted
##'
##' @return A list containing `seed`, the supplied seed, `initial_state`, the
##'   initial state of the RNG before the seed was set, and `kind`, the type of
##'   RNG used, as returned by [base::RNGkind()].
##'
##' @depend mgcv
##' @suggest brms
##' @importFrom dplyr mutate group_by summarize 
`abr_breaks.hgam` <- function(abr_fitobj,
                            breaktype = NULL,
                            stepsize = NULL,
                            delta = NULL, 
                            nsims = NULL,
                            transform = NULL,
                            aggregate = NULL,
                            threshold = NULL){
  
  
  t_range <- range(abr_fitobj$data$time)
  
  if(is.null(stepsize)){
    if(is.null(abr_fitobj$break_args$stepsize)){
      #this is currently pretty ugly; it would be nicer with pipes but
      #I need to figure out whether it's reasonable to import pipes into the
      #function
      stepsize <- group_by(abr_fitobj$data, variable)
      stepsize <- summarize(stepsize, tdiff = min(diff(time)))
      stepsize <- min(stepsize$tdiff)
      if(stepsize ==0){
        stepsize <- (t_range[2] -t_range[1]) /100
      }
    } else{
      stepsize <- abr_fitobj$break_args$stepsize
    }
  }
  
  if(is.null(nsims)){
    if(is.null(abr_fitobj$break_args$nsims)){
      nsims  <- 100
    } else nsims <- abr_fitobj$break_args$nsims
  }
  
  
  if(is.null(delta)){
    if(is.null(abr_fitobj$break_args$delta)){
      delta  <- 1e-6
    } else delta <- abr_fitobj$break_args$delta
  }
  
  
  if(is.null(transform)){
    if(is.null(abr_fitobj$break_args$transform)){
      transform  <- identity
    } else transform <- abr_fitobj$break_args$transform
  }
  
  
  
  if(is.null(threshold)){
    if(is.null(abr_fitobj$break_args$threshold)){
      threshold  <- 0
    } else threshold <- abr_fitobj$break_args$threshold
  }
  
  type <- abr_fitobj$type
  
  if(grepl("g", type,fixed = TRUE)){
    variable <- unique(abr_fitobj$data$variable)
  } else{
    variable <- NA
  }
  
  #Creating synthetic prediction data to 
  pred_at <- expand.grid(time = seq(t_range[1],t_range[2], by= stepsize),
                           variable = variable,stringsAsFactors = FALSE)
  
  pred_before <- pred_at
  pred_before$time <- pred_before$time - delta
  
  pred_after <- pred_at
  pred_after$time <- pred_after$time + delta
  
  lp_at <- mgcv::predict.gam(abr_fitobj$model, 
                             newdata = pred_at,
                             type = "lpmatrix")
  lp_before <- mgcv::predict.gam(abr_fitobj$model, 
                                 newdata = pred_before,
                                 type = "lpmatrix")
  lp_after <- mgcv::predict.gam(abr_fitobj$model, 
                                newdata = pred_after,
                                type = "lpmatrix")
  
  #taking advantage of the linearity of the derivative to just
  #calculate one lp matrix. Note that this won't work for nonlinear 
  #transformations of the data, but we can post-multiply the 
  #derivative by the derivative of the transformation w.r.t. 
  #x to get the new derivative. 
  lp_deriv <- calc_1st_deriv(lp_before, lp_after,delta)
  
  vcv <- abr_fitobj$model$Vc
  coef_means <- abr_fitobj$model$coefficients
  
  coef_sims <- t(mgcv::rmvn(n = nsims,mu = coef_means, V = vcv))
  deriv_fit <- lp_deriv%*%coef_sims
  
  if(length(threshold)==2){
    under_threshold <- rowMeans(deriv_fit< (threshold[1]))
    over_threshold  <- rowmeans(deriv_fit> (threshold[2]))
  }else {
    under_threshold <- rowMeans(deriv_fit< -abs(threshold))
    over_threshold  <- rowMeans(deriv_fit> abs(threshold))
  }
  #currently assuming two-sided intervals; need to add an option to allow for
  #one-sided intervals 
  pred_at$pos_change <- over_threshold
  pred_at$neg_change <- under_threshold
  pred_at$change_prob <- pmin(1-over_threshold, 1-under_threshold)
  
  
  pred_at
}
  

  
  
  
  