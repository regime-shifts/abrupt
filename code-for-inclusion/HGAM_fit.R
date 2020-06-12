library(lazyeval)


fit_hgam_model <- function(data, 
                           time_var, 
                           group_var,
                           obs_var,
                           family = Gamma(link="log"),
                           formula = NULL,
                           k_temporal = 10,
                           metrics = c("per-capita mean",
                                       "mean",
                                       "per-capita sd",
                                       "instantaneous fisher",
                                       "per-capita instantaneous fisher"),
                           n_posterior_sims=500,
                           CI_alpha = 0.05){
  
  #check if all variables are in the data: 
  if(!all(c(time_var,group_var,obs_var)%in%names(data))) stop("all named variables should be in the included data frame.")
  
  #check if group_var is a factor:
  if(!is.factor(data[,group_var])) stop("group_var should be a factor.")
  
  
  if(is.null(formula)){
    formula <-  glue("{obs_var} ~ s({time_var}, {group_var}, k = {k_temporal},bs='fs')")
  }
  
  hgam_model <- mgcv::gam(formula = as.formula(formula), family = family,data= data,
                          method= "REML")
  out <- list()
  out$hgam_model <- hgam_model
  
  out$metrics <- calc_hgam_metrics.gam(hgam_model,
                                   data =data,
                                   time_var = time_var, 
                                   group_var = group_var, 
                                   delta = 0.01,
                                   metrics = metrics,
                                   return_val = "summary",
                                   n_posterior_sims = n_posterior_sims,
                                   CI_alpha = CI_alpha )
  
  class(out) = c("abrupt", "abrupt_hgam")
  out
  
}



calc_deriv1 = function(before, after,delta) {
  (after-before)/(2*delta)
  
}

calc_deriv2 <-  function(at, before, after,delta) {
  (before + after-2*at)/delta^2
}

calc_fisher_instant <- function(deriv1,deriv2){
  sum(deriv1*deriv2)^2/sum(deriv1^2)^3
}



calc_hgam_metrics.default <- function(object, ...){
  UseMethod("calc_hgam_metrics", object)
}

calc_hgam_metrics.abrupt <- function(object, ...){
  NULL
  #still to be defined
}


calc_hgam_metrics.gam <- function(object, 
                              data, 
                              time_var,
                              group_var, 
                              delta=0.01, 
                              metrics = c("per-capita mean",
                                         "mean",
                                         "per-capita sd",
                                         "sd",
                                         "instantaneous fisher",
                                         "per-capita instantaneous fisher"),
                              return_val = c("summary", "sims"),
                              n_posterior_sims=500,
                              CI_alpha = 0.05){
  if(!all(metrics%in% c("per-capita mean","mean","per-capita sd","sd", "instantaneous fisher","per-capita instantaneous fisher"))) stop('All metrics should be one of "per-capita mean","mean","per-capita sd", "instantaneous fisher", or "per-capita instantaneous fisher"')
  
  return_val <- match.arg(return_val)
  
  
  #Randomly generate coefficients from the model posterior
  random_coefs <-  t(mgcv::rmvn(n_posterior_sims, 
                              mu = coef(object),
                              V = vcov(object)))
  
  data_after <- data
  data_after[,time_var] <- data[,time_var] + delta
  
  data_before <- data
  data_before[,time_var] <- data[,time_var] - delta
  
  pred        <- predict(object, data, type="lpmatrix") %*% random_coefs
  pred_before <- predict(object, data_before, type="lpmatrix") %*% random_coefs
  pred_after  <- predict(object, data_after, type="lpmatrix") %*% random_coefs
  
  simplified_data <- data[,c(time_var, group_var)]
  
  out <- list()
  
  for(i in metrics){
    
    if(length(grep("per-capita",i))>0 & object$family$link !="log") warning("per-capita indices for HGAM models will perform best for families using a log-link function. Non-log-link functions may suffer from numerical issues due to rounding of large numbers, or result in NaN values for indices.")
    if(i %in% c("mean", "instantaneous fisher", "sd")){
      pred_adj <- object$family$linkinv(pred)
      pred_adj_before <- object$family$linkinv(pred_before)
      pred_adj_after <- object$family$linkinv(pred_after)
    } else if(object$family$link !="log"){
      warning("per-capita indices for HGAM models will perform best for families using a log-link function. Non-log-link functions may suffer from numerical issues due to rounding of large numbers, or result in NaN values for indices.")
      pred_adj <- log(object$family$linkinv(pred))
      pred_adj_before <- log(object$family$linkinv(pred_before))
      pred_adj_after <- log(object$family$linkinv(pred_after))
    } else{
      pred_adj <- pred
      pred_adj_before <- pred_before
      pred_adj_after <-  pred_after
    }
    deriv1 <- calc_deriv1(pred_adj_before, pred_adj_after,delta)
    deriv1 <- as.data.frame(deriv1)
    colnames(deriv1) <- paste0("sim", 1:n_posterior_sims)
    deriv1 <- dplyr::bind_cols(simplified_data, deriv1)
    deriv1 <- tidyr::gather(deriv1,key = "sim", 
                            value = "deriv1",
                            tidyselect::starts_with("sim"))
    
    deriv2 <- calc_deriv2(pred_adj,pred_adj_before, pred_adj_after,delta)
    deriv2 <- as.data.frame(deriv2)
    colnames(deriv2) <- paste0("sim", 1:n_posterior_sims)
    deriv2 <- dplyr::bind_cols(simplified_data, deriv2)
    deriv2 <- tidyr::gather(deriv2,key = "sim", 
                            value = "deriv2",
                            tidyselect::starts_with("sim"))
    
    
    current_metric <- dplyr::left_join(deriv1, deriv2)
    current_metric <- ungroup(current_metric)
    current_metric <- dplyr::group_by(current_metric,
                                      .data[["sim"]],
                                      .data[[time_var]]
                                      )
    
    
    if(stringr::str_detect(i, "sd")){
      current_metric <- dplyr::summarize(current_metric, 
                                         metric = sd(.data[["deriv1"]]))
    } else if(stringr::str_detect(i, "mean")){
      current_metric <- dplyr::summarize(current_metric, 
                                         metric = mean(.data[["deriv1"]]))
    } else {
      current_metric <- dplyr::summarize(current_metric, 
                                         metric = calc_fisher_instant(.data[["deriv1"]], .data[["deriv2"]]))
    }
    
    if(return_val =="summary"){
      current_metric <- dplyr::ungroup(current_metric)
      current_metric <- dplyr::group_by(current_metric, 
                                        .data[[time_var]])
      current_metric <- dplyr::summarize(current_metric,
                                         metric_median = median(.data[["metric"]]),
                                         metric_lower = quantile(.data[["metric"]],
                                                                 probs = CI_alpha/2),
                                         metric_upper = quantile(.data[["metric"]],
                                                                 probs = 1-CI_alpha/2))
    }
    
    out[[i]] = current_metric
    
    
  }
  
  out
  
  
}


