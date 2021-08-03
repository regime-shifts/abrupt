##' Function that fits the Ricker model
##'
##' @param data ToDo
##'
##' @importFrom minpack.lm nlsLM
##' @importFrom stats AIC
`rickerfit` <- function (data){
    ## create an initial estimate of k to aide model convergence
    kest <- mean(data$Nt)

    ## supress warnings about failed convergence in oddball fits
    ##  - these models won't be favoured by AIC anyway
    op <- options(warn=-1)
    on.exit(options(op)) # reverse this on exit of function
    
    ## fit the model
    ricker.model <- tryCatch(nlsLM(Nt1 ~ Nt * exp(r * (1 - Nt/k)),
                                   start = list(r = 1.5, k = kest), data = data),
                             error = function(e) NULL)
    ## What outputs do we need from each run? AIC, r and k, and their resepective
    ## errors. Want to create a vector with this information in it so we can use
    ## this information later
    if(is.list(ricker.model)) {
        output <- c(AIC(ricker.model), #AIC
                    summary(ricker.model)$coefficients[1,1], # r
                    summary(ricker.model)$coefficients[1,2], # se for r
                    summary(ricker.model)$coefficients[2,1], # k
                    summary(ricker.model)$coefficients[2,2]) # se for k
    } else {
        ## if the model fails to converge, give it an an arbitrarily high but finite AIC
        output < -c(100000000, 0, 0, 0, 0)
    }
    
    output # return
}
