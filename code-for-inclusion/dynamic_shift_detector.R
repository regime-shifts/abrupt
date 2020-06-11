#script to analyse time series data to determine break points

#assume data is in form of data frame where first column is year, second is abundance (Nt)

#create function that makes an additional response variable Nt1
#where Nt1 is the abundance in the next year

addNt1<-function(data){
  #create empty vector to feed 'next year' abundance into
  Nt1 = c()
  #pull out values from the abundance column-
  #for each year, go to the NEXT year's abundance, put that in the vector
  for (i in 1:(length(data[,2])-1)) {
    Nt1 = c(Nt1,(data[,2])[i+1])
  }
  #cut out last sampling year, because there is no Nt+1 for that year
  data<-data[which(data[,1]<max(data[,1])),]
  #add the vector of new values to the data fram
  data$Nt1<-Nt1
  #rename columns in data set so that they're consistent no matter 
  #what they're named before importing
  names(data)[1:2]<-c("year", "Nt")
  #spit out that new data frame
  return(data)
}


#lets also create a function that will calculate an AICc correction factor
#which can be added to the AIC values we compute- because we have small sample
#sizes we need to more heavily penalize models with more parameters

AICcorrection<-function(data, breaks){
  a<-3*(breaks+1) # a is the number of parameters estimated each fit (r, k, error)
  correction<-(2*a*(a+1))/(nrow(data)-a-1)
  return(correction)
}


#load minpack.lm
#this is my prefered nonlinear package- fewer convergence issues.
#we'll want to use its nlsLM function to fit the Ricker equation.
## library(minpack.lm)

#create function that fits the ricker model

## use the rickerfit() in ricker-model.R

## rickerfit<-function (data){
##   #create an initial estimate of k to aide model convergence
##   kest<-mean(data$Nt)
##   #supress warnings about failed convergence in oddball fits- these models won't be favoured by AIC anyway
##   options(warn=-1)
##   #fit the model
##   ricker.model<-tryCatch(nlsLM(Nt1~ Nt*exp(r*(1- Nt/k)), start=list(r=1.5, k=kest), data=data), error=function(e) NULL)
##   #What outputs do we need from each run? AIC, r and k, and their resepective errors.
##   #we'll want to create a vecor with this information in it so we can use this information later
##   if(is.list(ricker.model)){
##     output<-c(AIC(ricker.model), #AIC
##               summary(ricker.model)$coefficients[1,1], # r
##               summary(ricker.model)$coefficients[1,2], # se for r
##               summary(ricker.model)$coefficients[2,1], # k
##               summary(ricker.model)$coefficients[2,2]) # se for k
##   }else{
##     output<-c(100000000, 0, 0, 0, 0)#if the model fails to converge, give it an an arbitrarily high but finite AIC
##   }
  
##   options(warn=0)#turn warnings back on
##   return(output)
## }



#now we need to build a tool that will cut a time series up,
#fit the model, and spit out relevant parameters

#we probably want to write a function for each break point combo
#then smash them into a function that figures out what the best model is


# Okay, so now we want a generalized model so it can handle N break point cases
# this is a bit more complex because it'll require some recursion from within the function
# and because there will be N break cases, we won't be able to plan a 2 dimensional data frame 
# that can have columns for all break cases, as I'd originally hoped. I think the best way to do this is by creating a 
# list of break points and an associated list of AICs for the fit, and then use the sum function
# for the AICs when including them in the output, but include the break point list as a single
#element in the output data frame
# so we want a 2 column data frame with a list of the breaks, and the 
#resultant AICc as columns, respectively.


breaks<-list() #create empty LIST for storing breaks
fit<-list() #create empty LIST for storing associated AICs
out.frame<-data.frame(matrix(vector(), 0, 2,
                             dimnames=list(c(), c("Breaks", "AICs"))),
                      stringsAsFactors=F)

splitnfit<-function(data, breaks, fit, out.frame){ #need to include vectors for breaks and fit to re-feed into this function
  #first fit no break model, put it in the data frame
  fit.0<-rickerfit(data) #fit the model
  
  out<-cbind(list(max(data$year)), list(fit.0[1])) #output vector with no breaks
  out.frame<-rbind(out.frame, out)
  
  Break1<-min(data$year)+3 #create first breakpoint four years into the time series to avoid overfitting
  while(Break1<(max(data$year))){
    part1<-data[which(data$year<Break1),] #create subsets at the breakpoint
    part2<-data[which(data$year>(Break1-1)),]
    if(nrow(part1)>3 & nrow(part2)>3){ #constrain model to run only when 4 or more points are present
      fit1<-rickerfit(part1) #fit the model to part 1
      fit2<-rickerfit(part2) #fit the model to part 2
      breaks.1<-c(breaks, max(part1$year), max(part2$year)) #breaks for one break
      fit.1<-c(fit, fit1[1], fit2[1]) #fit for one break
      out[1]<-list(breaks.1)#create output vector of two lists
      out[2]<-list(fit.1)
      out.frame<-rbind(out.frame, out) #bind it to previous results
    }
    Break1<-Break1+1 #move the break to next year
  }
  #rename columns in output for some reason
  colnames(out.frame)<- c("Breaks", "AICs")
  return(out.frame)
}


#lets take our results frame from splitnfit and test each break combinaton for the
# the ability to be broken up more. This situation will occur when the last two years
# in the breaks list are more than 5 years apart

findbreakable<-function(data){ #create a function that finds if the last subset of the data is still breakable
  breakable<-c() #create vector to put results in
  for (i in 1:nrow(data)){ #for each row in the data frame
    if (length(unlist(data$Breaks[i]))>1){ #if the data has been subset
      breakvector<-unlist(data$Breaks[i]) #create a vector of the breaks
      L<-length(breakvector) # find out how long breakvector is
      difference<-breakvector[L]-breakvector[L-1] #find out how big the last subset is
      if(difference>5){ #we can break it down more if the last subset has more than 5 points in it
        breakable.i<-TRUE
      }else{
        breakable.i<-FALSE #don't break more if the subset is 5 or smaller
      }
    }else{
      breakable.i<-FALSE # don't break it down more if the data is frm the zero breaks model
    }
    breakable<-c(breakable, breakable.i)
  }
  return(breakable)
  
}


#create a function that uses findbreakable to apply splitntfit to the datasets that are still breakble
out.frame<-data.frame(matrix(vector(), 0, 2, #create an empty data frame we can add to later
                             dimnames=list(c(), c("Breaks", "AICs"))),
                      stringsAsFactors=F)

subsequentsplit<-function(fitdata, rawdata){ 
  keepers<-findbreakable(fitdata) #find subsets that are still breakable
  newfitdata<-fitdata[which(keepers==TRUE),] #create new data frame with only these data in it
  breaklist<-newfitdata$Breaks #pull out our two operational objects out of data frame
  fitlist<-newfitdata$AICs
  result<-data.frame(matrix(vector(), 0, 2,
                            dimnames=list(c(), c("Breaks", "AICs"))),
                     stringsAsFactors=F)
  if(nrow(newfitdata)>0){ #if there are subsets with still breakable data
    for(i in 1:nrow(newfitdata)){ #for each row in the new frame we need to break down
      breakvector<-unlist(breaklist[i]) #turn the list element back into a vector
      fitvector<-unlist(fitlist[i])
      breakvector<-breakvector[-length(breakvector)]#remove the last element from each
      fitvector<-fitvector[-length(fitvector)]
      cullpoint<-max(breakvector) #find point of last break
      testdata<-rawdata[which(rawdata$year>cullpoint),]
      out<-splitnfit(testdata, breakvector, fitvector, out.frame)
      out<-out[-1,] #remove the no-break fit, we don't need that
      result<-rbind(result, out)
    }
  }else{ #if there is no more room for breaks in any of the fits
    result<-result #leave result empty
  }
  return(result)
}


#okay, now let's put this all together into a function that will fit all the breaks

nbreaker<-function(data){
  breaks<-list() #create empty LIST for storing breaks
  fit<-list() #create empty LIST for storing associated AICs
  out.frame<-data.frame(matrix(vector(), 0, 2,
                               dimnames=list(c(), c("Breaks", "AICs"))),
                        stringsAsFactors=F)
  onebreak<-splitnfit(data, breaks, fit, out.frame)#get fits for zero and one break
  feed<-onebreak #create derrived data to feed into the fit
  out<-onebreak #prepare these data to be output
  
  keepers<-findbreakable(onebreak) #find subsets that are still breakable
  newfitdata<-onebreak[which(keepers==TRUE),] #create new data frame with only these data in it
  stillbreakable<-nrow(newfitdata)
  while(stillbreakable>0){ #if there is data that can still be broken up
    feed<-subsequentsplit(feed, data) #fit the subsequent split using data fed from last breaks
    out<-rbind(out, feed) #attach this break iteration to the data frame
    #see if there's more breaks to be had
    keepers<-findbreakable(feed) #find subsets that are still breakable
    newfitdata<-feed[which(keepers==TRUE),] #create new data frame with only these data in it
    stillbreakable<-nrow(newfitdata)
    
  }
  return(out)
}


#victory! now we need to extract the data we've produced from this data frame, sum up AICs and add
# the corrections, and identify the best models

#function to tally up AICs
AICtally<-function(data){ #create a function that adds up the AICs in the list
  fitdata<-nbreaker(data)
  AICtots<-c() #create vector to put results in
  N_AIC<-c() #create a vector to put counts of AICs in
  for (i in 1:nrow(fitdata)){ #for each row in the data frame
    AICsvector<-unlist(fitdata$AICs[i]) #create a vector of the AICs for the fit
    total<-sum(AICsvector)#add up the AICs
    N<-length(AICsvector)#count the AICs
    AICtots<-c(AICtots, total)# add the total for each cell to the vector
    N_AIC<-c(N_AIC, N) #add the count of AICs
  }
  out<-as.data.frame(cbind(AICtots, N_AIC)) #bind the outputs into a data frame
  colnames(out)<- c("AICtot", "Nfits") #name the columns
  
  #now we want to calculate the AICc correction for each
  out$Nbreaks<-out$Nfits-1
  out$AICc<-out$AICtot+AICcorrection(data, out$Nbreaks) #n breaks = n fits-1, add correction to AIC
  
  return(out)
  
}


#see if AICtally outputs will stick to nbreaker outputs- and create a function that does this
allfits<-function(data){
  out<-as.data.frame(cbind(nbreaker(data), AICtally(data))) #stick output from two functions into a df
  return(out)
}




#create a function that finds equivalent fits in n breakpoint data
#add functionality for switching between AIC and AICc, with AIC as default

equivalentfit<-function(data, criterion){
  breakset<-allfits(data) #generate matrix of fits by breakpoints
  if(missing(criterion)){ #set AIC as the default criterion
    criterion<-"AIC"
  }
  if(criterion=="AIC"){
    AICbest<-min(breakset$AICtot) #find best AIC in the set
    deltaAIC<-AICbest+2 # create rule for equivalent models
    out.frame<-breakset[which(breakset$AICtot<deltaAIC),] #cut out all data but equivalent models
  }
  if (criterion=="AICc"){
    AICbest<-min(breakset$AICc) #find best AIC in the set
    deltaAIC<-AICbest+2 # create rule for equivalent models
    out.frame<-breakset[which(breakset$AICc<deltaAIC),] #cut out all data but equivalent models
  }
  return(out.frame)
}



#but equivalent fit is one thing- if there's equivalent fit and one of the model set has
# fewer parameters, we obviously want to go with that. create a function that does that but this time for n breaks

bestfit<-function(data, criterion){
  if(missing(criterion)){ #set AIC as the default criterion
    criterion<-"AIC"
  }

  breakset<-equivalentfit(data, criterion) #get set of eqivalent models
  if(criterion=="AIC"){
    AICbest<-min(breakset$AICtot) #find best AIC in the set
    out.frame<-breakset[which(breakset$AICtot==AICbest),] #pull models with numerically lowest AIC
  }
  if (criterion=="AICc"){
    AICbest<-min(breakset$AICc) #find best AIC in the set
    out.frame<-breakset[which(breakset$AICc==AICbest),] #pull models with numerically lowest AIC
  }
  
  return(out.frame)
}


#cool. looks like we can now adapt the model fitting function for these n breakpoint data

bestmodel<-function(data, criterion){
  modelspecs<-bestfit(data, criterion) #get the particulars of the best model
  out.frame<-data.frame(matrix(vector(), 0, 7,
                               dimnames=list(c(), 
                                             c("Year1", "Year2", "AIC", "r", "rse", "k", "kse"))),
                        stringsAsFactors=F)#Create a place to put our data
  breakvector<-unlist(modelspecs$Breaks[1])#pull out a vector of the max year in each fit
  
  if (length(breakvector)<=1){ #if there's no breaks
    fit<-rickerfit(data)
    output<-c(min(data$year), max(data$year), fit) #fit whole data series + output results
    out.frame<-rbind(out.frame, output)
    
  } else {
    for (i in 1:(length(breakvector))){ #for all breakpoints, including the end of the time series, in order
      part1<-data[which(data$year<breakvector[i]+1),] #create subsets at the breakpoint
      part2<-data[which(data$year>breakvector[i]),]
      fit1<-rickerfit(part1) #fit first segment
      output<-c(min(part1$year), max(part1$year), fit1)#save results of fitting segment in vector
      out.frame<-rbind(out.frame, output)#put output for segment in a data frame
      data<-part2 #update data to cull out already fitted segments 
    }
    
  } 
  colnames(out.frame)<- c("Year1", "Year2", "AIC", "r", "rse", "k", "kse")
  return(out.frame)
}

#adapt this function to provide specs of any fit information, as output in the format of the bestfit function

modelspecification<-function(specs, data){
  modelspecs<-specs #get the particulars of the best model
  out.frame<-data.frame(matrix(vector(), 0, 7,
                               dimnames=list(c(), 
                                             c("Year1", "Year2", "AIC", "r", "rse", "k", "kse"))),
                        stringsAsFactors=F)#Create a place to put our data
  breakvector<-unlist(modelspecs$Breaks[1])#pull out a vector of the max year in each fit
  
  if (modelspecs$Nbreaks[1]==0){ #if there's no breaks
    fit<-rickerfit(data)
    output<-c(min(data$year), max(data$year), fit) #fit whole data series + output results
    out.frame<-rbind(out.frame, output)
    
  } else {
    for (i in 1:(length(breakvector))){ #for all breakpoints, including the end of the time series, in order
      part1<-data[which(data$year<breakvector[i]+1),] #create subsets at the breakpoint
      part2<-data[which(data$year>breakvector[i]),]
      fit1<-rickerfit(part1) #fit first segment
      output<-c(min(part1$year), max(part1$year), fit1)#save results of fitting segment in vector
      out.frame<-rbind(out.frame, output)#put output for segment in a data frame
      data<-part2 #update data to cull out already fitted segments 
    }
    
  } 
  colnames(out.frame)<- c("Year1", "Year2", "AIC", "r", "rse", "k", "kse")
  return(out.frame)
}


#looks like that works! Okay! put it all together like we did for the 2 break model

DSdetector<-function(data, criterion){ #use raw time series data
  #plot the data
  plot(data)
  data1<-addNt1(data)
  plot(data1$Nt, data1$Nt1)
  #give an output of all possible break point combinations tested
  writeLines(paste("Here are the break points for all models tested"))
  print(allfits(data1))
  #output models with equivalent performance
  writeLines(paste("Here is the set of best performing models"))
  print(equivalentfit(data1, criterion))
  #output model with best performance
  writeLines(paste("Here is the best model- the one with the lowest AIC/AICc, or fewest parameters, in case of a tie"))
  print(bestfit(data1, criterion))
  # output regression parameters of best model
  writeLines(paste("Here is the set of regression parameters"))
  print(bestmodel(data1, criterion))
}



#looks like we have a working model! Boom!

#Okay, but what about 'break weights'- how 'strong' are the given breaks? Let's take some inspiration
#from multimodel interference and calculate the weight of a given break that the model finds
#first, we need to calculate the model weights for all of the fits

allweights<-function(data, criterion){
  breakset<-allfits(data) #generate matrix of fits by breakpoints
  if(missing(criterion)){ #set AIC as the default criterion
    criterion<-"AIC"
  }
  if(criterion=="AIC"){
    AICbest<-min(breakset$AICtot) #find best AIC in the set
    deltaAIC<-breakset$AICtot-AICbest
    sumdeltaAIC<-sum(exp(-deltaAIC/2))
  }
  if (criterion=="AICc"){
    AICbest<-min(breakset$AICc) #find best AIC in the set
    deltaAIC<-breakset$AICc-AICbest
    sumdeltaAIC<-sum(exp(-deltaAIC/2))
  }
  breakset$modelweights<-(exp(-deltaAIC/2))/sumdeltaAIC
  
  out.frame<-breakset[which(breakset$modelweights>0.001),] #cut out all models with weight less than 0.01, just to speed things up
  return(out.frame)
}

#now we need a function that calculates parameter weights for the breaks
breakweights<-function(data, criterion){
  modelset<-allweights(data, criterion)#find all models with weights over 0.001
  breaksfound<-sort(unique(unlist(modelset$Breaks))) #get a list of unique breaks, sort it in numerical order
  breakweights<-c() #blank vector to put the breaks in
  for(i in 1:length(breaksfound)){
    weightcounter<-0
    for (j in 1:length(modelset$modelweights)){
      if(breaksfound[i] %in% unlist(modelset$Breaks[j])){
        weightcounter<-weightcounter+modelset$modelweights[j] #add model weight of parameter where that break appears
      }else{
        weightcounter<-weightcounter+0
      }
    }
    breakweights<-c(breakweights, weightcounter)
  }
  correctedweights<-breakweights/max(breakweights) #breakweights is based on best model set- correct to normalize by break at the end of time series
  out.frame<-as.data.frame(cbind(breaksfound, breakweights, correctedweights))
  return(out.frame)
}

