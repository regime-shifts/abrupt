#script for creating simulated data under a variety of parameters
#and then determining if breakweights function correctly determines these parameters

#get the regime shift detector functions into memory
source("dynamic_shift_detector.R")

#create a function that will make fake data based on specified parameters
#assume change, noise is given in percent (0-100) scale, as is change to r, k

fakedata<-function(startyear, Nyears, startPop, noise, startK, startR, breaks, changeK, changeR){
  if(missing(startyear)){ #set default values for all paremeters
    startyear<-1900
  }
  if(missing(Nyears)){
    Nyears<-20
  }
  #in order to set the default breaks, I need to know what the last year in the time range is
  #so let's create a vector of years going into the time series now
  year<-seq(startyear, (startyear+Nyears-1)) #-1 because time range is inclusive of year 0
  lastyear<-max(year)-1# max is second last year because we don't have Nt1 for last year
  
  if(missing(startPop)){#let's make popultion size default to 1000
    startPop<-1000
  }
  if(missing(noise)){#no noise by default
    noise<-0
  }
  if(missing(startK)){ #so start population can grow by default
    startK<-1500
  }
  if(missing(startR)){
    startR<-1.5
  }
  if(missing(breaks)){
    breaks<-list()#no break model, null list
  }
  if(missing(changeK)){ # by default, don't change K after a break
    changeK<-0
  }
  if(missing(changeR)){ #same with r
    changeR<-0
  }
  
  #create a vector for noise for each year- it will be random, normally distributed error
  noisevector<-c()# make an empty vector
  for (i in 1:(length(year))){
    randomnoise<-rnorm(1, mean=0, sd=noise)#choose a random value with a sd of our % noise
    instant.buzz<-1+(randomnoise)/100 #generate an instantaneous buzz :)
    noisevector<-c(noisevector, instant.buzz) #add that to the vector
  }
  
  #create a vector of when regime shifts will occur
  change<-c(FALSE)# make a vector with first value false- cannot have a change in first year
  for (i in 1:(length(year)-1)){
    if(any(breaks==year[i])){
      switch<-TRUE
    }else{
      switch<-FALSE
    }
        change<-c(change, switch) #add that to the vector
  }
  
  #create a vector of changes to k 
  k<-c(startK)# initiate vector with start value at k
  for (i in 1:length(year)-1){
    if (change[i+1]){
      changesetK<-c(changeK, -changeK)
      nextk<-k[i]*(100+(sample(changesetK, 1)))/100 #randomly chose an increase or decrease in % change
    } else{
      nextk<-k[i] # or if it's not a break year, don't change k
    }
    k<-c(k, nextk)
  }

  # #create a vector of changes to r
  r<-c(startR)# initiate vector with start value at r
  for (i in 1:length(year)-1){
    if (change[i+1]){
      changesetR<-c(changeR, -changeR)
      nextr<-r[i]*(100+(sample(changesetR, 1)))/100 #randomly chose an increase or decrease in % change
    } else{
      nextr<-r[i] # or if it's not a break year, don't change r
    }
    r<-c(r, nextr)
  }
  #calculate Nt vector
  Nt<-c(startPop) #create population vector with starting population as entry 1
  for(i in 1:length(year)){
    Nt1<-Nt[i]*exp(r[i]*(1- Nt[i]/k[i]))*noisevector[i]
    Nt<-c(Nt, Nt1)

  }
  #now we need to make the simulated data into a data frame which would look like
  #one fed into the analysis
  addyear<-max(year)+1
  year<-c(year, addyear)
  simdata<-as.data.frame(cbind(year, Nt))
  
  return(simdata)
}
fakedata(noise=5, changeK=25, changeR=25, breaks=list("1905", "1910"))

#now we need to create a function that will take the simulated data, find the best break combination
#and compare the ones it finds to the ones the data was built with

'%ni%' <- Negate('%in%')


weight.fake.shifts<-function(startyear, Nyears, startPop, noise, startK, 
                             startR, breaks, changeK, changeR, criterion){
  #create simulated data based on input parameters
  test<-fakedata(startyear, Nyears, startPop, noise, startK, startR, breaks, changeK, changeR)
  endbreak<-startyear+Nyears-1 #add the end break to the break list so it's more comparable to output
  breaksin<-c(unlist(breaks), endbreak) #and make it a vector
  nbreaksin<-length(breaks)  
  output<-breakweights(addNt1(test), criterion)
  options(warn=-1) #turn off warnings, we'll handle special cases directly
  rightbreaks<-output[which(output$breaksfound %in% breaksin),]
  wrongbreaks<-output[which(output$breaksfound %ni% breaksin),]
  wrongweight<-mean(as.numeric(wrongbreaks$correctedweights)) #mean weight of incorrect breaks
  wrongmax<-max(wrongbreaks$correctedweights)# maximum weight of incorrect break
  #if there are no wrong breaks, we need to set a zero weight 
  if (is.nan(wrongweight)){
    wrongweight<-0
    wrongmax<-0
  }
  #right breaks need two cases- for no breaks and any breaks scenarios
  if(nbreaksin>0){
    weights<-rightbreaks$correctedweights[1:(length(rightbreaks$correctedweights)-1)]
    rightweight<-mean(weights)
    rightmin<-min(weights) #minimum weight of correct break
    #mean weight of correct breaks minus end of series
  }else{
    rightweight<-1 #end of series break has a weight of 1 by definition
    rightmin<-1 #only break is at the end of the series
  }
  
  #output needed information
  testconditions<-unlist(c(Nyears, startPop, noise, nbreaksin, startK, 
                           startR, changeK, changeR, rightweight, wrongweight, rightmin, wrongmax))
  return(testconditions)
  
}


#create a function that compiles sucesses and failures for iterations of fitting the model
# on simulated data produced under given conditions

break.it.down.2<-function(startyear, Nyears, startPop, noise, 
                        startK, startR, breaks, changeK, changeR, nIter, criterion){
  out.frame<-data.frame(matrix(vector(), 0, 12,
                               dimnames=list(c(), 
                                             c("Nyears", "startPop", "noise", "nbreaksin",
                                               "startK", "startR", "changeK", "changeR", "rightweight", 
                                               "wrongweight", "rightmin", "wrongmax"))),
                        stringsAsFactors=FALSE)#Create a place to put our data
  for (i in 1:nIter){
    test<-weight.fake.shifts(startyear, Nyears, startPop, noise, startK, 
                             startR, breaks, changeK, changeR, criterion)
    out.frame<-rbind(out.frame, test)#put output for segment in a data frame
  }
  colnames(out.frame)<- c("Nyears", "startPop", "noise", "nbreaksin",
                          "startK", "startR", "changeK", "changeR", "rightweight", 
                          "wrongweight", "rightmin", "wrongmax")
  return(out.frame)
  
}



#okay, now that we've got it all working, it's time to build out the tests. To prevent the permutations
# of possible tests from going to infinity, let's create a 'base scenario' that we modify one parameter
# at a time, and let's choose 1,2,3,4 break point scenarios in which to test these

#choose base parameters

startyear<-1 #should not affect output at all
Nyears<-25 #processing time goes up considerably with length of time series, so make this the base scenario
startPop<-3000 # arbtrary start point, but r, K need to be chosen in reasonable scale with this
noise<-1 #base scenario should have very little %noise, but needs some so  there's a wee bit of error in the fit 
startK<-2000 #seems reasonable for a startpop of 1500
startR<-2 #also reasonable r
changeK<-50# start with big, easily detected shifts
changeR<-0 # as with changeK
nIter<-5 # keep this low while we build the code

# create some script that randomly chooses the breaks, given certain rules
# recall that the model assumes breaks cannot occur less than three years apart 
# or from the start or end of the time series because of overfitting issues


#create a function that generates a list of breaks randomly from the available set of breaks

breaklist<-function(possibleBreaks, howmany){ #we'll cap it at 3 breaks for the simulations
  if (howmany>3){ #no cheating, we're capping this at 3 breaks for the simulations
    howmany<-3
  }
  if (howmany<1){ #seriously, don't try to break this here
    howmany<-1
  }
  firstbreak<-sample(possibleBreaks, 1)
  eliminatedSecondBreaks<-seq(firstbreak-3, firstbreak+3)
  possibleSecondBreaks<-possibleBreaks[!is.element(possibleBreaks, eliminatedSecondBreaks)]
  secondbreak<-tryCatch(sample(possibleSecondBreaks, 1), error=function(e) NULL)
  eliminatedThirdBreaks<-tryCatch(seq(secondbreak-3, secondbreak+3), error=function(e) NULL)
  possibleThirdBreaks<-possibleSecondBreaks[!is.element(possibleSecondBreaks, eliminatedThirdBreaks)]
  thirdbreak<-tryCatch(sample(possibleThirdBreaks, 1), error=function(e) NULL)

  if (howmany==1){
    #for one break, this is simple
    breaks=sample(possibleBreaks, 1)
  }else if (howmany==2){
    #for two breaks
        breaks<-sort(c(firstbreak, secondbreak))
  }else if (howmany==3){
    #for three breaks, follow from 2
    breaks<-sort(c(firstbreak, secondbreak, thirdbreak))
  }
  return(breaks)
}


#create a function that uses break.it.down to test the function in four break point scenarios

iterate.breakitdown.2<-function(startyear, Nyears, startPop, 
                              noise, startK, startR, 
                              changeK, changeR, nIter, numLoops, criterion){
  #figure out possible breaks
  #minumum break must be four years in or later
  minbreak<-startyear+4
  #maximum break must be four years prior to the end of the series or before, plus we lose the last year
  maxbreak<-startyear+Nyears-5
  #create a sequence of all posible breaks
  possibleBreaks<-seq(minbreak, maxbreak)
  #Create a place to put our data
  results.matrix<-data.frame(matrix(vector(), 0, 12, 
                                    dimnames=list(c(), c("Nyears", "startPop", "noise", "nbreaksin",
                                                         "startK", "startR", "changeK", "changeR", "rightweight", 
                                                         "wrongweight", "rightmin", "wrongmax"))),
                             stringsAsFactors=F)
  
   while (numLoops>0){
    #we want to test each scenario with  0-3 breaks
    breaks0<-list() #empty list for no break scenario
    breaks1<-breaklist(possibleBreaks, 1)
    breaks2<-breaklist(possibleBreaks, 2)
    breaks3<-breaklist(possibleBreaks, 3)
    result.matrix0<-break.it.down.2(startyear=startyear, Nyears=Nyears, startPop=startPop, 
                                  noise=noise, startK=startK, startR=startR, 
                                  breaks=breaks0, changeK=changeK, changeR=changeR, nIter=nIter, criterion=criterion)
    result.matrix1<-break.it.down.2(startyear=startyear, Nyears=Nyears, startPop=startPop, 
                                  noise=noise, startK=startK, startR=startR, 
                                  breaks=breaks1, changeK=changeK, changeR=changeR, nIter=nIter, criterion=criterion)
    result.matrix2<-break.it.down.2(startyear=startyear, Nyears=Nyears, startPop=startPop, 
                                  noise=noise, startK=startK, startR=startR, 
                                  breaks=breaks2, changeK=changeK, changeR=changeR, nIter=nIter, criterion=criterion)
    result.matrix3<-break.it.down.2(startyear=startyear, Nyears=Nyears, startPop=startPop, 
                                  noise=noise, startK=startK, startR=startR, 
                                  breaks=breaks3, changeK=changeK, changeR=changeR, nIter=nIter, criterion=criterion)
    
    results.matrix<-rbind(results.matrix, result.matrix0, result.matrix1, result.matrix2, result.matrix3)
    numLoops<-numLoops-1
   }
  return(results.matrix)
}




##########################################################

#Okay, now we're ready to generate some data on how well the RS detector works
#rerun from this point and alter parts here to fiddle with simulations

#first, create a frame to put the data in as we change the scenarios
simulation.results<-data.frame(matrix(vector(), 0, 12, 
                                      dimnames=list(c(), c("Nyears", "startPop", "noise", "nbreaksin",
                                                           "startK", "startR", "changeK", "changeR", "rightweight", 
                                                           "wrongweight", "rightmin", "wrongmax"))),
                               stringsAsFactors=F)#Create a place to put our data
clearsims<-simulation.results
test.iter<-data.frame(matrix(vector(), 0, 12, 
                                      dimnames=list(c(), c("Nyears", "startPop", "noise", "nbreaksin",
                                                           "startK", "startR", "changeK", "changeR", "rightweight", 
                                                           "wrongweight", "rightmin", "wrongmax"))),
                               stringsAsFactors=F)#Create a place to put our data

#create base simulation
#we will be holding these values completely constant for comparisons' sake 
startyear<-1
startPop<-3000
nIter<-1
numLoops<-1
startK<-2000
criterion<-"AIC"

#we also want to keep track of how long this takes to run, so
# Start the clock!
ptm <- proc.time()

#things we want to vary
Nyearslist<-c(15,20,25,30)
noiselist<-c(1,2,5,10,15)
startRlist<-c(-0.5, 0.5, 1, 1.5, 2)
changeRlist<-c(0,10,25,50,75)
changeKlist<-c(0,10,25,50,75)

##############
#base scenario
#variables with = should be altered to see how results change
test.iter<-iterate.breakitdown.2(startyear=startyear, startPop=startPop,
                               Nyears=Nyearslist[2],
                               startK=startK, noise=noiselist[2], 
                               startR=startRlist[5], 
                               changeK=changeKlist[5], changeR=changeRlist[3], 
                               nIter, numLoops, criterion)


# Stop the clock
proc.time() - ptm


#okay, let's do this as one iteration on a complete set, repeated x times
####################################################################
### Start runnng here if it breaks

simnumber<-250
nIter<-1
numLoops<-1
simulation.results<-clearsims
criterion="AICc"


###### replace number before :simnuber with last sucessful sim number
for (f in 1:simnumber){
  ptm <- proc.time()
  
  #first number of years on base scenario
  for (i in 1:length(Nyearslist)){
    test.iter<-iterate.breakitdown.2(startyear=startyear, startPop=startPop,
                                   Nyears=Nyearslist[i],
                                   startK=startK, noise=noiselist[2], 
                                   startR=startRlist[5], 
                                   changeK=changeKlist[5], changeR=changeRlist[3], 
                                   nIter, numLoops, criterion)
    
    #add these results to the data frame
    simulation.results<-rbind(simulation.results, test.iter)
    writeLines(paste("finished", Nyearslist[i], " years"))
  }
  
  #### starting values of r
  
  for(q in 1:length(startRlist)){
    #next changeR on base scenario
    for (i in 1:length(changeRlist)){
      test.iter<-iterate.breakitdown.2(startyear=startyear, startPop=startPop,
                                     Nyears=Nyearslist[2],
                                     startK=startK, noise=noiselist[2], 
                                     startR=startRlist[q], 
                                     changeK=changeKlist[5], changeR=changeRlist[i], 
                                     nIter, numLoops, criterion)
      #add these results to the data frame
      writeLines(paste("finished changeR ", changeRlist[i]))
      simulation.results<-rbind(simulation.results, test.iter)
    }
    
    #next changeK on base scenario
    for (i in 1:length(changeKlist)){
      test.iter<-iterate.breakitdown.2(startyear=startyear, startPop=startPop,
                                     Nyears=Nyearslist[2],
                                     startK=startK, noise=noiselist[2], 
                                     startR=startRlist[q], 
                                     changeK=changeKlist[i], changeR=changeRlist[3], 
                                     nIter, numLoops, criterion)
      writeLines(paste("finished changeK ", changeKlist[i]))
      #add these results to the data frame
      simulation.results<-rbind(simulation.results, test.iter)
    }
    writeLines(paste("finished startR ", startRlist[q]))
    
  }
  
  
  #noise
  for (j in 1:length(noiselist)){
    #next changeR on base scenario
    for (i in 1:length(changeRlist)){
      test.iter<-iterate.breakitdown.2(startyear=startyear, startPop=startPop,
                                     Nyears=Nyearslist[2],
                                     startK=startK, noise=noiselist[j], 
                                     startR=startRlist[5], 
                                     changeK=changeKlist[5], changeR=changeRlist[i], 
                                     nIter, numLoops, criterion)
      #add these results to the data frame
      writeLines(paste("finished changeR ", changeRlist[i]))
      simulation.results<-rbind(simulation.results, test.iter)
    }
    
    #next changeK on base scenario
    for (i in 1:length(changeKlist)){
      test.iter<-iterate.breakitdown.2(startyear=startyear, startPop=startPop,
                                     Nyears=Nyearslist[2],
                                     startK=startK, noise=noiselist[j], 
                                     startR=startRlist[5], 
                                     changeK=changeKlist[i], changeR=changeRlist[3], 
                                     nIter, numLoops, criterion)
      writeLines(paste("finished changeK ", changeKlist[i]))
      #add these results to the data frame
      simulation.results<-rbind(simulation.results, test.iter)
    }
    
    writeLines(paste("finished noise ", noiselist[j]))
    #save the simulation results
    
  }
  
  write.csv(simulation.results, file=paste0("simresults/Break_weights_AICc/simresultsweightsAICc_", f,".csv"))
  simulation.results<-clearsims
  
  # Stop the clock
  proc.time() - ptm
  writeLines(paste(proc.time() - ptm))
}


