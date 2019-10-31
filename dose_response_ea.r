#---------------------------------------------------#
# Fits a dose-response curve using an evolutionary  #
# algorithm                                         #

#---------------------------------------------------#

## Simulation Controller ##
## Takes observations, x-axis values (e.g. concentrations), and EA parameters and returns a fit ##        
eadrm <- function(obs,xvals,model='h',ea.params=c(1000,200,20,10,100))
{         
  pop.size <- ea.params[1]
  stable.pop.size <- ea.params[2]
  num.tournaments <- ea.params[3]
  tournament.size <- ea.params[4]
  num.generations <- ea.params[5]

  ## Initial Hillslope Estimates ##
  Top <- min(obs)
  Bottom <- max(obs)
  num.conc <- length(xvals)
  EC50 <- EC50.start(obs,xvals)
  W <- 2.5
  f<-0.5
  SST <- sum((obs - mean(obs))^2)
  num.children <- (stable.pop.size / num.tournaments) - 1     ### default 200/20 - 1 = 9
  hs.params <- c(Top,Bottom,EC50,W)  ## the third hs.params is EC50!!!
  
  
  Top2<-max(obs)
  Bottom2<-max(obs)
  hs3.params <- c(Top2,EC50,W)
  hs5.params <- c(Top2,Bottom2,EC50,W,f)
  #print('Starting values')
  #print('Top, Bottom, EC50, W')
  #print(hs.params)

  
  ## Initital Exponential Estimates ##
  B = 1
  k = 2.5
  exp.params <- c(B,k)
  
  ## Initial Population ##
  population <- initial.population(pop.size,obs,xvals,hs.params,hs3.params, hs5.params,exp.params,model)     ### default pop size = 1000
 
  ## Variable used for convergence measuring ##
  fitness.history <- NULL
    ###  Being Generational Simulation ###
    allTimeBestFitness <- -Inf           #####*******************#####
    allTimeBestFit <- NULL               #####*******************#####                              
    for(generation in 1:num.generations)               ### default num.generations = 100
    {
      ###  Tournament Selection ###
      max.int  <- length(population)
      next.generation <- NULL
      tournament.winners <- list()
      numWinners <- 0
      for(tourn in 1:num.tournaments)                   ### default num.tournaments = 20
      {
        tournament <- list()
        max.fitness <- -Inf
        winner <- -Inf
        ### selects 'player' with highest r^2 from 10 random individuals from population
        for(player in 1:tournament.size)                ### default tournament.size = 10
        {
          player.index <- as.integer(runif(1)*max.int+1)
          current.player <- population[[player.index]]
          tournament[[player]] <- current.player
          current.fitness <- as.numeric(current.player[2])
          if(current.fitness > max.fitness) 
          { 
            max.fitness <- current.fitness
            winner <- player.index
          }          
        }
        tournament.winners[[numWinners <- numWinners+1]] <- population[[winner]]  ### list of 20 tournament winners
      }
      
      next.gen <- list()
      gen.count <- 0 
      c <- 0.2*((num.generations - generation)/num.generations)  ## Set mutation decary rate
      fitnesses <- NULL
      n <- length(tournament.winners)
      for(i in 1:n)
      {
      	fitnesses <- c(fitnesses, tournament.winners[[1]][2])
      }
      ranks <- n - rank(fitnesses) + 1
      for(i in 1:n)
      {
        current <- tournament.winners[[i]]
        ##Allow tournament winner to have children##
        child.count <- 1
        #while(child.count <= num.children)
        totalChildren <- 2*num.children*(ranks[i]/num.tournaments)
        while(child.count <= totalChildren)
        {
           current.params <- as.numeric(current[3:length(current)])
           child <- new.individual(current.params, c)
           pred.vals <- NULL
           if(current[1] == 'h') pred.vals <- predict.Hillslope.values(current.params,xvals)
           if(current[1] == 'h3') pred.vals <- predict.Hillslope3.values(current.params,xvals)
           if(current[1] == 'h5') pred.vals <- predict.Hillslope5.values(current.params,xvals)
           if(current[1] == 'e') pred.vals <- predict.exponential.values(current.params,xvals)
           child.fitness <- fitness(current[1],obs,pred.vals)
           next.gen[[gen.count <- gen.count+1]] <- c(current[1],child.fitness,child)
           child.count <- child.count + 1
         }  
         ### Add parent back in###
         next.gen[[gen.count <- gen.count+1]] <- current
      }
      
      a <- 0                                       ### 'a' is the sum of r^2 among tournament winners
      for(i in 1:length(tournament.winners))
      {
        a <- a + as.numeric(tournament.winners[[i]][2])
      }
      mean.fitness <-  a / length(tournament.winners)
      fitness.history <- rbind(fitness.history,mean.fitness)
      population <- next.gen
            
                  
     
      winner.max <- as.numeric(population[[1]][2])
      winner.row <- 1
      for(i in 2:length(population)){
      	thisFitness <- as.numeric(population[[i]][2])
        if(thisFitness > winner.max){
        	winner.max <- thisFitness
        	winner.row <- i
        }
      }


      winner <- population[[winner.row]]
      if(winner.max > allTimeBestFitness){
         allTimeBestFitness <- winner.max
         allTimeBestFit <- winner                                             
      }
       
    }  # end generation loop
    #print(population[[1]][2])
    
    params <- NULL
    if(winner[1] == "h")
    {
      params <- rbind(c(winner[3:6]))
      colnames(params) <- c("EMAX","EMIN","EC50","W")
    } else if(winner[1] == "h3")
      {
        params <- rbind(c(winner[3:5]))
        colnames(params) <- c("EMAX","EC50","W")
    }else if(winner[1] == "h5")
    {
      params <- rbind(c(winner[3:7]))
      colnames(params) <- c("EMAX","EMIN","EC50","W","f")
    }else if(winner[1] == "e")
    {
      params <- rbind(c(winner[3:4]))
      colnames(params) <- c("B","K")
    }
    
	return(list(Model=winner[1],R2=winner[2],params=params)) 
         
     ##### Down to here #########
     
}

## Creates an initial population of specified size seeded by initial paramter estimates ##
initial.population <- function(size,obs,xvals,hs.params,hs3.params, hs5.params,exp.params,model)
{
  population <- list()
  SST <- sum((obs - mean(obs))^2)
  count <- 0
  if(model == 'h')
  {  
    ##Generate initial population##
    temp.pop <- NULL
    for(i in 1:size)
    {
      individual <- new.initial.individual(hs.params)               ### an 'individual' is just a set of parameters
      predicted.vals <- predict.Hillslope.values(individual,xvals)
      R.squared <- fitness('h',obs,predicted.vals)
      new.member <- c('h',R.squared,individual)                     ### each member of population is a model, an r^2, and an 'individual'
      population[[count <- count+1]] <- new.member
    }     
  } else if(model == 'h3')
    {  
      ##Generate initial population##
      for(i in 1:size)
      {
        individual <- new.initial.individual(hs3.params)
        predicted.vals <- predict.Hillslope3.values(individual,xvals)
        R.squared <- fitness('h3',obs,predicted.vals)
        new.member <- c('h3',R.squared,individual)
        population[[count <- count+1]] <- new.member
      } 
    } else if(model == 'h5')
    {  
      ##Generate initial population##
      for(i in 1:size)
      {
        individual <- new.initial.individual(hs5.params)
        predicted.vals <- predict.Hillslope5.values(individual,xvals)
        R.squared <- fitness('h5',obs,predicted.vals)
        new.member <- c('h5',R.squared,individual)
        population[[count <- count+1]] <- new.member
      } 
    } else if(model == 'e')
    {  
      ##Generate initial population##
      for(i in 1:size)
      {
        individual <- new.initial.individual(exp.params)
        predicted.vals <- predict.exponential.values(individual,xvals)
        R.squared <- fitness('e',obs,predicted.vals)
        new.member <- c('e',R.squared,individual)
        population[[count <- count+1]] <- new.member
      } 
    } else if(model == 'b')
      {  
        ##Generate initial population##
        for(i in 1:(size/4))
        {
          individual <- new.initial.individual(hs.params)
          predicted.vals <- predict.Hillslope.values(individual,xvals)
          R.squared <- fitness('h',obs,predicted.vals)
          new.member <- c('h',R.squared,individual)
          population[[count <- count+1]] <- new.member
        } 
      for(i in 1:(size/4))
      {
        individual <- new.initial.individual(hs3.params)
        predicted.vals <- predict.Hillslope3.values(individual,xvals)
        R.squared <- fitness('h3',obs,predicted.vals)
        new.member <- c('h3',R.squared,individual)
        population[[count <- count+1]] <- new.member
      } 
      for(i in 1:(size/4))
      {
        individual <- new.initial.individual(hs5.params)
        predicted.vals <- predict.Hillslope5.values(individual,xvals)
        R.squared <- fitness('h5',obs,predicted.vals)
        new.member <- c('h5',R.squared,individual)
        population[[count <- count+1]] <- new.member
      } 
        for(i in 1:(size/4))
        {
          individual <- new.initial.individual(exp.params)
          predicted.vals <- predict.exponential.values(individual,xvals)
          R.squared <- fitness('e',obs,predicted.vals)
          new.member <- c('e',R.squared,individual)
          population[[count <- count+1]] <- new.member
        } 
      }         
  return(population)
}

## Individual in initial population. Allowed to Mutate up to 100% from parental seed ##
new.initial.individual <- function(params)
{
    newParams <- NULL
    for(i in 1:length(params))
    {
      ran <- as.integer(runif(1)*100+1)
      if(ran%%2 == 1) ran <- -ran
      mutated <- params[i] + (params[i] * (ran/100)) 
      newParams <- c(newParams, mutated)
    }

    ## make sign of exponential parameter negative to fit inhibitory responses ## 
    ran <- as.integer(runif(1)*100+1)
    if(ran%%2 == 1) newParams[length(newParams)] <- newParams[length(newParams)] * -1        
    return(newParams)
}

#### Begin utility funtions ####

## Predict values based on HillSlope model and xvals ##
predict.Hillslope.values <- function(params,xvals)
{
  predicted.yvals <- NULL
  for(x in xvals)
  {
    y.val <- params[1] - (params[1]-params[2])/(1+(x/params[3])^(params[4]))
    predicted.yvals <- c(predicted.yvals,y.val)
  }
  return(predicted.yvals)
}

##Predict values based on 3-p logistic model
predict.Hillslope3.values <- function(params,xvals)
{
  predicted.yvals <- NULL
  for(x in xvals)
  {
    y.val <- (params[1])/(1+exp(params[3]*(log(x)-log(params[2]))))
    predicted.yvals <- c(predicted.yvals,y.val)
  }
  return(predicted.yvals)
}

##Predict values based on 5-p logistic model

predict.Hillslope5.values <- function(params,xvals)
{
  predicted.yvals <- NULL
  for(x in xvals)
  {
    y.val <- params[2]+(params[1]-params[2])/((1+exp(params[4]*(log(x)-log(params[3]))))^params[5])
    predicted.yvals <- c(predicted.yvals,y.val)
  }
  return(predicted.yvals)
}


## Predict values based on exponential model and xvals ##
predict.exponential.values <- function(params,xvals)
{
  predicted.yvals <- NULL
  for(x in xvals)
  {
    y.val <- params[1]*exp(params[2]*x)
    predicted.yvals <- c(predicted.yvals,y.val)
  }
  return(predicted.yvals)
}


new.individual <- function(params, c)
{
	newParams <- params*(1 + runif(length(params), -c, c))
    return(newParams)
}


## Calculates fitness of indivdual ##
fitness <- function(model,obs,predicted.vals)
{
  n1<-length(obs)
  if(model == 'h') aiccminus<--(n1*log(sum( ((obs) - (predicted.vals))^2 )/(n1))+4*log(n1))
  if(model == 'h3') aiccminus<--(n1*log(sum( ((obs) - (predicted.vals))^2 )/(n1))+3*log(n1))
  if(model == 'h5') aiccminus<--(n1*log(sum( ((obs) - (predicted.vals))^2 )/(n1))+5*log(n1))
  if(model == 'e') aiccminus<--(n1*log(sum( ((obs) - (predicted.vals))^2 )/(n1))+2*log(n1))
  return(aiccminus)
}


## Get start value for EC50, which is the xval in the "middle" ##
EC50.start <- function(obs,xvals)
{   
  if(length(xvals) %% 2 == 0)
  {
      first <- length(xvals) / 2
      last  <- first + 1
      EC50  <- 10^( (log10(xvals[first]) + log10(xvals[last]))/2 )
  } else {
      EC50 <- xvals[round(length(xvals) / 2)] 
  }
 
  return(EC50)
}
