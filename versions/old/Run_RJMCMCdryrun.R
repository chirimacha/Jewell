## Set working directory to location of Jewell Code 
# setwd("home/data/Data/data")
setwd("/Users/snutman/Documents/CHAGAS_DATA/jewell/DryRun/")


#TODO FOR EACH CHAIN:
#set parameters
params <- list(
  banditarm=1,
  iterations=10,
  #beta is beta for each chain
  beta=0.1
)

#import libraries
library(doSNOW)
library(foreach)
library(doRNG)
library(lubridate)
library(PBSmapping)
library(plyr)
library(inline)
library(Rcpp)


#source Jewell MCMC code
source("RJMCMC.R")

#reset WD to data 
# setwd("/Users/snutman/Documents/CHAGAS_DATA/jewell/DryRun/")
setwd("/Users/snutman/Documents/CHAGAS_DATA/bandit/run_bandit/data/")

##set up timer
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
}

toc <- function()
{
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  toc - tic
}

#set seed
number <- sample(1:100000, 1)
print(paste("Seed: ",number,sep=""))
set.seed(number)

results<-run.mcmc(params$banditarm,params$iterations,params$beta)

write.csv(results, file=paste("Results","tiabaya",params$banditarm,"Seed",number,"Iterations",params$iterations,"Beta",params$beta,".csv",sep="_"),row.names=FALSE)

