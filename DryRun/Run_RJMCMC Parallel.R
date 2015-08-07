## Set working directory to location of Jewell Code/Data
setwd("home/data/Data/data")

#import libraries
library(doSNOW)
library(foreach)
library(doRNG)
library(lubridate)
library(PBSmapping)
library(plyr)
library(inline)
library(Rcpp)

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

#list packages to use with parallel processing
packages<-c("lubridate","PBSmapping","plyr","inline","Rcpp","doRNG")

#source Jewell MCMC code
source("RJMCMC.R")

#set parameters
params <- list(
  banditarm=1,
  iterations=10,
  #beta is beta for each chain
  beta=data.frame(0.1,0.3,0.6)
)

#make clusters

#TO DO: PUT IN PATH TO NAMES OF MACHINES
machines <-readLines("PATH")
machines <-rep(c("localhost", machines), each=1)
clust<-makeCluster(machines,type="SOCK")
registerdoSNOW(clust)
# 
# clu<-makeCluster(3)
# registerDoParallel(clu)

tic()
#set seed
number <- sample(1:100000, 1)
print(paste("Seed: ",number,sep=""))
set.seed(number)

  
#run bandit on each of the chains  
results.jewell <- foreach(i=params$beta,.packages=packages) %dorng% {
  run.mcmc(params$banditarm,params$iterations,i)
}
toc()

#turn off cluster
stopCluster(clust)

#Get results and average
v1 <-data.frame(results.jewell[1])
v2 <-data.frame(results.jewell[2])
v3 <-data.frame(results.jewell[3])

v1 <- v1[order(v1$unicode),]
v2 <- v2[order(v2$unicode),]
v3 <- v3[order(v3$unicode),]

test<-sum(v1$occult.prob,v2$occult.prob)
probmean <- ((as.numeric(as.character(v1$occult.prob)) + as.numeric(as.character(v2$occult.prob)) + as.numeric(as.character(v3$occult.prob)))/3)
v1$probmean <- probmean
v <- v1[order(v1$probmean, decreasing = TRUE),]
Ranking <- c(1:dim(v)[1])
averageRanking <- cbind(v, Ranking)
averageRanking$occult.prob <- NULL 

## New single csv should be in same folder as working directory 
allbetas=""
for (i in 1:length(params$beta)) {
  allbetas <-paste(allbetas,params$beta[i],sep="_")
}
write.csv(averageRanking, file=paste("AverageResults","tiabaya",params$banditarm,"seed",number,"iterations",params$iterations,"beta",allbetas,sep="_"))

##save results from each chain 
for (i in 1:length(results.jewell)) {
  write.csv(results.jewell[i],file=paste("Jewell_Results","tiabaya",params$banditarm,"seed",number,"iterations",params$iterations,"beta",params$beta[i],"chain",i,sep="_"))
}



