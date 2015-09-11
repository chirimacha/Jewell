setwd("~/Jewell/JSM")
#set seed
set.seed(1234)

#load libraries
library("lubridate")
library("PBSmapping")
library("plyr")
library("inline")
library("Rcpp")

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

###########################
###read data##############
#############################

tic()

#set working drive
setwd("/Users/EMWB/Jewell/Data")
#setwd("~/Desktop/Levy Lab")
#setwd("~/Users/e/Jewell/Data")
#setwd("/Users/mzlevy/Jewell/Data")



##################################################
#######generate data##############################
###data simulation###
#NOTE: this simulation assumes Markov properties
#we can later extend to non-Markov chains
#simulation statistic vectors initialized


#for (s in 1:S.sim){
totalnuminf=1
N=200
indicator=indicator2=0
while(sum(indicator)<6){
  
  id<-c(1:N) #simulate 100 houses
  
  #locations uniform on [0,1]x[0,1] grid
  location<-matrix(c(NA,NA),nrow=max(id),ncol=2) 
  for (i in 1:max(id)) {location[i,]<-c(runif(1),runif(1))}
  #plot(location[,1],location[,2])
  
  #calculate distances between houses
  distance<-matrix(NA,nrow=max(id),ncol=max(id))
  for (i in 1:max(id)){
    for (j in 1:max(id)){
      distance[i,j]=sqrt((location[i,1]-location[j,1])^2+(location[i,2]-location[j,2])^2)
    }
  }
  
  #simulate epidemic
  beta=truebeta=.2
  Rb=trueRb=1.09
  K=1000 #carrying capacity
  maxt=100
  tobs=99
  
  beverton.holt<-function(id,K,R,bugs,trueremovaltime,trueinfectiontime){
    for(t in trueinfectiontime:(trueremovaltime-1)){
      bugs[id,(t+1)]=ceiling(R*bugs[id,t]/(1+bugs[id,t]/(K/(R-1))))
    }
    return(bugs[id,])
  }
  
  S=I=R=matrix(NA, nrow=max(id),ncol=N*N) 
  bugs=bugs.Rb=matrix(0,nrow=N,ncol=maxt)
  
  #probability of infestation differs by hops (<.3) or jumps (>.3)
  threshold=ifelse(distance<.2,1,.01)
  
  #initial state vectors
  S[,1]<-rep(1,N)
  I[,1]<-rep(0,N)
  R[,1]<-rep(0,N)
  
  #initial infective
  initialinfective<-sample(1:N,1)
  S[initialinfective,1]=0
  I[initialinfective,1]=1
  
  
  #keeping counts of S, E, and I at all the time points of the simulation
  StoI=ItoR=rep(0,N)
  infectiontime=rep(0,N)
  removaltime=rep(0,N)
  
  #probability that i infects j
  h1<-function(t,r,I,i,j,beta){
    #n=bugs[i,t]
    n=1
    #N.prime=ifelse(t-I[i]>0,Rb^(t-I[i]),0)
    deriv <- r^(t-I[i])*log(r)-(K-1)*K*r^(t-I[i])*log(r)/(K+r^(t-I[i])-1)^2
    #deriv=K*n^2*log(r)/((K-n)*exp(-log(r)*(t-I[i]))+n)^2
    hazard=ifelse(deriv>0,1-(1-beta/N*threshold[i,j])^(deriv),0)
    return(hazard) 
  }
  
  #find initial infectives notification and recovery times
  infectiontime[initialinfective]<-1
  bugs[initialinfective,1]<-1
  removaltime[initialinfective]<-maxt
  bugs[initialinfective,]=rpois(maxt,beverton.holt(initialinfective,K,Rb,bugs,maxt,infectiontime[initialinfective]))
  bugs[initialinfective,]=ifelse(bugs[initialinfective,]==0,1,bugs[initialinfective,])
  bug.ind=matrix(NA,nrow=N,ncol=N)
  for(t in 2:(maxt-1)){   #looping through each day of the simulation
    for(i in 1:N){
      for(j in 1:N){
        bug.ind[i,j]=1-ifelse(t>=infectiontime[i]&infectiontime[i]>0,h1(t,Rb,infectiontime,i,j,beta),0)
      }}
    for(j in 1:N){
      ifelse(S[j,(t-1)]==1, StoI[j]<-rbinom(1,S[j,(t-1)],1-prod(bug.ind[,j])),StoI[j]<-0)
      if(StoI[j]==1){
        infectiontime[j]<-t
        bugs[j,t]=1
        removaltime[j]<-maxt
        bugstemp=beverton.holt(j,K,Rb,bugs,maxt,infectiontime[j])
        bugs[j,(t:maxt)]=rpois((maxt-t+1),bugstemp[t:maxt])
      }
    } 
    
    #Now update the daily counts
    S[,t]=S[,(t-1)]-StoI
    I[,t]=I[,(t-1)]+StoI-ItoR
    R[,t]=R[,(t-1)]
  }  
  
  
  S=S[,1:(maxt-1)]
  I=I[,1:(maxt-1)]
  R=R[,1:(maxt-1)]
  
  
  #look at distribution at a few times
  par(mfrow=c(2,2))
  t=1
  plot(location[,1],location[,2],col="gray",main = substitute(paste(Time, " = ", t),list(t=t)),xlab="X",ylab="Y")
  for (i in 1:N) if(I[i,t]==1) points(location[i,1],location[i,2],pch=19,col="red")
  t=round(maxt/3)
  plot(location[,1],location[,2],col="gray",xlab="X",ylab="Y",main = substitute(paste(Time, " = ", t),list(t=t)))
  for (i in 1:N) if(I[i,t]==1) points(location[i,1],location[i,2],pch=19,col="red")
  t=tobs
  plot(location[,1],location[,2],col="gray",xlab="X",ylab="Y",main = substitute(paste(Time, " = ", t),list(t=t)))
  for (i in 1:N) if(I[i,t]==1) points(location[i,1],location[i,2],pch=19,col="red")
  t=maxt-1
  plot(location[,1],location[,2],col="gray",xlab="X",ylab="Y",main = substitute(paste(Time, " = ", t),list(t=t)))
  for (i in 1:N) if(I[i,t]==1) points(location[i,1],location[i,2],pch=19,col="red",xlab="X",ylab="Y",main = substitute(paste(Time, " = ", t),list(t=t)))
  
  
  #count total number infected
  totalnuminf=sum(S[,1])-sum(S[,maxt-1])
  
  check3=ifelse(infectiontime<=tobs&infectiontime!=0,bugs[,tobs], Inf)
  indicator=length(check3[which(check3>0&check3<Inf)])
  indicator2=length(infectiontime[which(infectiontime>tobs)])
}
#make dataset
data=cbind(infectiontime, removaltime,bugs)
truebugs=bugs


#delete a few observations
#update data to reflect only OBSERVED data
#temp<-which(infectiontime!=0&infectiontime!=1)
#delete.number=floor(1/3*length(temp))
#true.occult<-sample(temp,2)
#infectiontime[true.occult]<-0
#bugs[true.occult,]=0
predprobs <- rep((totalnuminf+1)/N,N)

########################################################
#########MCMC algorithm###################################
#######################################################
tic()   #begin timer


##Jewell MCMC
M <- 100000 #length of simulation
m <- 1 #first iteration

tobs <- rep(maxt,N)

#define number of houses
occult.sum.new<-rep(0,N)
inspected <- rep(0,N)
sum.insp.temp <- bugs[,99]

#determine random observation times
for (i in 1:length(sum.insp.temp)){
  if(sum.insp.temp[i]>0){
    tobs[i] <- round(runif(1,min = infectiontime[i], max = maxt))
    inspected[i]=1
  }
}

sum.insp <- rep(0,N)
for(i in 1:N) sum.insp[i] <- bugs[i,tobs[i]]

infectiontime<-rep(Inf,N)
check3<- rep(Inf,N) #initialize data vector
T_b <- 0.2 #threshold for bug infectiousness
jumpprob <- .02 #probability of jump vs. hop
bugs <- matrix(0,nrow=N,ncol=maxt) #initialize but matrix
maxbugs <- max(sum.insp) #find most observed bugs in data
initialinfective <- which(sum.insp==maxbugs) #set this house as initialinfective
id=1:N #generate ids
K=1000 #carrying capacity
tuning <- 0.01 #tuning parameter for RJ


check3<-ifelse(sum.insp>0,sum.insp,Inf) #replace with observed bug counts
I=ifelse(check3!=Inf&check3!=Inf,50,Inf) #set initial values for infection times

trueremovaltime=ifelse(check3<Inf,maxt+1,Inf) #set recovery times 
detectiontime=ifelse(check3>0&check3<Inf,tobs,Inf) #set detection time vector

#initialize parameters 
beta=Rb=rep(0,M)
beta[1]=.6
betastar=.6
betastar.I=matrix(0,nrow=N,ncol=N)
betastar.sum=rep(0,N)
S=H.mat=matrix(0,nrow=N,ncol=N)
U=rep(0,N)
accept.beta <- rep(0,M)
accept.Rb <- rep(0,M)
accept.I <- rep(0,M)
thresholdsum <- apply(threshold,1,sum)

#keep track of occult infestations
occult=rep(0,N)

#initialize infection times
I[initialinfective]=1
Istar=I

#define indicators for infected houses
N_N=which(check3!=Inf)
N_I=which(I!=Inf)
#take out initial infective
N_N=N_N[which(N_N!=initialinfective)]
N_I=N_I[which(N_I!=initialinfective)]
infectedhouses=rep(0,N)
infectedhouses[N_N]=1
infectedhousesI=rep(0,N)
infectedhousesI[N_I]=1

#intialize bug mean vector
lambda_t=rep(0,N)
lambda_t[1]=1

#keep track of rankings
rankings1 <- rankings2 <- matrix(0,nrow=N,ncol=M/100)

#initialize Rb parameter
Rb[1]=1.2
############################################
###############functions####################


`%notin%` <- function(x,y) !(x %in% y) 

sampleWithoutSurprises <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } else {
    return(sample(x,1))
  }
}


#BH function to update bug counts given matrix
beverton.holt<-function(id,K,R,bugs,maxt,trueinfectiontime){
  for(t in trueinfectiontime:(maxt-1)){
    bugs[id,(t+1)]=ceiling(R*bugs[id,t]/(1+bugs[id,t]/(K/(R-1))))
  }
  return(bugs[id,])
}

#BH function to update infection times
beverton.holt.I<-function(update,K,R,check3,tobs,bugs){
  tobs1=tobs
  time=ifelse(update %in% N_N, -log((K-check3[update])/(check3[update]*K-check3[update]))/log(R),-log((K-bugs[update,tobs1])/(bugs[update,tobs1]*K-bugs[update,tobs1]))/log(R))
  I[update]=tobs1-time
  return(I[update])
}


#BH function to udpate bug counts given vector
beverton.holt.update<-function(K,R,bugs,trueremovaltime,trueinfectiontime){
  for(t in trueinfectiontime:(min(maxt-1,trueremovaltime))){
    bugs[(t+1)]=ceiling(R*bugs[t]/(1+bugs[t]/(K/(R-1))))
  }
  return(bugs)
}


#poisson piece of likelihood for matrix
f_D<-function(i,bugs,I,check3,Rb){
  den=0
  lambda_t=beverton.holt.update(K,Rb,lambda_t,maxt,1)
  if(check3[i]<Inf){
    adjusted.bugs=bugs[i,which(bugs[i,]!=0)]
    den=den+dpois(adjusted.bugs,lambda_t[1:length(adjusted.bugs)],log=TRUE)
  }
  
  return(den)
}

#poisson piece of likelihood for vector
f_D.update<-function(i,bugs,I,check3,Rb){
  den=0
  lambda_t=beverton.holt.update(K,Rb,lambda_t,maxt,1)
  if(check3[i]<Inf){
    adjusted.bugs=bugsstar[which(bugsstar!=0)]
    den=den+dpois(adjusted.bugs,lambda_t[1:length(adjusted.bugs)],log=TRUE)
  }
  
  return(den)
}


#' First piece of likelihood
#'
#' @param I Numeric vector of times.
#' @param beta Float hazard rate for infection.
#' @param initialinfective Cardinal of initial infected house.
#' @param r Growth rate, a float.
#' @param K Distance factor?, a float.
#' @param N Number of houses.
#' @param N_I Those houses that were surveiled? A list of indices.
#' @param threshold NxN matrix of float cutoff distances.
#' @return Array of likelihoods, as numeric vector.

firstpiece <- function(I, beta, initialinfective, r, K, N, N_I, threshold) {
  H.mat <- matrix(1, nrow=N, ncol=N)
  beta.sum <- rep(1, N)
  for (j in 1:N) {
    for (i in 1:N) {
      if (i %in% N_I | i==initialinfective) {
        if (I[i]<I[j] & I[j]<Inf) {
          t<-I[j]-I[i]
          deriv <- r^t*log(r)-(K-1)*K*r^t*log(r)/(K+r^t-1)^2
          H.mat[i,j] <- ifelse(deriv>0, (1-beta*threshold[i,j])^(deriv), 1)
        }
      }
    }
    beta.sum[j] <- 1- prod(H.mat[,j])
  }
  ifelse(is.na(beta.sum), 0, beta.sum)
}

first.include <- '
#include <set>
#include <cmath>
'

firstpiece.wrap <- cxxfunction(signature(IS="numeric", betaS="float",initialinfectiveS="int", rS="float", KS="float", NS="int", N_IS="numeric",thresholdS="numeric"), plugin="Rcpp", incl=first.include, body='
                               Rcpp::NumericVector I(IS);
                               double beta=Rcpp::as<double>(betaS);
                               int initialinfective=Rcpp::as<int>(initialinfectiveS);
                               double r=Rcpp::as<double>(rS);
                               double K=Rcpp::as<double>(KS);
                               int N=Rcpp::as<int>(NS);
                               Rcpp::NumericVector N_I(N_IS);
                               Rcpp::NumericVector threshold(thresholdS);
                               std::set<int> infecteds(N_I.begin(), N_I.end());
                               typedef std::set<int>::const_iterator infiter;
                               infecteds.insert(initialinfective);
                               Rcpp::NumericVector H(N, 1.0);
                               for (int j=0; j<N; ++j) {
                               double total=1;
                               for (infiter it=infecteds.begin(); it!=infecteds.end(); ++it) {
                               int i=*it - 1; // -1 to convert to 0-based indexing.
                               if (I[i]<I[j] && std::isfinite(I[j])) {
                               double deriv=std::pow(r,(I[j]-I[i]))*std::log(r)-(K-1)*K*std::pow(r,(I[j]-I[i]))*std::log(r)/std::pow((K+std::pow(r,(I[j]-I[i]))-1),2);
                               if (deriv>0) {
                               double add=std::pow(1-beta*threshold[i*N+j], deriv);
                               //std::cout << "i "<<i<<" j " << j<<" I i "<<I[i]<<" I[j] "<< I[j]
                               //<<" deriv " << deriv << " hmat " << add <<std::endl;
                               total*=add;
                               }
                               }
                               }
                               H[j]=1-total;
                               }
                               return H;
                               ')

H <- function(i, j, t, r, I, beta, K, threshold) {
  h <- rep(0,t)
  for(a in 1:t){
    deriv <- r^a*log(r)-(K-1)*K*r^a*log(r)/(K+r^a-1)^2
    h[a] <- ifelse(deriv>0, 1-(1-beta*threshold[i,j])^(deriv), 0)
  }
  sum(h)
}

secondpiece <- function(I, trueremovaltime, beta, r, K, N, maxt, threshold) {
  S1 <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    if (I[i]!=Inf) {
      for (j in 1:N) {
        t <- min(maxt, I[j],trueremovaltime[i]) - min(I[i], I[j])
        if (t>0) {
          S1[i,j] <- H(i, j, t, r, I, beta, K, threshold)
        }
      }
    }
  }
  S1<-ifelse(S1=="NaN", 0, S1)  
  sum(S1)
}


#' Second piece of likelihood
#'
#' @param I Numeric vector of times.
#' @param beta Float hazard rate for infection.
#' @param r Growth rate, a float.
#' @param K Distance factor?, a float.
#' @param N Number of houses.
#' @param maxt A maximum time, a float.
#' @param threshold NxN matrix of float cutoff distances.
#' @return Single sum, a float.

secondpiece.wrap <- function(I, beta, r, K, N, maxt, threshold,thresholdsum) {
  S1 <- matrix(0, nrow=N, ncol=N)
  for (i in 1:N) {
    if (I[i]!=Inf) {
      for (j in 1:N) {
        t <- min(maxt, I[j]) - min(I[i], I[j])
        if (t>0) {
          S1[i,j] <- H(i, j, t, r, I, beta, K, threshold)
        }
      }
    }
  }
  S1<-ifelse(S1=="NaN", 0, S1)  
  sum(S1)/N
}

secondpiece.wrap <- cxxfunction(signature(IS="numeric", trueremovaltimeS="numeric",betaS="float",rS="float", KS="float", NS="int", maxtS="float",thresholdS="numeric",thresholdsumS="numeric"), plugin="Rcpp", incl=first.include, body='
                                Rcpp::NumericVector I(IS);
                                Rcpp::NumericVector trueremovaltime(trueremovaltimeS);
                                double beta=Rcpp::as<double>(betaS);
                                double r=Rcpp::as<double>(rS);
                                double K=Rcpp::as<double>(KS);
                                int N=Rcpp::as<int>(NS);
                                double maxt=Rcpp::as<double>(maxtS);
                                Rcpp::NumericVector threshold(thresholdS);
                                Rcpp::NumericVector thresholdsum(thresholdsumS);
                                
                                double total=0;
                                double add=0;
                                for (int i=0; i<N; ++i) {
                                if (std::isfinite(I[i])) {
                                for (int j=0; j<N; ++j) {
                                double t=std::min(std::min(trueremovaltime[i],maxt), I[j]) - std::min(I[i], I[j]);
                                if (t>0) {
                                for (int a=1; a<=t;++a){
                                double deriv=std::pow(r,a)*std::log(r)-(K-1)*K*std::pow(r,a)*std::log(r)/std::pow((K+std::pow(r,a)-1),2);
                                double result=1-std::pow(1-beta*threshold[i*N+j], deriv);
                                add+= result/N;
                                if (std::isfinite(add)) {
                                total+=add;
                                }
                                }
                                }
                                }
                                }
                                }
                                return Rcpp::wrap(add);
                                ')

##############################
#########for loop begins#######
##############################
bugsize=NULL
infectiontime=I
id = 1:N

#find initial infectives notification and recovery times
infectiontime[initialinfective]<-1
bugs[initialinfective,1]<-1
bugs[initialinfective,]=rpois(maxt,beverton.holt(initialinfective,K,Rb[1],bugs,maxt,infectiontime[initialinfective]))

#initialize bug counts
for (i in which(I!=Inf)){
  bugs[i,I[i]]=1
  bugstemp=beverton.holt(i,K[1],Rb[1],bugs,maxt,I[i])
  bugs[i,(I[i]:min(trueremovaltime[i],maxt))]=rpois((min(maxt,trueremovaltime[i])-I[i]+1),bugstemp[I[i]:min(maxt,trueremovaltime[i])])
  bugs[i,tobs[i]]=ifelse(check3[i]<Inf,check3[i],bugs[i,tobs[i]])
}

for (m in 2:M){
   
  ################################
  ######update Rb##################
  ####################################
  
  Rbstar=rnorm(1,Rb[m-1],.02)
  if(Rbstar<1) Rbstar=1+(1-Rbstar)
  Q=Qstar=rep(NA,length(I[which(I!=Inf)]))
  
  for (i in which(I!=Inf)){
    Qstar[which(I!=Inf)==i]=f_D(i,bugs,Istar,check3,Rbstar)[length(f_D(i,bugs,Istar,check3,Rbstar))]
    Q[which(I!=Inf)==i]=f_D(i,bugs,I,check3,Rb[m-1])[length(f_D(i,bugs,I,check3,Rb[m-1]))]
  }
  thirdpieceloglike <- sum(Qstar[which(Qstar!="NA")])-sum(Q[which(Q!="NA")])
  
  logfirstpieceIstar <- log(firstpiece.wrap(Istar, beta[m-1], initialinfective, Rbstar, K, N, N_I, threshold))
  logfirstpieceIstar <- ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
  loglikestar <- sum(logfirstpieceIstar)-secondpiece.wrap(I, trueremovaltime,beta[m-1], Rbstar, K, N, maxt, threshold,thresholdsum)
  
  logfirstpiece<-log(firstpiece.wrap(I, beta[m-1], initialinfective, Rb[m-1], K, N, N_I, threshold))
  logfirstpiece <- ifelse(logfirstpiece=="-Inf",0,logfirstpiece)
  loglike<- sum(logfirstpiece)-secondpiece.wrap(I, trueremovaltime,beta[m-1], Rb[m-1], K, N, maxt, threshold,thresholdsum)
  
  priors <- dgamma(Rbstar,shape=trueRb,scale=1,log=TRUE)-dgamma(Rb[m-1],shape=trueRb,scale=1,log=TRUE)
  
  Rbloglike <- loglikestar+thirdpieceloglike-loglike+priors
  #Metropolis step
  mstep.Rb=min(0,Rbloglike)
  #if(mstep.Rb=="NaN") mstep.Rb=0
  R=log(runif(1))
  if(R<mstep.Rb){
    Rb[m]=Rbstar
    loglike<-loglikestar
    Q<-Qstar
    accept.Rb[m] <- 1
  }else{
    Rb[m]=Rb[m-1]
    accept.Rb[m] <- 0}
  
  ###############
  ##update beta##
  ###############
  a.beta <- 10*beta[m-1]
  b.beta <- 10-beta[m-1]*10
  betastar=rbeta(1 , a.beta , b.beta)
  #if(betastar>1){betastar <- 1-(betastar-1)}
  logfirstpiecestar<-log(firstpiece.wrap(I, betastar, initialinfective, Rb[m], K, N, N_I, threshold))
  logfirstpiecestar=ifelse(logfirstpiecestar=="-Inf",0,logfirstpiecestar)
  loglikestar=sum(logfirstpiecestar)-secondpiece.wrap(I, trueremovaltime, betastar, Rb[m], K, N, maxt, threshold,thresholdsum)
  betapriors <- dbeta(betastar,shape1=truebeta*10,shape2=10-truebeta*10,log=TRUE)-dbeta(beta[m-1],shape1=truebeta*10,shape2=10-truebeta*10,log=TRUE)
  mstep.beta=min(1,exp(loglikestar-loglike+betapriors-dbeta(betastar, a.beta, b.beta, log=TRUE)+dbeta(beta[m-1], a.beta, b.beta,log=TRUE)))
  if(mstep.beta=="NaN") mstep.beta=1
  R=runif(1)
  if(R<mstep.beta){
    beta[m]=betastar
    loglike <- loglikestar
    accept.beta[m]=1
  }else{
    beta[m]=beta[m-1]
    accept.beta[m]=0
  }
  
  
  ########################################################
  #####decide whether to update I, add I, or delete I#####
  ########################################################
  add.del.move<-sample(c("add","del","move"),1)
  
  if(add.del.move=="move"){
  ############
  ##update I##
  ############
  
  ##pick a house to update the time out infected houses
  update=sample(N_I,1,replace=TRUE)
  if(bugs[update,min(maxt,trueremovaltime[update])]==0) bugs[update,min(maxt,trueremovaltime[update])]=1
  Istar[update] <- sample(c(2:(maxt-2)),1,replace=TRUE)
  bugsstar=rep(0,maxt)
  bugsstar[Istar[update]]=1
  bugsstar=beverton.holt.update(K,Rb[m],bugsstar,maxt,Istar[update])
  bugsstar[tobs[update]]=check3[update]
  logfirstpieceI<-log(firstpiece.wrap(I, beta[m], initialinfective, Rb[m], K, N, N_I, threshold))
  logfirstpieceI=ifelse(logfirstpieceI=="-Inf",0,logfirstpieceI)
  loglike=sum(logfirstpieceI)-secondpiece.wrap(I, trueremovaltime,beta[m], Rb[m], K, N, maxt, threshold,thresholdsum)
  logfirstpieceIstar<-log(firstpiece.wrap(Istar, beta[m], initialinfective, Rb[m], K, N, N_I, threshold))
  logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
  loglikestar=sum(logfirstpieceIstar)-secondpiece.wrap(Istar, trueremovaltime,beta[m], Rb[m], K, N, maxt, threshold,thresholdsum)
  Q=sum(f_D.update(update,bugsstar,Istar,check3,Rb[m]))-sum(f_D(update,bugs,I,check3,Rb[m]))
  
  #Metroplis step; decide whether to accept new time#
  mstep.I=min(1,exp(loglikestar-loglike+Q))
  if(mstep.I=="NaN") mstep.I=1
  
  R=runif(1)
  if(R<mstep.I){
    I<-Istar
    bugs[update,]=bugsstar
    loglike<-loglikestar
    accept.I[m] <- 1
  }else{
    Istar<-I
    accept.I[m] <- 0
  } 
  
  }else if(add.del.move=="add"){ 
    
    ###########
    ###add I###
    ##########
    
    addinf<-which(I==Inf)
    if(length(addinf)>1){
      update=sample(addinf,1)
      Istar[update]=floor(runif(1,min=2,max=maxt-2))
      trueremovaltime[update]=maxt+1
      detectiontime[update]=maxt
      tobs[update]=maxt
      bugsstar=rep(0,maxt)
      bugsstar[Istar[update]]=1
      bugsstar=beverton.holt.update(K,Rb[m],bugsstar,trueremovaltime[update],Istar[update])
      bugstest=bugs
      bugstest[update,]=bugsstar[1:maxt]
      check3[update]<-ifelse(inspected[update]==1,0,bugsstar[tobs])
      check3[update]<-ifelse(check3[update]>K,K,check3[update])
      logfirstpieceI<-log(firstpiece.wrap(I, beta[m], initialinfective, Rb[m], K, N, N_I, threshold))
      logfirstpieceI <- ifelse(logfirstpieceI=="-Inf",0,logfirstpieceI)
      loglike <- sum(logfirstpieceI)-secondpiece.wrap(I, trueremovaltime, beta[m], Rb[m], K, N, maxt, threshold,thresholdsum)
      logfirstpieceIstar <- log(firstpiece.wrap(Istar, beta[m], initialinfective, Rb[m], K, N, N_I, threshold))
      logfirstpieceIstar <- ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
      loglikestar <- sum(logfirstpieceIstar)-secondpiece.wrap(Istar, trueremovaltime, beta[m], Rb[m], K, N, maxt, threshold,thresholdsum)
      
      alpha.p <- predprobs[update]
      beta.p <- 1-alpha.p
      probifadded <- (sum(occult[update])+1)/m
      probifnotadded <- sum(occult[update]+.000000000000001)/m
      extra.piece=(length(addinf))/(length(N_I)-length(N_N)+1)*dbeta(probifadded, alpha.p, beta.p)*maxt #/dbeta(probifnotadded, alpha.p, beta.p)
      
      #metropolis hastings step for adding an infection
      mstep.I=min(1,exp(loglikestar-loglike)*extra.piece)
      if(mstep.I=="NaN") mstep.I=1
      R=runif(1)
      if(R<mstep.I){
        I<-Istar
        N_I<-c(N_I,update)
        infectedhousesI[N_I]=1
        bugs[update,]=bugsstar
        loglike <-loglikestar
        Q<-Qstar
        accept.I[m] <- 2
      }else{
        Istar<-I
        trueremovaltime[update]=tobs[update]=maxt
        detectiontime[update] <- Inf
        check3[update]<-Inf
        accept.I[m] <- 3}
    }
  }else{
    
    ###############
    ####delete I###
    ###############
    
    if(length(N_I)>length(N_N)){ 
      
      #pick which house to delete
      update=sampleWithoutSurprises(N_I[!(N_I %in% N_N)])
      Istar[update]=Inf
      check3[update]=Inf
      tobs[update]=maxt
      trueremovaltime[update]=detectiontime[update]=Inf
      bugstest <- bugs
      bugstest[update,]=rep(0,maxt)
      logfirstpieceI<-log(firstpiece.wrap(I, beta[m], initialinfective, Rb[m], K, N, N_I, threshold))
      logfirstpieceI=ifelse(logfirstpieceI=="-Inf",0,logfirstpieceI)
      loglike <- sum(logfirstpieceI)-secondpiece.wrap(I, trueremovaltime, beta[m], Rb[m], K, N, maxt, threshold,thresholdsum)
      logfirstpieceIstar<-log(firstpiece.wrap(Istar, beta[m], initialinfective, Rb[m], K, N, N_I, threshold))
      logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
      loglikestar <- sum(logfirstpieceIstar)-secondpiece.wrap(Istar, trueremovaltime, beta[m], Rb[m], K, N, maxt, threshold,thresholdsum)
      alpha.p <- predprobs[update]
      beta.p <- 1-alpha.p
      probifdeleted <- (sum(occult[update])-1)/m
      probifnotdeleted <- sum(occult[update])/m
      extra.piece=(length(N_I)-length(N_N))/(length(addinf))/dbeta(probifnotdeleted, alpha.p, beta.p)/maxt #/dbeta(probifnotdeleted, alpha.p, beta.p)
      #decide whether to accept new I
      mstep.I=min(1,exp(loglikestar-loglike)*extra.piece)
      if(mstep.I=="NaN") mstep.I=1
      R=runif(1)
      if(R<mstep.I){
        I<-Istar
        N_I<-N_I[which(N_I!=update)]
        accept.I[m] <- 5
        loglike <- loglikestar
        infectedhousesI[update]=0
        Q<-Qstar
        bugs<-bugstest
      }else{
        Istar<-I
        trueremovaltime[update]=maxt+1
        tobs[update]=maxt
        detectiontime[update]=Inf
        accept.I[m]=6
        check3[update]=bugs[update,tobs[update]]}
      }
  }
  occult[N_I[!(N_I %in% N_N)]]=occult[N_I[!(N_I %in% N_N)]]+1
  occult.prob<- occult/m
  rank <- 1:N
  occult.prob.ids <- cbind(id, occult.prob, location[,1], location[,2])
  occult.prob.ids <- occult.prob.ids[order(occult.prob, decreasing = TRUE),]
  occult.prob.ids <- cbind(occult.prob.ids, rank)
  
   if(m%%100==0) {print(m)
                 print(beta[m])
                 print(Rb[m])
                 rankings1[,(m/100)]<-occult.prob.ids[,1]
                 rankings2[,(m/100)]<-occult.prob.ids[,2]
  colfunc = gray.colors(length(unique(occult.prob.ids[,2])),start=1,end=0)[as.factor(occult.prob.ids[,2])]
  par(mfrow=c(1,1))
  par(mar=c(1, 1, 1, 1), xpd=TRUE)
  plot(occult.prob.ids[,3], occult.prob.ids[,4],col = colfunc,pch=16,cex=occult.prob.ids[,2]*50)
  #points(location[true.occult,1], location[true.occult,2],col = "blue")
  for (i in 1:N) if(sum.insp[i]>0) points(location[i,1],location[i,2],pch=18,col="firebrick4",cex=1)
  }
}

toc()

for (i in 1:N) if(sum.insp[i]>0) points(location[i,1],location[i,2],pch=18,col="firebrick",cex=2)
top <- occult.prob.ids[1:10,]
points(location[top,1],location[top,2],pch=18,col="gold",cex=1)

par(mfrow=c(1,1))
par(mar=c(3,3,3,3), xpd=TRUE)
plot(beta[1:m],type="l",main="beta")
abline(h=0.3,col="red")
median(beta[1000:m])
par(mar=c(3,3,3,3), xpd=TRUE)
plot(Rb[1:m-1],type="l",main="Rb")
abline(h=1.14,col="red")
median(Rb[1000:m])
par(mar=c(3,4.2,3,3), xpd=TRUE)
plot(1:(M/1000),rep(-100,(M/1000)),,ylim=c(200,1),xlab="",ylab="Ranking",axes=FALSE)
axis(2, at=c(1,50,100,150,200),cex=0.75)
mtext("Iteration", 1, line=1)
for(i in 1:200){
  for(j in 1:(M/100)){
 if(j%%10==0) points(j/10,(which(rankings1[,j]==i)),type="o",col=i,pch=18,cex=.5)
}
}

