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
N=500
indicator=indicator2=0
thresholdsum<-apply(threshold,1,sum)
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


check3<-ifelse(sum.insp>0,sum.insp,Inf) #replace with observed bug counts
I=ifelse(check3!=Inf&check3!=Inf,50,Inf) #set initial values for infection times

trueremovaltime=ifelse(check3<Inf,maxt+1,Inf) #set recovery times 
detectiontime=ifelse(check3>0&check3<Inf,tobs,Inf) #set detection time vector

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

#second piece of likelihood
secondpiece.update<-function(i,trueremovaltime,detectiontime,I,beta,r){
  S1=rep(0,N)
  for (j in 1:N){
    t=min(maxt,I[i],trueremovaltime[i])-min(I[j],I[i])
    if(t>0) S1[j]<-H(i, j, t, r, I, beta, K, threshold)
  }
  S1<-ifelse(S1=="NaN",0,S1)  
  sumofS=sum(S1)
  return(sumofS)
}

firstpiece.update<-function(update,I,beta,initialinfective,r){
  beta.I=H.mat=H.mat1=matrix(0,nrow=N,ncol=N)
  beta.sum=rep(0,N)
  j=update
  t = tobs[update]
  for (i in 1:N) if(i %in% N_I | i==initialinfective) {
    if(I[i]<I[j]&I[j]<trueremovaltime[i]) H.mat[i,j]=ht(t, r, I, i, j, beta, K, threshold)
  }
  beta.sum=sum(H.mat[,j])
  tempthres1=N
  beta.sum<-ifelse(is.na(beta.sum),0,beta.sum)
  
  return(beta.sum)
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


mleest<-function(par, data, threshold){
  beta=par[1]
  Rb=par[2]
  thresholdsum <- apply(threshold,1,sum)
  trueinfectiontime=data[,1]
  trueremovaltime=data[,2]
  bugs=data[,(3:dim(data)[2])]
  trueinfectiontime=ifelse(trueinfectiontime==0,Inf,trueinfectiontime)
  trueremovaltime=ifelse(trueremovaltime==0,Inf,trueremovaltime)
  initialinfective=which(trueinfectiontime==1)
  lambda_t=rep(0,N)
  lambda_t[1]=1
  Q=NULL
  for (i in which(trueinfectiontime!=Inf)){
    Q[i]=sum(f_D(i,bugs,trueinfectiontime,check3,Rb))
  }
  thirdpieceloglike=sum(Q[which(Q!="NA")])
  logfirstpiece<-log(firstpiece.wrap(trueinfectiontime, beta, initialinfective, Rb, K, N, N_I, threshold))
  logfirstpiece=ifelse(logfirstpiece=="-Inf"|logfirstpiece=="NaN",0,logfirstpiece)
  loglike=sum(logfirstpiece)-secondpiece.wrap(trueinfectiontime, trueremovaltime, beta, Rb, K, N, maxt, threshold,thresholdsum)
  likelihood<- loglike+thirdpieceloglike
  return(-likelihood)}

optim(par=c(.1,1.4),mleest, data=data,threshold=threshold,method="Nelder-Mead") 



likelihood=NULL
i=seq(0,1,by=0.001)
for(beta in i){
  Rb=1.1
  thresholdsum <- apply(threshold,1,sum)
  trueinfectiontime=data[,1]
  trueremovaltime=data[,2]
  bugs=data[,(3:dim(data)[2])]
  trueinfectiontime=ifelse(trueinfectiontime==0,Inf,trueinfectiontime)
  trueremovaltime=ifelse(trueremovaltime==0,Inf,trueremovaltime)
  initialinfective=which(trueinfectiontime==1)
  for (i in which(trueinfectiontime!=Inf)){
    Q[i]=sum(f_D(i,bugs,trueinfectiontime,check3,Rb))
  }
  thirdpieceloglike=sum(Q[which(Q!="NA")])
  logfirstpiece<-log(firstpiece.wrap(trueinfectiontime, beta, initialinfective, Rb, K, N, N_I, threshold))
  logfirstpiece=ifelse(logfirstpiece=="-Inf",0,logfirstpiece)
  loglike=sum(logfirstpiece)-secondpiece.wrap(trueinfectiontime, trueremovaltime, beta, Rb, K, N, maxt, threshold,thresholdsum)/N
  likelihood[beta*100]<- loglike+thirdpieceloglike
  
}
plot(likelihood,type="l")
which(likelihood==max(likelihood))

