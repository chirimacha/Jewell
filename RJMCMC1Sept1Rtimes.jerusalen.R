#set working drive
#setwd("/home/ebillig/Jewell_data")
setwd("/Users/EMWB/Jewell/Data/jerusalen")
#setwd("~/Desktop/Levy Lab")
#setwd("~/Users/e/Jewell/Data")
#setwd("/Users/mzlevy/Jewell/Data")
#setwd("C:/Users/ebillig/Documents/Jewell/Data/jerusalen")

#set seed
set.seed(9754)

#run function
#vary Rbstart between 1.05 and 1.4
Rbstart=1.4

#vary betastart between 0 and 1
betastart=0.3

#how long
totaliterations=10000

#run.mcmc <- function(totaliterations,Rbstart, betastart){

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

#read in data
jerusalen <- read.csv("bdjerusalen.csv")
datedata <- read.csv("tbl_vivienda.csv")
priors <- read.csv("Corentins_Predictions_Jun-24-2015_07-13-06.csv")
names(priors)[names(priors)=="UNICODE"] <- "unicode"

dataset<- merge(jerusalen, datedata, by="unicode",all.x=TRUE)
dataset <- merge(dataset, priors, by="unicode",all.x=TRUE)


##################################################
#######define parameters##############################
#########################################

########set up times#######

#Import dates as separate columns for month, day and year in that respective order
date <- function(m,d,y){
  #Convert it into one string
  right.date <- paste(m,d,y,sep = "/", collapse = NULL)
  #Read it as a date in the right format for R
  
  new.dates <- as.Date(right.date, "%m/%d/%Y")
  
  return(new.dates)
}
#outputs dates in the correct format that R uses
dataset$date <- as.Date(dataset$FECHA, "%m/%d/%Y")



############################################################
############# UPDATED CHANGES AND NEW FOR LOOP #############
############################################################
#if more than one observation for a house, pick most recent

#identify unique unicodes
unicode<-as.character(dataset$unicode)
unique.unicodes <- unique(unicode)
dates <- dataset$date

#find repeated unicodes
repeated.unicodes <- unicode[which(duplicated(unicode) == TRUE)]
#setting an empty vector for the for loop
unique.dates <- c(1:length(unique.unicodes)*NA)

#for houses with more than on observation, only take into account the one with latest date

for (i in 1:length(unique.unicodes)){
  
  u <- which(unicode==unique.unicodes[i])
  fecha <- dataset$date[u]
  
  if(is.na(fecha[1]) ==FALSE) {
    maxf <- max(fecha, na.rm = TRUE)
    v <- which(dataset$date == maxf)
    unique.dates[i] <- intersect(u,v)
  }
  else {
    unique.dates[i] <- max(u)
  }
}


#new dataframe with single unicodes 
unique.data <-dataset[unique.dates,]
dataset <- unique.data


earliest <- sort(dataset$date)[1]
latest <- sort(dataset$date)[length(dataset$date[which(!is.na(dataset$date))])]
timetest <- (dataset$date - earliest)/90
initialtime <- as.Date("01/01/04", "%m/%d/%Y")
today <- as.Date("10/01/15", "%m/%d/%Y")
timefrombeginning <- round((earliest - initialtime)/90)
tobs <- ceiling(timetest) + timefrombeginning
maxt <- round((today - latest)/90)+max(tobs[which(!is.na(tobs))])
maxt <- as.numeric(maxt)

#set removal times
trueremovaltimetest <- rep(Inf, length(dataset$unicode))

#identify which houses weren't inspected
inspected <- ifelse(is.na(tobs),0,1)

#replace NAs with max time
tobs = ifelse(is.na(tobs), maxt,tobs)

#sum inspecciones in districts 4,5,6
sum.insp <- ifelse(is.na(dataset$capt_upch),0,dataset$capt_upch)



#Replace NA values for prior probability with median value
#find median value of those that are not NA


#replace NAs with this value
predprobs <- rep(0.02,length(sum.insp))

#get unicodes as strings
unicode<-as.character(dataset$unicode)


########################################################
#########MCMC algorithm###################################
#######################################################
tic()   #begin timer


##Jewell MCMC
M <- totaliterations
m <- 1 #first iteration


#define number of houses
N <- dim(dataset)[1]
occult.sum.new<-rep(0,N)
infectiontime<-rep(Inf,N)

#calculate distances between houses
distance<-matrix(NA,nrow=N,ncol=N)
for (i in 1:N){
  for (j in 1:N){
    distance[i,j]=sqrt((dataset$easting[i]-dataset$easting[j])^2+(dataset$northing[i]-dataset$northing[j])^2)
  }
}


check3<- rep(Inf,N) #initialize data vector
T_b <- 50 #threshold for bug infectiousness
jumpprob <- 0.01 #probability of jump vs. hop
lambda <- 0.300 #block factor from streets paper (Barbu, 2013)
delta <- 9.00 #sd from streets paper (Barbu, 2013)
bugs <- matrix(0,nrow=N,ncol=maxt) #initialize but matrix
maxbugs <- max(sum.insp) #find most observed bugs in data
initialinfective <- which(sum.insp==maxbugs) #set this house as initialinfective
id=1:N #generate ids
K=1000 #carrying capacity

# #option1: define threshold as block
# noblock <- which(is.na(dataset$uniblock))
# dataset$UNICODE[noblock]
# dataset$uniblock[noblock[1:2]]<-dataset$uniblock[127]
# dataset$uniblock[noblock[3]]<-dataset$uniblock[130]
# dataset$uniblock[noblock[4]]<-dataset$uniblock[138]
# dataset$uniblock[noblock[5]]<-dataset$uniblock[159]
# dataset$uniblock[noblock[6]]<-dataset$uniblock[190]

# thresholdblocks<-matrix(0,nrow=N,ncol=N)
# for(i in 1:N){
#   for(j in 1:N){
#     thresholdblocks[i,j] <- ifelse(dataset$uniblock[i]==dataset$uniblock[j], 1 , lambda)
#   }
# }
# 
# threshold1<-matrix(0,nrow=N,ncol=N)
# for(i in 1:N){
#   for(j in 1:N){
#     threshold1[i,j] <- thresholdblocks[i,j]*exp(-distance[i,j]/delta)
#   }
# }


#option2: define threshold by radius of T_b
#probability of infestation differs by hops (<T_b m) or jumps (>T_b m)
threshold2 <- ifelse(distance<T_b, 1 , jumpprob)


#decide which one you are using; this is what to change
threshold <- threshold2
#threshold <- threshold1

thresholdsum <- apply(threshold,1,sum)

check3<-ifelse(sum.insp>0,sum.insp,Inf) #replace with observed bug counts
I=ifelse(check3!=Inf&check3!=Inf,2,Inf) #set initial values for infection times

detectiontime=ifelse(check3>0&check3<Inf,tobs,Inf) #set detection time vector

#initialize parameters 
beta=Rb=rep(0,M)
Rb[1]=Rbstart
beta[1]=betastart
betastar=.3
truebeta=0.5
trueRb = 1.1
betastar.I=matrix(0,nrow=N,ncol=N)
betastar.sum=rep(0,N)
S=H.mat=matrix(0,nrow=N,ncol=N)
U=rep(0,N)
accept.beta=accept.I=accept.Rb=rep(0,M)

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

#those treated before X date are counted as susceptible
trueremovaltime <- ceiling(ifelse(trueremovaltimetest<=43&I==Inf,Inf,trueremovaltimetest))
trueremovaltime <- ifelse(trueremovaltime<tobs&check3<Inf,Inf,trueremovaltime)

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
tic()
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
  a.beta <- 100*beta[m-1]
  b.beta <- 100-beta[m-1]*100
  #betastar=rbeta(1 , a.beta , b.beta)
  betastar=abs(rnorm(1,beta[m-1],0.1))
  if(betastar>1){betastar <- 1-(betastar-1)}
  logfirstpiecestar<-log(firstpiece.wrap(I, betastar, initialinfective, Rb[m], K, N, N_I, threshold))
  logfirstpiecestar=ifelse(logfirstpiecestar=="-Inf",0,logfirstpiecestar)
  loglikestar=sum(logfirstpiecestar)-secondpiece.wrap(I, trueremovaltime, betastar, Rb[m], K, N, maxt, threshold,thresholdsum)
  betapriors <- dbeta(betastar,shape1=truebeta*10,shape2=10-truebeta*10,log=TRUE)-dbeta(beta[m-1],shape1=truebeta*10,shape2=10-truebeta*10,log=TRUE)
  mstep.beta=min(1,exp(loglikestar-loglike+betapriors)) #-dbeta(betastar, a.beta, b.beta, log=TRUE)+dbeta(beta[m-1], a.beta, b.beta,log=TRUE)
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
    if(inspected[update]==1) Istar[update] <- ifelse(sum.insp[update]==0,sample(c((tobs[update]+1):(maxt-2)),1,replace=TRUE),sample(c(2:(tobs[update]-1)),1,replace=TRUE))
    if(inspected[update]==0) Istar[update] <- sample(c(2:(maxt-2)),1,replace=TRUE)
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
    
    addinf<-which(I==Inf&(inspected==0 | (inspected == 1 & (tobs<(maxt-2) ))))
    if(length(addinf)>1){
      update=sample(addinf,1,replace=TRUE)
      Istar[update] <- ifelse(inspected[update]==0,floor(runif(1,min=2,max=maxt-2)),floor(runif(1,min=tobs[update],max=maxt-2)))
      N_Istar <- c(N_I, update)
      trueremovaltime[update]=maxt+1
      detectiontime[update]=maxt
      if(is.na(tobs[update])) tobs[update]=maxt
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
      logfirstpieceIstar <- log(firstpiece.wrap(Istar, beta[m], initialinfective, Rb[m], K, N, N_Istar, threshold))
      logfirstpieceIstar <- ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
      loglikestar <- sum(logfirstpieceIstar)-secondpiece.wrap(Istar, trueremovaltime, beta[m], Rb[m], K, N, maxt, threshold,thresholdsum)
      
      alpha.p <- 10*predprobs[update]
      beta.p <- 10-alpha.p
      probifadded <- (sum(occult[update])+1)/m
      probifnotadded <- sum(occult[update])/m
      extra.piece=(length(addinf))/(length(N_I)-length(N_N)+1)*dunif(1, min=1, max=length(addinf))/dunif(1,min=0,max=(length(N_I)-length(N_N)+1)) #*dbeta(probifadded, alpha.p, beta.p)
      
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
        trueremovaltime[update]=maxt
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
      addinf<-which(I==Inf&(inspected==0 | (inspected == 1 & (tobs<maxt-2))))
      update=sampleWithoutSurprises(N_I[!(N_I %in% N_N)])
      Istar[update] <- Inf
      check3[update] <- Inf
      #tobs[update] <- maxt
      N_Istar<-N_I[which(N_I!=update)]
      trueremovaltime[update] = detectiontime[update] <- Inf
      bugstest <- bugs
      bugstest[update,] <- rep(0,maxt)
      logfirstpieceI <- log(firstpiece.wrap(I, beta[m], initialinfective, Rb[m], K, N, N_I, threshold))
      logfirstpieceI <- ifelse(logfirstpieceI=="-Inf",0,logfirstpieceI)
      loglike <- sum(logfirstpieceI)-secondpiece.wrap(I, trueremovaltime, beta[m], Rb[m], K, N, maxt, threshold,thresholdsum)
      logfirstpieceIstar<-log(firstpiece.wrap(Istar, beta[m], initialinfective, Rb[m], K, N, N_Istar, threshold))
      logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
      loglikestar <- sum(logfirstpieceIstar)-secondpiece.wrap(Istar, trueremovaltime, beta[m], Rb[m], K, N, maxt, threshold,thresholdsum)
      alpha.p <- 10*predprobs[update]
      beta.p <- 10-alpha.p
      probifdeleted <- (sum(occult[update]))/m
      probifnotdeleted <- sum(occult[update]+1)/m
      extra.piece <- (length(N_I)-length(N_N))/(length(addinf)+1)*dunif(1,min=0,max=(length(N_I)-length(N_N)))/dunif(1, min=1, max=length(addinf)+1)  #/dbeta(probifnotdeleted, alpha.p, beta.p)
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
        #tobs[update]=maxt
        detectiontime[update]=Inf
        accept.I[m]=6
        check3[update]=bugs[update,tobs[update]]}
    }
  }
  occult[N_I[!(N_I %in% N_N)]]=occult[N_I[!(N_I %in% N_N)]]+1
  #occult.sum <- apply(occult,1,sum)
  occult.prob<- occult/m
  occult.prob.ids <- data.frame(id, occult.prob, dataset$easting, dataset$northing, unicode)
  occult.prob.ids.ordered <- occult.prob.ids[order(occult.prob, decreasing = TRUE),]
  if(m%%10000==0) {
  	print(m)
  	write.csv(occult.prob.ids.ordered, file=paste("/home/ebillig/Jewell_data/Oct5/Rb",Rbstart,"beta",betastart,"ResultsSept13.csv", sep=""))
    save.image(paste("/home/ebillig/Jewell_data/Oct5/Rb",Rbstart,"beta",betastart,"ResultsSept13.Rdata", sep=""))}
if(m%%100==0){
  print(m)
  print(N_I)
  print(beta[m])
  print(Rb[m])
  colfunc = gray.colors(length(unique(as.numeric(occult.prob.ids.ordered[,2]))),start=1,end=0)[as.factor(occult.prob.ids.ordered[,2])]
  plot(as.numeric(occult.prob.ids.ordered[,3]), as.numeric(occult.prob.ids.ordered[,4]),col = colfunc,pch=16,cex=as.numeric(occult.prob.ids.ordered[,2])*50) #as.numeric(Results1[,2])*2000)
  points(occult.prob.ids.ordered[1:10,3], occult.prob.ids.ordered[1:10,4],col = "gold",pch=18) 
  for (i in 1:N) if(sum.insp[i]>0) points(dataset$easting[i],dataset$northing[i],pch=18,col="firebrick3")
}
  if(m%%M==0) print(toc())
}
toc()
#return(occult.prob.ids.ordered)

#}
# 
# #run code
# Results <- run.mcmc(totaliterations,Rbstart,betastart)
# 
# #record results
# write.csv(Results, file=paste("Rb",Rbstart,"beta",betastart,"ResultsSept13.csv", sep=""))
# #save.image("ResultsSept7.Rdata")
# 
# 
# #plot results
# colfunc = gray.colors(length(unique(as.numeric(Results[,2]))),start=1,end=0)[as.factor(Results[,2])]
# plot(as.numeric(Results[,3]), as.numeric(Results[,4]),col = colfunc,pch=16,cex=as.numeric(Results[,2])*50000) #as.numeric(Results1[,2])*2000)
# for (i in 1:N) if(sum.insp[i]>0) points(dataset$X[i],dataset$Y[i],pch=18,col="firebrick3")
# legend("topleft",c("Known Infested House", "Top Probability of Infestation"),pch=c(18,18),col=c("firebrick3","gold"),bty="n")
# points(dataset$X,dataset$Y,col=threshold1[N_N,],cex=threshold1[N_N,])
# 
#   if(m%%100==0){
#     print(N_I)
#     print(beta[m])
#     print(Rb[m])
#   colfunc = gray.colors(length(unique(as.numeric(occult.prob.ids.ordered[,2]))),start=1,end=0)[as.factor(occult.prob.ids.ordered[,2])]
#   plot(as.numeric(occult.prob.ids.ordered[,3]), as.numeric(occult.prob.ids.ordered[,4]),col = colfunc,pch=16,cex=as.numeric(occult.prob.ids.ordered[,2])*10) #as.numeric(Results1[,2])*2000)
#   points(occult.prob.ids.ordered[1:10,3], occult.prob.ids.ordered[1:10,4],col = "gold",pch=18) 
#   for (i in 1:N) if(sum.insp[i]>0) points(dataset$X[i],dataset$Y[i],pch=18,col="firebrick3")
   points(dataset$easting,dataset$northing,col="green",cex=1-inspected)
#   points(dataset$X,dataset$Y,col="green",cex=1-inspected)}
# 
# 
