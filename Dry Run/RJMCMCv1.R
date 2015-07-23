#set working drive
setwd("/home/ebillig/Jewell_data")
setwd("/Users/EMWB/Jewell/Data")
#setwd("~/Desktop/Levy Lab")
#setwd("~/Users/e/Jewell/Data")
#setwd("/Users/mzlevy/Jewell/Data")



run.mcmc <- function(banditarm, totaliterations){

#set seed
set.seed(123456)

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
inspecciones = read.csv("inspecciones.csv")
data = read.csv(file=paste("tiabaya",banditarm,"_waddlhouses.csv",sep=""))
priors = read.csv("Corentins_Predictions_Jun-24-2015_07-13-06.csv")

#drop X variable from priors dataset
priors$X=NULL

#capitalize unicode in data
data <- rename(data,c("unicode"="UNICODE"))
#merge data
insp.data <- merge(inspecciones,data, by="UNICODE",all.y=TRUE)
dataset <- merge(insp.data, priors, by="UNICODE",all.x=TRUE)

#Replace NA values with zeros for sum
dataset$PD_TCAP_TOT <- ifelse(is.na(dataset$PD_TCAP_TOT),0,dataset$PD_TCAP_TOT)
dataset$IN_TCAP_TOT <- ifelse(is.na(dataset$IN_TCAP_TOT),0,dataset$IN_TCAP_TOT)

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
dataset$date <- date(dataset$MES,dataset$DIA,dataset$ANIO)


############################################################
############# UPDATED CHANGES AND NEW FOR LOOP #############
############################################################
#if more than one observation for a house, pick most recent

#identify unique unicodes
unicode<-as.character(dataset$UNICODE)
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


#if more than one observation for a house, pick most recent
dataset$date <- date(dataset$MES,dataset$DIA,dataset$ANIO)

#new dataframe with single unicodes 
unique.data <-dataset[unique.dates,]

#rename dataset
dataset <- unique.data

earliest <- sort(dataset$date)[1]
latest <- sort(dataset$date)[length(dataset$date[which(!is.na(dataset$date))])]
timetest <- (dataset$date - earliest)/90
initialtime <- date(12, 31, 2004)
today <- date(7, 22, 2015)
timefrombeginning <- round((earliest - initialtime)/90)
tobs <- ceiling(timetest) + timefrombeginning
maxt <- round((today - latest)/90)+max(tobs[which(!is.na(tobs))])
maxt <- as.numeric(maxt)

#identify which houses weren't inspected
inspected <- ifelse(is.na(tobs),0,1)

#replace NAs with max time
tobs = ifelse(is.na(tobs), maxt,tobs)

#sum inspecciones in districts 4,5,6
sum.insp <- dataset$PD_TCAP_TOT + dataset$IN_TCAP_TOT


#Replace NA values for prior probability with median value
#find median value of those that are not NA
median.pred.prob <- median(dataset$predicteddensity[which(!is.na(dataset$predicteddensity))])

#replace NAs with this value
predprobs <- ifelse(is.na(dataset$predicteddensity), median.pred.prob, dataset$predicteddensity)

#get unicodes as strings
unicode<-as.character(dataset$UNICODE)


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
    distance[i,j]=sqrt((dataset$X[i]-dataset$X[j])^2+(dataset$Y[i]-dataset$Y[j])^2)
  }
}


check3<- rep(Inf,N) #initialize data vector
T_b <- 30 #threshold for bug infectiousness
jumpprob <- .01 #probability of jump vs. hop
bugs <- matrix(0,nrow=N,ncol=maxt) #initialize but matrix
maxbugs <- max(sum.insp) #find most observed bugs in data
initialinfective <- which(sum.insp==maxbugs) #set this house as initialinfective
id=1:N #generate ids
K=1000 #carrying capacity
tuning <- 0.01 #tuning parameter for RJ

#option1: define threshold as block

#first fill in NAs
#dataset$polygon[128]=191
#dataset$polygon[129]=191
#dataset$polygon[131]=191
#dataset$polygon[139]=195
#dataset$polygon[160]=256
#dataset$polygon[191]=255

#threshold1<-matrix(0,nrow=N,ncol=N)
#for(i in 1:N){
#  for(j in 1:N){
#    threshold1[i,j] <- ifelse(dataset$polygon[i]==i.v.gps$polygon[j], 1 , jumpprob)
#  }
#}


#option2: define threshold by radius of T_b
#probability of infestation differs by hops (<T_b m) or jumps (>T_b m)
threshold2 <- ifelse(distance<T_b, 1 , jumpprob)


#decide which one you are using; this is what to change
threshold <- threshold2
#threshold <- threshold1

check3<-ifelse(sum.insp>0,sum.insp,Inf) #replace with observed bug counts
I=ifelse(check3!=Inf&check3!=Inf,2,Inf) #set initial values for infection times

trueremovaltime=ifelse(check3<Inf,maxt+1,Inf) #set recovery times 
detectiontime=ifelse(check3>0&check3<Inf,tobs,Inf) #set detection time vector

#initialize parameters 
beta=Rb=rep(0,M)
beta[1]=.3
betastar=.3
betastar.I=matrix(0,nrow=N,ncol=N)
betastar.sum=rep(0,N)
S=H.mat=matrix(0,nrow=N,ncol=N)
U=rep(0,N)
accept.beta=accept.Iadd=accept.Idel=rep(0,M)

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

#initialize Rb parameter
Rb[1]=1.11
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
beverton.holt.I<-function(update,K,R,check3,tobs){
  tobs1=tobs
  time=ifelse(update %in% N_N, -log((K-check3[update])/(check3[update]*K-check3[update]))/log(R[m-1]),-log((K-bugs[update,tobs1])/(bugs[update,tobs1]*K-bugs[update,tobs1]))/log(R[m-1]))
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


#hazard function
ht <- function(t, r, I, i, j, beta, K, threshold) {
  n <- 1
  deriv <- K*n^2*log(r)/((K-n)*exp(-log(r)*(t-I[i]))+n)^2
  ifelse(deriv>0, 1-(1-beta*threshold[i,j])^(deriv), 0)
}

#cumulative hazard function
H<-function(i,j,t,Rb,I,beta){
  r=log(Rb)
  Ht=t*(1-(1-beta*threshold[i,j])^(r/K))-
    ((K-1)*r^2*t^2*log(1-beta*threshold[i,j])*(1-beta*threshold[i,j])^(r/K))/(K^2)-
    t^3*((K-1)*r^3*log(1-beta*threshold[i,j])*(1-beta*threshold[i,j])^(r/K)*(2*(K-1)*r*
                                                                               log(1-beta*threshold[i,j])+K*(2*K-3)))/(3*K^4)
  return(Ht)}

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
  H.mat <- matrix(0, nrow=N, ncol=N)
  beta.sum <- rep(0, N)
  for (j in 1:N) {
    for (i in 1:N) {
      if (i %in% N_I | i==initialinfective) {
        if (I[i]<I[j] & I[j]<Inf) {
          n <- 1
          deriv <- K*n^2*log(r)/((K-n)*exp(-log(r)*(I[j]-I[i]))+n)^2
          H.mat[i,j] <- ifelse(deriv>0, 1-(1-beta*threshold[i,j])^(deriv), 0)
          # print(paste("i", i, "j", j, "I[i]", I[i], "I[j]", I[j],
          #   "deriv", deriv, "hmat", H.mat[i,j]))
        }
      }
    }
    beta.sum[j] <- sum(H.mat[,j])
  }
  ifelse(is.na(beta.sum), 0, beta.sum)
}

first.include <- '
#include <set>
#include <cmath>
'

firstpiece.wrap <- cxxfunction(signature(IS="numeric", betaS="float",
                                         initialinfectiveS="int", rS="float", KS="float", NS="int", N_IS="numeric",
                                         thresholdS="numeric"), plugin="Rcpp", incl=first.include, body='
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
                               Rcpp::NumericVector H(N, 0.0);
                               const int n=1;
                               for (int j=0; j<N; ++j) {
                               double total=0;
                               for (infiter it=infecteds.begin(); it!=infecteds.end(); ++it) {
                               int i=*it - 1; // -1 to convert to 0-based indexing.
                               if (I[i]<I[j] && std::isfinite(I[j])) {
                               double deriv=K*std::pow(n, 2)*std::log(r)/
                               std::pow(
                               (K-n)*std::exp(
                               -std::log(r)*(I[j]-I[i])
                               )+n,
                               2
                               );
                               if (deriv>0) {
                               double add=1-std::pow(1-beta*threshold[i*N+j], deriv);
                               //std::cout << "i "<<i<<" j " << j<<" I i "<<I[i]<<" I[j] "<< I[j]
                               //  <<" deriv " << deriv << " hmat " << add <<std::endl;
                               total+=add;
                               }
                               }
                               }
                               H[j]=total;
                               }
                               return H;
                               ')


H <- function(i, j, t, Rb, I, beta, K, threshold) {
  r <- log(Rb)
  obt <- (1-beta*threshold[i,j])^(r/K)
  t*(1-obt)-((K-1)*r^2*t^2*log(1-beta*threshold[i,j])*obt)/(K^2)-
    t^3*((K-1)*r^3*log(1-beta*threshold[i,j])*obt*(2*(K-1)*r*
                                                     log(1-beta*threshold[i,j])+K*(2*K-3)))/(3*K^4)
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
secondpiece <- function(I, beta, r, K, N, maxt, threshold) {
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
  sum(S1)
}


secondpiece.wrap <- cxxfunction(signature(IS="numeric", betaS="float",
                                          rS="float", KS="float", NS="int", maxtS="float",
                                          thresholdS="numeric"), plugin="Rcpp", incl=first.include, body='
                                Rcpp::NumericVector I(IS);
                                double beta=Rcpp::as<double>(betaS);
                                double logr=std::log(Rcpp::as<double>(rS));
                                double K=Rcpp::as<double>(KS);
                                int N=Rcpp::as<int>(NS);
                                double maxt=Rcpp::as<double>(maxtS);
                                Rcpp::NumericVector threshold(thresholdS);
                                
                                double total=0;
                                for (int i=0; i<N; ++i) {
                                if (std::isfinite(I[i])) {
                                for (int j=0; j<N; ++j) {
                                double t=std::min(maxt, I[j]) - std::min(I[i], I[j]);
                                if (t>0) {
                                double obt=std::pow(1-beta*threshold[i*N+j], logr/K);
                                double add=t*(1-obt)-((K-1)*logr*logr*t*t*std::log(1-beta*threshold[i*N+j])*obt)/(K*K)-
                                std::pow(t, 3)*((K-1)*std::pow(logr, 3)*std::log(1-beta*threshold[i*N+j])*obt*(2*(K-1)*logr*
                                std::log(1-beta*threshold[i*N+j])+K*(2*K-3)))/(3*std::pow(K, 4));
                                if (std::isfinite(add)) {
                                total+=add;
                                }
                                }
                                }
                                }
                                }
                                return Rcpp::wrap(total);
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
  
  ###############
  ##update beta##
  ###############
  
  betastar=abs(rnorm(1,beta[m-1],.01))
  if(betastar > 1) betastar <- 1-(betastar-1)
  logfirstpiecestar<-log(firstpiece.wrap(I, betastar, initialinfective, Rb[m-1], K, N, N_I, threshold) )
  logfirstpiecestar=ifelse(logfirstpiecestar=="-Inf",0,logfirstpiecestar)
  logfirstpiece<-log(firstpiece.wrap(I, beta[m-1], initialinfective, Rb[m-1], K, N, N_I, threshold))
  logfirstpiece=ifelse(logfirstpiece=="-Inf",0,logfirstpiece)
  dbetastar=sum(logfirstpiecestar)-secondpiece.wrap(I, betastar, Rb[m-1], K, N, maxt, threshold)+dunif(betastar,min=.001,max=.99,log=TRUE)
  dbeta=sum(logfirstpiece)-secondpiece.wrap(I, beta[m-1], Rb[m-1], K, N, maxt, threshold)+dunif(beta[m-1],min=.001,max=.99,log=TRUE)
  mstep.beta=min(1,exp(sum(dbetastar)-sum(dbeta)))
  if(mstep.beta=="NaN") mstep.beta=1
  R=runif(1)
  if(R<mstep.beta){
    beta[m]=betastar
    #accept.beta[m]=1
  }else{
    beta[m]=beta[m-1]
    #accept.beta[m]=0
  }
  
  
  #######################################################
  #############update times##############################
  #######################################################
  
  
  ############
  ##update I##
  ############
  
  ##pick a house to update the time out infected houses
  update=sample(N_I,1)
  if(bugs[update,min(maxt,trueremovaltime[update])]==0) bugs[update,min(maxt,trueremovaltime[update])]=1
  Istar[update]=ifelse(check3[update]>0,round(rnorm(1,beverton.holt.I(update,K,Rb,check3,tobs[update]),1)),floor(rnorm(1,beverton.holt.I(update,K,Rb,bugs[,min(maxt,trueremovaltime[update])], min(maxt,trueremovaltime[update])),1)))
  if(Istar[update]=="NaN") {bugs[update,min(maxt,trueremovaltime[update])]=1000
                            Istar[update]=ifelse(check3[update]>0,round(beverton.holt.I(update,K,Rb,check3,tobs[update])),floor(beverton.holt.I(update,K,Rb,bugs[,min(maxt,trueremovaltime[update])], min(maxt,trueremovaltime[update]))))
  }
  if(Istar[update]>=maxt) Istar[update]=maxt-1
  Istar[update]=ifelse(Istar[update]<1,2,Istar[update])
  bugsstar=rep(0,maxt)
  bugsstar[Istar[update]]=1
  bugsstar=beverton.holt.update(K,Rb[m-1],bugsstar,maxt,Istar[update])
  bugsstar[tobs[update]]=check3[update]
  logfirstpieceI<-log(firstpiece.wrap(I, beta[m], initialinfective, Rb[m-1], K, N, N_I, threshold))
  logfirstpieceI=ifelse(logfirstpieceI=="-Inf",0,logfirstpieceI)
  loglike.I=sum(logfirstpieceI)-secondpiece.wrap(I, beta[m], Rb[m-1], K, N, maxt, threshold)
  logfirstpieceIstar<-log(firstpiece.wrap(Istar, beta[m], initialinfective, Rb[m-1], K, N, N_I, threshold) )
  logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
  loglike.Istar=sum(logfirstpieceIstar)-secondpiece.wrap(Istar, beta[m], Rb[m-1], K, N, maxt, threshold)
  Q=sum(f_D.update(update,bugsstar,Istar,check3,Rb[m-1]))-sum(f_D(update,bugs,I,check3,Rb[m-1]))
  
  
  #Metroplis step; decide whether to accept new time#
  
  mstep.I=min(1,exp(loglike.Istar-loglike.I+Q))
  if(mstep.I=="NaN") mstep.I=1
  
  R=runif(1)
  if(R<mstep.I){
    I<-Istar
    bugs[update,]=bugsstar
    loglike.I<-loglike.Istar
  }else{
    Istar<-I
  }
  
  ################################
  ######update Rb##################
  ####################################
  
  Rbstar=rnorm(1,Rb[m-1],.1)
  if(Rbstar <=1) Rbstar <- 1+(1-Rbstar)
  Q=Qstar=rep(NA,length(I!=Inf))
  
  for (i in which(I!=Inf)){
    Qstar[i]=sum(f_D(i,bugs,I,check3,Rbstar))
    Q[i]=sum(f_D(i,bugs,I,check3,Rb[m-1]))
  }
  thirdpieceloglike=sum(Qstar[which(Qstar!="NA")])-sum(Q[which(Q!="NA")])+dunif(Rbstar,min=1,max=2,log=TRUE)-dunif(Rb[m-1],min=1,max=2,log=TRUE)
  logfirstpieceIstar<-log(firstpiece.wrap(I, beta[m], initialinfective, Rbstar, K, N, N_I, threshold) )
  logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
  loglike.Istar=sum(logfirstpieceIstar)-secondpiece.wrap(I, beta[m], Rbstar, K, N, maxt, threshold)
  Rbloglike=loglike.Istar+thirdpieceloglike-loglike.I
  
  #Metropolis step
  mstep.Rb=min(0,Rbloglike)
  #if(mstep.Rb=="NaN") mstep.Rb=0
  R=log(runif(1))
  if(R<mstep.Rb){
    Rb[m]=Rbstar
    loglike.I<-loglike.Istar
    Q<-Qstar
  }else{
    Rb[m]=Rb[m-1]}
  
  ########################################################
  #####decide whether to update I, add I, or delete I#####
  ########################################################
  
  add.del<-sample(c("add","del"),1)
  if(add.del=="add"){ 
    
    ###########
    ###add I###
    ##########
    
    addinf<-which((I==Inf&tobs<maxt-2&inspected==1&dataset$add.house==0)|(I==Inf&inspected==0&dataset$add.house==0))
    if(length(addinf)>1){
      update=sample(addinf,1)
      Istar[update]=ifelse(inspected[update]==1,floor(runif(1,min=tobs[update],max=maxt-1)),floor(runif(1,min=2,max=maxt-2)))
      trueremovaltime[update]=maxt+1
      detectiontime[update]=maxt
      tobs[update]=maxt
      bugsstar=rep(0,maxt)
      bugsstar[Istar[update]]=1
      bugsstar=beverton.holt.update(K,Rb[m],bugsstar,trueremovaltime[update],Istar[update])
      bugstest=bugs
      bugstest[update,]=bugsstar[1:maxt]
      check3[update]<-ifelse(inspected[update]==1,0,bugsstar[tobs])
      logfirstpieceIstar<-log(firstpiece.update(update,Istar,beta[m], initialinfective,Rb[m]))
      logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
      loglike=sum(logfirstpieceIstar)-secondpiece.update(update,trueremovaltime,detectiontime,Istar,beta[m],Rb[m])
      if(loglike==0) loglike=-Inf
      alpha.p <- predprobs[update]
      beta.p <- 1-alpha.p
      probifadded <- (sum(occult[update])+1)/m
      probifnotadded <- sum(occult[update]+.000000000000001)/m
      extra.piece=(length(addinf))/(maxt*(length(N_I)-length(N_N)+1))*dbeta(probifadded, alpha.p, beta.p) #/dbeta(probifnotadded, alpha.p, beta.p)
      
      #metropolis hastings step for adding an infection
      mstep.I=min(1,exp(loglike)*extra.piece)
      if(mstep.I=="NaN") mstep.I=1
      R=runif(1)
      if(R<mstep.I){
        I<-Istar
        N_I<-c(N_I,update)
        infectedhousesI[N_I]=1
        #accept.Iadd[m]=1
        bugs[update,]=bugsstar
        loglike.I<-loglike.Istar
        Q<-Qstar
      }else{
        Istar<-I
        #accept.Iadd[m]=2
        trueremovaltime[update]=tobs[update]=maxt
        detectiontime[update] <- Inf
        check3[update]<-Inf}
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
      bugsstest<-bugs
      bugstest[update,]=0
      logfirstpieceIstar<-log(firstpiece.update(update,I,beta[m], initialinfective,Rb[m]))
      logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
      loglike.Istar=-sum(logfirstpieceIstar)+secondpiece.update(update,trueremovaltime,detectiontime,I,beta[m],Rb[m])
      loglike=loglike.Istar 
      alpha.p <- predprobs[update]
      beta.p <- 1-alpha.p
      probifdeleted <- (sum(occult[update])-1)/m
      probifnotdeleted <- sum(occult[update])/m
      extra.piece=maxt*(length(N_I)-length(N_N))/(length(addinf))*dbeta(probifdeleted, alpha.p, beta.p) #/dbeta(probifnotdeleted, alpha.p, beta.p)
      #decide whether to accept new I
      mstep.I=min(1,exp(loglike.Istar)*extra.piece)
      if(mstep.I=="NaN") mstep.I=1
      R=runif(1)
      if(R<mstep.I){
        I<-Istar
        N_I<-N_I[which(N_I!=update)]
        #accept.Idel[m]=1
        loglike.I<-loglike.Istar
        infectedhousesI[update]=0
        Q<-Qstar
        bugs<-bugstest
      }else{
        Istar<-I
        trueremovaltime[update]=maxt+1
        tobs[update]=maxt
        detectiontime[update]=Inf
        #accept.Idel[m]=2
        check3[update]=bugs[update,tobs[update]]}
      
      
    }
  }
  occult[N_I[!(N_I %in% N_N)]]=occult[N_I[!(N_I %in% N_N)]]+1
  #occult.sum <- apply(occult,1,sum)
  occult.prob<- occult/m
  occult.prob.ids <- cbind(id, occult.prob, dataset$X, dataset$Y)
  occult.prob.ids.ordered <- cbind(occult.prob.ids,unicode)
  occult.prob.ids.ordered <- occult.prob.ids.ordered[order(occult.prob, decreasing = TRUE),]
  }
toc()
return(occult.prob.ids.ordered)

}

