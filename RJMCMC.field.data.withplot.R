set.seed(123456)


library("lubridate")
library("PBSmapping")
library("plyr")
library("inline")
library("Rcpp")
###read data###

#set working drive
setwd("/Users/EMWB/Jewell/Data")
#setwd("~/Desktop/Levy Lab")
#setwd("~/Users/e/Jewell/Data")
#setwd("/Users/mzlevy/Jewell/Data")

#set up tiabaya_gps dataset
tiabaya.gps = read.csv("Tiabaya_GPS.csv")

getUTM<-function(id, x,y){
  coord <- data.frame(x,y)
  colnames(coord) <- c("X","Y")
  attributes(coord)$projection <- "LL"
  utm.coord <- convUL(coord)
  utm.coord$X<-utm.coord$X*1000
  utm.coord$Y<-utm.coord$Y*1000
  return(utm.coord)
}
tiabaya.test <- getUTM(x=tiabaya.gps$LONGITUDE,y=tiabaya.gps$LATITUDE)

#link unicodes with coordinates
tiabaya.gps <- cbind(tiabaya.gps$UNICODE,tiabaya.test)
tiabaya.gps <- rename(tiabaya.gps,c("tiabaya.gps$UNICODE" = "UNICODE"))


#read in data
inspecciones = read.csv("inspecciones.csv")
vig = read.csv("byHouse_fullEID.csv")


#merge data
i.v <- merge(inspecciones,vig, by="UNICODE",all=TRUE)
i.v.gps <- merge(i.v,tiabaya.gps,by="UNICODE")

#subset data to only include localities 4,5,6
i.v.gps.456 <- i.v.gps[which(i.v.gps$L.y==4 | i.v.gps$L.x==4 | i.v.gps$L.y==5 | i.v.gps$L.x==5 |i.v.gps$L.y==6 | i.v.gps$L.x==6),]


#Replace NA values with zeros for sum
i.v.gps.456$PD_TCAP_TOT <- ifelse(is.na(i.v.gps.456$PD_TCAP_TOT),0,i.v.gps.456$PD_TCAP_TOT)
i.v.gps.456$IN_TCAP_TOT <- ifelse(is.na(i.v.gps.456$IN_TCAP_TOT),0,i.v.gps.456$IN_TCAP_TOT)

#sum inspecciones in districts 4,5,6
sum.insp <- i.v.gps.456$PD_TCAP_TOT + i.v.gps.456$IN_TCAP_TOT

#plot new dataset
plot(i.v.gps.456$X.y,i.v.gps.456$Y.y,col="slategray3",xaxt="n",yaxt="n",pch=18,cex=.4)
points(i.v.gps.456$X.y[which(i.v.gps.456$L.y==4)],i.v.gps.456$Y.y[which(i.v.gps.456$L.y==4)],pch=19,cex=.3,col="green")
points(i.v.gps.456$X.y[which(i.v.gps.456$L.y==5)],i.v.gps.456$Y.y[which(i.v.gps.456$L.y==5)],pch=19,cex=.3,col="purple")
points(i.v.gps.456$X.y[which(i.v.gps.456$L.y==6)],i.v.gps.456$Y.y[which(i.v.gps.456$L.y==6)],pch=19,cex=.3,col="orange")
for (i in 1:length(sum.insp)) points(i.v.gps.456$X.y[i][which(sum.insp[i]>0)],i.v.gps.456$Y.y[i][which(sum.insp[i]>0)],pch=18,cex=sum.insp[i]*.1+1,col="red")


##set up timer & basic functions
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

`%notin%` <- function(x,y) !(x %in% y) 

sampleWithoutSurprises <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } else {
    return(sample(x,1))
  }
}


##################################################
#######generate data##############################

N=dim(i.v.gps.456)[1]

###data simulation###
#NOTE: this simulation assumes Markov properties
#we can later extend to non-Markov chains
#simulation statistic vectors initialized
S.sim=1 #number of simulations
p=1:50/100
occult.sum.new<-rep(0,N)
infectiontime<-rep(Inf,N)
true.occult=total=neg=rep(NA,S.sim)
total.prob=true.pos=true.neg=matrix(NA,nrow=S.sim,ncol=length(p))
beta.sim=0
Rb.sim=0



#set up time scale
year<-i.v.gps.456$ANIO-2004

  #calculate distances between houses
  distance<-matrix(NA,nrow=N,ncol=N)
  for (i in 1:N){
    for (j in 1:N){
      distance[i,j]=sqrt((i.v.gps.456$X.y[i]-i.v.gps.456$X.y[j])^2+(i.v.gps.456$Y.y[i]-i.v.gps.456$Y.y[j])^2)
    }
  }
  
  
  K=1000 #carrying capacity
  maxt=15
  tobs=rep(maxt-1,N)
  
  beverton.holt<-function(id,K,R,bugs,trueremovaltime,trueinfectiontime){
    for(t in trueinfectiontime:(trueremovaltime-1)){
      bugs[id,(t+1)]=ceiling(R*bugs[id,t]/(1+bugs[id,t]/(K/(R-1))))
    }
    return(bugs[id,])
  }
  

T_b=30 #threshold for bug infectiousness
maxt=15 #by months, total of 6 years
tobs = rep(maxt,N) #observing after 4 years
tobs = ifelse(is.na(year), maxt,year)
bugs=matrix(0,nrow=N,ncol=maxt)
maxbugs<-max(sum.insp)
initialinfective<-which(sum.insp==maxbugs)
K=1000
Rb=1.1
id=1:N


beverton.holt<-function(id,K,R,bugs,trueremovaltime,trueinfectiontime){
  for(t in trueinfectiontime:(trueremovaltime-1)){
    bugs[id,(t+1)]=ceiling(R*bugs[id,t]/(1+bugs[id,t]/(K/(R-1))))
  }
  return(bugs[id,])
}

#probability of infestation differs by hops (<30m) or jumps (>30m)
threshold <- ifelse(distance<50,1,.01)

#find initial infectives notification and recovery times
infectiontime[initialinfective]<-1
bugs[initialinfective,1]<-1
bugs[initialinfective,]=rpois(maxt,beverton.holt(initialinfective,K,Rb,bugs,maxt,infectiontime[initialinfective]))

  
########################################################
#########MCMC algorithm###################################
#######################################################
tic()   #begin timer


##Jewell MCMC
M=4 #length of simulation
m=1 #first iteration
check3=rep(Inf,N)

check3<-ifelse(sum.insp>0,sum.insp,Inf)
I=ifelse(check3!=Inf&check3!=Inf,2,Inf)

trueremovaltime=ifelse(check3<Inf,101,Inf)
detectiontime=ifelse(check3>0&check3<Inf,tobs,Inf)
for(i in 1:N) tobs[i] <- ifelse(I[i]==Inf, maxt, year[i])


#initialize parameters 
beta=Rb=rep(0,M)
beta[1]=.3
betastar=.3
betastar.I=matrix(0,nrow=N,ncol=N)
betastar.sum=rep(0,N)
S=H.mat=matrix(0,nrow=N,ncol=N)
U=rep(0,N)
accept.beta=accept.Iadd=accept.Idel=rep(0,M)
K=1000

#keep track of occult infestations
occult=matrix(0,nrow=N,ncol=M)

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



#initialize Rb parameter
Rb[1]=1.11

#BH function to update infection times
beverton.holt.I<-function(update,K,R,check3,tobs){
  tobs1=tobs
  time=ifelse(update %in% N_N, -log((K-check3[update])/(check3[update]*K-check3[update]))/log(R[m-1]),-log((K-bugs[update,tobs1])/(bugs[update,tobs1]*K-bugs[update,tobs1]))/log(R[m-1]))
  I[update]=tobs1-time
  return(I[update])
}

#BH function to update bug counts given matrix
beverton.holt<-function(id,K,R,bugs,trueremovaltime,trueinfectiontime){
  for(t in trueinfectiontime:(maxt-1)){
    bugs[id,(t+1)]=ceiling(R*bugs[id,t]/(1+bugs[id,t]/(K/(R-1))))
  }
  return(bugs[id,])
}

#BH function to udpate bug counts given vector
beverton.holt.update<-function(K,R,bugs,trueremovaltime,trueinfectiontime){
  for(t in trueinfectiontime:(min(maxt-1,trueremovaltime))){
    bugs[(t+1)]=ceiling(R*bugs[t]/(1+bugs[t]/(K/(R-1))))
  }
  return(bugs)
}

#initialize bug counts
for (i in which(I!=Inf)){
  bugs[i,I[i]]=1
  bugstemp=beverton.holt(i,K[1],Rb[1],bugs,trueremovaltime,I[i])
  bugs[i,(I[i]:min(trueremovaltime[i],maxt))]=rpois((min(maxt,trueremovaltime[i])-I[i]+1),bugstemp[I[i]:min(maxt,trueremovaltime[i])])
  bugs[i,tobs[i]]=ifelse(check3[i]<Inf,check3[i],bugs[i,tobs[i]])
}


#first piece of likelihood
firstpiece<-function(I,beta,initialinfective,r){
  beta.I=H.mat=matrix(0,nrow=N,ncol=N)
  beta.sum=rep(0,N)
  for (j in 1:N) {
    for (i in 1:N) if(i %in% N_I | i==initialinfective) {
      if(I[i]<I[j]&I[j]<trueremovaltime[i]) H.mat[i,j]=ht(t, r, I, i, j, beta, K, threshold)
    }
    beta.sum[j]=sum(H.mat[,j])
  }
  beta.sum<-ifelse(is.na(beta.sum),0,beta.sum)
  
  return(beta.sum)
}


#second piece of likelihood
secondpiece<-function(trueremovaltime,detectiontime,I,beta,r){
  S1=matrix(0,nrow=N,ncol=N)
  for (i in 1:N){
    if(I[i]!=Inf) {
      for (j in 1:N){
        t=min(maxt,I[j],trueremovaltime[i])-min(I[i],I[j])
        if(t>0) S1[i,j]<-H(i, j, t, r, I, beta, K, threshold) #/sum(H.mat1[,j])
      }}
  }
  S1<-ifelse(S1=="NaN",0,S1)  
  doublesumofS=sum(S1)
  return(doublesumofS)
}


ht<-function(t, r, I, i, j, beta, K, threshold){
  n=1
  deriv=K*n^2*log(r)/((K-n)*exp(-log(r)*(t-I[i]))+n)^2
  hazard=ifelse(deriv>0,1-(1-beta*threshold[i,j])^(deriv),0)
  return(hazard) 
}


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
  for (i in 1:N) if(i %in% N_I | i==initialinfective) {
    if(I[i]<I[j]&I[j]<trueremovaltime[i]) H.mat[i,j]=ht(t, r, I, i, j, beta, K, threshold)
  }
  beta.sum=sum(H.mat[,j])
  tempthres1=N
  beta.sum<-ifelse(is.na(beta.sum),0,beta.sum)
  
  return(beta.sum)
}


#intialize bug mean vector
lambda_t=rep(0,N)
lambda_t[1]=1

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




ht <- function(t, r, I, i, j, beta, K, threshold) {
  n <- 1
  deriv <- K*n^2*log(r)/((K-n)*exp(-log(r)*(t-I[i]))+n)^2
  ifelse(deriv>0, 1-(1-beta*threshold[i,j])^(deriv), 0)
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

for (m in 2:M){
  #occult.sum=apply(occult,1,sum)
  #occult.sum.new=occult.sum/m
  #colfunc = gray.colors(length(unique(occult.sum.new)),start=1,end=0)[as.factor(occult.sum.new)]
  #plot(i.v.gps.456$X.y,i.v.gps.456$Y.y,col="gray",xlab="",ylab="")
  #for (i in 1:N) points(i.v.gps.456$X.y[i],i.v.gps.456$Y.y[i],pch=16,col=colfunc[i],xlab="X",ylab="Y",main = substitute(paste(Time, " = ", t),list(t=t)))
 #for (i in 1:N) if(sum.insp[i]>0) points(i.v.gps.456$X.y[i],i.v.gps.456$Y.y[i],pch=18,col="firebrick4") #,cex=check3[i]*.04+1)
  #legend("topright", inset=c(-0.8,0),c("Observed Infested","True Occult Infestation"),bty="n",col=c("firebrick4","firebrick1"),pch=c(18,1))
  
  ###############
  ##update beta##
  ###############
  
  betastar=abs(rnorm(1,beta[m-1],.01))
  logfirstpiecestar<-log(firstpiece.wrap(I, betastar, initialinfective, Rb[m-1], K, N, N_I, threshold) )
  logfirstpiecestar=ifelse(logfirstpiecestar=="-Inf",0,logfirstpiecestar)
  logfirstpiece<-log(firstpiece.wrap(I, beta[m-1], initialinfective, Rb[m-1], K, N, N_I, threshold))
  logfirstpiece=ifelse(logfirstpiece=="-Inf",0,logfirstpiece)
  dbetastar=sum(logfirstpiecestar)-secondpiece.wrap(I, betastar, Rb[m-1], K, N, maxt, threshold)+dunif(betastar,min=.001,max=.5,log=TRUE)
  dbeta=sum(logfirstpiece)-secondpiece.wrap(I, beta[m-1], Rb[m-1], K, N, maxt, threshold)+dunif(beta[m-1],min=.001,max=.5,log=TRUE)
  mstep.beta=min(1,exp(sum(dbetastar)-sum(dbeta)))
  if(mstep.beta=="NaN") mstep.beta=1
  R=runif(1)
  if(R<mstep.beta){
    beta[m]=betastar
    accept.beta[m]=1
  }else{
    beta[m]=beta[m-1]
    accept.beta[m]=0
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
  
  Rbstar=rnorm(1,Rb[m-1],.01)
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
    
    addinf<-which(I==Inf)
    if(length(addinf)>1){
      update=sample(addinf,1)
      points(i.v.gps.456$X.y[update],i.v.gps.456$Y.y[update],col="pink",pch=18)
      Istar[update]=floor(runif(1,min=2,max=maxt-2))
      trueremovaltime[update]=101
      detectiontime[update]=maxt
      tobs[update]=maxt
      bugsstar=rep(0,maxt)
      bugsstar[Istar[update]]=1
      bugsstar=beverton.holt.update(K,Rb[m],bugsstar,trueremovaltime[update],Istar[update])
      bugstest=bugs
      bugstest[update,]=bugsstar[1:maxt]
      check3[update]<-bugsstar[tobs[update]]
      logfirstpieceIstar<-log(firstpiece.update(update,Istar,beta[m], initialinfective,Rb[m]))
      logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
      loglike=sum(logfirstpieceIstar)-secondpiece.update(update,trueremovaltime,detectiontime,Istar,beta[m],Rb[m])
      if(loglike==0) loglike=-Inf
      extra.piece=(N-length(N_I)-1)/(length(N_I)-length(N_N)+1)*.00015 #*exp(sum(f_D.update(update,bugsstar,Istar,check3,Rb[m])))
      
      #metropolis hastings step for adding an infection
      mstep.I=min(1,exp(loglike)*extra.piece)
      if(mstep.I=="NaN") mstep.I=1
      R=runif(1)
      if(R<mstep.I){
        I<-Istar
        N_I<-c(N_I,update)
        infectedhousesI[N_I]=1
        accept.Iadd[m]=1
        bugs[update,]=bugsstar
        loglike.I<-loglike.Istar
        Q<-Qstar
        points(i.v.gps.456$X.y[update],i.v.gps.456$Y.y[update],col="green",pch=18)
      }else{
        Istar<-I
        accept.Iadd[m]=2
        trueremovaltime[update]=detectiontime[update]=tobs[update]=maxt
        check3[update]<-Inf}
    }
  }else{
    
    ###############
    ####delete I###
    ###############
    
    if(length(N_I)>length(N_N)){ 
      
      #pick which house to delete
      update=sampleWithoutSurprises(N_I[!(N_I %in% N_N)])
      points(i.v.gps.456$X.y[update],i.v.gps.456$Y.y[update],col="blue",pch=18)
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
      extra.piece=(length(N_I)-length(N_N)+1)/(N-length(N_I)-1)/.00015 #exp(sum(f_D.update(update,bugsstar,Istar,check3,Rb[m])))
      #decide whether to accept new I
      mstep.I=min(1,exp(loglike.Istar)*extra.piece)
      if(mstep.I=="NaN") mstep.I=1
      R=runif(1)
      if(R<mstep.I){
        I<-Istar
        N_I<-N_I[which(N_I!=update)]
        accept.Idel[m]=1
        loglike.I<-loglike.Istar
        infectedhousesI[update]=0
        Q<-Qstar
        bugs<-bugstest
        points(i.v.gps.456$X.y[update],i.v.gps.456$Y.y[update],col="purple",pch=18)
      }else{
        Istar<-I
        trueremovaltime[update]=101
        detectiontime[update]=tobs[update]=maxt
        accept.Idel[m]=2
        check3[update]=bugs[update,tobs[update]]}
      
      
    }
  }
 occult[N_I[!(N_I %in% N_N)],m]=1
 occult.sum <- apply(occult,1,sum)
 occult.prob<- occult.sum/m
 occult.prob.ids <- cbind(id, occult.prob, i.v.gps.456$X.y, i.v.gps.456$Y.y)
 
 occult.prob.ids <- occult.prob.ids[order(occult.prob, decreasing = TRUE),]
 colfunc2 = gray.colors(length(unique(occult.prob.ids[,2])),start=1,end=0)[as.factor(occult.prob.ids[,2])]
 par(mfrow=c(1,1))
 par(mar=c(1, 1, 1, 1), xpd=TRUE)
 plot(occult.prob.ids[,3], occult.prob.ids[,4],col = colfunc2,pch=16,cex=occult.prob.ids[,2]*20)
 top <- occult.prob.ids[1:10,]
 points(top[,3], top[,4],col = "blue")
 for (i in 1:N) if(sum.insp[i]>0) points(i.v.gps.456$X.y[i],i.v.gps.456$Y.y[i],pch=18,col="firebrick4",cex=.5) #,cex=check3[i]*.04+1)
 
 
  if(m%%1==0) {print(I) 
               print(m)
               print(N_I)
               print(Rb[m])
               print(beta[m])}

  #readline(prompt="Press [enter] to continue")
}
