###not yet working###
set.seed(123456789)

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
###data simulation###
#NOTE: this simulation assumes Markov properties
#we can later extend to non-Markov chains
#simulation statistic vectors initialized
S.sim=1 #number of simulations
p=1:50/100
true.occult=total=neg=rep(NA,S.sim)
total.prob=true.pos=true.neg=matrix(NA,nrow=S.sim,ncol=length(p))
beta.sim=0
Rb.sim=0

#for (s in 1:S.sim){
totalnuminf=1
N=100
indicator=indicator2=0
while(sum(indicator)<5 | indicator2<2){
  
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
  beta=.02
  Rb=1.2
  K=1000 #carrying capacity
  maxt=100
  tobs=75
  
  beverton.holt<-function(id,K,R,bugs,trueremovaltime,trueinfectiontime){
    for(t in trueinfectiontime:(trueremovaltime-1)){
      bugs[id,(t+1)]=ceiling(R*bugs[id,t]/(1+bugs[id,t]/(K/(R-1))))
    }
    return(bugs[id,])
  }
  
  S=I=R=matrix(NA, nrow=max(id),ncol=N*N) 
  bugs=bugs.Rb=matrix(0,nrow=N,ncol=maxt)
  
  #probability of infestation differs by hops (<.3) or jumps (>.3)
  threshold=ifelse(distance<.1,1,.02)
  
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
  h<-function(t,r,I,i,j,beta){
    #n=bugs[i,t]
    n=1
    #N.prime=ifelse(t-I[i]>0,Rb^(t-I[i]),0)
    deriv=K*n^2*log(r)/((K-n)*exp(-log(r)*(t-I[i]))+n)^2
    hazard=ifelse(deriv>0,1-(1-beta*threshold[i,j])^(deriv),0)
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
      bug.ind[i,j]=ifelse(t>=infectiontime[i]&infectiontime[i]>0,h(t,Rb,infectiontime,i,j,beta),0)
      }}
    for(j in 1:N){
      ifelse(S[j,(t-1)]==1, StoI[j]<-rbinom(1,S[j,(t-1)],sum(bug.ind[,j])/N),StoI[j]<-0)
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



########################################################
#########MCMC algorithm###################################
#######################################################
tic()   #begin timer


##Jewell MCMC
M=100000 #length of simulation
tobs=75 #when observation occurs
m=1 #first iteration

#update data to reflect only OBSERVED data
check3=ifelse(infectiontime<=tobs&infectiontime!=0,bugs[,tobs], Inf)
trueremovaltime=ifelse(check3<Inf,101,Inf)
detectiontime=ifelse(check3>0&check3<Inf,tobs,Inf)
initialinfective=which(infectiontime==1)


#initialize parameters 
#threshold=ifelse(distance<.2,1,.01)
beta=Rb=rep(0,M)
beta[1]=.1
betastar=.1
betastar.I=matrix(0,nrow=N,ncol=N)
betastar.sum=rep(0,N)
I=rep(0,N)
S=H.mat=matrix(0,nrow=N,ncol=N)
U=rep(0,N)
accept.beta=accept.Iadd=accept.Idel=rep(0,M)
K=1000

#keep track of occult infestations
occult=matrix(0,nrow=N,ncol=M)

#initialize infection times
I=ifelse(check3!=Inf&check3!=0,2,Inf)
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
Rb[1]=1.15

#BH function to update infection times
beverton.holt.I<-function(update,K,R,check3,tobs){
  time=-log((K-check3[update])/(check3[update]*K-check3[update]))/log(R[m-1])
  I[update]=tobs-time
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
  for(t in trueinfectiontime:(maxt-1)){
    bugs[(t+1)]=ceiling(R*bugs[t]/(1+bugs[t]/(K/(R-1))))
  }
  return(bugs)
}

#initialize bug counts
for (i in which(I!=Inf)){
  bugs[i,I[i]]=1
  bugstemp=beverton.holt(i,K[1],Rb[1],bugs,maxt,I[i])
  bugs[i,(I[i]:maxt)]=rpois((maxt-I[i]+1),bugstemp[I[i]:maxt])
  bugs[i,tobs]=ifelse(check3[i]<Inf,check3[i],bugs[i,tobs])
}


#first piece of likelihood
firstpiece<-function(I,beta,initialinfective,r){
  beta.I=H.mat=H.mat1=matrix(0,nrow=N,ncol=N)
  beta.sum=rep(0,N)
  for (j in 1:N) {
    for (i in 1:N) if(i %in% N_I | i==initialinfective) {
      if(I[i]<I[j]&I[j]<Inf) H.mat[i,j]=ht(t,r,I,i,j,beta)
    }}
  tempthres1=N
  for (j in 1:N) if(j %in% N_I | j==initialinfective) {
    for (i in 1:N){
      if(I[i]<I[j] & I[j]<=maxt) {beta.I[i,j]=H.mat[i,j]/sum(H.mat[,j])}
    }
    beta.sum[j]=sum(beta.I[j,])
  }
  
  beta.sum<-ifelse(is.na(beta.sum),0,beta.sum)
  
  return(beta.sum)
}


#second piece of likelihood
secondpiece<-function(trueremovaltime,detectiontime,I,beta,r){
  S1=S2=H.mat1=H.mat2=matrix(0,nrow=N,ncol=N)
  for (j in 1:N) {
    for (i in 1:N){
      t=min(maxt,I[j])-min(I[i],I[j])
      H.mat1[i,j]=ifelse(t>0,H(i,j,t,r,I,beta),0)
    }}
  tempthres1=H.mat1
  for (i in 1:N){
    if(I[i]!=Inf) {
      for (j in 1:N){
        t=min(maxt,I[j])-min(I[i],I[j])
        if(t>0) S1[i,j]<-H(i,j,t,r,I,beta)/sum(H.mat1[i,])
      }}
  }
  S1<-ifelse(S1=="NaN",0,S1)  
  doublesumofS=sum(S1)
  return(doublesumofS)
}


ht<-function(t,r,I,i,j,beta){
  n=1
  deriv=K*n^2*log(r)/((K-n)*exp(-log(r)*(t-I[i]))+n)^2
  hazard=ifelse(deriv>0,1-(1-beta*threshold[i,j])^(deriv),0)
  return(hazard) 
}


H<-function(i,j,tt,Rb,I,beta){
  r=log(Rb)
  t=tt-I[i]
  Ht=t*(1-(1-beta*threshold[i,j])^(r/K))-
    ((K-1)*r^2*t^2*log(1-beta*threshold[i,j])*(1-beta*threshold[i,j])^(r/K))/(K^2)-
    t^3*((K-1)*r^3*log(1-beta*threshold[i,j])*(1-beta*threshold[i,j])^(r/K)*(2*(K-1)*r*
                log(1-beta*threshold[i,j])+K*(2*K-3)))/(3*K^4)
  return(Ht)}

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

##############################
#########for loop begins#######
##############################

for (m in 18586:M){
  
  ###############
  ##update beta##
  ###############
  
  betastar=abs(rnorm(1,beta[m-1],.01))
  logfirstpiecestar<-log(firstpiece(I,betastar, initialinfective,Rb[m-1]))
  logfirstpiecestar=ifelse(logfirstpiecestar=="-Inf",0,logfirstpiecestar)
  logfirstpiece<-log(firstpiece(I,beta[m-1], initialinfective,Rb[m-1]))
  logfirstpiece=ifelse(logfirstpiece=="-Inf",0,logfirstpiece)
  dbetastar=sum(logfirstpiecestar)-secondpiece(trueremovaltime,detectiontime,I,betastar,Rb[m-1])+dgamma(betastar,.1,scale=.1,log=TRUE)
  dbeta=sum(logfirstpiece)-secondpiece(trueremovaltime,detectiontime,I,beta[m-1],Rb[m-1])+dgamma(beta[m-1],.1,scale=.1,log=TRUE)
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
  if(bugs[update,maxt]==0) bugs[update,maxt]=1
  Istar[update]=ifelse(check3[update]>0,round(rnorm(1,beverton.holt.I(update,K,Rb,check3,tobs),1)),floor(rnorm(1,beverton.holt.I(update,K,Rb,bugs[,maxt], maxt),1)))
  if(Istar[update]=="NaN") {bugs[update,maxt]=1000
                            Istar[update]=ifelse(check3[update]>0,round(beverton.holt.I(update,K,Rb,check3,tobs)),floor(beverton.holt.I(update,K,Rb,bugs[,maxt], maxt)))
  }
  if(Istar[update]==maxt) Istar[update]=maxt-1
  Istar[update]=ifelse(Istar[update]<1,2,Istar[update])
  bugsstar=rep(0,maxt)
  bugsstar[Istar[update]]=1
  bugsstar=beverton.holt.update(K,Rb[m-1],bugsstar,maxt,Istar[update])
  bugsstar[tobs]=check3[update]
  logfirstpieceI<-log(firstpiece(I,beta[m], initialinfective,Rb[m-1]))
  logfirstpieceI=ifelse(logfirstpieceI=="-Inf",0,logfirstpieceI)
  loglike.I=sum(logfirstpieceI)-secondpiece(trueremovaltime,detectiontime,I,beta[m],Rb[m-1])
  logfirstpieceIstar<-log(firstpiece(Istar,beta[m], initialinfective,Rb[m-1]))
  logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
  loglike.Istar=sum(logfirstpieceIstar)-secondpiece(trueremovaltime,detectiontime,Istar,beta[m],Rb[m-1])
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
  
  Rbstar=rnorm(1,Rb[m-1],.003)
  Q=Qstar=rep(NA,length(I!=Inf))
  
  for (i in which(I!=Inf)){
    Qstar[i]=sum(f_D(i,bugs,I,check3,Rbstar))
    Q[i]=sum(f_D(i,bugs,I,check3,Rb[m-1]))
  }
  thirdpieceloglike=sum(Qstar[which(Qstar!="NA")])-sum(Q[which(Q!="NA")])+dgamma(Rbstar,shape=1.1,scale=1,log=TRUE)-dgamma(Rb[m-1],shape=1.1,scale=1,log=TRUE)
  logfirstpieceIstar<-log(firstpiece(I,beta[m], initialinfective,Rbstar))
  logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
  loglike.Istar=sum(logfirstpieceIstar)-secondpiece(trueremovaltime,detectiontime,I,beta[m],Rbstar)
  Rbloglike=loglike.Istar+thirdpieceloglike-loglike.I
  
  #Metropolis step
  mstep.Rb=min(0,Rbloglike)
  if(mstep.Rb=="NaN") mstep.Rb=0
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
      Istar[update]=floor(runif(1,min=tobs,max=maxt-2))
      trueremovaltime[update]=101
      check3[update]=0
      bugsstar=rep(0,maxt)
      bugsstar[Istar[update]]=1
      bugsstar=beverton.holt.update(K,Rb[m],bugsstar,trueremovaltime[update],Istar[update])
      bugstest=bugs
      bugstest[update,]=bugsstar[1:maxt]
      #update log likelihood with new infection
      for (i in which(I!=Inf)){
        Qstar[i]=sum(f_D(i,bugstest,Istar,check3,Rb[m]))
      }
      thirdpieceloglike=sum(Qstar[which(Qstar!="NA")])-sum(Q[which(Q!="NA")])
      logfirstpieceIstar<-log(firstpiece(Istar,beta[m], initialinfective,Rb[m]))
      logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
      loglike.Istar=sum(logfirstpieceIstar)-secondpiece(trueremovaltime,detectiontime,Istar,beta[m],Rb[m])
      loglike=loglike.Istar+thirdpieceloglike-loglike.I
      extra.piece=(N-length(N_I)-1)/(length(N_I)-length(N_N)+1)*dunif(Istar[update],min=tobs,max=maxt-2) #*exp(sum(f_D.update(update,bugsstar,Istar,check3,Rb[m])))

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
      }else{
        Istar<-I
        accept.Iadd[m]=0
        trueremovaltime[update]=Inf
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
      for (i in which(I!=Inf)){
        Qstar[i]=sum(f_D(i,bugs,Istar,check3,Rb[m]))
      }
      thirdpieceloglike=sum(Qstar[which(Qstar!="NA")])-sum(Q[which(Q!="NA")])
      logfirstpieceIstar<-log(firstpiece(Istar,beta[m], initialinfective,Rb[m]))
      logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
      loglike.Istar=sum(logfirstpieceIstar)-secondpiece(trueremovaltime,detectiontime,Istar,beta[m],Rb[m])+thirdpieceloglike
      extra.piece=(length(N_I)-length(N_N))/(N-length(N_I))*dunif(I[update],min=tobs,max=maxt-2) #*exp(sum(f_D.update(update,bugsstar,Istar,check3,Rb[m])))
      #decide whether to accept new I
      mstep.I=min(1,exp(loglike.Istar-loglike.I)*extra.piece)
      if(mstep.I=="NaN") mstep.I=1
      R=runif(1)
      if(R<mstep.I){
        I<-Istar
        N_I<-N_I[which(N_I!=update)]
        accept.Idel[m]=1
        loglike.I<-loglike.Istar
        infectedhousesI[update]=0
        trueremovaltime[update]=Inf
        Q<-Qstar
      }else{
        Istar<-I
        accept.Idel[m]=0
        check3[update]=0}
      
      
    }
  }
  if(m%%1==0) {print(I) 
               print(m)
               print(N_I)
               print(Rb[m])}
  occult[N_I[!(N_I %in% N_N)],m]=1
  #readline(prompt="Press [enter] to continue")
}

occult.sum=apply(occult,1,sum)
occult.sum.new=occult.sum/m

#simulation statistics
total[s]=N-length(N_N)-1
true.occult[s]=length(infectiontime[which(infectiontime>tobs)])
neg[s]=length(infectiontime[which(infectiontime==0)])

p=1:10/100
for(i in seq(along=p)){
  total.prob[s,i]=length(occult.sum.new[which(occult.sum.new>p[i])])
  true.pos[s,i]=length(occult.sum.new[which(occult.sum.new>p[i]&infectiontime>tobs)])
  true.neg[s,i]=length(occult.sum.new[which(occult.sum.new<=p[i]&infectiontime<=tobs)])
}
beta.sim=c(beta.sim,beta)
Rb.sim=c(Rb.sim, Rb)
print(s)
#}
toc()
bugsize=bugs[,100]/100
colfunc = gray.colors(length(unique(occult.sum.new)),start=1,end=0)[as.factor(occult.sum.new)]
par(mfrow=c(1,1))
par(mar=c(4.1, 3.1, 3.1, 13), xpd=TRUE)
plot(location[,1],location[,2],col="gray",xlab="",ylab="")
for (i in 1:N) points(location[i,1],location[i,2],pch=19,col=colfunc[i],xlab="X",ylab="Y",main = substitute(paste(Time, " = ", t),list(t=t)))
for (i in 1:N) if(infectiontime[i]>=tobs) points(location[i,1],location[i,2],pch=1,col="firebrick1")
for (i in 1:N) if(infectiontime[i]>=1&infectiontime[i]<tobs) points(location[i,1],location[i,2],pch=18,col="firebrick4",cex=bugsize[i]/3)
legend("topright", inset=c(-0.7,0),c("Observed Infested","True Occult Infestation"),bty="n",col=c("firebrick4","firebrick1"),pch=c(18,1))

