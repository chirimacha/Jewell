###data simulation###
#NOTE: this simulation assumes Markov properties
#we can later extend to non-Markov chains

id<-c(1:500) #simulate 500 houses

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
beta=.5
#gamma=0.15
shape=10
scale=.7
r=1 #notification parameter

#stochastic.seasonal.sir<-function(id, distance, beta, gamma){
N=max(id)
S=I=Not=R=matrix(NA, nrow=max(id),ncol=10000) 

#probability of infestation differs by hops (<.3) or jumps (>.3)
threshold=ifelse(distance<.3,1,.1)
W=beta*threshold
  
  #initial state vectors
  S[,1]<-rep(1,N)
  I[,1]<-rep(0,N)
  Not[,1]<-rep(0,N)
  R[,1]<-rep(0,N)
  
  #initial infective
  initialinfective<-sample(1:N,1)
  S[initialinfective,1]=0
  I[initialinfective,1]=1

  
  #keeping counts of S, E, and I at all the time points of the simulation
StoI=ItoN=NtoR=rep(0,N)
infectiontime=rep(0,N)
notificationtime=rep(0,N)
removaltime=rep(0,N)


#find initial infectives notification and recovery times
  infectiontime[initialinfective]<-1
  notificationtime[initialinfective]<-round(rexp(1,r))+infectiontime[initialinfective]
  removaltime[initialinfective]<-notificationtime[initialinfective]+ceiling(rexp(1,.5))


  t=2
 while(sum(Not[,(t-1)]+I[,(t-1)])>0){   #looping through each day of the simulation
    for(i in 1:N){
      ifelse(S[i,(t-1)]==1, StoI[i]<-rbinom(1,S[i,(t-1)],beta*sum(threshold[i,]/sum(threshold[i,])*(I[,(t-1)]+Not[,(t-1)]))),StoI[i]<-0)
      if(StoI[i]==1){
        infectiontime[i]<-t
        notificationtime[i]<-round(rexp(1,r))+infectiontime[i]
        removaltime[i]<-notificationtime[i]+ceiling(rexp(1,.5))
      }
      ifelse(t==notificationtime[i], ItoN[i]<-1,ItoN[i]<-0)
      ifelse(t==removaltime[i], NtoR[i]<-1, NtoR[i]<-0)
    } 
  
      
    #Now update the daily counts
    S[,t]=S[,(t-1)]-StoI
    I[,t]=I[,(t-1)]+StoI-ItoN
    Not[,t]=Not[,(t-1)]+ItoN-NtoR
    R[,t]=R[,(t-1)]+NtoR
    
    t=t+1
  }  
maxt=t-1
maxt
  S=S[,1:(maxt)]
  I=I[,1:(maxt)]
  Not=Not[,1:(maxt)]
  R=R[,1:(maxt)]
#  return(list(S,I,R))
#}

#look at distribution at a few times
par(mfrow=c(2,2))
t=1
plot(location[,1],location[,2],col="gray",main = substitute(paste(Time, " = ", t),list(t=t)),xlab="X",ylab="Y")
for (i in 1:N) if(I[i,t]==1 | Not[i,t]==1) points(location[i,1],location[i,2],pch=19,col="red")
t=round(maxt/3)
plot(location[,1],location[,2],col="gray",xlab="X",ylab="Y",main = substitute(paste(Time, " = ", t),list(t=t)))

for (i in 1:N) if(I[i,t]==1 | Not[i,t]==1) points(location[i,1],location[i,2],pch=19,col="red")
t=round(maxt/3*2)
plot(location[,1],location[,2],col="gray",xlab="X",ylab="Y",main = substitute(paste(Time, " = ", t),list(t=t)))

for (i in 1:N) if(I[i,t]==1 | Not[i,t]==1) points(location[i,1],location[i,2],pch=19,col="red")
t=maxt
plot(location[,1],location[,2],col="gray",xlab="X",ylab="Y",main = substitute(paste(Time, " = ", t),list(t=t)))

for (i in 1:N) if(I[i,t]==1 | Not[i,t]==1) points(location[i,1],location[i,2],pch=19,col="red",xlab="X",ylab="Y",main = substitute(paste(Time, " = ", t),list(t=t)))


#count total number infected
totalnuminf=sum(S[,1])-sum(S[,maxt])

#make dataset
data=cbind(infectiontime, notificationtime, removaltime)
