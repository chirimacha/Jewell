###not yet working###


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


tic()
##Jewell MCMC
M=10000 #length of simulation
tobs=maxt #when observation occurs

#for simulation, mark true values
trueinfectiontime=data[,1]
truenotificationtime=data[,2]
trueremovaltime=data[,3]
trueinfectiontime=ifelse(trueinfectiontime==0,Inf,trueinfectiontime)
truenotificationtime=ifelse(truenotificationtime==0,Inf,truenotificationtime)
trueremovaltime=ifelse(trueremovaltime==0,Inf,trueremovaltime)
initialinfective=which(trueinfectiontime==1)


#initialize parameters 
threshold=ifelse(distance<.3,1,.1)
beta=rep(0,M)
beta[1]=1
betastar=1
betastar.I=matrix(0,nrow=N,ncol=N)
betastar.sum=rep(0,N)
I=rep(0,N)
W=beta*threshold
S=matrix(0,nrow=N,ncol=N)
U=rep(0,N)
accept.beta=accept.r=rep(0,M)

#define indicators for infected houses
temp=which(truenotificationtime!=Inf)
#take out initial infective
temp=temp[which(temp!=initialinfective)]
infectedhouses=rep(0,N)
infectedhouses[temp]=1


mu = .25 #non-centering parameter

#initialize infection times
I=ifelse(truenotificationtime!=Inf,2,Inf)
I[initialinfective]=1
Istar=I
D=truenotificationtime-I
D=ifelse(D=="NaN",0,D)

firstpiece<-function(I,beta,initialinfective){
  beta.I=matrix(0,nrow=N,ncol=N)
  beta.sum=NULL
  for (j in 1:N) {
    for (i in 1:N){
      if(I[i]<I[j] & I[j]<=trueremovaltime[i]) {beta.I[i,j]=beta*threshold[i,j]/sum(threshold[j,])}
    }
    beta.sum[j]=sum(beta.I[,j])
  }
  return(beta.sum)
}



secondpiece<-function(tobs,trueremovaltime,truenotificationtime,I,beta){
  S1=S2=matrix(0,nrow=N,ncol=N)
  for (i in 1:N){
    if(I[i]!=Inf) {
      for (j in 1:N){
        S1[i,j]<-(min(tobs,truenotificationtime[i],I[j])-min(I[i],I[j]))*beta*threshold[i,j]/sum(threshold[i,])
      }}
  }
  for (i in 1:N){
    if(truenotificationtime[i]!=Inf) {
      for (j in 1:N){
        S2[i,j]<-(min(tobs,trueremovaltime[i],I[j])-min(truenotificationtime[i], I[j]))*beta*threshold[i,j]/sum(threshold[i,])
      }}
  }
  S1<-ifelse(S1=="NaN",0,S1)  
  S2<-ifelse(S2=="NaN",0,S2) 
  doublesumofS=sum(S1)+sum(S2)
  return(doublesumofS)
}
#initialize r
r=rep(0,M)
r[1]=1
m=1
W=beta[1]*threshold


for (m in 2:M){
 # for(p in 1:5){
    logfirstpieceI<-log(firstpiece(I,beta[m-1],initialinfective))
    logfirstpieceI=ifelse(logfirstpieceI=="-Inf",0,logfirstpieceI)
    loglike.I=sum(logfirstpieceI)-secondpiece(tobs,trueremovaltime,truenotificationtime,I,beta[m-1])+sum(dexp(D,r[m-1],log=TRUE))
    
##update I
update=sample(temp,1)
Istar[update]=truenotificationtime[update]-rexp(1,r[m-1])
Dstar=truenotificationtime-Istar
Dstar=ifelse(Dstar=="NaN",0,Dstar)

logfirstpieceIstar<-log(firstpiece(Istar,beta[m-1],initialinfective))
logfirstpieceIstar=ifelse(logfirstpieceIstar=="-Inf",0,logfirstpieceIstar)
loglike.Istar=sum(logfirstpieceIstar)-secondpiece(tobs,trueremovaltime,truenotificationtime,Istar,beta[m-1])+sum(dexp(Dstar,r[m-1],log=TRUE))
  
#decide whether to accept new I
mstep.I=min(1,exp(loglike.Istar-loglike.I))
if(mstep.I=="NaN") mstep.I=1

R=runif(1)
if(R<mstep.I){
  I<-Istar
  D<-Dstar
}else{Istar<-I}
#}

##update r
rstar=abs(rnorm(1,r[m-1],.1))

drstar<-dexp(rstar,D,log=TRUE)
drstar<-ifelse(drstar=="-Inf",0,drstar)
dr<-dexp(r[m-1],D,log=TRUE)
dr<-ifelse(dr=="-Inf",0,dr)
#decide whether to accept new r
mstep.r=min(1,exp(sum(drstar)-sum(dr)))
R=runif(1)
if(R<mstep.r){
  r[m]<-rstar
  accept.r[m]=1
}else{
  r[m]=r[m-1]
  accept.r[m]=0
}


##update Z, and hence, A and B
z=rbinom(totalnuminf,1,mu)
#get which indices are in group A and B
A=sample(which(infectedhouses==1),sum(z))
U[A]=r[m]*D[A]
B=ifelse(U==0&infectedhouses==1,1,0)

#update r given z (this will also update I if in group A)
rstar=abs(rnorm(1,r[m-1],.1))
drstar<-dexp(rstar,D,log=TRUE)
drstar<-ifelse(drstar=="-Inf",0,drstar)
dr<-dexp(r[m-1],D,log=TRUE)
dr<-ifelse(dr=="-Inf",0,dr)

#decide whether to accept new r using partial non-centering
#group A is non-centered; these infection times are also updated in this step
mstep.r2=min(1,exp(sum(drstar[B==1])-sum(dr[B==1])))
R=runif(1)
if(R<mstep.r2){
  r[m]=rstar
  I[A]= truenotificationtime[A]-1/rstar*U[A]
  D[A]=truenotificationtime[A]-I[A]
}


##update theta=c(beta)
betastar=abs(rnorm(1,beta[m-1],.1))
Wstar=betastar*threshold
W=beta[m-1]*threshold
logfirstpiecestar<-log(firstpiece(I,betastar,initialinfective))
logfirstpiecestar=ifelse(logfirstpiecestar=="-Inf",0,logfirstpiecestar)
logfirstpiece<-log(firstpiece(I,beta[m-1],initialinfective))
logfirstpiece=ifelse(logfirstpiece=="-Inf",0,logfirstpiece)
dbetastar=sum(logfirstpiecestar)-secondpiece(tobs,trueremovaltime,truenotificationtime,I,betastar)
dbeta=sum(logfirstpiece)-secondpiece(tobs,trueremovaltime,truenotificationtime,I,beta[m-1])
mstep.beta=min(1,exp(sum(dbetastar)-sum(dbeta)))
if(mstep.beta=="NaN") mstep.beta=1
R=runif(1)
if(R<mstep.beta){
  beta[m]=betastar
  accept.beta[m]=1
  W=betastar*threshold
}else{
  beta[m]=beta[m-1]
  accept.beta[m]=0
}
print(m)
print(c(accept.beta[m], accept.r[m],beta[m],r[m]))
if(m%%10==0) print(I)
}

toc()
