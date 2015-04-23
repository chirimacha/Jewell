###for debugging purposes###
#assumes all known infection,notification,removal times
#just trying to estimate beta and r parameters


firstpiece<-function(I,beta,initialinfective){
  beta.I=matrix(0,nrow=N,ncol=N)
  beta.sum=rep(0,)
  for (j in 1:N) {
    for (i in 1:N){
      if(I[i]<I[j] & I[j]<=trueremovaltime[i]) {beta.I[i,j]=beta*threshold[i,j]/sum(threshold[j,])}
    }
    beta.sum[j]=sum(beta.I[,j])
  }
  return(beta.sum)
}



secondpiece<-function(trueremovaltime,truenotificationtime,I,beta){
  S1=S2=matrix(0,nrow=N,ncol=N)
  for (i in 1:N){
    if(trueinfectiontime[i]!=Inf) {
      for (j in 1:N){
        S1[i,j]<-(min(truenotificationtime[i],I[j])-min(I[i],I[j]))*beta*threshold[i,j]/sum(threshold[i,])
      }}
  }
  for (i in 1:N){
    if(truenotificationtime[i]!=Inf) {
      for (j in 1:N){
        S2[i,j]<-(min(trueremovaltime[i],I[j])-min(truenotificationtime[i], I[j]))*beta*threshold[i,j]/sum(threshold[i,])
      }}
  }
  S1<-ifelse(S1=="NaN",0,S1)  
  S2<-ifelse(S2=="NaN",0,S2) 
  doublesumofS=sum(S1)+sum(S2)
  return(doublesumofS)
}

mleest<-function(par, data, threshold){
  beta=par[1]
  r=par[2]
  trueinfectiontime=data[,1]
  truenotificationtime=data[,2]
  trueremovaltime=data[,3]
  trueinfectiontime=ifelse(trueinfectiontime==0,Inf,trueinfectiontime)
  truenotificationtime=ifelse(truenotificationtime==0,Inf,truenotificationtime)
  trueremovaltime=ifelse(trueremovaltime==0,Inf,trueremovaltime)
  D=truenotificationtime-trueinfectiontime
  D=D[which(D!="NaN")]
  piece1<-log(firstpiece(trueinfectiontime, beta,initialinfective))
  piece1[initialinfective]<-0
  piece1<-ifelse(trueinfectiontime!=Inf,piece1,0)
  piece1<-ifelse(piece1=="NaN",0,piece1)
  likelihood<- sum(piece1)-secondpiece(trueremovaltime, truenotificationtime, trueinfectiontime, beta)+sum(dexp(D, r, log=TRUE))
  return(-likelihood)}

optim(par=c(1,1),mleest, data=data,method="Nelder-Mead")

#optimize(mleest, data=data, interval=c(0,2))

