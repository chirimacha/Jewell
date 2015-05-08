#plots manzanas by  prediction from Barbu 2015

setwd(paste(Sys.getenv("SPATIAL_UNCERTAINTY"), "/code", sep="")) 

d<-read.csv("byManz2015-03-10__18-42-57763.csv")

X<-d$X
Y<-d$Y
DenPredicted<-d$DenPredicted

#Mapa completa
plot(X,Y,cex=DenPredicted)

#Mapa sin puntos lejanos
sel<-which(d$X>222000)
plot(X[sel],Y[sel],cex=DenPredicted[sel]
     
#Mapa por Districto
mapD<-function(D) {
  sel=which(d$D==D)
  plot(X[sel],Y[sel],cex=DenPredicted[sel]*5)
  }

#tiabaya
mapD(24)
#Sachaca
mapD(18)




