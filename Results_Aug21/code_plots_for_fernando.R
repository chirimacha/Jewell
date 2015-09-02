
#after running the main code to create dataset
plot(dataset$X,dataset$Y,cex=10*dataset$predicteddensity)
points(dataset$X,dataset$Y,cex=log(sum.insp+1),col="red",pch=19)

setwd("~/Jewell/Results_Aug21/")

Results2<- read.csv("Results2Aug21.csv")
Results21<- read.csv("Results21Aug21.csv")
Results22<- read.csv("Results22Aug21.csv")
Results23<- read.csv("Results23Aug21.csv")

#set whichever dataset you want to plot as Results1
Results1 <- Results21
colfunc = gray.colors(length(unique(Results1[,3])),start=1,end=0)[as.factor(Results1[,3])]

#plot priors
plot(dataset$X,dataset$Y,cex=15*dataset$predicteddensity)

#plot denuncia scaled by size
points(dataset$X,dataset$Y,cex=log(sum.insp+1),col="red",pch=16)

#plot posteriors
plot(Results1[,4], Results1[,5],col = colfunc,pch=16,cex=1.2)
points(dataset$X[sum.insp>=1],dataset$Y[sum.insp>=1], pch="x",col="black")


#points(dataset$X,dataset$Y,cex=log(sum.insp+1),col="red",pch=16)


points(dataset$X,dataset$Y,cex=log(sum.insp+1),col="red",pch=16)
