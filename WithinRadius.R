unicode[which(sum.insp>0)]
dataset$id <- id 
idswithbugs <-id[which(sum.insp>0)]
threshold2 <- ifelse(distance<T_b, 1 , jumpprob)
abrevthreshold <- threshold2[idswithbugs,]
houseswithinradius <- NULL
allhouses <- NULL
for(i in 1:length(idswithbugs)){
  #find houses within radius 
  houseswithinradius <- which(threshold2[i,]==1) 
  #remove the house itself
  houseswithinradius <- houseswithinradius[houseswithinradius!=i]
  allhouses <- c(allhouses, houseswithinradius)
}
uninfested <- allhouses 
for(i in uninfested) {
  inspected[i]=1
}
allhouses <- write.table(allhouses, file = "/Users/patrickemedom/Desktop/Levy_Lab/Jewell/allhouses.aug24")



locatedidx<-c()
locatedidy<-c()
for(i in allhouses){
  locatedidx <- c(locatedidx, dataset[i, c(69)])
  locatedidy <- c(locatedidy, dataset[i, c(70)])
  #homes <- c(homes, which(dataset$id==i))
}
locationid <- data.frame(allhouses, locatedidx, locatedidy)

plot((locationid$locatedidx), (locationid$locatedidy), pch=16, axes=FALSE, xlab="",ylab="")
for (i in 1:N) if(sum.insp[i]>0) points(dataset$X[i],dataset$Y[i],pch=18,col="firebrick3",cex=1.5)
legend("bottomleft",c("Known Infested House"),pch=c(18,18),col=c("firebrick3"),bty="n", cex = 0.75)
