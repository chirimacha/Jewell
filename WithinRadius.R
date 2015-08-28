#define parameter (we wont actually run algorithm)
totaliterations<- 2

#set your working directory to where data is
setwd("/Users/EMWB/Jewell/Data")


#before running this code, run lines 12-334 of RJMCMC1Aug10Rtimes.smallerdata.R
unicode[which(sum.insp>0)]
dataset$id <- id 
idswithbugs <-id[which(sum.insp>0)]
threshold2 <- ifelse(distance<T_b, 1 , jumpprob)
abrevthreshold <- threshold2[idswithbugs,]
houseswithinradius <- NULL
allhouses <- NULL
for(i in 1:length(idswithbugs)){
  #find houses within radius 
  houseswithinradius <- which(abrevthreshold[i,]==1) 
  allhouses<- c(allhouses,houseswithinradius)

}
uninfested <- allhouses 
for(i in uninfested) {
  inspected[i]=1
}
allhouses <- write.table(allhouses, file = "/Users/patrickemedom/Desktop/Levy_Lab/Jewell/allhouses.aug24")


plot(dataset$X,dataset$Y)
points(dataset$X[allhouses],dataset$Y[allhouses],col="red")
for (i in 1:N) if(sum.insp[i]>0) points(dataset$X[i],dataset$Y[i],pch=18,col="firebrick3",cex=1.5)
legend("bottomleft",c("Known Infested House"),pch=c(18,18),col=c("firebrick3"),bty="n", cex = 0.75)
