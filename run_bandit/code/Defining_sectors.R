'
######################################################
Antes de correr la funciona "ClusterData" 
Necesitamos la base (o cualquier otro distrito) en este formato:
(no necesitamos eso nombres)
Las functionas CleanAggDataDist y CleanAggDataLoc puede hacerlo 
por el districto y la localidad
/
DIVIDE.VAR  uni   #/uni #/divide.var
  7         701   30     149
  7         702   22     149
  7         703   37     149
  7         704   23     149
  7         705   4      149
  7         706   25     149

Divide.var: el lugar que quiere dividir (distrito,localidad,etc.)
uni: el uni de los "bloques" (localidad, mananzas, etc.)


/
########################################################'


### Install library
library(classInt)


#Set directories to spatial data path
if(Sys.getenv("SPATIAL_DATA") != "") {
  setwd(paste(Sys.getenv("SPATIAL_DATA"), sep=""))  
} else if(Sys.getenv("USER") == "snutman") {
  setwd("/Users/snutman/Documents/CHAGAS_DATA/spatial_surveillence/Spatial_data")
} else if(Sys.getenv("USER") == "USER") {
  setwd("PATH_TO_SPATIAL_DATA")
}

#prepare for output path: 

#TODO: create a variable with the output path or create an environment variable SPATIAL_UNCERTAINTY

if(Sys.getenv("SPATIAL_UNCERTAINTY") != "") {
  output <-Sys.getenv("SPATIAL_UNCERTAINTY")
} else if(Sys.getenv("USER") == "sgutfraind") {
  output <- "/Users/sgutfraind/academic/chagas/bandit"
} else if(Sys.getenv("USER") == "mlevy") {
  #TODO: check your USER string (your username on your computer) and add appropriate directory here with 
  output<-"PATH_TO_BANDIT"
} 

#Function for summary of data by locality
CleanAggDataDist <- function(dataset) {
  
  block<-data.frame(table(dataset[,c("dist","loc")]))
  colnames(block)=c("dist","loc","locsize")
  district<-data.frame(table(dataset[c("dist")]))
  colnames(district)=c("dist","distsize")
  block<-merge(block,district,by="dist")
  block$dist<-as.numeric(as.character(block$dist))
  block$loc <- sub("[A-Z]",".5",as.character(block$loc),ignore.case=TRUE)
  block$loc<-as.numeric(block$loc)
  block$uniloc<-block$dist*1000+block$loc
  block<-block[,c(1,5,3,4)]
  return(block)
}

#Function of summary of data by Manzana
CleanAggDataLoc <- function(dataset) {
  block<-data.frame(table(dataset[,c("dist","loc","block")]))
  colnames(block)=c("dist","loc","block","blocksize")
  locality<-data.frame(table(dataset[c("dist","loc")]))
  colnames(locality)=c("dist","loc","locsize")
  block<-merge(block,locality,by=c("dist","loc"))
  block$dist<-as.numeric(as.character(block$dist))
  block$block<-as.numeric(as.character(block$block))
  block$loc <- sub("[A-Z]",".5",as.character(block$loc),ignore.case=TRUE)
  block$loc<-as.numeric(block$loc)
  block$uniblock<-as.double(block$dist*1000+10*block$loc+block$block/10000)
  #block<-block[,c(2,6,3,4)]
  return(block)  
}

ClusterData <- function(dataset,divide.var,size.var,cutpoints) {
  
  clusterLabel <- 0
  blockAssignments <- c()
  
  #FIND WHATEVER YOU ARE DIVIDING (E.G. DISTRICTS, LOCALITIES)
  divide<-dataset[,divide.var]
  
  #calculate clusters
  dataset$Nclusters <- (findInterval(dataset[,size.var], cutpoints))


  
  loc.comp<-  as.list(unique(divide))
  times<-length(loc.comp)

  #BEGIN LOOP: THIS WILL RUN FOR EACH THING YOU ARE DIVIDING
  for (i in 1:times) {
    currentLoc <- as.numeric(loc.comp[i])
    print(currentLoc)
    
    pull<-which(dataset[,divide.var]==currentLoc)
  
    tmpTable <- dataset[pull,]
   
    numClusts <- tmpTable[1,5]
    
 
    if(numClusts == 1){
      outRows <- as.numeric(as.character(tmpTable[,2])) 
      outRows <- cbind(outRows,rep(clusterLabel + 1,length(tmpTable[,1])))
      clusterLabel <- clusterLabel + 1
    } else {
      numBlocks <- length(tmpTable[,1])
      inData<-cbind(as.numeric(as.character(tmpTable[,2])),as.numeric(tmpTable[,3]))
      #sort to make sure it is in the right order 
      inData<-inData[order(inData[,1]),]
    
      as.numeric(tmpTable[,2])
      # Create a vector of counts for each block of the block ID
      # equal to the amount of houses on the block
      # WILL NEED TO BE CHANGED FOR DATA STRUCTURE
      myV <- c()
      for (j in 1:length(inData[,2])){
        for (k in 1:inData[j,2]){
          myV <- append(myV,inData[j,1])
        }
      }
      # Creates Intervals based on quantile interval selection
      # Should return an almost optimal selection of clusters with
      # as close to an equal amount of houses per cluster without
      # splitting the block
      c1 <- classIntervals(myV,numClusts,style = "quantile")
      t1 <- c1$brks
        
      inData<-cbind(inData,rep(0,numBlocks))
      #t1
     
        
        ### OPTIONAL
        # Based on data structure, add another column with cluster
        for(j in 1:length(inData[,1])){
            for(k in 1:numClusts){
                if (inData[j,1] <= t1[k+1]){
                    inData[j,3] <- clusterLabel + k
                    break
                }
            }
        }
        clusterLabel <- clusterLabel + numClusts
        outRows <- cbind(inData[,1],inData[,3])
    }
    blockAssignments <- rbind(blockAssignments,outRows)
  }
  colnames(blockAssignments) <-c("uniID", "clusterID")
  
  return(data.frame(blockAssignments))
}



#SET CUTPOINTS
cutpoints <- c(1000,2000,3000,4000,5000,6000,7000,8000,9000,
               10000,11000,12000,13000,14000,15000,16000,17000,
               18000,19000,20000,21000,22000,23000,24000,25000,
               26000,27000,28000,29000,30000,31000,32000,33000,
               34000,35000,36000,37000,38000,39000,40000,41000,
               42000,43000,44000)

#Run for Tiabaya, based on locality 
tiabaya <- read.csv("output/Tiabaya_uniblock.csv")
tiabaya.loc <- CleanAggDataDist(tiabaya) #organize data at district/locality level
#Tiabaya has a problem where locality 14 is actually close to 1. So we artificially make it less than 1
change.loc <-which(tiabaya.loc$uniloc==24014)
tiabaya.loc$uniloc[change.loc]=24000.75
#Same thing with locality 12
change.loc <-which(tiabaya.loc$uniloc==24012)
tiabaya.loc$uniloc[change.loc]=24000.8

#Additionally localities 3,3.5,4 are close to 10-13 whereas 6,7,8,9 are not close to 10
change.loc <-which(tiabaya.loc$uniloc==24003)
tiabaya.loc$uniloc[change.loc]=24013.8
change.loc <-which(tiabaya.loc$uniloc==24003.5)
tiabaya.loc$uniloc[change.loc]=24013.9
change.loc <-which(tiabaya.loc$uniloc==24004)
tiabaya.loc$uniloc[change.loc]=24014.1

tiabaya.clusters <- ClusterData(tiabaya.loc,"dist","distsize",cutpoints=cutpoints)  #cluster by locality
# View(tiabaya.clusters)
#remake district/locality 
tiabaya.clusters$dist <-floor(tiabaya.clusters$uniID/1000)
#change localities  back
change.loc <-which(tiabaya.clusters$uniID==24000.75)
tiabaya.clusters$uniID[change.loc]=24014
change.loc <-which(tiabaya.clusters$uniID==24000.8)
tiabaya.clusters$uniID[change.loc]=24012
change.loc <-which(tiabaya.clusters$uniID==24013.8)
tiabaya.clusters$uniID[change.loc]=24003
change.loc <-which(tiabaya.clusters$uniID==24013.9)
tiabaya.clusters$uniID[change.loc]=24003.5
change.loc <-which(tiabaya.clusters$uniID==24014.1)
tiabaya.clusters$uniID[change.loc]=24004

tiabaya.clusters$loc <-as.character(tiabaya.clusters$uniID-(tiabaya.clusters$dist*1000))
tiabaya.clusters$loc <- sub(".5","A",as.character(tiabaya.clusters$loc),ignore.case=TRUE)
tiabaya.clusters$uniID <-NULL

# Locality 9 is still grouped with 10+. This is bad because it is far away. Switch them to cluster 2. 
change.sectors <-which(tiabaya.clusters$loc==9)
tiabaya.clusters$clusterID[change.sectors]=2

#Check data

#merge cluster info
tiabaya<-merge(tiabaya,tiabaya.clusters,by=c("dist","loc"))

#cluster sizes 
table(tiabaya[c("clusterID")])
# 1    2    3    4 
# 1160 1206  476  663 

ZPLOT<-function(dataset) {
  max_lon=max(dataset$lon)
  min_lon=min(dataset$lon)
  max_lat=max(dataset$lat)
  min_lat=min(dataset$lat)
  plot(dataset$lon,dataset$lat,ylim=c(min_lat,max_lat),xlim=c(min_lon,max_lon),col=dataset$clusterID,
       xlab="Longitude",ylab="Latitude")
  #identify(dataset$lon,dataset$lat, labels=dataset$unicode)
}
ZPLOT(tiabaya)

#create sector names
tiabaya.clusters$cluster.name<-paste("tiabaya_",tiabaya.clusters$clusterID,sep="")

#export sector list
#change working directory to bandit
setwd(output)
write.csv(tiabaya.clusters,"run_bandit/data/sectores_tiabaya.csv",row.names=FALSE)







# 
# ### To view the block assignments, look at the blockAssignmetns matrix
# 
# blockAssignments1 <- as.data.frame(blockAssignments)
# names(blockAssignments1) <- c("uni_block", "clusterID")
# 
# Clusters <-merge(block, blockAssignments1, by= "uni_block")
# Clusters_info <- merge(Clusters, block_percent, by.x = "uni_block", by.y = "uniMZ", all.x = T)
# 
# write.csv(Cluster_data, file = "Cluster_data.csv")
# 
# ##
# clusterSizes <- aggregate(blocksize ~ clusterID, data=Clusters, FUN=sum)
# hist(clusterSizes[,2])
# sortedSizes <- clusterSizes[order(clusterSizes$blocksize),]
# ##
# 
# Clusters_collapsed <- aggregate(cbind(loc, locsize, Nclusters, Cinfested, Csprayed, Cclosed,
#                                       Creluctant, Cuninhabited) ~ clusterID, data = Cluster_data, FUN = mean)
# write.csv(Clusters_collapsed, file = "Clusters_collapsed.csv")
