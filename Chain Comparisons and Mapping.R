'
################################
This code:
  -Pulls in results from a set of Jewell MCMC Chains
  -Creates tables showing relative rankings based on reference chain 
  -Maps the top N houses for each chain

CREATOR: SN 
######################################
'

# Set working directory to location of Jewell Results
setwd("/Users/snutman/Downloads/resultsfromamazoncluster/")


TimeNow <- function() {
  x <- format(Sys.time(), "%b-%d-%Y_%H-%M-%S")
  return(x)
}


files <- list.files(pattern="*Results*")

names=vector()
#read data
for (i in 1:length(files)) {
  name=paste("chain",i,sep="_")
  #note this data has to be in standard CSV form
  data <-data.frame(read.csv(files[i]))
  data$ranking <- rank(-data$occult.prob,ties.method = "random")
  data <-data[order(data$occult.prob,decreasing=TRUE),]
  data$unicode <-as.character(data$unicode)
  #chain assignments based on order in the "files" list
  data$chain <-i
  assign(name,data,envir = .GlobalEnv)
  names[i]<-name
}

#pull tiabaya 1 to get lon/lat
tiabaya <-read.csv("/Users/snutman/Documents/CHAGAS_DATA/bandit/run_bandit/data/tiabaya1_waddlhouses.csv")
tiabaya <-data.frame(unicode=tiabaya$unicode,lon=tiabaya$lon,lat=tiabaya$lat)

SubsetData <- function (names,houses) {
  
  data.subset <-get(names[1],envir=.GlobalEnv)
  data.subset<-data.subset[1:houses,]
  table<-data.frame(unicode=data.subset$unicode,Ranking1=data.subset$ranking)
  
  for (i in 2:length(names)) {
    temp <- get(names[i],envir=.GlobalEnv)
    sel <- which(temp$unicode %in% data.subset$unicode)
    temp <- temp[sel,]
    temp.table <- data.frame(unicode=temp$unicode,ranking=temp$ranking)
    table <- merge(table,temp.table,by="unicode",all=TRUE)
  }
  colname <-unlist(list("unicode",names))
  colnames(table) <-colname
  return(table)
}
   
SubsetMapping <- function (names,houses) {
  mapping.data <-data.frame()
  for (i in 1:length(names)) {
    data <- get(names[i],envir=.GlobalEnv)
    data <- data[1:houses,]
    data$chain=i
    mapping.data <- rbind(mapping.data,data)
  }
  mapping.data$chain <-as.factor(mapping.data$chain)
  return(mapping.data)
}

MapPoints <- function(dataset,group) {
  max_lon=max(dataset$lon)
  min_lon=min(dataset$lon)
  max_lat=max(dataset$lat)
  min_lat=min(dataset$lat)
  plot(dataset$lon,dataset$lat,ylim=c(min_lat,max_lat),xlim=c(min_lon,max_lon),col=group,pch=as.integer(group), xlab="Longitude",ylab="Latitude")
}

MapPoints.alltia <- function(dataset,group) {
  max_lon=max(tiabaya$lon)
  min_lon=min(tiabaya$lon)
  max_lat=max(tiabaya$lat)
  min_lat=min(tiabaya$lat)
  plot(dataset$lon,dataset$lat,ylim=c(min_lat,max_lat),xlim=c(min_lon,max_lon),col=group,pch=as.integer(group), xlab="Longitude",ylab="Latitude")
}

table <- SubsetData(names,100)
map.data <- SubsetMapping(names,15)
map.data <-merge(map.data,tiabaya,by="unicode")

MapPoints(map.data,map.data$chain)

MapPoints.alltia(map.data,map.data$chain)


