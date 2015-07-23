'
############################################
BORDER ZONE BLOCKING - FOR TIABAYA

This code:
  -Pulls in generalized border zone blocking code 
  -Conducts border zone blocking for Tiabaya clusters

CREATOR: SN
#################################

'
#set path to spatial_surveillence 
if(Sys.getenv("SPATIAL_DATA") != "") {
  setwd(paste(Sys.getenv("SPATIAL_DATA"), sep=""))  
} else if(Sys.getenv("USER") == "snutman") {
  setwd("/Users/snutman/Documents/CHAGAS_DATA/spatial_surveillence/Spatial_data")
} else if(Sys.getenv("USER") == "USER") {
  setwd("PATH_TO_SPATIAL_DATA")
}

#set bandit path

if(Sys.getenv("SPATIAL_UNCERTAINTY") != "") {
  bandit.path <-Sys.getenv("SPATIAL_UNCERTAINTY")
} else if(Sys.getenv("USER") == "sgutfraind") {
  bandit.path <- "/Users/sgutfraind/academic/chagas/bandit"
} else if(Sys.getenv("USER") == "mlevy") {
  #TODO: check your USER string (your username on your computer) and add appropriate directory here with 
  bandit.path<-"PATH_TO_BANDIT"
} 

#timing
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


#pull in border code 
source("code/border_zone_additions.R")

open.file<-function(path) {
  dataset<-read.csv(path)
  dataset$X<-NULL
  dataset$Y<-NULL
  UTM <- getUTM(x=dataset$lon,y=dataset$lat)
  dataset <- cbind(dataset,UTM)
  return(dataset)
}

#pull tiabaya  
tiabaya <- open.file("output/Tiabaya_uniblock.csv")

#pull surroudning districts: sachaca, uchumayo, socabaya, hunter
sachaca<- open.file("output/Sachaca_uniblock.csv")
uchumayo<- open.file("output/Uchumayo_uniblock.csv")
hunter<- open.file("output/Hunter_uniblock.csv")

#pull bandit data 

sectores<-read.csv(paste(bandit.path,"/run_bandit/data/sectores_tiabaya.csv",sep=""))
sectores$cluster.name<-NULL
#merge cluster info
tiabaya<-merge(tiabaya,sectores,by=c("dist","loc"))


AddLocality <- function (target.district,other.districts,sector) {
  target.data<-target.district[target.district$clusterID==sector,]
  target.add <-target.district[target.district$clusterID!=sector,]
  target.add$clusterID<-NULL
  total.add <-rbind(target.add,other.districts)
  add.points <-find.additional.points(target=target.data,additional.points = total.add)
  return(add.points)
}

#make datasets for additional points and run

#sector1 is in the south. It can't touch sachaca
  other.districts<-rbind(hunter,uchumayo)
  tic()
  tiabaya1 <- AddLocality(target.district=tiabaya,other.districts=other.districts,sector=1)
  toc()
#sector 2 is in the east. it can't touch uchumayo 
  tic()
  other.districts<-rbind(hunter,sachaca)
  tiabaya2 <- AddLocality(target.district=tiabaya,other.districts=other.districts,sector=2)
  toc()
#sector 3 can touch all three surrounding districts 
  tic()
  other.districts<-rbind(hunter,sachaca,uchumayo)
  tiabaya3 <- AddLocality(target.district=tiabaya,other.districts=other.districts,sector=3)
  toc()
#function to add houses
AddExtraHouses <- function(target.data,add.houses,sector) {
  target.data<-target.data[target.data$clusterID==sector,]
  target.data$mark<-0
  add.houses$clusterID <-NA
  target.data<-rbind(target.data,add.houses)
  names(target.data)[names(target.data) == 'mark'] <- 'add.house'
  return(target.data)
}

#now add the districts and make the final dataset

tiabaya1.fin <-AddExtraHouses(tiabaya,tiabaya1,1)
tiabaya2.fin <-AddExtraHouses(tiabaya,tiabaya2,2)
tiabaya3.fin <-AddExtraHouses(tiabaya,tiabaya3,3)

write.csv(tiabaya1.fin,paste(bandit.path,"/run_bandit/data/tiabaya1_waddlhouses",".csv",sep=""),row.names = FALSE)
write.csv(tiabaya2.fin,paste(bandit.path,"/run_bandit/data/tiabaya2_waddlhouses",".csv",sep=""),row.names = FALSE)
write.csv(tiabaya3.fin,paste(bandit.path,"/run_bandit/data/tiabaya3_waddlhouses",".csv",sep=""),row.names = FALSE)



  