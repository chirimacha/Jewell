#set seed
set.seed(1234)

#load libraries
library("lubridate")
library("PBSmapping")
library("plyr")
library("inline")
library("Rcpp")

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

###########################
###read data##############
#############################

tic()

#set working drive
setwd("/Users/EMWB/Jewell/Data")
#setwd("~/Desktop/Levy Lab")
#setwd("~/Users/e/Jewell/Data")
#setwd("/Users/mzlevy/Jewell/Data")

#set up tiabaya_gps dataset
tiabaya.gps = read.csv("Tiabaya_GPS.csv")

getUTM<-function(id, x,y){
  coord <- data.frame(x,y)
  colnames(coord) <- c("X","Y")
  attributes(coord)$projection <- "LL"
  utm.coord <- convUL(coord)
  utm.coord$X<-utm.coord$X*1000
  utm.coord$Y<-utm.coord$Y*1000
  return(utm.coord)
}
tiabaya.test <- getUTM(x=tiabaya.gps$LONGITUDE,y=tiabaya.gps$LATITUDE)

#link unicodes with coordinates
tiabaya.gps <- cbind(tiabaya.gps$UNICODE,tiabaya.test)
tiabaya.gps <- rename(tiabaya.gps,c("tiabaya.gps$UNICODE" = "UNICODE"))


#read in data
inspecciones = read.csv("inspecciones.csv")
vig = read.csv("byHouse_fullEID.csv")
priors = read.csv("Corentins_Predictions_Jun-24-2015_07-13-06.csv")


#merge data
i.v <- merge(inspecciones,vig, by="UNICODE",all=TRUE)
i.v.gps <- merge(i.v,tiabaya.gps,by="UNICODE")
i.v.gps <- merge(i.v.gps, priors, by="UNICODE",all.x=TRUE)

#subset data to only include localities 4,5,6
i.v.gps.456 <- i.v.gps[which(i.v.gps$L.y==4 | i.v.gps$L.x==4 | i.v.gps$L.y==5 | i.v.gps$L.x==5 |i.v.gps$L.y==6 | i.v.gps$L.x==6),]


#Replace NA values with zeros for sum
i.v.gps.456$PD_TCAP_TOT <- ifelse(is.na(i.v.gps.456$PD_TCAP_TOT),0,i.v.gps.456$PD_TCAP_TOT)
i.v.gps.456$IN_TCAP_TOT <- ifelse(is.na(i.v.gps.456$IN_TCAP_TOT),0,i.v.gps.456$IN_TCAP_TOT)

#Replace NA values for prior probability with median value
#find median value of those that are not NA
median.pred.prob <- median(i.v.gps.456$predicteddensity[which(!is.na(i.v.gps.456$predicteddensity))])

#replace NAs with this value
predprobs <- ifelse(is.na(i.v.gps.456$predicteddensity), median.pred.prob, i.v.gps.456$predicteddensity)

#sum inspecciones in districts 4,5,6
sum.insp <- i.v.gps.456$PD_TCAP_TOT + i.v.gps.456$IN_TCAP_TOT

#get unicodes as strings
unicode<-as.character(i.v.gps.456$UNICODE)

#plot new dataset
plot(i.v.gps.456$X.y,i.v.gps.456$Y.y,col="slategray3",xaxt="n",yaxt="n",pch=18,cex=.4)
points(i.v.gps.456$X.y[which(i.v.gps.456$L.y==4)],i.v.gps.456$Y.y[which(i.v.gps.456$L.y==4)],pch=19,cex=.3,col="green")
points(i.v.gps.456$X.y[which(i.v.gps.456$L.y==5)],i.v.gps.456$Y.y[which(i.v.gps.456$L.y==5)],pch=19,cex=.3,col="purple")
points(i.v.gps.456$X.y[which(i.v.gps.456$L.y==6)],i.v.gps.456$Y.y[which(i.v.gps.456$L.y==6)],pch=19,cex=.3,col="orange")
for (i in 1:length(sum.insp)) points(i.v.gps.456$X.y[i][which(sum.insp[i]>0)],i.v.gps.456$Y.y[i][which(sum.insp[i]>0)],pch=18,cex=sum.insp[i]*.1+1,col="red")

##################################################
#######define parameters##############################
#########################################

########set up times#######

#Import dates as separate columns for month, day and year in that respective order
date <- function(m,d,y){
  #Convert it into one string
  right.date <- paste(m,d,y,sep = "/", collapse = NULL)
  #Read it as a date in the right format for R
  
  new.dates <- as.Date(right.date, "%m/%d/%Y")
  
  return(new.dates)
}
#outputs dates in the correct format that R uses
test <- date(i.v.gps.456$MES,i.v.gps.456$DIA,i.v.gps.456$ANIO)

############################################################
#############FEEL FREE TO CHANGE STARTING HERE##############
############################################################
#if more than one observation for a house, pick most recent

#identify unique unicodes
unique.unicodes <- unique(unicode)

#make dataset with just unicodes and dates
unicodes.dates <- cbind(unicode, test)

#for houses with more than on observation, delete all observations except latest inspection
#this part doesn't work yet...
for (i in 1:length(unique.unicodes)){
  house <- unicode[which(unicode==unique.unicodes[i])] 
  numberobs <- length(house)
  if(numberobs>1){
    houseobs <- unicodes.dates[order(unicode,test),] 
    tokeep <- by(unicodes.dates,unicode,tail,n=1)
    tokeepd <- do.call("rbind", as.list(tokeep))
    i.v.gps.456[which(i.v.gps.456$UNICODE==house),]
  }
}
