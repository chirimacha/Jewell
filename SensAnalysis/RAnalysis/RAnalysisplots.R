library(ggplot2)
library(grid)
source("http://peterhaschke.com/Code/multiplot.R")

## Read the results of the different chains
setwd("~/Desktop/Levy_Lab/Jewell_output/RAnalysis/output")
v1.rb1 <- data.frame(read.table("1000000ResultsAug10v1Rb1"))
v1.rb11 <- data.frame(read.table("1000000ResultsAug10v1Rb11"))
v1.rb12 <- data.frame(read.table("1000000ResultsAug10v1Rb12"))
v2.rb1 <- data.frame(read.table("1000000ResultsAug10v2Rb1"))
v2.rb11 <- data.frame(read.table("1000000ResultsAug10v2Rb11"))
v2.rb12 <- data.frame(read.table("1000000ResultsAug10v2Rb12"))
v3.rb1 <- data.frame(read.table("1000000ResultsAug10v3Rb1"))
v3.rb11 <- data.frame(read.table("1000000ResultsAug10v3Rb11"))
v3.rb12 <- data.frame(read.table("1000000ResultsAug10v3Rb12"))

## Since the probabilities are ordered by decreasing already, bind the rankings
 # to the data frames
Ranking <- c(1:dim(v1.rb1)[1])
v1.rb1 <- cbind(v1.rb1, Ranking)
v1.rb11 <- cbind(v1.rb11, Ranking)
v1.rb12 <- cbind(v1.rb12, Ranking)
v2.rb1 <- cbind(v2.rb1, Ranking)
v2.rb11 <- cbind(v2.rb11, Ranking)
v2.rb12 <- cbind(v2.rb12, Ranking)
v3.rb1 <- cbind(v3.rb1, Ranking)
v3.rb11 <- cbind(v3.rb11, Ranking)
v3.rb12 <- cbind(v3.rb12, Ranking)

## merge each of the three chains of dataset by unicode into one
rb1 <- merge(v1.rb1, v2.rb1, by = "unicode")
rb1 <- merge(rb1, v3.rb1, by = "unicode")
rb1 <- rb1[order(rb1$Ranking.x),]
rb11 <- merge(v1.rb11, v2.rb11, by = "unicode")
rb11 <- merge (rb11, v3.rb11, by = "unicode")
rb11 <- rb11[order(rb11$Ranking.x),]
rb12 <- merge(v1.rb12, v2.rb12, by = "unicode")
rb12 <- merge (rb12, v3.rb12, by = "unicode")
rb12 <- rb12[order(rb12$Ranking.x),]

## Graph the rankings between the chains
P1 <- ggplot(rb1, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("Chain 1 at r = 1.00005") + ylab("Chain 2 at r = 1.00005") + ggtitle("Rb Analysis")
print(P1)

P2 <- ggplot(rb1, aes(Ranking.x, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 1 at r = 1.00005") + ylab("Chain 3 at r = 1.00005") + ggtitle("Rb Analysis")
print(P2)

P3 <- ggplot(rb1, aes(Ranking.y, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 2 at r = 1.00005") + ylab("Chain 3 at r = 1.00005") + ggtitle("Rb Analysis")
print(P3)

P4 <- ggplot(rb11, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("Chain 1 at r = 1.1") + ylab("Chain 2 at r = 1.1") + ggtitle("Rb Analysis")
print(P4)

P5 <- ggplot(rb11, aes(Ranking.x, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 1 at r = 1.1") + ylab("Chain 3 at r = 1.1") + ggtitle("Rb Analysis")
print(P5)

P6 <- ggplot(rb11, aes(Ranking.y, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 2 at r = 1.1") + ylab("Chain 3 at r = 1.1") + ggtitle("Rb Analysis")
print(P6)

P7 <- ggplot(rb12, aes(Ranking.y, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 2 at r = 1.2") + ylab("Chain 3 at r = 1.2") + ggtitle("Rb Analysis")
print(P7)

P8 <- ggplot(rb12, aes(Ranking.x, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 1 at r = 1.2") + ylab("Chain 3 at r = 1.2") + ggtitle("Rb Analysis")
print(P8)

P9 <- ggplot(rb12, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("Chain 1 at r = 1.2") + ylab("Chain 2 at r = 1.2") + ggtitle("Rb Analysis")
print(P9)

## merge one chain of each r value by unicode into one 
rb.v1 <- merge(v1.rb1, v1.rb11, by = "unicode")
rb.v1 <- merge(rb.v1, v1.rb12, by = "unicode")
rb.v2 <- merge(v2.rb1, v2.rb11, by = "unicode")
rb.v2 <- merge(rb.v2, v2.rb12, by = "unicode")
rb.v3 <- merge(v3.rb1, v3.rb11, by = "unicode")
rb.v3 <- merge(rb.v3, v3.rb12, by = "unicode")
rb.v1 <- rb.v1[order(rb.v1$Ranking.x),]
rb.v2 <- rb.v2[order(rb.v2$Ranking.x),]
rb.v3 <- rb.v3[order(rb.v3$Ranking.x),]

##Plots between different chains 
P10 <- ggplot(rb.v1, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("Chain 1 at r = 1.00005") + ylab("Chain 1 at r = 1.11") + ggtitle("Rb Analysis")
print(P10)

P11 <- ggplot(rb.v1, aes(Ranking.x, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 1 at r = 1.00005") + ylab("Chain 1 at r = 1.12") + ggtitle("Rb Analysis")
print(P11)

P12 <- ggplot(rb.v1, aes(Ranking.y, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 1 at r = 1.11") + ylab("Chain 1 at r = 1.12") + ggtitle("Rb Analysis")
print(P12)

P13 <- ggplot(rb.v2, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("Chain 2 at r = 1.00005") + ylab("Chain 2 at r = 1.11") + ggtitle("Rb Analysis")
print(P13)

P14 <- ggplot(rb.v2, aes(Ranking.x, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 2 at r = 1.00005") + ylab("Chain 2 at r = 1.12") + ggtitle("Rb Analysis")
print(P14)

P15 <- ggplot(rb.v1, aes(Ranking.y, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 2 at r = 1.11") + ylab("Chain 2 at r = 1.12") + ggtitle("Rb Analysis")
print(P15)

P16 <- ggplot(rb.v3, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("Chain 3 at r = 1.00005") + ylab("Chain 3 at r = 1.11") + ggtitle("Rb Analysis")
print(P16)

P17 <- ggplot(rb.v3, aes(Ranking.x, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 3 at r = 1.00005") + ylab("Chain 3 at r = 1.12") + ggtitle("Rb Analysis")
print(P17)

P18 <- ggplot(rb.v3, aes(Ranking.y, Ranking)) + geom_point(alpha = 5/10) +
  xlab("Chain 3 at r = 1.11") + ylab("Chain 3 at r = 1.12") + ggtitle("Rb Analysis")
print(P18)

## Create a table containing rankings of the top houses
#  Function for the top rankings of r = 1.00005
toprankings.rb1<-function(top)
{
  j <- rb1[,c("id.x", "unicode", "Ranking.x", "Ranking.y", "Ranking")]
  togo <- data.frame(j[0:top,])
  names(togo) <- c("Id", "Unicode", "Chain 1", "Chain 2", "Chain 3")
  rownames(togo) <- NULL
  return(togo)
}

# Function for the top rankings of r = 1.1
toprankings.rb11<-function(top)
{
  j <- rb11[,c("id.x", "unicode", "Ranking.x", "Ranking.y", "Ranking")]
  togo <- data.frame(j[0:top,])
  names(togo) <- c("Id", "Unicode", "Chain 1", "Chain 2", "Chain 3")
  rownames(togo) <- NULL
  return(togo)
}
# Function for the top rankings of r= 1.2
toprankings.rb12<-function(top)
{
  j <- rb12[,c("id.x", "unicode", "Ranking.x", "Ranking.y", "Ranking")]
  togo <- data.frame(j[0:top,])
  names(togo) <- c("Id", "Unicode", "Chain 1", "Chain 2", "Chain 3")
  rownames(togo) <- NULL
  return(togo)
}

write.csv(toprankings.rb1(15), file="toprankings.rb1.csv",row.names=FALSE)
write.csv(toprankings.rb11(15), file="toprankings.rb11.csv",row.names=FALSE)
write.csv(toprankings.rb12(15), file="toprankings.rb12.csv",row.names=FALSE)

## Maps of the Top Houses 

#set working drive to /Jewell/Data 
setwd("/Users/patrickemedom/Desktop/Levy_Lab/Jewell/Data")

library("lubridate")
library("PBSmapping")
library("plyr")
library("inline")
library("Rcpp")

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
blockdata = read.csv("Tiabaya_uniblock.csv")
blockdata <- subset(blockdata, select = c(5,6)) #keep only unicodes and blocks
blockdata<- rename(blockdata, c("unicode" = "UNICODE"))


#merge data
i.v <- merge(inspecciones,vig, by="UNICODE",all=TRUE)
i.v.gps <- merge(i.v,tiabaya.gps,by="UNICODE")
i.v.gps <- merge(i.v.gps, priors, by="UNICODE",all.x=TRUE)
i.v.gps <- merge(i.v.gps, blockdata, by = "UNICODE", all=TRUE)

#subset data to only include localities 4,5,6
i.v.gps.456 <- i.v.gps[which(i.v.gps$L.y==4 | i.v.gps$L.x==4 | i.v.gps$L.y==5 | i.v.gps$L.x==5 |i.v.gps$L.y==6 | i.v.gps$L.x==6),]


#Replace NA values with zeros for sum
i.v.gps.456$PD_TCAP_TOT <- ifelse(is.na(i.v.gps.456$PD_TCAP_TOT),0,i.v.gps.456$PD_TCAP_TOT)
i.v.gps.456$IN_TCAP_TOT <- ifelse(is.na(i.v.gps.456$IN_TCAP_TOT),0,i.v.gps.456$IN_TCAP_TOT)

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
i.v.gps.456$date <- date(i.v.gps.456$MES,i.v.gps.456$DIA,i.v.gps.456$ANIO)


############################################################
############# UPDATED CHANGES AND NEW FOR LOOP #############
############################################################
#if more than one observation for a house, pick most recent

#identify unique unicodes
unicode<-as.character(i.v.gps.456$UNICODE)
unique.unicodes <- unique(unicode)
dates <- i.v.gps.456$date

#find repeated unicodes
repeated.unicodes <- unicode[which(duplicated(unicode) == TRUE)]

#setting an empty vector for the for loop
unique.dates <- c(1:length(unique.unicodes)*NA)

#for houses with more than on observation, only take into account the one with latest date

for (i in 1:length(unique.unicodes)){
  
  u <- which(unicode==unique.unicodes[i])
  fecha <- i.v.gps.456$date[u]
  
  if(is.na(fecha[1]) ==FALSE) {
    maxf <- max(fecha, na.rm = TRUE)
    v <- which(i.v.gps.456$date == maxf)
    unique.dates[i] <- intersect(u,v)
  }
  else {
    unique.dates[i] <- max(u)
  }
}

#new dataframe with single unicodes 
unique.data <-i.v.gps.456[unique.dates,]


#if more than one observation for a house, pick most recent
i.v.gps.456$date <- date(i.v.gps.456$MES,i.v.gps.456$DIA,i.v.gps.456$ANIO)

#new dataframe with single unicodes 
unique.data <-i.v.gps.456[unique.dates,]

#rename dataset
i.v.gps.456 <- unique.data

earliest <- sort(i.v.gps.456$date)[1]
latest <- sort(i.v.gps.456$date)[length(i.v.gps.456$date[which(!is.na(i.v.gps.456$date))])]
timetest <- (i.v.gps.456$date - earliest)/90
initialtime <- date(12, 31, 2004)
today <- date(7, 22, 2015)
timefrombeginning <- round((earliest - initialtime)/90)
tobs <- ceiling(timetest) + timefrombeginning
maxt <- round((today - latest)/90)+max(tobs[which(!is.na(tobs))])
maxt <- as.numeric(maxt)

#identify which houses weren't inspected
inspected <- ifelse(is.na(tobs),0,1)

#replace NAs with max time
tobs = ifelse(is.na(tobs), maxt,tobs)

#sum inspecciones in districts 4,5,6
sum.insp <- i.v.gps.456$PD_TCAP_TOT + i.v.gps.456$IN_TCAP_TOT

#base r = 1.00005 with top houses of varied chains
par(mar=c(1,1,1,1), xpd=TRUE)
colfunc = gray.colors(length(unique(v1.rb1[,2])),start=1,end=0)[as.factor(v1.rb1[,2])]
plot((v1.rb1$X), (v1.rb1$X.1), col = colfunc, pch=16, axes=FALSE, xlab="",ylab="")
points((v1.rb1$X[1:10]), (v1.rb1$X.1[1:10]), col = "green",pch=20, cex = 1.2)
points((v2.rb1$X[1:10]), (v2.rb1$X.1[1:10]), col = "blue",pch=14)
points((v3.rb1$X[1:10]), (v3.rb1$X.1[1:10]), col = "gold",pch=7, cex=1.5)
for (i in 1:867) if(sum.insp[i]>0) points(i.v.gps.456$X.y[i],i.v.gps.456$Y.y[i],pch=18,col="firebrick3",cex=1.5)
legend("topleft",c("Known Infested House", "Top Probability of Infestation at Chain 1", "Top Probability of Infestation at Chain 2",
                   "Top Probability of Infestation at Chain 3"), cex=0.65 ,pch=c(18,20,14,7),col=c("firebrick3","green", "blue", "gold"),bty="n")

#base r = 1.11 with top houses of varied chains
par(mar=c(1,1,1,1), xpd=TRUE)
colfunc = gray.colors(length(unique(v1.rb11[,2])),start=1,end=0)[as.factor(v1.rb11[,2])]
plot((v1.rb11$X), (v1.rb11$X.1), col = colfunc, pch=16, axes=FALSE, xlab="",ylab="")
points((v1.rb11$X[1:10]), (v1.rb11$X.1[1:10]), col = "green",pch=20, cex = 1.2)
points((v2.rb11$X[1:10]), (v2.rb11$X.1[1:10]), col = "blue",pch=14)
points((v3.rb11$X[1:10]), (v3.rb11$X.1[1:10]), col = "gold",pch=7, cex=1.5)
for (i in 1:867) if(sum.insp[i]>0) points(i.v.gps.456$X.y[i],i.v.gps.456$Y.y[i],pch=18,col="firebrick3",cex=1.5)
legend("topleft",c("Known Infested House", "Top Probability of Infestation at Chain 1", "Top Probability of Infestation at Chain 2",
                   "Top Probability of Infestation at Chain 3"), cex=0.65 ,pch=c(18,20,14,7),col=c("firebrick3","green", "blue", "gold"),bty="n")

#base r = 1.12 with top houses of varied chains
par(mar=c(1,1,1,1), xpd=TRUE)
colfunc = gray.colors(length(unique(v1.rb12[,2])),start=1,end=0)[as.factor(v1.rb12[,2])]
plot((v1.rb12$X), (v1.rb12$X.1), col = colfunc, pch=16, axes=FALSE, xlab="",ylab="")
points((v1.rb12$X[1:10]), (v1.rb12$X.1[1:10]), col = "green",pch=20, cex = 1.2)
points((v2.rb12$X[1:10]), (v2.rb12$X.1[1:10]), col = "blue",pch=14)
points((v3.rb12$X[1:10]), (v3.rb12$X.1[1:10]), col = "gold",pch=7, cex=1.5)
for (i in 1:867) if(sum.insp[i]>0) points(i.v.gps.456$X.y[i],i.v.gps.456$Y.y[i],pch=18,col="firebrick3",cex=1.5)
legend("topleft",c("Known Infested House", "Top Probability of Infestation at Chain 1", "Top Probability of Infestation at Chain 2",
                   "Top Probability of Infestation at Chain 3"), cex=0.65 ,pch=c(18,20,14,7),col=c("firebrick3","green", "blue", "gold"),bty="n")

toplocation <- function(d) 
  {
par(mar=c(1,1,1,1), xpd=TRUE)
colfunc = gray.colors(length(unique(d[,2])),start=1,end=0)[as.factor(d[,2])]
plot((d$X), (d$X.1), col = colfunc, pch=16, axes=FALSE, xlab="",ylab="")
points((d$X[1:10]), (d$X.1[1:10]), col = "gold",pch=18)
for (i in 1:867) if(sum.insp[i]>0) points(i.v.gps.456$X.y[i],i.v.gps.456$Y.y[i],pch=18,col="firebrick3",cex=1.5)
legend("topleft",c("Known Infested House", "Top Probability of Infestation"),pch=c(18,18),col=c("firebrick3","gold"),bty="n")
}