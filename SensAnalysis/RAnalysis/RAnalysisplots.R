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
rb1 <- merge (rb1, v3.rb1, by = "unicode")
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

## Create a table containing rankings of the top houses
#  Function for the top rankings of r = 1.00005
toprankings.rb1<-function(top)
{
  j <- rb1[,c("id.x","Ranking.x", "Ranking.y", "Ranking")]
  togo <- data.frame(j[0:top,])
  names(togo) <- c("Id", "Chain 1", "Chain 2", "Chain 3")
  rownames(togo) <- NULL
  return(togo)
}

# Function for the top rankings of r = 1.1
toprankings.rb11<-function(top)
{
  j <- rb11[,c("id.x","Ranking.x", "Ranking.y", "Ranking")]
  togo <- data.frame(j[0:top,])
  names(togo) <- c("Id", "Chain 1", "Chain 2", "Chain 3")
  rownames(togo) <- NULL
  return(togo)
}
# Function for the top rankings of r= 1.2
toprankings.rb12<-function(top)
{
  j <- rb12[,c("id.x","Ranking.x", "Ranking.y", "Ranking")]
  togo <- data.frame(j[0:top,])
  names(togo) <- c("Id", "Chain 1", "Chain 2", "Chain 3")
  rownames(togo) <- NULL
  return(togo)
}

write.csv(toprankings.rb1(15), file="toprankings.rb1.csv",row.names=FALSE)
write.csv(toprankings.rb11(15), file="toprankings.rb11.csv",row.names=FALSE)
write.csv(toprankings.rb12(15), file="toprankings.rb12.csv",row.names=FALSE)