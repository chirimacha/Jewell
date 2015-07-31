library(ggplot2)
library(grid)
source("http://peterhaschke.com/Code/multiplot.R")

setwd("~/Desktop/Levy_Lab/Jewell/SensAnalysis/K.Sens.Results")
v1.100 <- data.frame(read.table("100.1000000ResultsKv1"))
v2.100 <- data.frame(read.table("100.1000000ResultsKv2"))
v3.100 <- data.frame(read.table("100.1000000ResultsKv3"))
v1.500 <- data.frame(read.table("500.1000000ResultsKv1"))
v2.500 <- data.frame(read.table("500.1000000ResultsKv2"))
v3.500 <- data.frame(read.table("500.1000000ResultsKv3"))
v1.1000 <- data.frame(read.table("1000.1000000Resultsv1"))
v2.1000 <- data.frame(read.table("1000.1000000Resultsv2"))
v3.1000 <- data.frame(read.table("1000.1000000Resultsv3"))

v1.100 <- v1.100[order(v1.100$unicode),]
v2.100 <- v2.100[order(v2.100$unicode),]
v3.100 <- v3.100[order(v3.100$unicode),]

probmean <- ((v1.100$occult.prob + v2.100$occult.prob + v3.100$occult.prob)/3)

v1.100$probmean <- probmean
v.100 <- v1.100[order(v1.100$probmean, decreasing = TRUE),]

Ranking <- c(1:dim(v.100)[1])
v.100 <- cbind(v.100, Ranking)

v1.500 <- v1.500[order(v1.500$unicode),]
v2.500 <- v2.500[order(v2.500$unicode),]
v3.500 <- v3.500[order(v3.500$unicode),]

probmean2 <- ((v1.500$occult.prob + v2.500$occult.prob + v3.500$occult.prob)/3)
v1.500$probmean2 <- probmean2

v.500 <- v1.500[order(v1.500$probmean2, decreasing = TRUE),]
v.500 <- cbind(v.500, Ranking)

vdata.1.5 <- merge(v.100, v.500, by = "unicode")
vdata.1.5 <- vdata.1.5[order(vdata.1.5$Ranking.x),]

#plot(vdata.1.5$Ranking.x, vdata.1.5$Ranking.y) 

v1.1000 <- v1.1000[order(v1.1000$unicode),]
v2.1000 <- v2.1000[order(v2.1000$unicode),]
v3.1000 <- v3.1000[order(v3.1000$unicode),]

probmean3 <- ((v1.1000$occult.prob + v2.1000$occult.prob + v3.1000$occult.prob)/3)
v1.1000$probmean3 <- probmean3

v.1000 <- v1.1000[order(v1.1000$probmean3, decreasing = TRUE),]
v.1000 <- cbind(v.1000, Ranking)

vdata.1.1 <- merge(v.1000, v.100, by = "unicode")
vdata.1.1 <- vdata.1.1[order(vdata.1.1$Ranking.x),]

#plot.1000.100 <- plot(vdata.1.1$Ranking.x, vdata.1.1$Ranking.y)

vdata.10.5 <- merge(v.1000, v.500, by = "unicode")
vdata.10.5 <- vdata.10.5[order(vdata.10.5$Ranking.x),]

#plot(vdata.10.5$Ranking.x, vdata.10.5$Ranking.y)

P1 <- ggplot(vdata.1.5, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) + 
  xlab("K100") + ylab("K500") + ggtitle("Sensitivity Analysis of Carrying Capacity at 100 and 500")
print(P1)
P2 <- ggplot(vdata.1.1, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("K1000") + ylab("K100") + ggtitle("Sensitivity Analysis of Carrying Capacity at 100 and 1000")
print(P2)
P3 <- ggplot(vdata.10.5, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("K1000") + ylab("K500")+ ggtitle("Sensitivity Analysis of Carrying Capacity at 500 and 1000")
print(P3)
#multiplot(P1, P2, P3, cols = 3) 


