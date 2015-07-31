## Set working directory to location of output csv
setwd("~/Desktop/Levy_Lab/Jewell/SensAnalysis/K.Sens.Results")

## In v1, v2, and v3 insert the name of the csv and run all
v1 <- data.frame(read.table("ResultsKv1"))
v2 <- data.frame(read.table("ResultsKv2"))
v3 <- data.frame(read.table("ResultsKv3"))
v1 <- v1[order(v1$unicode),]
v2 <- v2[order(v2$unicode),]
v3 <- v3[order(v3$unicode),]
probmean <- ((v1$occult.prob + v2$occult.prob + v3$occult.prob)/3)
v1$probmean <- probmean
v <- v1[order(v1$probmean, decreasing = TRUE),]
Ranking <- c(1:dim(v)[1])
averageRanking <- cbind(v, Ranking)
averageRanking$occult.prob <- NULL 

## New single csv should be in same folder as working directory 
write.csv(averageRanking, file="AverageResults",row.names=FALSE)
