setwd("~/Desktop/Levy_Lab/Jewell/SensAnalysis/K.Sens.Results")
v1 <- data.frame(read.table("100.1000000ResultsKv1"))
v2 <- data.frame(read.table("100.1000000ResultsKv2"))
v3 <- data.frame(read.table("100.1000000ResultsKv3"))
v1 <- v1[order(v1$unicode),]
v2 <- v2[order(v2$unicode),]
v3 <- v3[order(v3$unicode),]
probmean <- ((v1$occult.prob + v2$occult.prob + v3$occult.prob)/3)
v1$probmean <- probmean
v <- v1[order(v1$probmean, decreasing = TRUE),]
Ranking <- c(1:dim(v)[1])
averageRanking <- cbind(v, Ranking)
averageRanking$occult.prob <- NULL 
write.csv(averageRanking, file="AverageResults",row.names=FALSE)
