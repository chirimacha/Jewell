library(ggplot2)

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

##### K = 100 #####
v1.100 <- v1.100[order(v1.100$occult.prob, decreasing = TRUE),]
v2.100 <- v2.100[order(v2.100$occult.prob, decreasing = TRUE),]
v3.100 <- v3.100[order(v3.100$occult.prob, decreasing = TRUE),]

Ranking <- c(1:dim(v1.100)[1])
v1.100 <- cbind(v1.100, Ranking)
v2.100 <- cbind(v2.100, Ranking)
v3.100 <- cbind(v3.100, Ranking)

vdata.1.2.100 <- merge(v1.100, v2.100, by = "unicode")
vdata.1.2.100 <- vdata.1.2.100[order(vdata.1.2.100$Ranking.x),]
vdata.1.3.100 <- merge(v1.100, v3.100, by = "unicode")
vdata.1.3.100 <- vdata.1.3.100[order(vdata.1.3.100$Ranking.x),]
vdata.2.3.100 <- merge(v2.100, v3.100, by = "unicode")
vdata.2.3.100 <- vdata.2.3.100[order(vdata.2.3.100$Ranking.x),]

P1 <- ggplot(vdata.1.2.100, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("V1") + ylab("V2") + ggtitle("Sensitivity Analysis of Carrying Capacity at 100")
print(P1)
P2 <- ggplot(vdata.1.3.100, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("V1") + ylab("V3") + ggtitle("Sensitivity Analysis of Carrying Capacity at 100")
print(P2)
P3 <- ggplot(vdata.2.3.100, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("V2") + ylab("V3") + ggtitle("Sensitivity Analysis of Carrying Capacity at 100")
print(P3)

##### k = 500 #####
v1.500 <- v1.500[order(v1.500$occult.prob, decreasing = TRUE),]
v2.500 <- v2.500[order(v2.500$occult.prob, decreasing = TRUE),]
v3.500 <- v3.500[order(v3.500$occult.prob, decreasing = TRUE),]

Ranking <- c(1:dim(v1.500)[1])
v1.500 <- cbind(v1.500, Ranking)
v2.500 <- cbind(v2.500, Ranking)
v3.500 <- cbind(v3.500, Ranking)

vdata.1.2.500 <- merge(v1.500, v2.500, by = "unicode")
vdata.1.2.500 <- vdata.1.2.500[order(vdata.1.2.500$Ranking.x),]
vdata.1.3.500 <- merge(v1.500, v3.500, by = "unicode")
vdata.1.3.500 <- vdata.1.3.500[order(vdata.1.3.500$Ranking.x),]
vdata.2.3.500 <- merge(v2.500, v3.500, by = "unicode")
vdata.2.3.500 <- vdata.2.3.500[order(vdata.2.3.500$Ranking.x),]

P4 <- ggplot(vdata.1.2.500, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("V1") + ylab("V2") + ggtitle("Sensitivity Analysis of Carrying Capacity at 500")
print(P4)
P5 <- ggplot(vdata.1.3.500, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("V1") + ylab("V3") + ggtitle("Sensitivity Analysis of Carrying Capacity at 500")
print(P5)
P6 <- ggplot(vdata.2.3.500, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) + 
  xlab("V2") + ylab("V3") + ggtitle("Sensitivity Analysis of Carrying Capacity at 1000")
print(P6)

##### k = 1000 #####
v1.1000 <- v1.1000[order(v1.1000$occult.prob, decreasing = TRUE),]
v2.1000 <- v2.1000[order(v2.1000$occult.prob, decreasing = TRUE),]
v3.1000 <- v3.1000[order(v3.1000$occult.prob, decreasing = TRUE),]

Ranking <- c(1:dim(v1.1000)[1])
v1.1000 <- cbind(v1.1000, Ranking)
v2.1000 <- cbind(v2.1000, Ranking)
v3.1000 <- cbind(v3.1000, Ranking)

vdata.1.2.1000 <- merge(v1.1000, v2.1000, by = "unicode")
vdata.1.2.1000 <- vdata.1.2.1000[order(vdata.1.2.1000$Ranking.x),]
vdata.1.3.1000 <- merge(v1.1000, v3.1000, by = "unicode")
vdata.1.3.1000 <- vdata.1.3.1000[order(vdata.1.3.1000$Ranking.x),]
vdata.2.3.1000 <- merge(v2.1000, v3.1000, by = "unicode")
vdata.2.3.1000 <- vdata.2.3.1000[order(vdata.2.3.1000$Ranking.x),]

P7 <- ggplot(vdata.1.2.1000, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("V1") + ylab("V2") + ggtitle("Sensitivity Analysis of Carrying Capacity at 1000")
print(P7) 
P8 <- ggplot(vdata.1.3.1000, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) + 
  xlab("V1") + ylab("V3") + ggtitle("Sensitivity Analysis of Carrying Capacity at 1000")
print(P8)
P9 <- ggplot(vdata.2.3.1000, aes(Ranking.x, Ranking.y)) + geom_point(alpha = 5/10) +
  xlab("V2") + ylab("V3") + ggtitle("Sensitivity Analysis of Carrying Capacity at 1000")
print(P9) 
