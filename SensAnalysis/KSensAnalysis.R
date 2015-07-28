v1.100 <- data.frame(read.table("100.1000000ResultsKv1"))
v2.100 <- data.frame(read.table("100.1000000ResultsKv2"))
v3.100 <- data.frame(read.table("100.1000000ResultsKv3"))

v1.100 <- v1.100[order(v1.100$unicode),]
v2.100 <- v2.100[order(v2.100$unicode),]
v3.100 <- v3.100[order(v3.100$unicode),]

probmean <- ((v1.100$occult.prob + v2.100$occult.prob + v3.100$occult.prob)/3)

v1.100$probmean <- probmean
v.100 <- v1.100[order(v1.100$probmean, decreasing = TRUE),]

Ranking <- c(1:dim(v.100)[1])
v.100 <- cbind(v.100, Ranking)

v1.500 <- data.frame(read.table("500.1000000ResultsKv1"))
v2.500 <- data.frame(read.table("500.1000000ResultsKv2"))
v3.500 <- data.frame(read.table("500.1000000ResultsKv3"))

v1.500 <- v1.500[order(v1.500$unicode),]
v2.500 <- v2.500[order(v2.500$unicode),]
v3.500 <- v3.500[order(v3.500$unicode),]

probmean2 <- ((v1.500$occult.prob + v2.500$occult.prob + v3.500$occult.prob)/3)
v1.500$probmean2 <- probmean2

v.500 <- v1.500[order(v1.500$probmean2, decreasing = TRUE),]
v.500 <- cbind(v.500, Ranking)

vdata <- merge(v.100, v.500, by = "unicode")
vdata <- vdata[order(vdata$Ranking.x),]

plot(vdata$Ranking.x, vdata$Ranking.y)
