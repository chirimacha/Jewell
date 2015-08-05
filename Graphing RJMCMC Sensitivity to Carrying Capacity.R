# Set working directory to location of Jewell Data for Sensitivy Analyses

setwd("/Users/snutman/Documents/CHAGAS_DATA/jewell/SensAnalysis/K.Sens.Results/")



TimeNow <- function() {
  x <- format(Sys.time(), "%b-%d-%Y_%H-%M-%S")
  return(x)
}


files <- list.files(pattern="*Results*")

#read data
for (i in 1:length(files)) {
  test<-unlist(strsplit(files[i],".1"))
  name=gsub("0","",test[2],ignore.case=TRUE)
  name=paste("data",test[1],name,sep="_")
  data <-data.frame(read.table(files[i]))
  data <-data[order(data$unicode),]
  data$unicode <-as.character(data$unicode)
  assign(name,data,envir = .GlobalEnv)
}


#Average Chains
probmean100 <- ((data_100_ResultsK$occult.prob + data_100_ResultsKv2$occult.prob + data_100_ResultsKv3$occult.prob)/3)
probmean500 <- ((data_500_ResultsK$occult.prob + data_500_ResultsKv2$occult.prob + data_500_ResultsKv3$occult.prob)/3)
probmean1000 <- ((data_1000_Results$occult.prob + data_1000_Resultsv2$occult.prob + data_1000_Resultsv3$occult.prob)/3)


MakeRankings<-function(datain,results) {
  datain$probmean<-results
  v <- datain[order(datain$probmean, decreasing = TRUE),]
  Ranking <- c(1:dim(v)[1])
  averageRanking <- cbind(v, Ranking)
return(averageRanking)
}
 
data100 <-MakeRankings(data_100_ResultsK,probmean100)
data500<- MakeRankings(data_500_ResultsK,probmean500)
data1000<-MakeRankings(data_1000_Results,probmean1000)

#Take the top 100 in the 1000 carrying capacity dataset
data1000_subset <-data1000[1:100,]

#For the other two datasets find the unicodes in the subset 1000-carrying capacity dataset
sel <- which(data500$unicode %in% data1000_subset$unicode)
data500_subset <-data500[sel,]
sel <- which(data100$unicode %in% data1000_subset$unicode)
data100_subset <-data100[sel,]


#For graphing purposes we want to separate the top 20 out from the rest 
data1000_subset$ranknum<-as.numeric(data1000_subset$Ranking)
sel <-which(data1000_subset$Ranking>20)
data1000_subset$ranknum[sel]<-100


#Create coloring variables for graphing and merge with 500 and 100 datasets
formerge<-as.data.frame(cbind(as.character(data1000_subset$unicode),data1000_subset$ranknum))
colnames(formerge)<-c("unicode","ranknum")
formerge$unicode<-as.character(formerge$unicode)
formerge$ranknum<-as.numeric(as.character(formerge$ranknum))

data500_subset<-merge(data500_subset,formerge,by="unicode")
data100_subset<-merge(data100_subset,formerge,by="unicode")

#Put in dummy xcoordinates for graphing
data1000_subset$Xcoord<-1000
data500_subset$Xcoord<-500
data100_subset$Xcoord<-100

#Find max ranking in data for Y axis
maxranking=max(max(data1000_subset$Ranking),max(data500_subset$Ranking),max(data100_subset$Ranking))

#Create master graphing dataset
test<-rbind(data100_subset,data1000_subset,data500_subset)
test$rankcolor<-as.factor(test$ranknum)
test$Xcoord<-as.factor(test$Xcoord)

#plot
x<-qplot(test$Xcoord, test$Ranking, ylim=rev(c(0,maxranking)),cex=2,xlab="Carrying Capacity",ylab="Ranking")+
    theme(legend.position="none")+aes(color=test$rankcolor)+scale_color_manual(values = c("red", "#FF4900",
    "#FF9200","#FFDB00","#DBFF00","#92FF00","#49FF00","green","#00FF49","#00FF92","#00FFDB","#00DBFF","#0092FF","#0049FF",
      "blue","#4900FF","#9200FF","#DB00FF","#FF00DB","#FF0092","black")) +labs(title="Ranking Sensitivity to Carrying Capacity")


#export
pdf(paste("carrying_capacity_sensitivities_",TimeNow(),".pdf",sep=""))
print(x)
dev.off()


