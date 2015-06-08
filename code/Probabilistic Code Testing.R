library(reshape2)
library(ggplot2)
library(pROC)

setwd(paste(Sys.getenv("SPATIAL_UNCERTAINTY"), "/code", sep=""))

params_block <- list(
  sim_prob_shape1=0.8,
  sim_prob_shape2=0.8,
  denuncias_path="../data/surveillance/vigilancia_denuncias.csv",
  bases_tiabaya_casa_path="../data/bases_tiabaya/Tiabaya_Points_Casa-blocks.csv",
  byManz_fullEID_path="../data/corentinEID/byManz_fullEID.csv",
  byHouse_fullEID_path="../data/corentinEID/byHouse_fullEID.csv"
)

z_analyze_manzanas <- function(params) {
  #copied from Corentin
  #byManz$Den1<-byManz$id_manz %in% byHouse$id_manz[byHouse$nbDen1]
  byManz <- read.csv(params[["byManz_fullEID_path"]])
  
  #byHouse <- read.csv(params[["byHouse_fullEID_path"]])
  # A<-glm((Den1>0)~log(nbPos+1)*log(Unspray+1)*ageAP+log(Total+1)+as.factor(D),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"aic:",A$aic,"\n")
  manz.model <-glm(Den1~log(nbPos+1)*log(Unspray+1)*ageAP+log(Total+1)+as.factor(D),data=byManz,family=binomial());
  print(manz.model)
  print("Odds Ratios:")
  #browser()
  cat("N:",dim(byManz)[1],"aic:",manz.model$aic,"\n")
  
  byManz$DenPredicted <- mean(manz.model$fitted.values)  #WARNING:  imputation
  prob_denuncias_predicted <- manz.model$fitted.values #predict.glm(manz.model, type="response")
  byManz[names(prob_denuncias_predicted),]$DenPredicted <- prob_denuncias_predicted
  #y imput
  y_imput <- manz.model$y 
  
  #byManz[names(y_imput),]$y2 <- y_imput
  
  byManz$LogDenPredicted <- -log(byManz$DenPredicted)



  #ADD IN RANKINGS - Highest density gets highest rank
  byManz$DenPredRank <- rank(byManz$DenPredicted, ties.method="first")
  
  #write.csv(byManz, time_stamp("output/byManz", ".csv"))
  
  assign("byManz",byManz, envir=.GlobalEnv)
  assign("y_imput",y_imput, envir=.GlobalEnv)
  assign("prob_denuncias_predicted",prob_denuncias_predicted, envir=.GlobalEnv)
  
  
  


  #fitted<-assign(manz.model$fitted.values)
}  

z_analyze_manzanas(params=params_block)


byManz$yimput

byManz[names(y_imput),]$yimput<-y_imput

test<-data.frame(y_imput,names(y_imput)=)
byManz2<-merge(byManz,test,by="row.names")


###test 

  
  #diagnostics
  roccurve<-roc(manz.model$y~manz.model$fitted.values)
  auc(roccurve)
  
  #table(manz.model$y, ifelse(manz.model$fitted.values>0.14, 1, 0))
  #print(exp(coef(manz.model)))
  #print(exp(cbind(OR = coef(manz.model), confint(manz.model))))
  #print(summaryOR(manz.model),digits=3);
  
  nullmod <- glm(manz.model$y~1, family="binomial")
  mcfadden_pr2 <- 1-logLik(manz.model)/logLik(nullmod)
  print(mcfadden_pr2)
  # fine and probably easier to explain
}