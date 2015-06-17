library(reshape2)
library(ggplot2)
library(pROC)
library(xlsx)
library(stringr)
library(lme4)
library(plyr)

timeNow <- format(Sys.time(), "%b-%d-%Y_%H-%M-%S")

setwd(paste(Sys.getenv("SPATIAL_UNCERTAINTY"), "/code", sep=""))

#Set parameters
params_block <- list(
  sim_prob_shape1=0.8,
  sim_prob_shape2=0.8,
  denuncias_path="../data/surveillance/vigilancia_denuncias.csv",
  bases_tiabaya_casa_path="../data/bases_tiabaya/Tiabaya_Points_Casa-blocks.csv",
  byManz_fullEID_path="../data/corentinEID/byManz_fullEID.csv",
  byHouse_fullEID_path="../data/corentinEID/byHouse_fullEID.csv",
  inspect_path="../data/2013_surveillance_data/inspecciones.csv",
  rociado_path="../data/2013_surveillance_data/rociado.xlsx",
  block_mapping="../data/2013_surveillance_data/tiabaya.csv",
  district_restrict=24
)

time_stamp <- function(prefix="", suffix="", outmsg=TRUE) {
  #generates a unique time stamp for every data generated
  t <- format(Sys.time(), "%Y-%m-%d__%H-%M-%S");
  s <- as.integer(runif(1, max=1000))
  filename <- paste(prefix, t, s, suffix, sep="")
  if (outmsg) {
    print(filename)
  }
  return(filename)
}

#READ IN HOUSE DATA
byHouse <- read.csv(params_block$byHouse_fullEID_path)
#READ IN TIABAYA DATA
tiabaya <- read.csv(params_block$block_mapping)

#/tests for missing data/
byHouse_sel <-byHouse[ which(byHouse$D==24),]
tiabaya$UNICODE<-str_trim(paste(tiabaya$prov,tiabaya$dist,tiabaya$loc,tiabaya$viv, sep="."))
#Compare EID surveillence to Tiabaya Mapping Database
count(byHouse_sel$UNICODE %in% tiabaya$UNICODE) # 3 in EID not in Tiabaya
count(tiabaya$UNICODE %in% byHouse$UNICODE) # 557 in Tiabaya not in EID

#Bring in Denuncias Data
den<-read.csv(params_block$denuncias_path)
den<-merge(den,tiabaya,by="UNICODE", all.x=TRUE)
#Bring in inspections data
inspect<-read.csv(params_block$inspect_path)


#COPIED FROM CORENTIN
summaryOR<-function(model){
  if(class(model)[1] %in% c("mer","lme4","glmerMod","bglmerMod","lmerMod")){
    require("lme4")
    # using Wald Zs
    if(class(model)[1] == "mer"){
      coefs<-fixef(summary(model))
      tableCoefs <- summary(model)@coefs
    }else{
      coefs<-fixef(model)
      tableCoefs <- summary(model)[["coefficients"]]
    }
    se<-sqrt(diag(vcov(model)))
    
    lower<-exp(coefs + qnorm(.025)*se)
    upper<-exp(coefs + qnorm(.975)*se)
    B<-as.data.frame(cbind(exp(coefs),lower,upper))
    names(B)<-c("OR","2.5%","97.5%")
    
    pvalcol<-grep("Pr",colnames(tableCoefs))
    ztcol<-grep("value",colnames(tableCoefs))
    if(length(pvalcol==1)){
      pvalues<- tableCoefs[,pvalcol]
    }else if(length(ztcol)==1){
      pvalues<-2*pnorm(-abs(tableCoefs[,ztcol]))
    }else{
      cat("can't get the pvalue for mer model\n")
    }
    
  }else{
    coefs<-coef(model)
    B<-as.data.frame(cbind(exp(coefs),exp(confint.default(model))),digits=3)
    names(B)[1]<-"OR"
    pvalcol<-grep("Pr",colnames(coef(summary(model))))
    pvalues<-coef(summary(model))[,pvalcol]
  }
  # EDITED OUT DUE TO LACK OF STAR.ON.PVALUES FUNCTION, problem counting Na
  #add NAs
  #if(count(is.na(coefs))){
    #NAtoAdd<-which(is.na(model$coefficients))
    #for(i in 1:length(NAtoAdd)){
      #pvalues<-insert(NA,pvalues,NAtoAdd[i])
    #}
  #}
  #signifs<-star.on.pvalues(pvalues)
  #B$signifs<-signifs
  return(B)
}


Z_Plot_Dist<-function(district) {
  sel=which(byHouse$D==district)
  sel_den=which(den$D==district)
  max_lon=max(byHouse$lon[sel],den$lon[sel_den],na.rm=TRUE)
  min_lon=min(byHouse$lon[sel],den$lon[sel_den],na.rm=TRUE)
  max_lat=max(byHouse$lat[sel],den$lat[sel_den],na.rm=TRUE)
  min_lat=min(byHouse$lat[sel],den$lat[sel_den],na.rm=TRUE)
  plot(byHouse$lon[sel],byHouse$lat[sel],cex=PRED[sel]*10,ylim=c(min_lat,max_lat),xlim=c(min_lon,max_lon),col="GREEN",xlab="Latitude",ylab="Longitude")
  par(new=TRUE)
  plot(den$lon[sel_den],den$lat[sel_den],ylim=c(min_lat,max_lat),xlim=c(min_lon,max_lon),col="RED",xlab="",ylab="")
  title(main="Predicted Density vs. Denuncias",sub=paste("District:",district,sep= " "))
  }


Z_Plot_Loc<-function(district,local) {
  sel=which(byHouse$D==district & byHouse$L==local)
  sel_den=which(den$D==district & den$L==local)
  max_lon=max(byHouse$lon[sel],den$lon[sel_den],na.rm=TRUE)
  min_lon=min(byHouse$lon[sel],den$lon[sel_den],na.rm=TRUE)
  max_lat=max(byHouse$lat[sel],den$lat[sel_den],na.rm=TRUE)
  min_lat=min(byHouse$lat[sel],den$lat[sel_den],na.rm=TRUE)
  plot(byHouse$lon[sel],byHouse$lat[sel],cex=PRED[sel]*10,ylim=c(min_lat,max_lat),xlim=c(min_lon,max_lon),col="GREEN",xlab="Latitude",ylab="Longitude")
  par(new=TRUE)
  plot(den$lon[sel_den],den$lat[sel_den],ylim=c(min_lat,max_lat),xlim=c(min_lon,max_lon),col="RED",xlab="",ylab="")
  title(main="Predicted Density vs. Denuncias",sub=paste("District:",district,"Locality:",local,sep= " "))
  }

#COPIED FROM CORENTIN (Data name is byHouse not dat)
    A<-glmer(Insp1Pos>0 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=byHouse,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(byHouse)[1],"(",count(byHouse$Insp1Pos>0),"+) aic:",AIC(A),"\n")
    #                                  OR     2.5%    97.5% signifs
    # (Intercept)                2.11e-05 3.77e-06 1.18e-04     ***
    # infFAstatus > 0TRUE        3.21e+00 2.02e+00 5.10e+00     ***
    # rocFAstatus == 0TRUE       2.15e+01 4.02e+00 1.15e+02     ***
    # ageAP                      1.55e+00 1.19e+00 2.02e+00      **
    # rocFAstatus == 0TRUE:ageAP 7.28e-01 5.47e-01 9.68e-01       *
    # N: 56491 ( 116 +) aic: 1415.484 
    
    #PREDICT INFESTATION
    PRED<-exp(predict(A))
    #summary(PRED)
#  Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 1.759e-05 1.500e-04 3.544e-04 1.826e-03 7.961e-04 2.031e-01 


pdf(paste("output/Pred_density_vs_denuncias_",timeNow,".pdf",sep=""))
Z_Plot_Dist(24)
for (n in 1:14) {
  Z_Plot_Loc(24,n)
}
dev.off()

    
  
  
  
 

