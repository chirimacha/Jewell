# install.packages("Rlib/lme4_1.1-6.tar.gz",lib="Rlib/")

source("spatcontrol/spatcontrol.R",chdir=TRUE)
library(lme4)
pour<-function(num) paste("(",round(100*num,2),"%)",sep="")
nameSimul <- "glmVigilance"
#===================
# import data
#===================
load("byHouseVig.img")
load("rocvig.img")
load("commonVig.img")
load("commonVigByYearManz.img")
load("byManzVig.img")
load("byLocRociado.img")
den<-read.csv("DENUNCIAS_30-01-2013.csv",sep=",")

#===================
# Complementary columns
#===================

byHouse$nbDen1<-byHouse$UNICODE %in% den$UNICODE[den$TYPO_DENUNCIA==1]
byHouse$Pos1 <- byHouse$UNICODE %in% den$UNICODE[den$numero.registro %in% commonVig$NRO_DENUNCIA[commonVig$TYPO_DENUNCIA==1 & commonVig$PosVig] & den$numero.registro>0]
byHouse$Insp1<-byHouse$UNICODE %in% commonVig$UNICODE[commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1]
byHouse$Insp1Pos<-byHouse$UNICODE %in% commonVig$UNICODE[commonVig$INSP_POSITIVA==1 & commonVig$TYPO_DENUNCIA==1]
den<-merge(den,byHouse[,c("UNICODE","id_manz")],all.x=TRUE)


#===================
# Restrict sample
#===================
#-------------------
# Columns needed for restriction of the sample
#-------------------
rocvig$codeLoc<-with(rocvig,paste(P,D,L,sep="."))
commonVig$codeLoc<-with(commonVig,paste("1",D,L,sep="."))
commonVigByYearManz$codeLoc<-with(commonVigByYearManz,paste("1",D,L,sep="."))
byManz$codeLoc<-with(byManz,paste(P,D,L,sep="."))
byHouse$codeLoc<-with(byHouse,paste(P,D,L,sep="."))
den$codeLoc<-with(den,paste(P,D,L,sep="."))

# -------------------
# Define sample
#-------------------
# in locality sprayed before 2012
LocFocus<-byLoc$codeLoc[which(byLoc$Ispray>0 & byLoc$IIspray>0 &byLoc$YearLastII<12)]

# in city block with at least one house sprayed in I and one house 
# sprayed in II
ManzFocus<-byManz$id_manz[byManz$nbSprayI>0 & byManz$nbSprayII>0]

# only households that existed in AFII
HouseFocus<-byHouse$UNICODE[!is.na(byHouse$AFIIstatus)]

DFocus<-unique(byHouse$D[byHouse$codeLoc %in% LocFocus])
#-------------------
# Flag sample
#-------------------

# both Loc and manz restricted
byHouse$focus<-byHouse$codeLoc %in% LocFocus & byHouse$id_manz %in% ManzFocus & byHouse$UNICODE %in% HouseFocus
rocvig$focus<-rocvig$codeLoc %in% LocFocus & rocvig$id_manz %in% ManzFocus & rocvig$UNICODE %in% HouseFocus
commonVig$focus<-commonVig$codeLoc %in% LocFocus & commonVig$id_manz %in% ManzFocus & commonVig$UNICODE %in% HouseFocus
commonVigByYearManz$focus<-commonVigByYearManz$codeLoc %in% LocFocus & commonVigByYearManz$id_manz %in% ManzFocus 
byManz$focus<-byManz$codeLoc %in% LocFocus & byManz$id_manz %in% ManzFocus  
den$focus<-den$id_manz %in% ManzFocus

# # only Loc and manz restricted
# byHouse$focus<-byHouse$codeLoc %in% LocFocus
# rocvig$focus<-rocvig$codeLoc %in% LocFocus
# commonVig$focus<-commonVig$codeLoc %in% LocFocus
# commonVigByYearManz$focus<-commonVigByYearManz$codeLoc %in% LocFocus
# byManz$focus<-byManz$codeLoc %in% LocFocus

#-------------------
# Describe sample restriction 
#-------------------
# plot
par(mfrow=c(2,2))
with(byHouse,plot_reel(X,Y,Iinfested,base=0,main="By Household"))
with(byManz,plot_reel(X,Y,nbPos>0,base=0,main="By city-block"))
distsH<-cumsum(seq(0,90,7))
hAutoCor<-which(!is.na(byHouse$X)&!is.na(byHouse$Y) & !is.na(byHouse$id_manz) & (byHouse$Ispray | byHouse$IIspray) & byHouse$D %in% c(14,10))
# # sample autocorrelation
# mats_neighH<-with(byHouse[hAutoCor,],gen.mats.neigh(distsH,X,Y,group=id_manz))
# SMI<-structured.moransI(mats_neighH,byHouse[hAutoCor,"infFAstatus"]>0,nb_rep_sign=30)
# plot.structured.moransI(SMI)

# limit the reports to the district actually in surveillance
totRep<-count(den$TYPO_DENUNCIA==1 & den$D %in% DFocus)
totRepLoc<-count(den$TYPO_DENUNCIA==1 & den$codeLoc %in% LocFocus)
totRepFoc<-count(den$TYPO_DENUNCIA==1 & den$focus)
repTot<-c(totRep,totRepLoc,totRepFoc)

totRepH<-length(unique(den$UNICODE[den$TYPO_DENUNCIA==1 & den$D %in% DFocus]))
totRepHtargLoc<-length(unique(den$UNICODE[den$TYPO_DENUNCIA==1 & den$codeLoc %in% LocFocus]))
totRepHtarg<-count(unique(den$UNICODE[den$TYPO_DENUNCIA==1]) %in% byHouse$UNICODE[byHouse$focus])

# cat("Since 2009",totRep,"community reports from",totRepH,"households",totRepHtarg,pour(totRepHtarg/totRepH),"of which in city-blocks targeted in attack phase \n")

totInsp1<-count(commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1 & commonVig$D %in% DFocus)
totInsp1LocFoc<-count(commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1 & commonVig$codeLoc %in% LocFocus)
totInsp1Foc<-count(commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1 & commonVig$focus)
totInsp1s<-c(totInsp1,totInsp1LocFoc,totInsp1Foc)

totInsp1Pos<-count(commonVig$INSP_POSITIVA==1 & commonVig$TYPO_DENUNCIA==1 & commonVig$D %in% DFocus)
totInsp1PosLocFoc<-count(commonVig$INSP_POSITIVA==1 & commonVig$TYPO_DENUNCIA==1 & commonVig$codeLoc %in% LocFocus)
totInsp1PosFoc<-count(commonVig$INSP_POSITIVA==1 & commonVig$TYPO_DENUNCIA==1 & commonVig$focus)
totInsp1sPos<-c(totInsp1Pos,totInsp1PosLocFoc,totInsp1PosFoc)



totInsp1H<-length(unique(commonVig$UNICODE[commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1]))

# cat("We performed",totInsp1,"inspections in",totInsp1H,"households \n")  


households<-c(count(byHouse$D %in% DFocus & !is.na(byHouse$AFIstatus)),
	      count(byHouse$codeLoc %in% LocFocus & !is.na(byHouse$AFIstatus)),
	      length(unique(byHouse$UNICODE[byHouse$focus])))
cbs<-c(length(unique(byHouse$id_manz[byHouse$D %in% DFocus & !is.na(byHouse$AFIstatus)])),
       length(unique(byHouse$id_manz[byHouse$codeLoc %in% LocFocus & !is.na(byHouse$AFIstatus)])),
       length(unique(byHouse$id_manz[byHouse$focus & !is.na(byHouse$AFIstatus)])))
target<-c("District","Locality","City Block")
repH<-c(totRepH,totRepHtargLoc,totRepHtarg)

# fisher tests on reports
ORloc<-with(byHouse[byHouse$D %in% DFocus,],fisher.test(codeLoc %in% LocFocus,nbDen1))
ORcb<-with(byHouse[byHouse$codeLoc %in% LocFocus,],fisher.test(focus,nbDen1))

ORtarget<-c(NA,ORloc$estimate,ORcb$estimate)
names(ORtarget)<-c()
ORtarget<-signif(as.numeric(ORtarget),3)
ORpvalue<-signif(as.numeric(c("",ORloc$p.value,ORcb$p.value)),3)
ORstars<-c("",star.on.pvalues(ORloc$p.value),star.on.pvalues(ORcb$p.value))

# overall: focus a risk factor in the district? No.
library(lme4)
# A<-glmer(nbDen1~ focus +(1|D),data=byHouse[byHouse$D %in% DFocus,],family=binomial())
# Generalized linear mixed model fit by the Laplace approximation 
# Formula: nbDen1 ~ focus + (1 | D) 
#    Data: byHouse[byHouse$D %in% DFocus, ] 
#   AIC  BIC logLik deviance
#    3597 3624  -1795     3591
#   Random effects:
#    Groups Name        Variance Std.Dev.
#    D      (Intercept) 0.98107  0.99049 
#    Number of obs: 80039, groups: D, 8
# 
#    Fixed effects:
#                Estimate Std. Error z value Pr(>|z|)    
#    (Intercept)  -5.8056     0.3774 -15.384   <2e-16 ***
#    focusTRUE     0.1169     0.1536   0.761    0.447    
#    ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
#    Correlation of Fixed Effects:
#              (Intr)
#    focusTRUE -0.294

byHouse$Ltarget<-byHouse$codeLoc %in% LocFocus

# glmer(nbDen1~ Ltarget +(1|D),data=byHouse[byHouse$D %in% DFocus,],family=binomial())
# Generalized linear mixed model fit by the Laplace approximation 
# Formula: nbDen1 ~ Ltarget + (1 | D) 
#    Data: byHouse[byHouse$D %in% DFocus, ] 
#   AIC  BIC logLik deviance
#    3596 3624  -1795     3590
#   Random effects:
#    Groups Name        Variance Std.Dev.
#    D      (Intercept) 0.97427  0.98705 
#    Number of obs: 80039, groups: D, 8
# 
#    Fixed effects:
#                Estimate Std. Error z value Pr(>|z|)    
#    (Intercept)  -5.9217     0.3997 -14.817   <2e-16 ***
#    LtargetTRUE   0.2348     0.2021   1.162    0.245    
#    ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
#    Correlation of Fixed Effects:
#                (Intr)
#    LtargetTRUE -0.437

## restriction on households inspected/positive
A<-glmer(Insp1Pos ~ focus+(1|Ltarget),data=byHouse[byHouse$Insp1>0,],family=binomial())
# Generalized linear mixed model fit by the Laplace approximation 
# Formula: Insp1Pos ~ focus + (1 | Ltarget) 
#    Data: byHouse[byHouse$Insp1 > 0, ] 
#    AIC   BIC logLik deviance
#     786.3 800.4 -390.2    780.3
#    Random effects:
#     Groups  Name        Variance Std.Dev.
#     Ltarget (Intercept)  0        0      
#     Number of obs: 799, groups: Ltarget, 2
# 
#     Fixed effects:
#                 Estimate Std. Error z value Pr(>|z|)    
#     (Intercept) -1.39303    0.18368  -7.584 3.35e-14 ***
#     focusTRUE   -0.06201    0.21065  -0.294    0.768    
#     ---
#     Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
#     Correlation of Fixed Effects:
#               (Intr)
#     focusTRUE -0.872
#     > summaryOR(A)
#                        OR      2.5%     97.5% signifs
#     (Intercept) 0.2483221 0.1732473 0.3559298     ***
#     focusTRUE   0.9398757 0.6219691 1.4202736        

Insps<-with(byHouse,c(count(Insp1>0 & D %in% DFocus),count(byHouse$Insp1>0 & byHouse$codeLoc %in% LocFocus),
		      count(Insp1>0 & focus)))
InspsPos<-with(byHouse,c(count(Insp1Pos>0 & D %in% DFocus),count(byHouse$Insp1Pos>0 & byHouse$codeLoc %in% LocFocus),
			 count(Insp1Pos>0 & focus)))

print(as.data.frame(cbind(target,cbs,households,repTot,repH,ORtarget,ORpvalue,ORstars,totInsp1s,Insps,totInsp1sPos,InspsPos),digits=3))

aimForDen<-function(denNum){
count(commonVig$NRO_DENUNCIA == denNum)
}
den$nbInspPlan<-sapply(den$numero.registro,aimForDen)
inspForDen<-function(denNum){
count(commonVig$NRO_DENUNCIA == denNum & commonVig$INSP_COMPLETA ==1)
}
den$nbInsp<-sapply(den$numero.registro,inspForDen)
infForDen<-function(denNum){
count(commonVig$NRO_DENUNCIA == denNum & commonVig$INSP_POSITIVA ==1)
}
den$nbInspPos<-sapply(den$numero.registro,infForDen)

# observed number of infested arround a report type 1
den$FECHA<-as.Date(den$FECHA)
den$TypeInsp<-0
den$TypeInsp[den$FECHA>as.Date("2009-04-01")]<-1
den$TypeInsp[den$FECHA>as.Date("2012-03-01")]<-2
with(den[den$TYPO_DENUNCIA==1 & den$codeLoc %in% LocFocus,],table(nbInspPlan,TypeInsp))
with(den[den$TYPO_DENUNCIA==1 & den$codeLoc %in% LocFocus,],table(nbInsp,TypeInsp))

# # problems: without observations, 0 or one inspected
# # for each plot the city block as well as reports and inspected
# probDen<-den[which(den$TypeInsp==1 & den$TYPO_DENUNCIA==1 & den$codeLoc %in% LocFocus & den$OBSERVACIONES=="NULL"),"numero.registro"]
# # probDen<-den[which(den$TypeInsp==2 & den$TYPO_DENUNCIA==1 & den$codeLoc %in% LocFocus &den$nbInspPlan %in% c(0,1,2) & den$OBSERVACIONES=="NULL"),"numero.registro"]
# par(mfrow=c(4,6))
# i<-1
# for(regden in probDen){
# 	idDenDen<-which(den$numero.registro==regden)
# 	uniDen<-den$UNICODE[idDenDen]
# 	idDenBH<-which(byHouse$UNICODE==uniDen)
# 	yearDen<-den$ANIO[idDenDen]
# 	id_manz<-byHouse$id_manz[idDenBH]
# 	with(byHouse[byHouse$id_manz==id_manz,],plot(X,Y,main=paste(den$FECHA[idDenDen],den$OBSERVACIONES[idDenDen]),asp=1))
# 	with(byHouse,text(X[idDenBH],Y[idDenBH],"R"))
# 	with(byHouse[byHouse$id_manz==id_manz & byHouse$UNICODE %in% commonVig$UNICODE[commonVig$NRO_DENUNCIA==regden],],lines(X,Y,col=4,type="p"))
# 	i<-i+1
# 	if(i%%24==1){
# 		dev.new()
# 		par(mfrow=c(4,6))
# 	}
# }

# plot.loc.arequipa("^1.13.*")
# with(byHouseAll[byHouseAll$D == 13,],lines(X,Y,pch=".",type="p",cex=1,asp=1,col="black"))
# # with(byHouseAll[byHouseAll$D %in% DFocus,],lines(X,Y,pch=".",type="p",cex=1,asp=1,col="black"))
# with(byHouseAll[byHouseAll$D == 13 & byHouseAll$codeLoc %in% LocFocus,],lines(X,Y,pch=".",type="p",asp=1,col="blue"))
# with(byHouseAll[byHouseAll$D == 13 & byHouseAll$focus,],lines(X,Y,pch=".",type="p",asp=1,col="yellow"))


#--------------------
# Restrict to sample
#--------------------
byHouseAll<-byHouse
byHouse<-byHouse[byHouse$focus,]
rocvig<-rocvig[rocvig$focus,]
commonVig<-commonVig[commonVig$focus,]
commonVigByYearManz<-commonVigByYearManz[commonVigByYearManz$focus,]
byManz<-byManz[byManz$focus,]

cat("#===================
# By Year data table
#===================
")

AFI<-with(byHouse,table(Ispray,YearLastAFII,Iinfested))
AFI<-rbind(AFI[2,,2],AFI[2,,1],AFI[1,,1])
rownames(AFI)<-c("AFI+","AFI-","noAFI")

AFII<-with(byHouse,table(IIspray,YearLastAFII,IIinfested))
AFII<-rbind(AFII[2,,2],AFII[2,,1],AFII[1,,1])
# rownames(AFII)<-c("+","-","null")
rownames(AFII)<-c("AFII+","AFII-","noAFII")


surv<-with(byHouse,table(nbInsp>0|nbRoc>0,YearLastAFII,nbInspPos>0|nbRocPos>0))
surv<-rbind(surv[2,,2],surv[2,,1],surv[1,,1])
rownames(surv)<-c("S+","S-","noS")
# rownames(surv)<-c("+","-","null")

byYearTot<-rbind(AFI,AFII,surv)
colnames(byYearTot)<-as.integer(colnames(byYearTot))+2000
print(byYearTot)

# library(odfWeave)
# file.in<-"~/Documents/recherche/Penn/ArticleBaseVigilance/TestODFweave.odt"
# file.out<-"~/Documents/recherche/Penn/ArticleBaseVigilance/TestODFweave_out.odt"
# byYearTot<-as.data.frame(byYearTot)
# odfWeave(file.in, file.out)
# # tables are not processed correctly
# stop()

cat("#===================
# By Year and AF infestation table
#===================\n")
    # status: positive:0
    #         negative:1
    #         NA      :2
# both inspections and rociado
# byHouse$Sstatus<-as.integer(byHouse$nbInsp>0 | byHouse$nbRoc>0)
# byHouse$Sstatus[byHouse$nbRocPos>0 | byHouse$nbInspPos>0]<-2
# byHouse$Sstatus<-2-byHouse$Sstatus

# only inspections 
byHouse$Sstatus<-as.integer(byHouse$Insp1>0)
byHouse$Sstatus[byHouse$Insp1Pos>0]<-2
byHouse$Sstatus<-2-byHouse$Sstatus

byHouse$AFstatus<-as.integer(byHouse$rocFAstatus>0)
byHouse$AFstatus[byHouse$infFAstatus>0]<-2
byHouse$AFstatus<- 2-byHouse$AFstatus

SbyYearInfAF<-with(byHouse,table(AFstatus,YearLastAFII,Sstatus))
SbyYearInfAF<-rbind(SbyYearInfAF[,,1],SbyYearInfAF[,,2],SbyYearInfAF[,,3])
colnames(SbyYearInfAF)<-as.integer(colnames(SbyYearInfAF))+2000
SbyYearInfAF<- as.data.frame(SbyYearInfAF,row.names=seq(1,9))
SbyYearInfAF$Total<-apply(SbyYearInfAF,1,sum)

Surv<-rep(" ",9)
Surv[2]<- "+"
Surv[5]<- "-"
Surv[8]<- "∅"

AFtreat<-rep(c("+","-","∅"),3)

SbyYearInfAF<-cbind(Surv,AFtreat,SbyYearInfAF)
print(SbyYearInfAF,quote=TRUE)


cat("#===================
# By Year and AF infestation table (chronological order)
#===================\n")

SbyYearInfAF<-with(byHouse,table(Sstatus,YearLastAFII,AFstatus))
SbyYearInfAF<-rbind(SbyYearInfAF[,,1],SbyYearInfAF[,,2],SbyYearInfAF[,,3])
colnames(SbyYearInfAF)<-as.integer(colnames(SbyYearInfAF))+2000
SbyYearInfAF<- as.data.frame(SbyYearInfAF,row.names=seq(1,9))
SbyYearInfAF$Total<-apply(SbyYearInfAF,1,sum)

AFtreat<-rep(" ",9)
AFtreat[2]<- "+"
AFtreat[5]<- "-"
AFtreat[8]<- "∅"

Surv<-rep(c("+","-","∅"),3)

SbyYearInfAF<-cbind(AFtreat,Surv,SbyYearInfAF)
print(SbyYearInfAF,quote=TRUE)



cat("#===================
# By Year and AF treatments table
#===================
")
SbyYearTreatAF<-with(byHouse,table(rocFAstatus,YearLastAFII,Sstatus))
SbyYearTreatAF<-rbind(SbyYearTreatAF[,,1],SbyYearTreatAF[,,2],SbyYearTreatAF[,,3])
colnames(SbyYearTreatAF)<-as.integer(colnames(SbyYearTreatAF))+2000
SbyYearTreatAF<- as.data.frame(SbyYearTreatAF,row.names=seq(1,9))
SbyYearTreatAF$Total<-apply(SbyYearTreatAF,1,sum)

Surv<-rep(" ",9)
Surv[2]<- "+"
Surv[5]<- "-"
Surv[8]<- "∅"

AFtreat<-rep(c("0","1","2"),3)

SbyYearTreatAF<-cbind(Surv,AFtreat,SbyYearTreatAF)
print(SbyYearTreatAF)
# based on byHouse
cat("A total of",count(byHouse$Sstatus!=2),"households were inspected in surveillance and",
    count(byHouse$Sstatus==0),pour(count(byHouse$Sstatus==0)/count(byHouse$Sstatus!=2)),
    "were infested, of which",with(byHouse,count(Sstatus==0 & Ispray & IIspray)),
    pour(with(byHouse,count(Sstatus==0 & Ispray & IIspray))/count(byHouse$Sstatus==0)),
    "were sprayed twice",with(byHouse,count(Sstatus==0 & ((Ispray & !IIspray) | (!Ispray & IIspray)))),
    pour(with(byHouse,count(Sstatus==0 & ((Ispray & !IIspray) | (!Ispray & IIspray)))/count(Sstatus==0))),
    "were sprayed once and",with(byHouse,count(Sstatus==0 & (!Ispray & !IIspray))),
    pour(with(byHouse,count(Sstatus==0 & (!Ispray & !IIspray))/count(Sstatus==0))),"never sprayed\n")

# based on commonVig
dat<-commonVig[(commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1 & commonVig$id_manz %in% byManz$id_manz[byManz$nbSpray>0]),]
cat("A total of",count(dat$INSP_COMPLETA==1),"inspections were conducted in",length(unique(dat$UNICODE[dat$INSP_COMPLETA==1])),"households and",count(dat$INSP_POSITIVA==1),pour(count(dat$INSP_POSITIVA==1)/count(dat$INSP_COMPLETA==1)),"\n")
A<-with(byHouse[byHouse$Sstatus==0 & byHouse$rocFAstatus==2,],table(infFAstatus==0))
cat("In addition among households infested in surveillance sprayed twice in attack-phase",A[2],pour(A[2]/sum(A)),"had no history of infestation in attack-phase\n")

cat("#===================
# By Year and AF most detailed
#===================
")
SbyYearHistAF<-with(byHouse,table(AFIstatus,YearLastAFII,AFIIstatus,Sstatus))
SbyYearHistAF<-as.data.frame(rbind(SbyYearHistAF[,,1,1],SbyYearHistAF[,,2,1],SbyYearHistAF[,,3,1],SbyYearHistAF[,,1,2],SbyYearHistAF[,,2,2],SbyYearHistAF[,,3,2],SbyYearHistAF[,,1,3],SbyYearHistAF[,,2,3],SbyYearHistAF[,,3,3]),row.names=seq(1,dim(SbyYearHistAF)[1]))
names(SbyYearHistAF)<-as.integer(colnames(SbyYearHistAF))+2000
Total<-apply(SbyYearHistAF,1,sum)
SbyYearHistAF<-cbind(SbyYearHistAF,Total)
Total<-apply(SbyYearHistAF,2,sum)
SbyYearHistAF<-rbind(SbyYearHistAF,Total)

surv<-rep(" ",28)
surv[5]<-"+"
surv[14]<-"-"
surv[23]<-"∅"
AFII<-rep(" ",9)
AFII[2]<-"+"
AFII[5]<-"-"
AFII[8]<- "∅"
AFII<-c(rep(AFII,3)," ")
AFI<-c(rep(c("+","-","∅"),9),"Total")
SbyYearHistAF<-cbind(surv,AFII,AFI,SbyYearHistAF)
print(SbyYearHistAF,quote=TRUE)
write.csv(SbyYearHistAF,"SbyYearInfAF.csv")

#===================
# Analysis
#===================
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
	# add NAs
	if(count(is.na(coefs))){
		NAtoAdd<-which(is.na(model$coefficients))
		for(i in 1:length(NAtoAdd)){
			pvalues<-insert(NA,pvalues,NAtoAdd[i])
		}
	}
	signifs<-star.on.pvalues(pvalues)
	B$signifs<-signifs
	return(B)
}

# basic OR 
fisher.test(with(byHouse,table(rocFAstatus==0,nbInspPos>0)))
# 	Fisher's Exact Test for Count Data
# data:  with(byHouse, table(rocFAstatus == 0, nbInspPos > 0)) 
# p-value = 0.1007
# alternative hypothesis: true odds ratio is not equal to 1 
# 95 percent confidence interval:
#  0.9149669 2.0659882 
# sample estimates:
# odds ratio 
#   1.391934 

fisher.test(with(commonVig[commonVig$INSP_COMPLETA>0,],table(rocFAstatus==0,INSP_POSITIVA>0)))
# 	Fisher's Exact Test for Count Data
# data:  
# p-value = 1.728e-05
# alternative hypothesis: true odds ratio is not equal to 1 
# 95 percent confidence interval:
#  1.620101 3.730011 
# sample estimates:
# odds ratio 
#   2.470177 

fisher.test(with(byHouse[byHouse$Sstatus<2,],table(rocFAstatus==0,Sstatus==0)))

# # on within surveillance (commonVig) by House
# cat("----------------
#     byHouse within surveillance, all commonVig
# ------------------\n")
# A<-glm(PosVig ~ (infFAstatus>0)*(rocFAstatus==0) * (ageAP),data=commonVig[commonVig$VisitVig==1,],family=binomial());
# #=> very conclusive
# print(summaryOR(A),digits=2)
# 
# # within surveillance treatment only (rocvig) byHouse
# rocvig$YearLastAFII<-byLoc$YearLastII[rocvig$codeLoc%in% byLoc$codeLoc]
# rocvig$ageAP<-rocvig$ANIO-rocvig$YearLastAFII
# A<-glm(infestedRV ~ (infFAstatus>0)*(rocFAstatus==0) * (ageAP),data=rocvig[rocvig$complete==1,],family=binomial());
# print(summaryOR(A),digits=2)
# #=> absolutely not conclusive
# 
# # on within surveillance inspections only (commonVig) by House
# A<-glm(PosVig ~ (infFAstatus>0)*(rocFAstatus==0) * (ageAP)+as.factor(D),data=commonVig[commonVig$INSP_COMPLETA==1,],family=binomial())
# print(summaryOR(A),digits=2)
# #=> very conclusive, but infestation in AF does mater here.
# cat("----------------
#     byManz within surveillance, all commonVig
# ------------------\n")
# 
# A<-glm((PosVig>0)~(nbPos>0)+log(Unspray+1)*ageAPatVig+as.factor(D),data=commonVigByYearManz[commonVigByYearManz$VisitVig>0,],family=binomial())
# print(summaryOR(A),digits=2)
# #=> serious problem here, cannot get former results
# 
# cat("----------------
#     byHouse within all
# ------------------\n")
# 
# A<-glm((nbRocPos>0) ~ (rocFAstatus==0) * (ageAP)+(infFAstatus>0)+as.factor(D),data=byHouse[is.finite(byHouse$YearLastAFII) & !is.na(byHouse$YearLastAFII),],family=binomial());
# print(summaryOR(A),digits=2)
# 
# # or similar but less
# A<-glm((nbInspPos>0) ~ (rocFAstatus==0) * (ageAP)*(infFAstatus>0)+as.factor(D),data=byHouse,family=binomial());
# print(summaryOR(A),digits=2)
# #=> ok
# 
# cat("----------------
#     byManz within all
# ------------------\n")
# A<-glm((nbRoc>0)~log(nbPos+1)+(nbPos>0)+log(Unspray+1)*ageAP+as.factor(D),data=byManz,family=binomial())
# print(summaryOR(A),digits=2)
# #=> very good, can keep that.
# 
# cat("----------------
#     Denuncias type 1 byHouse within all
# ------------------\n")
# 
# A<-glm((nbDen1>0) ~ (rocFAstatus==0) * (ageAP)+(infFAstatus>0)+as.factor(D),data=byHouse[is.finite(byHouse$YearLastAFII) & !is.na(byHouse$YearLastAFII),],family=binomial());
# print(summaryOR(A),digits=2)
# #=> loosing the significance for the interaction
# 
# cat("----------------
#     denuncias byManz within all
# ------------------\n")
# A<-glm((DenH>0)~(nbPos>0)+log(Unspray+1)*ageAP+as.factor(D),data=byManz,family=binomial());print(summaryOR(A),digits=2);A$aic
# # The best
# 
# A<-glm((DenH>0)~log(nbPos+1)*log(Unspray+1)*ageAP+as.factor(D),data=byManz,family=binomial());print(summaryOR(A),digits=2);A$aic
# print(summaryOR(A),digits=2)
# # best for explanation in methods
# 
# cat("----------------
#     from denuncias type 1 (insp byHouse within all)/(insp by house within surv) / (den by Manz within all)
# ------------------\n")

dat<-byHouse
# A<-glm(Insp1Pos ~ (infFAstatus>0)*(rocFAstatus==0)*(ageAP)+as.factor(D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"aic:",A$aic,"\n")

# on within surveillance inspections from Den1 only (commonVig) by House
dat<-commonVig[commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1,]
# A<-glm(INSP_POSITIVA==1 ~ (infFAstatus>0)*(rocFAstatus==0)*(ageAP)+as.factor(D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"aic:",A$aic,"\n")
#=> not conclusive, but infestation in AF does mater here.

# byManz$Den1<-byManz$id_manz %in% byHouse$id_manz[byHouse$nbDen1]
# A<-glm((Den1>0)~(nbPos>0)+log(Unspray+1)*ageAP+as.factor(D),data=byManz,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"aic:",A$aic)
# # great!

byManz$Den1<-byManz$id_manz %in% byHouse$id_manz[byHouse$nbDen1]
dat<-byManz
# A<-glm((Den1>0)~log(nbPos+1)*log(Unspray+1)*ageAP+log(Total+1)+as.factor(D),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"aic:",A$aic,"\n")
# fine and probably easier to explain

# Can also be delt with using real random effect (mixed models)
# similar results but takes forever
# this one doesn't work
# dat<-byHouse
# A<-glmer(Insp1Pos ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# the infestation in inspection is linked with a log of things
# controlling by codeLoc: too detailed, loosing signif
# controlling by D and id_manz a bit too much
# controlling for D ok
# controlling for id_manz a bit too much 

# A<-glmer(nbDen1 ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# the reports in themselves are not that directly linked with
# infestation or anything in attack-phase at individual level

### P(report~ infestation*age*treatment+cityBlockRandomEffect)
# dat<-byHouse
# among all households localities treated in AP
# > A<-glmer(nbDen1 ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|id_manz),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# Warning message:
# In mer_finalize(ans) : false convergence (8)
#                                 OR    2.5%   97.5% signifs
# (Intercept)                0.00006 1.5e-05 2.3e-04     ***
# infFAstatus > 0TRUE        1.34237 1.9e-01 9.6e+00        
# ageAP                      1.34574 1.1e+00 1.7e+00      **
# rocFAstatus == 0TRUE       1.82684 3.3e-01 1.0e+01        
# infFAstatus > 0TRUE:ageAP  1.10687 8.2e-01 1.5e+00        
# ageAP:rocFAstatus == 0TRUE 0.96402 7.2e-01 1.3e+00        
# N: 67218 aic: 2671.937 

# among households in city-blocks treated in AP
dat<-byHouse[which(byHouse$id_manz %in% byManz$id_manz[byManz$nbSpray>0]),]
## city block random effect
# A<-glmer(nbDen1 ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|id_manz),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# No convergence 

cat("-----------------
    Infestation among all households
-----------------\n")

dat<-byHouse
# A<-glmer(Insp1Pos>0 ~ (infFAstatus>0)+(rocFAstatus==0)+(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(",count(dat$Insp1Pos>0),"+) aic:",AIC(A),"\n")
#                           OR    2.5%   97.5% signifs
# (Intercept)          4.6e-05 9.3e-06 0.00023     ***
# infFAstatus > 0TRUE  3.3e+00 2.0e+00 5.30698     ***
# rocFAstatus == 0TRUE 3.7e+00 2.1e+00 6.45201     ***
# ageAP                1.3e+00 1.0e+00 1.74005       *
# N: 56491 ( 116 +) aic: 1418.837 

dat<-byHouse
write.csv(byHouse,file="byHouse_fullEID.csv")
A<-glmer(Insp1Pos>0 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(",count(dat$Insp1Pos>0),"+) aic:",AIC(A),"\n")
#                                  OR     2.5%    97.5% signifs
# (Intercept)                2.11e-05 3.77e-06 1.18e-04     ***
# infFAstatus > 0TRUE        3.21e+00 2.02e+00 5.10e+00     ***
# rocFAstatus == 0TRUE       2.15e+01 4.02e+00 1.15e+02     ***
# ageAP                      1.55e+00 1.19e+00 2.02e+00      **
# rocFAstatus == 0TRUE:ageAP 7.28e-01 5.47e-01 9.68e-01       *
# N: 56491 ( 116 +) aic: 1415.484 
MainModel <- A
MainModelNull <- glmer(Insp1Pos>0 ~ 1+(1|codeLoc),data=dat,family=binomial())
MainModelNullNoRe <- glm(Insp1Pos>0 ~ 1,data=dat,family=binomial())

NagelkerkeFromModels <- function(FullModel,NullModel,NullModelNoRe){
    n<-dim(NullModelNoRe$data)[1]
    dev <- deviance(FullModel)
    nullDev <- deviance(NullModel)
    R2toNullRe <- (1 - exp((dev- nullDev)/n))/(1 - exp(-nullDev/n))
    cat("Nagelkerke R2 to Null with RE",R2toNullRe,"\n")

    nullNoReDev <- deviance(MainModelNullNoRe)
    R2toNullNoRe <- (1 - exp((dev- nullNoReDev)/n))/(1 - exp(-nullNoReDev/n))
    cat("Nagelkerke R2 to Null no RE",R2toNullNoRe,"\n")
}
NagelkerkeFromModels(MainModel,MainModelNull,MainModelNullNoRe)

# normal glm
A<-glm(Insp1Pos>0 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(",count(dat$Insp1Pos>0),"+) aic:",AIC(A),"\n")
#                                  OR    2.5 %   97.5 % signifs
# (Intercept)                4.35e-05 1.46e-05 1.30e-04     ***
# infFAstatus > 0TRUE        3.43e+00 2.24e+00 5.25e+00     ***
# rocFAstatus == 0TRUE       4.04e+01 9.74e+00 1.67e+02     ***
# ageAP                      1.77e+00 1.51e+00 2.07e+00     ***
# rocFAstatus == 0TRUE:ageAP 6.43e-01 5.04e-01 8.21e-01     ***
# N: 56491 ( 116 +) aic: 1570.09 

# using bayesian allows to avoid the warning (But might be hidden)
# by district
library(blme)
A<-bglmer(Insp1Pos>0 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(",count(dat$Insp1Pos>0),"+) aic:",AIC(A),"\n")
#                                  OR     2.5%    97.5% signifs
# (Intercept)                2.89e-05 3.96e-06 2.11e-04     ***
# infFAstatus > 0TRUE        3.56e+00 2.32e+00 5.47e+00     ***
# rocFAstatus == 0TRUE       5.00e+01 9.83e+00 2.55e+02     ***
# ageAP                      1.84e+00 1.34e+00 2.52e+00     ***
# rocFAstatus == 0TRUE:ageAP 6.18e-01 4.71e-01 8.12e-01     ***
# N: 56491 ( 116 +) aic: 1524.488 

# by locality
A<-bglmer(Insp1Pos>0 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(",count(dat$Insp1Pos>0),"+) aic:",AIC(A),"\n")
#                                  OR     2.5%    97.5% signifs
# (Intercept)                1.98e-05 2.86e-06 1.37e-04     ***
# infFAstatus > 0TRUE        3.20e+00 1.96e+00 5.24e+00     ***
# rocFAstatus == 0TRUE       2.09e+01 3.26e+00 1.34e+02      **
# ageAP                      1.54e+00 1.14e+00 2.09e+00      **
# rocFAstatus == 0TRUE:ageAP 7.31e-01 5.35e-01 1.00e+00       *
# N: 56491 ( 116 +) aic: 1415.547



# # difference in error "proportional reduction in error" by the model compared to the null model
# 1-(sum((dat$Insp1Pos-predict(Ablgmer_byDis,type="response"))^2)/sum((dat$Insp1Pos-mean(dat$Insp1Pos))^2))





# correct accounting for age would be:
# A<-glmer(Insp1Pos>0 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP)+(1|codeLoc)+offset(log(pmin(ageAP,4))),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(",count(dat$Insp1Pos>0),"+) aic:",AIC(A),"\n")
#                                 OR     2.5%    97.5% signifs
# (Intercept)                9.35e-06 1.27e-06 6.89e-05     ***
# infFAstatus > 0TRUE        3.19e+00 1.95e+00 5.22e+00     ***
# rocFAstatus == 0TRUE       2.65e+01 3.82e+00 1.85e+02     ***
# ageAP                      1.40e+00 1.02e+00 1.92e+00       *
# rocFAstatus == 0TRUE:ageAP 7.01e-01 5.05e-01 9.72e-01       *
# N: 56491 ( 116 +) aic: 1419.693 
#=> aic no better but conclusions robust

# or maybe:
#                                      OR     2.5%    97.5% signifs
# (Intercept)                     1.98e-05 0.000002 1.96e-04     ***
# infFAstatus > 0TRUE             3.20e+00 1.955014 5.23e+00     ***
# rocFAstatus == 0TRUE            1.43e+01 1.845750 1.11e+02       *
# log(ageAP)                      1.87e+00 0.506987 6.92e+00        
# rocFAstatus == 0TRUE:log(ageAP) 4.22e-01 0.128438 1.39e+00        
# N: 56491 ( 116 +) aic: 1424.851 
#=> really bad



# # on glm the following works fine to predict at 0 age:
# A<-glm(Insp1Pos>0 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(",count(dat$Insp1Pos>0),"+) aic:",AIC(A),"\n")
# datInit<-dat
# datInit$ageAP<-0
# 
# sum(predict(A,datInit,type="response")[datInit$Ispray & datInit$IIspray!=0])
# sum(predict(A,datInit,type="response")[datInit$rocFAstatus!=0])

# # on the glmer need to do it by hand inspired by:
# library(lme4)
# library(ggplot2) # Plotting
# data("Orthodont",package="MEMSS")
# fm1 <- lmer(
# 	        formula = distance ~ age*Sex + (age|Subject)
# 		    , data = Orthodont
# 		)
# newdat <- expand.grid(
# 		          age=c(8,10,12,14)
# 			      , Sex=c("Male","Female")
# 			      , distance = 0
# 			      )
# mm <- model.matrix(terms(fm1),newdat)
# newdat$distance <- mm %*% fixef(fm1)
# pvar1 <- diag(mm %*% tcrossprod(vcov(fm1),mm))
# tvar1 <- pvar1+VarCorr(fm1)$Subject[1]  ## must be adapted for more complex models
# newdat <- data.frame(
# 		         newdat
# 			     , plo = newdat$distance-2*sqrt(pvar1)
# 			     , phi = newdat$distance+2*sqrt(pvar1)
# 			         , tlo = newdat$distance-2*sqrt(tvar1)
# 			         , thi = newdat$distance+2*sqrt(tvar1)
# 				 )
# #plot confidence
# g0 <- ggplot(newdat, aes(x=age, y=distance, colour=Sex))+geom_point()
# g0 + geom_errorbar(aes(ymin = plo, ymax = phi))+
#     opts(title="CI based on fixed-effects uncertainty ONLY")
# #plot prediction
# g0 + geom_errorbar(aes(ymin = tlo, ymax = thi))+
#     opts(title="CI based on FE uncertainty + RE variance")

# A<-glmer(Insp1Pos>0 ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(",count(dat$Insp1Pos>0),"+) aic:",AIC(A),"\n")
#                                 OR    2.5%   97.5% signifs
# (Intercept)                2.2e-05 2.8e-06 1.8e-04     ***
# infFAstatus > 0TRUE        2.6e+00 1.6e-01 4.4e+01        
# ageAP                      1.5e+00 1.1e+00 2.1e+00       *
# rocFAstatus == 0TRUE       2.0e+01 2.6e+00 1.6e+02      **
# infFAstatus > 0TRUE:ageAP  1.0e+00 6.8e-01 1.6e+00        
# ageAP:rocFAstatus == 0TRUE 7.3e-01 5.2e-01 1.0e+00       .
# N: 56491 ( 116 +) aic: 1417.459


## locality random effect
# A<-glmer(nbDen1 ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
#                                OR    2.5%   97.5% signifs
# (Intercept)                0.00015 4.6e-05 0.00051     ***
# infFAstatus > 0TRUE        1.26171 2.4e-01 6.74610        
# ageAP                      1.48068 1.2e+00 1.80059     ***
# rocFAstatus == 0TRUE       2.43587 6.1e-01 9.72917        
# infFAstatus > 0TRUE:ageAP  1.08924 8.5e-01 1.40054        
# ageAP:rocFAstatus == 0TRUE 0.89874 7.1e-01 1.14011        
# N: 59279 aic: 2697.632 


# ### P(insp+ arround reports ~ infestation * not treated * age + (1|NRO_DEN)
# ### only households in treated city blocks
# dat<-commonVig[(commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1 & commonVig$id_manz %in% byManz$id_manz[byManz$Total != byManz$Unspray]),]
# A<-glmer(INSP_POSITIVA==1 ~ (infFAstatus>0)*(ageAP)+as.factor(rocFAstatus)*(ageAP)+(1|D)+(1|NRO_DENUNCIA),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(",count(dat$INSP_POSITIVA==1),"+) aic:",AIC(A),"\n")
# A<-glmer(INSP_POSITIVA==1 ~ (infFAstatus>0)+(rocFAstatus==0)+(1|NRO_DENUNCIA),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(",count(dat$INSP_POSITIVA==1),"+) aic:",AIC(A),"\n")

# probably controlling for the denuncia number make even more sens
# A<-glmer(INSP_POSITIVA==1 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP)+(1|NRO_DENUNCIA),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# similar but stronger impact of non-spray

dat<-byManz
# A<-glmer(Den1~(nbPos>0)*ageAP+log(Unspray+1)*ageAP+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# best but asks for more explanations
# A<-glmer(Den1~(nbPos>0)*log(Unspray+1)*ageAP+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# good but still the log to explain
# A<-glmer(Den1~(nbPos>0)*(Unspray)*ageAP+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# not good because everything gets significance
# much less good but ok 
# A<-glmer(Den1~log(nbPos+1)*log(Unspray+1)*ageAP+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# good and little explanations, 1|codeLoc loose everything

# A<-glmer(Den1~(nbPos>0)*log(Unspray+1)*(ageAP)+(1|codeLoc)+Total,data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")

# A<-glmer(Den1~(nbPos)*(Unspray)*ageAP+Total+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# really not good



# current selection
cat("-----------------
    confirmed reports among all city blocks
-----------------\n")
dat<-byManz
# A<-glmer(Den1 ~log(nbPos+1)*ageAP+log(Unspray+1)*ageAP+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$Den1),") aic:",AIC(A),"\n")
#                              OR     2.5%    97.5% signifs
# (Intercept)            0.000232 2.87e-05  0.00187     ***
# log(nbPos + 1)         1.790870 7.82e-01  4.10166        
# ageAP                  1.846637 1.34e+00  2.53555     ***
# log(Unspray + 1)       5.601234 2.34e+00 13.38435     ***
# log(nbPos + 1):ageAP   0.999767 8.74e-01  1.14401        
# ageAP:log(Unspray + 1) 0.777630 6.70e-01  0.90254     ***
# N: 3727 (+ 164 ) aic: 1050.025 

# A<-glmer(Den1 & nbInspPos>0 ~log(nbPos+1)+ageAP+log(Unspray+1)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$Den1 & dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                        OR     2.5%   97.5% signifs
# (Intercept)      0.000354 4.67e-05 0.00268     ***
# log(nbPos + 1)   2.093098 1.42e+00 3.09108     ***
# ageAP            1.286001 9.62e-01 1.71951       .
# log(Unspray + 1) 1.689805 1.05e+00 2.71440       *
# N: 3727 (+ 80 ) aic: 585.5192 

dat<-byManz
# A<-glmer(Den1 & nbInspPos>0 ~(nbPos>0)+log(Unspray+1)*ageAP+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$Den1 & dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                              OR     2.5%   97.5% signifs
# (Intercept)            9.04e-05 5.83e-06  0.0014     ***
# nbPos > 0TRUE          3.36e+00 1.62e+00  6.9276      **
# log(Unspray + 1)       4.81e+00 1.42e+00 16.3465       *
# ageAP                  1.69e+00 1.13e+00  2.5195       *
# log(Unspray + 1):ageAP 8.16e-01 6.62e-01  1.0049       .
# N: 3727 (+ 80 ) aic: 585.5598 

dat<-byManz
write.csv(byManz,"byManz_fullEID.csv")
A<-glmer(Den1 & nbInspPos>0 ~log(nbPos+1)+log(Unspray+1)*ageAP+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$Den1 & dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                              OR     2.5%    97.5% signifs
# (Intercept)            8.22e-05 4.54e-06  0.00149     ***
# log(nbPos + 1)         2.08e+00 1.42e+00  3.04420     ***
# log(Unspray + 1)       4.85e+00 1.34e+00 17.48516       *
# ageAP                  1.69e+00 1.10e+00  2.58137       *
# log(Unspray + 1):ageAP 8.21e-01 6.61e-01  1.02070       .
# N: 3727 (+ 80 ) aic: 583.8216 
n<-dim(dat)[1]
FullModel <- A
byManzModel <-A
byManzModelNull <- glmer(Den1 & nbInspPos>0 ~ 1+(1|codeLoc),data=dat,family=binomial())
byManzModelNullNoRe <- glm(Den1 & nbInspPos>0 ~ 1,data=dat,family=binomial())
NagelkerkeFromModels(byManzModel,byManzModelNull,byManzModelNullNoRe)

# test that using fit adapted to small number makes same results  than normal glm
# library("logistf")
# ltf <- logistf(Insp1Pos/ones ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP),data=dat,family=binomial())


# correct accounting for time allowed for reports would be:
# A<-glmer(Den1 & nbInspPos>0 ~log(nbPos+1)+log(Unspray+1)*ageAP+(1|codeLoc)+offset(log(pmin(ageAP,4))),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$Den1 & dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                              OR     2.5%    97.5% signifs
# (Intercept)            2.99e-05 1.42e-06 6.31e-04     ***
# log(nbPos + 1)         2.09e+00 1.42e+00 3.09e+00     ***
# log(Unspray + 1)       5.96e+00 1.57e+00 2.27e+01      **
# ageAP                  1.56e+00 9.93e-01 2.45e+00       .
# log(Unspray + 1):ageAP 7.91e-01 6.30e-01 9.93e-01       *
# N: 3727 (+ 80 ) aic: 588.5868 
#=> but not a better AIC


# dat<-byManz
# A<-glmer(Den1 & nbInspPos>0 ~log(nbPos+1)*ageAP+log(Unspray+1)*ageAP+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$Den1 & dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                              OR     2.5%    97.5% signifs
# (Intercept)            0.000164 7.03e-06  0.00381     ***
# log(nbPos + 1)         1.077921 3.34e-01  3.48147        
# ageAP                  1.478850 9.14e-01  2.39340        
# log(Unspray + 1)       4.501639 1.22e+00 16.66419       *
# log(nbPos + 1):ageAP   1.124597 9.27e-01  1.36503        
# ageAP:log(Unspray + 1) 0.832307 6.67e-01  1.03838        
# N: 3727 (+ 80 ) aic: 584.0742 

# A<-glmer(Den1 & nbInspPos>0 ~log(nbPos+1)*ageAP+log(Unspray+1)*ageAP+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$Den1 & dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                              OR     2.5%    97.5% signifs
# (Intercept)            0.000353 2.69e-05  0.00462     ***
# log(nbPos + 1)         1.382500 5.56e-01  3.43781        
# ageAP                  1.635513 1.09e+00  2.45114       *
# log(Unspray + 1)       5.404811 2.17e+00 13.45349     ***
# log(nbPos + 1):ageAP   1.043871 9.01e-01  1.20879        
# ageAP:log(Unspray + 1) 0.784354 6.70e-01  0.91798      **
# N: 3727 (+ 80 ) aic: 610.2778 

cat("-----------------
    Infestation among households inspected
-----------------\n")
dat<-commonVig[(commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1 & commonVig$id_manz %in% byManz$id_manz[byManz$nbSpray>0]),]
dat$FA<-as.numeric(dat$infFAstatus>0)
dat$FA[dat$FA==1]<-"Infested"
dat$FA[dat$FA=="0"]<-"Non-Infested"
dat$FA[dat$rocFAstatus==0]<-"Unsprayed"
dat$FA<-factor(dat$FA,levels=c("Non-Infested","Infested","Unsprayed"))
# A<-glmer(INSP_POSITIVA==1 ~ FA+(ageAP)+(1|NRO_DENUNCIA),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$INSP_POSITIVA==1),") aic:",AIC(A),"\n")
#                 OR   2.5% 97.5% signifs
# (Intercept) 0.0825 0.0401  0.17     ***
# FAInfested  1.9819 1.1993  3.28      **
# FAUnsprayed 5.1177 2.8131  9.31     ***
# ageAP       1.0408 0.9172  1.18        
# N: 775 (+ 130 ) aic: 650.2869 

# A<-glmer(INSP_POSITIVA==1 ~ (infFAstatus>0)+(rocFAstatus==0)+(ageAP)+(1|NRO_DENUNCIA),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$INSP_POSITIVA==1),") aic:",AIC(A),"\n")
#                          OR   2.5% 97.5% signifs
# (Intercept)          0.0825 0.0401  0.17     ***
# infFAstatus > 0TRUE  1.9819 1.1993  3.28      **
# rocFAstatus == 0TRUE 5.1177 2.8131  9.31     ***
# ageAP                1.0408 0.9172  1.18        
# N: 775 (+ 130 ) aic: 650.2869 

# dat<-byHouse[byHouse$Insp1>0,]
# A<-glmer(Insp1Pos>0 ~ (infFAstatus>0)+(rocFAstatus==0)+(ageAP|L)+(1|id_manz),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$Insp1Pos==1),") aic:",AIC(A),"\n")
# > A<-glmer(INSP_POSITIVA==1 ~ (infFAstatus>0)+(rocFAstatus==0)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$INSP_POSITIVA==1),") aic:",AIC(A),"\n")
#                         OR   2.5% 97.5% signifs
# (Intercept)          0.108 0.0706 0.167     ***
# infFAstatus > 0TRUE  1.950 1.1654 3.264       *
# rocFAstatus == 0TRUE 3.796 2.1918 6.574     ***
# N: 636 (+ 102 ) aic: 543.5821 

# A<-glmer(INSP_POSITIVA==1 ~ (infFAstatus>0)+(rocFAstatus==0)+(1|NRO_DENUNCIA),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$INSP_POSITIVA==1),") aic:",AIC(A),"\n")                         OR   2.5% 97.5% signifs
# (Intercept)          0.0963 0.0654 0.142     ***
# infFAstatus > 0TRUE  1.9454 1.1018 3.435       *
# rocFAstatus == 0TRUE 3.8646 2.1597 6.916     ***
# N: 636 (+ 102 ) aic: 543.2284

# > A<-glmer(INSP_POSITIVA==1 ~ (infFAstatus>0)+(rocFAstatus==0)+(ageAP|L),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$INSP_POSITIVA==1),") aic:",AIC(A),"\n")                        OR   2.5%  97.5% signifs
# (Intercept)          0.082 0.0538  0.125     ***
# infFAstatus > 0TRUE  2.241 1.2972  3.873      **
# rocFAstatus == 0TRUE 5.778 3.1610 10.561     ***
# N: 636 (+ 102 ) aic: 516.1753 

write.csv(commonVig,"commonVig_fullEID.csv")
dat<-commonVig[(commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1),]
# in addition keep only the first time inspected
dat<-dat[order(dat$ANIO_INSP),]
dat<-dat[match(unique(dat$UNICODE),dat$UNICODE),]
byHouse$Insp1posfirst<- (byHouse$UNICODE %in% dat[dat$INSP_POSITIVA==1,"UNICODE"])

A<-glmer(INSP_POSITIVA==1 ~ (infFAstatus>0)+(rocFAstatus==0)+(ageAP|L)+(1|NRO_DENUNCIA),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$INSP_POSITIVA==1),") aic:",AIC(A),"\n")
#                          OR   2.5%  97.5% signifs
# (Intercept)          0.0789 0.0512  0.122     ***
# infFAstatus > 0TRUE  2.2362 1.2802  3.906      **
# rocFAstatus == 0TRUE 5.8833 3.1801 10.884     ***
# N: 636 (+ 102 ) aic: 518.0342 
byHouseInInsp<-A
byHouseInInspNull<- glmer(INSP_POSITIVA==1 ~ 1+(1|NRO_DENUNCIA),data=dat,family=binomial())
byHouseInInspNullNoRE <- glm(INSP_POSITIVA==1 ~ 1,data=dat,family=binomial())

NagelkerkeFromModels(byHouseInInsp,byHouseInInspNull,byHouseInInspNullNoRE)

# dat<-commonVig[(commonVig$INSP_COMPLETA==1 & commonVig$TYPO_DENUNCIA==1 & commonVig$id_manz %in% byManz$id_manz[byManz$nbSpray>0]),]
# A<-glmer(INSP_POSITIVA==1 ~ (infFAstatus>0)*ageAP+(rocFAstatus==0)*(ageAP)+(1|NRO_DENUNCIA),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$INSP_POSITIVA==1),") aic:",AIC(A),"\n")
#                                 OR   2.5% 97.5% signifs
# (Intercept)                 0.0548 0.0188  0.16     ***
# infFAstatus > 0TRUE         2.8941 0.6381 13.13        
# ageAP                       1.1246 0.9256  1.37        
# rocFAstatus == 0TRUE       13.3781 3.3554 53.34     ***
# infFAstatus > 0TRUE:ageAP   0.9283 0.7023  1.23        
# ageAP:rocFAstatus == 0TRUE  0.7856 0.5782  1.07        
# N: 775 (+ 130 ) aic: 651.8858 



#=====================================
## reports in all households 
## Nota: non-participating may not be very likely to report
#=====================================
dat<-byHouse
#---------------------------------------
# reports (may not be exact but very close to exact)
#---------------------------------------

# A<-glmer(nbDen1 ~ (infFAstatus>0)+(ageAP)+(rocFAstatus==0)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"aic:",AIC(A),"\n")
# Warning message:
# In mer_finalize(ans) : false convergence (8)
#                           OR    2.5%   97.5% signifs
# (Intercept)          0.00025 7.3e-05 0.00083     ***
# infFAstatus > 0TRUE  2.51753 1.9e+00 3.38822     ***
# ageAP                1.55417 1.3e+00 1.88691     ***
# rocFAstatus == 0TRUE 1.40041 9.5e-01 2.05746       .
# N: 59279 aic: 2829.81 

# A<-glmer(nbDen1 ~ (infFAstatus>0)+(ageAP)+(rocFAstatus==0)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1),") aic:",AIC(A),"\n")
#                           OR    2.5%   97.5% signifs
# (Intercept)          0.00016 5.6e-05 0.00047     ***
# infFAstatus > 0TRUE  2.36536 1.7e+00 3.29082     ***
# ageAP                1.46761 1.2e+00 1.74348     ***
# rocFAstatus == 0TRUE 1.38879 9.0e-01 2.15147        
# N: 59279 (+ 234 ) aic: 2642.73 

# A<-glmer(nbDen1 ~ (infFAstatus>0)+(ageAP)*(rocFAstatus==0)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1),") aic:",AIC(A),"\n")
#                                 OR    2.5%   97.5% signifs
# (Intercept)                0.00014 4.2e-05 4.3e-04     ***
# infFAstatus > 0TRUE        2.19093 1.6e+00 3.0e+00     ***
# ageAP                      1.50932 1.3e+00 1.8e+00     ***
# rocFAstatus == 0TRUE       2.75168 7.2e-01 1.1e+01        
# ageAP:rocFAstatus == 0TRUE 0.88191 7.0e-01 1.1e+00        
# N: 59279 (+ 243 ) aic: 2696.146 

# A<-glmer(nbDen1 ~ (infFAstatus>0)+(ageAP)*(rocFAstatus==0)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1),") aic:",AIC(A),"\n")
#                                 OR    2.5%   97.5% signifs
# (Intercept)                0.00015 4.1e-05 5.4e-04     ***
# infFAstatus > 0TRUE        2.33110 1.7e+00 3.1e+00     ***
# ageAP                      1.71837 1.4e+00 2.1e+00     ***
# rocFAstatus == 0TRUE       5.34430 1.7e+00 1.7e+01      **
# ageAP:rocFAstatus == 0TRUE 0.79429 6.5e-01 9.7e-01       *
# N: 59279 (+ 243 ) aic: 2916.21 

# A<-glmer(nbDen1 ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1),") aic:",AIC(A),"\n")
#                                OR    2.5%   97.5% signifs
# (Intercept)                0.00015 4.6e-05 0.00051     ***
# infFAstatus > 0TRUE        1.26171 2.4e-01 6.74610        
# ageAP                      1.48068 1.2e+00 1.80059     ***
# rocFAstatus == 0TRUE       2.43587 6.1e-01 9.72917        
# infFAstatus > 0TRUE:ageAP  1.08924 8.5e-01 1.40054        
# ageAP:rocFAstatus == 0TRUE 0.89874 7.1e-01 1.14011        
# N: 59279 aic: 2697.632 

cat("---------------------------------------
	Reports confirmed (same household)
---------------------------------------\n")
# A<-glmer(nbDen1 & nbInspPos>0 ~ (infFAstatus>0)+(rocFAstatus==0)+(ageAP)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1& dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                           OR    2.5%   97.5% signifs
# (Intercept)          4.1e-05 6.8e-06 0.00025     ***
# infFAstatus > 0TRUE  4.0e+00 2.4e+00 6.77523     ***
# rocFAstatus == 0TRUE 2.6e+00 1.4e+00 4.83641      **
# ageAP                1.6e+00 1.2e+00 2.19253     ***
# N: 59279 (+ 74 ) aic: 1051.978 

# A<-glmer(nbDen1 & nbInspPos>0 ~ (infFAstatus>0)+(rocFAstatus==0)+(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1& dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                           OR    2.5%   97.5% signifs
# (Intercept)          1.8e-05 2.8e-06 0.00012     ***
# infFAstatus > 0TRUE  3.9e+00 2.1e+00 7.18984     ***
# rocFAstatus == 0TRUE 2.9e+00 1.3e+00 6.24420      **
# ageAP                1.5e+00 1.1e+00 2.02948      **
# N: 59279 (+ 74 ) aic: 995.5366 

# A<-glmer(nbDen1 & nbInspPos>0 ~ (infFAstatus>0)+(rocFAstatus==0)+(ageAP)+(1|id_manz),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1& dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                           OR    2.5%  97.5% signifs
# (Intercept)          2.2e-07 2.4e-10 0.0002     ***
# infFAstatus > 0TRUE  3.2e+00 1.6e+00 6.5356     ***
# rocFAstatus == 0TRUE 3.3e+00 1.4e+00 7.7012      **
# ageAP                1.6e+00 5.5e-01 4.4550        
# N: 59279 (+ 74 ) aic: 841.5942 

# A<-glmer(nbDen1 & nbInspPos>0 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1& dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                                 OR    2.5%   97.5% signifs
# (Intercept)                1.6e-05 1.9e-06 1.3e-04     ***
# infFAstatus > 0TRUE        4.0e+00 2.4e+00 6.7e+00     ***
# rocFAstatus == 0TRUE       2.1e+01 2.6e+00 1.7e+02      **
# ageAP                      1.9e+00 1.4e+00 2.7e+00     ***
# rocFAstatus == 0TRUE:ageAP 7.0e-01 4.9e-01 9.9e-01       *
# N: 59279 (+ 74 ) aic: 1049.732 

# A<-glmer(nbDen1 & nbInspPos>0 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1& dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                                 OR    2.5%   97.5% signifs
# (Intercept)                1.2e-05 1.3e-06 1.1e-04     ***
# infFAstatus > 0TRUE        3.9e+00 2.1e+00 7.2e+00     ***
# rocFAstatus == 0TRUE       9.7e+00 7.3e-01 1.3e+02       .
# ageAP                      1.6e+00 1.2e+00 2.3e+00      **
# rocFAstatus == 0TRUE:ageAP 8.1e-01 5.3e-01 1.2e+00        
# N: 59279 (+ 74 ) aic: 996.1918 

# A<-glmer(nbDen1 & nbInspPos>0 ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1& dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                                 OR    2.5%   97.5% signifs
# (Intercept)                2.0e-05 1.9e-06 2.1e-04     ***
# infFAstatus > 0TRUE        2.4e+00 1.4e-01 4.2e+01        
# ageAP                      1.9e+00 1.3e+00 2.7e+00      **
# rocFAstatus == 0TRUE       1.7e+01 1.6e+00 1.8e+02       *
# infFAstatus > 0TRUE:ageAP  1.1e+00 7.1e-01 1.6e+00        
# ageAP:rocFAstatus == 0TRUE 7.2e-01 4.9e-01 1.1e+00       .
# N: 59279 (+ 74 ) aic: 1051.608 

# A<-glmer(nbDen1 & nbInspPos>0 ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$nbDen1& dat$nbInspPos>0),") aic:",AIC(A),"\n")
#                                 OR    2.5%   97.5% signifs
# (Intercept)                1.6e-05 1.4e-06 1.8e-04     ***
# infFAstatus > 0TRUE        1.6e+00 4.8e-02 5.3e+01        
# ageAP                      1.6e+00 1.1e+00 2.3e+00       *
# rocFAstatus == 0TRUE       7.3e+00 4.6e-01 1.1e+02        
# infFAstatus > 0TRUE:ageAP  1.1e+00 6.9e-01 1.9e+00        
# ageAP:rocFAstatus == 0TRUE 8.4e-01 5.4e-01 1.3e+00        
# N: 59279 (+ 74 ) aic: 997.8161 

cat("---------------------------------------
	Household with community report leading to some infestation
---------------------------------------\n")
# A<-glmer(Pos1 ~ (infFAstatus>0)+(rocFAstatus==0)+(ageAP)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$Pos1),") aic:",AIC(A),"\n")
#                           OR    2.5%   97.5% signifs
# (Intercept)          0.00016 3.4e-05 0.00077     ***
# infFAstatus > 0TRUE  4.16622 2.6e+00 6.66030     ***
# rocFAstatus == 0TRUE 2.35538 1.3e+00 4.17352      **
# ageAP                1.33789 1.0e+00 1.71922       *
# N: 59279 (+ 90 ) aic: 1260.179 

# A<-glmer(Pos1 ~ (infFAstatus>0)+(rocFAstatus==0)+(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$Pos1),") aic:",AIC(A),"\n")
#                           OR    2.5%   97.5% signifs
# (Intercept)          0.00006 1.5e-05 0.00024     ***
# infFAstatus > 0TRUE  4.19600 2.4e+00 7.22249     ***
# rocFAstatus == 0TRUE 2.48203 1.3e+00 4.92541      **
# ageAP                1.36970 1.1e+00 1.70946      **
# N: 59279 (+ 90 ) aic: 1203.814 

# A<-glmer(Pos1 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$Pos1),") aic:",AIC(A),"\n")
#                                 OR    2.5%   97.5% signifs
# (Intercept)                5.1e-05 1.1e-05 2.4e-04     ***
# infFAstatus > 0TRUE        4.2e+00 2.4e+00 7.2e+00     ***
# rocFAstatus == 0TRUE       4.1e+00 6.2e-01 2.7e+01        
# ageAP                      1.4e+00 1.1e+00 1.8e+00      **
# rocFAstatus == 0TRUE:ageAP 9.1e-01 6.6e-01 1.3e+00        
# N: 59279 (+ 90 ) aic: 1205.417 

dat<-byHouse
A<-glmer(Pos1 ~ (infFAstatus>0)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$Pos1),") aic:",AIC(A),"\n")
#                                 OR     2.5%    97.5% signifs
# (Intercept)                0.00005 9.87e-06 2.53e-04     ***
# infFAstatus > 0TRUE        4.06656 2.34e+00 7.06e+00     ***
# rocFAstatus == 0TRUE       4.24006 6.26e-01 2.87e+01        
# ageAP                      1.38947 1.07e+00 1.80e+00       *
# rocFAstatus == 0TRUE:ageAP 0.91505 6.59e-01 1.27e+00        
# N: 56550 (+ 91 ) aic: 1159.98 


# A<-glmer(Pos1 ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|D),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$Pos1),") aic:",AIC(A),"\n")
#                                 OR    2.5%  97.5% signifs
# (Intercept)                0.00018 3.1e-05  0.001     ***
# infFAstatus > 0TRUE        0.52128 4.7e-02  5.834        
# ageAP                      1.33160 1.0e+00  1.767       *
# rocFAstatus == 0TRUE       3.89542 7.1e-01 21.482        
# infFAstatus > 0TRUE:ageAP  1.37397 9.6e-01  1.966       .
# ageAP:rocFAstatus == 0TRUE 0.90036 6.7e-01  1.212        
# N: 59279 (+ 90 ) aic: 1260.179 

# A<-glmer(Pos1 ~ (infFAstatus>0)*(ageAP)+(rocFAstatus==0)*(ageAP)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=2);cat("N:",dim(dat)[1],"(+",count(dat$Pos1),") aic:",AIC(A),"\n")
#                                 OR    2.5%   97.5% signifs
# (Intercept)                9.5e-05 1.9e-05 4.8e-04     ***
# infFAstatus > 0TRUE        4.3e-01 2.5e-02 7.5e+00        
# ageAP                      1.3e+00 9.8e-01 1.7e+00       .
# rocFAstatus == 0TRUE       2.2e+00 3.3e-01 1.5e+01        
# infFAstatus > 0TRUE:ageAP  1.4e+00 9.3e-01 2.2e+00        
# ageAP:rocFAstatus == 0TRUE 1.0e+00 7.2e-01 1.4e+00        
# N: 59279 (+ 90 ) aic: 1260.179 

#-------------------------
# accounting separatly for infestation in second phase
# => do not bother
#-------------------------

dat<-byHouse
# A<-glm(Insp1Pos>0 ~ as.factor(2-rocFAstatus)+Iinfested+IIinfested+(ageAP)y=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(",count(d"+) aic:",AIC(A),"\n")
#     			         OR    2.5 %   97.5 % signifs
# (Intercept)                 0.000117 5.19e-05 2.63e-04     ***
# as.factor(2 - rocFAstatus)1 1.329207 7.70e-01 2.30e+00        
# as.factor(2 - rocFAstatus)2 3.634544 2.22e+00 5.96e+00     ***
# IinfestedTRUE               3.247986 2.08e+00 5.08e+00     ***
# IIinfestedTRUE              4.891326 1.51e+00 1.58e+01      **
# ageAP                       1.507499 1.34e+00 1.69e+00     ***
# N: 56491 ( 116 +) aic: 1584.539 


# ------------------------------
# likelihood to be inspected
# ------------------------------
# A<-glmer(nbInsp>0 ~ (rocFAstatus==0)*ageAP+(infFAstatus>0)+(1|codeLoc),data=dat,family=binomial());print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(",count(dat$nbInsp>0),"+) aic:",AIC(A),"\n")
#                                  OR     2.5%   97.5% signifs
# (Intercept)                0.000474 0.000161 0.00139     ***
# rocFAstatus == 0TRUE       2.283233 1.235992 4.21779      **
# ageAP                      1.400197 1.167098 1.67985     ***
# infFAstatus > 0TRUE        1.799014 1.488900 2.17372     ***
# rocFAstatus == 0TRUE:ageAP 0.867172 0.771809 0.97432       *
# N: 56491 ( 116 +) aic: 6891.04 

cat("---------------------------------------
	what about spray in surveillance?
---------------------------------------\n")
write.csv(rocvig,"rocvig_fullEID.csv")
dat<-rocvig[rocvig$complete==1,]
A<-lmer(infestedRV~ (rocFAstatus==0) + (infFAstatus>0)+(1|codeLoc),data=dat);
print(summaryOR(A),digits=3);cat("N:",dim(dat)[1],"(+",count(dat$infestedRV==1),") aic:",AIC(A),"\n")

# -------------------------------
# Comparison with estimated residual OR
# -------------------------------
fisher.test(with(byHouse,table(nbInspPos>0,rocFAstatus>0)))
fisher.test(with(byHouse[byHouse$nbInsp>0&byHouse$Ispray & byHouse$Iinfested>0,],table(nbInspPos>0,IIspray)))
fisher.test(with(byHouse,table(nbInspPos>0|nbRocPos>0,rocFAstatus==0)))

# =====================================
# Correlation participation aux deux traitements
# =====================================
dat<-byHouse
A<-glmer(IIspray ~ Ispray + (1|D),data=dat);summaryOR(A)
A<-glmer(IIspray ~ Ispray + (1|codeLoc),data=dat);summaryOR(A)

# =====================================
# correlation of initial infestation and openning
# =====================================
byHouse$firstObs<-NA
byHouse$firstObs[which(byHouse$IIspray)]<-byHouse$IIinfested[which(byHouse$IIspray)]
byHouse$firstObs[which(byHouse$Ispray)]<-byHouse$Iinfested[which(byHouse$Ispray)]
dat<-byHouse[which(!is.na(byHouse$firstObs)),]

A<-glmer(firstObs ~ Ispray + (1|codeLoc),data=dat);summaryOR(A)
A<-glm(firstObs ~ rocFAstatus==1 ,data=dat);summaryOR(A)
A<-glmer(firstObs ~ (rocFAstatus==1) + (1|codeLoc),data=dat);summaryOR(A)
A<-glm(rocFAstatus==2 ~ firstObs ,data=dat);summaryOR(A)
A<-glmer(rocFAstatus==2 ~ firstObs + (1|codeLoc),data=dat);summaryOR(A)
A<-glmer(firstObs ~ (rocFAstatus==2) + (1|codeLoc),data=dat);summaryOR(A)

dat<-byHouse[which( byHouse$IIspray) ,]
A<-glm(firstObs ~ (rocFAstatus==2) ,data=dat,family=binomial);summaryOR(A)
A<-glmer(firstObs ~ (rocFAstatus==2) + (1|codeLoc) ,data=dat,family=binomial);summaryOR(A)

dat<-byHouse[which( byHouse$Ispray) ,]
A<-glm(firstObs ~ rocFAstatus==2 ,data=dat,family=binomial);summaryOR(A)
A<-glmer(firstObs ~ (rocFAstatus==2) +(1|codeLoc) ,data=dat,family=binomial);summaryOR(A)

dat<-byHouse[which(!is.na(byHouse$firstObs)&(!byHouse$Ispray | !byHouse$IIspray)) ,]
A<-glm(firstObs ~ Ispray==1 ,data=dat,family=binomial);summaryOR(A)
A<-glmer(firstObs ~ (Ispray==1) +(1|codeLoc) ,data=dat,family=binomial);summaryOR(A)

save.image(file="glmVigilance.img")

