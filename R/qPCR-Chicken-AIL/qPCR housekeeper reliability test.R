# qPCR housekeeper reliability test
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Shijie Lyu
# last modified June, 2016
# first written June, 2016

setwd("D:/Chicken/Rdata/qPCR") 
DataAll <- read.table("RawData/qPCR-housekeeper-test.txt", header=TRUE,sep="\t",na.strings = "",check.names = FALSE)

housekeepers <- c("HPRT","HMBS","ACTB")
NHIgroup <- c(432,430,437,427,426)
WL77group <- c(102,88,78,82,90,94)

par(mfrow=c(3,2))
for (HP in housekeepers){
  TempHP <- DataAll[which(DataAll[,"Target_Name"]==HP),]
  TempHPNHI <- TempHP[which(TempHP[,"Sample_Name"] %in% NHIgroup),]
  TempHPWL77 <- TempHP[which(TempHP[,"Sample_Name"] %in% WL77group),]
  plot(c(0,8),c(0,36),main=paste("Ct-mean distribution of",HP),xlab="Individuls",ylab="Ct mean",t="n")
  points(TempHPNHI[,"Ct_Mean"],col="blue",pch=19)
  points(TempHPWL77[,"Ct_Mean"],col="red",pch=19)
  legend("topright",c("NHI","WL77"),pch=c(19,19),col=c("blue","red"))
  boxplot(TempHPNHI[,"Ct_Mean"],TempHPWL77[,"Ct_Mean"],main=paste("Ct-mean comparison of",HP),ylab="Ct mean",xlab="Groups")
  axis(1,c("NHI","WL77"),at=1:2)
  cat(HP,"all chickens' variance is",var(TempHP[,"Ct_Mean"]),"\n")
  cat(paste0(HP,"-NHI"),"mean is",mean(TempHPNHI[,"Ct_Mean"]),"and variance is",var(TempHPNHI[,"Ct_Mean"]),"\n")
  cat(paste0(HP,"-WL77"),"mean is",mean(TempHPWL77[,"Ct_Mean"]),"and variance is",var(TempHPWL77[,"Ct_Mean"]),"\n")
  cat(paste0(HP,"-NHI"),"vs",paste0(HP,"-WL77"),"ttest pvalue is",t.test(TempHPNHI[,"Ct_Mean"],TempHPWL77[,"Ct_Mean"])$p.value,"\n")
}




HPRT <- DataAll[which(DataAll[,"Target_Name"]=="HPRT"),]

HPRTNHI <- HPRT[which(HPRT[,"Sample_Name"] %in% NHIgroup),]
HPRTWL77 <- HPRT[which(HPRT[,"Sample_Name"] %in% WL77group),]

#variance among individuls in one group
mean(HPRTNHI[,"Ct_Mean"])
var(HPRTNHI[,"Ct_Mean"])
mean(HPRTWL77[,"Ct_Mean"])
var(HPRTWL77[,"Ct_Mean"])

#differences between two groups
mean(HPRTNHI[,"Ct_Mean"]) - mean(HPRTWL77[,"Ct_Mean"])
t.test(HPRTNHI[,"Ct_Mean"],HPRTWL77[,"Ct_Mean"])