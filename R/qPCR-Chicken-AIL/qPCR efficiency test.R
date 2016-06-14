# qPCR efficiency test for 5 pairs of primers 
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Shijie Lyu
# last modified June, 2016
# first written June, 2016

setwd("D:/Chicken/Rdata/qPCR") 
ctMean <- read.table("RawData/qPCR-efficieny-test-3.txt", header=TRUE,sep="\t",na.strings = "",check.names = FALSE)

#ctMean <- ctMean[-c(seq(1,30,6)),]  # Remove the 0.01ng loading amount


Genes <- c("LCORL","NCAPG","HPRT","HMBS","ACTB")  
Housekeepers <- c("HPRT","HMBS","ACTB")

#CountInAmount <- function(gene){
 # CountInAmount <- length(which(diff(ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"])<0))+1
 # return(CountInAmount)
#}

slope <- function(gene){
  CountInAmount <- length(which(diff(ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"])<0))+1
  ct <- ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"][1:CountInAmount]
  amount <- ctMean[which(ctMean[,"Target_Name"]==gene),"Amount"][1:CountInAmount]
  return(slope = round(as.numeric(lm(ct~log10(amount))[[1]][2]),3))
}

par(mfrow=c(2,3))                                                                                                                          # plot the original ctmean for every amount cDNA
for(gene in Genes){
  plot(ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"],t="l",xlab="cDNA amount",ylab="CT value",xaxt="n",main=gene)
  axis(1, at=seq(1,6),round(log10(ctMean[which(ctMean[,"Target_Name"]==gene),"Amount"]),2))
  points(ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"])
  CountInAmount <- length(which(diff(ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"])<0))+1
  points(ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"][1:CountInAmount],col="red",pch=19)
  legend("topright", legend = c(paste0("Slope = ",slope(gene)),paste0("Efficiency =",round((10^(-1/slope(gene))-1)*100)),"%"))
}


DeltaCtSlope <- function(gene,housekeeper){   
  CountInAmount <- min(length(which(diff(ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"])<0))+1,length(which(diff(ctMean[which(ctMean[,"Target_Name"]==housekeeper),"CT_mean"])<0))+1)
  diffCt <- ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"]-ctMean[which(ctMean[,"Target_Name"]==housekeeper),"CT_mean"]
  diffCt <- diffCt[which(!is.na(ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"]))][1:CountInAmount]
  amount <- ctMean[which(ctMean[,"Target_Name"]==gene),"Amount"][which(!is.na(ctMean[which(ctMean[,"Target_Name"]==gene),"CT_mean"]))][1:CountInAmount]
  return(DeltaCtSlope = round(as.numeric(lm(diffCt~log10(amount))[[1]][2]),3))
}


for(gene in Genes){
   for (housekeeper in Housekeepers){
     cat(paste("Delta Ct Slope",gene,"vs",housekeeper,"is",DeltaCtSlope(gene,housekeeper)),"\n")
   }
}

