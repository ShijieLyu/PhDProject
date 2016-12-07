# Expression data analysis for 3 tissues at 13 weeks
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Shijie Lyu
# last modified August, 2016
# first written August, 2016

setwd("D:/Chicken/Rdata/qPCR") 
ExpData <- read.table("RawData/Expression of 13 week chicken.txt", header=TRUE,sep="\t",check.names = FALSE)

Tissues <- as.character(unique(ExpData[,"Tissue"]))
Breed <- as.character(unique(ExpData[,"Breed"]))
Genes <- c("LCORL","NCAPG")

ExpData <- cbind(ExpData, LCORL_RelativeExpression = 2^-(ExpData[,"LCORL_Ct"] - ExpData[,"GeoMean_Hosekeeper"]), 
                          NCAPG_RelativeExpression = 2^-(ExpData[,"NCAPG_Ct"] - ExpData[,"GeoMean_Hosekeeper"]))  ### Calculate the relative expression using 2^-Delat Ct

ExpAnalysis <- vector("list", length(Tissues)) 
names(ExpAnalysis) <- Tissues
for (Tissue in Tissues){
  TissueMatrix <- matrix(NA,2,6)
  rownames(TissueMatrix) <- Breed
  colnames(TissueMatrix) <- c("LCORL_FC_WL77:NHI", "LCORL_SD", "LCORL_Pvalue","NCAPG_FC_WL77:NHI", "NCAPG_SD","NCAPG_Pvalue")
  TemTissue <- ExpData[which(ExpData[,"Tissue"] == Tissue),]
  BM_NHI_LCORL <- TemTissue[which(TemTissue[,"Breed"]=="NHI"),][,"LCORL_RelativeExpression"]
  BM_WL77_LCORL <- TemTissue[which(TemTissue[,"Breed"]=="WL77"),][,"LCORL_RelativeExpression"]
  BM_NHI_NCAPG <- TemTissue[which(TemTissue[,"Breed"]=="NHI"),][,"NCAPG_RelativeExpression"]
  BM_WL77_NCAPG <- TemTissue[which(TemTissue[,"Breed"]=="WL77"),][,"NCAPG_RelativeExpression"]
  TissueMatrix["NHI","LCORL_FC_WL77:NHI"] <- 1
  TissueMatrix["WL77","LCORL_FC_WL77:NHI"] <- round(mean(BM_WL77_LCORL)/mean(BM_NHI_LCORL),3)
  TissueMatrix["NHI","LCORL_SD"] <- sd(BM_NHI_LCORL)
  TissueMatrix["WL77","LCORL_SD"] <- sd(BM_WL77_LCORL)
  TissueMatrix["NHI","LCORL_Pvalue"] <- round(t.test(BM_NHI_LCORL,BM_WL77_LCORL, alternative = "two.sided",var.equal = TRUE)$p.value,3)
  #cat(Tissue, "NHI_LCORL",shapiro.test(BM_NHI_LCORL)$p.value,"\n")
  #cat(Tissue, "WL77_LCORL",shapiro.test(BM_WL77_LCORL)$p.value,"\n")
  #cat(Tissue,"LCORL-NHI vs WL77 Homogeneity",var.test(BM_NHI_LCORL,BM_WL77_LCORL)$p.value,"\n")
  
  TissueMatrix["NHI","NCAPG_FC_WL77:NHI"] <- 1
  TissueMatrix["WL77","NCAPG_FC_WL77:NHI"] <- round(mean(BM_WL77_NCAPG)/mean(BM_NHI_NCAPG),3)
  TissueMatrix["NHI","NCAPG_SD"] <- sd(BM_NHI_NCAPG)
  TissueMatrix["WL77","NCAPG_SD"] <- sd(BM_WL77_NCAPG)
  TissueMatrix["NHI","NCAPG_Pvalue"] <- round(t.test(BM_NHI_NCAPG,BM_WL77_NCAPG, alternative = "two.sided",var.equal = TRUE)$p.value,3)
  #cat(Tissue, "NHI_NCAPG",shapiro.test(BM_NHI_NCAPG)$p.value,"\n")
  #cat(Tissue, "WL77_NCAPG",shapiro.test(BM_WL77_NCAPG)$p.value,"\n")
  #cat(Tissue,"NCAPG-NHI vs WL77 Homogeneity",var.test(BM_NHI_NCAPG,BM_WL77_NCAPG)$p.value,"\n")
  
  ExpAnalysis[[Tissue]] <- TissueMatrix
}









