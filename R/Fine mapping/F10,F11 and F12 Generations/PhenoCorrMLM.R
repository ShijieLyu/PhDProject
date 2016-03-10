# Correct the phenotypic data for the chicken (All 3 Generations)
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/FineMapping")
QTLdataAll <- read.table("Analysis/QTLdataAll.txt", sep = "\t",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("character", 9), rep("factor",6)))

library(lme4)
GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

AllCorPheno <- vector("list", length(GrowthTraits)) 
names(AllCorPheno) <- GrowthTraits

for (Phenotype in GrowthTraits){
  idx <- which(!is.na(QTLdataAll[,Phenotype]))
  onlyEnv <- lmer(QTLdataAll[,Phenotype][idx] ~ QTLdataAll[,"Batch"][idx] + (1|QTLdataAll[,"Parents"][idx]), REML=FALSE) 
  Intercept <- coef(onlyEnv)$`QTLdataAll[, "Parents"][idx]`["(Intercept)"]
  corPheno <- NULL
  for (x in 1:length(idx)){
    corPhenoind <- round(resid(onlyEnv)[x] + Intercept[which(rownames(Intercept)==QTLdataAll[idx[x],"Parents"]),])
    corPheno <- rbind(corPheno,corPhenoind)
  }
  rownames(corPheno) <- QTLdataAll[idx,"Parents"]
  colnames(corPheno) <- paste0("Corrected",Phenotype)
  AllCorPheno[[Phenotype]] <- corPheno
}

### QTL effect size on BW20
BW20ID <- which(!is.na(QTLdataAll[,"Gew_20Wo"]))
BW20 <- QTLdataAll[BW20ID,"Gew_20Wo"]
BW20Geno <- QTLdataAll[BW20ID,"rs14490774"]
mean(BW20[which(BW20Geno=="CC")])-mean(BW20[which(BW20Geno=="TT")])

### Plot for BW10, BW20, BW510, BW1520

BW10ID <- which(!is.na(QTLdataAll[,"Gew_10Wo"]))
BW10Geno <- QTLdataAll[BW10ID,"rs14490774"]
NBW10Geno <- factor(BW10Geno, levels=levels(BW10Geno)<-c("TT","CT","CC")) 
BW10 <- QTLdataAll[BW10ID,"Gew_10Wo"]
BW10df <- data.frame(as.numeric(BW10),as.character(BW10Geno))

BW20ID <- which(!is.na(QTLdataAll[,"Gew_20Wo"]))
BW20Geno <- QTLdataAll[BW20ID,"rs14490774"]
NBW20Geno <- factor(BW20Geno, levels=levels(BW20Geno)<-c("TT","CT","CC")) 
BW20 <- QTLdataAll[BW20ID,"Gew_20Wo"]
BW20df <- data.frame(as.numeric(BW20),as.character(BW20Geno))

BWG510ID <- which(!is.na(QTLdataAll[,"BWG510"]))
BWG510Geno <- QTLdataAll[BWG510ID,"rs14490774"]
NBWG510Geno <- factor(BWG510Geno, levels=levels(BWG510Geno)<-c("TT","CT","CC")) 
BWG510 <- QTLdataAll[BWG510ID,"BWG510"]
BWG510df <- data.frame(as.numeric(BWG510),as.character(BWG510Geno))

BWG1520ID <- which(!is.na(QTLdataAll[,"BWG1520"]))
BWG1520Geno <- QTLdataAll[BWG1520ID,"rs14490774"]
NBWG1520Geno <- factor(BWG1520Geno, levels=levels(BWG1520Geno)<-c("TT","CT","CC")) 
BWG1520 <- QTLdataAll[BWG1520ID,"BWG1520"]
BWG1520df <- data.frame(as.numeric(BWG1520),as.character(BWG1520Geno))

Trans <- function(CWdf){
  Geno <- as.character(CWdf[,2])
  Geno[which(Geno == "TT")] <- 1
  Geno[which(Geno == "CT")] <- 2
  Geno[which(Geno == "CC")] <- 3
  return(as.numeric(Geno))
}

tiff("Analysis/4 growth traits phenotypes.tif", res = 300,width = 2400, height = 2400, compression = "lzw")
par(mfrow=c(2,2))
boxplot(BW10~NBW10Geno, ylab= "Body weight at 10 weeks (g)", xlab = "Genotypes", sub = "rs14490774", cex.lab=1.3) 
abline(lm(BW10df[,1]~Trans(BW10df)),col="red")
legend("bottomright", legend = paste0("P = ",signif(summary(lm(BW10df[,1]~Trans(BW10df)))$coefficients[2,4],3)),bty = "n")

boxplot(BW20~NBW20Geno, ylab= "Body weight at 20 weeks (g)", xlab = "Genotypes", sub = "rs14490774", cex.lab=1.3) 
abline(lm(BW20df[,1]~Trans(BW20df)),col="red")
legend("bottomright", legend = paste0("P = ",signif(summary(lm(BW20df[,1]~Trans(BW20df)))$coefficients[2,4],3)),bty = "n")

boxplot(BWG510~NBWG510Geno, ylab= "Body weight gain from 5 to 10 weeks (g)", xlab = "Genotypes", sub = "rs14490774", cex.lab=1.3) 
abline(lm(BWG510df[,1]~Trans(BWG510df)),col="red")
legend("bottomright", legend = paste0("P = ",signif(summary(lm(BWG510df[,1]~Trans(BWG510df)))$coefficients[2,4],3)),bty = "n")

boxplot(BWG1520~NBWG1520Geno, ylab= "BWG1520 (g)", xlab = "Genotypes", sub = "rs14490774", cex.lab=1.3) 
abline(lm(BWG1520df[,1]~Trans(BWG1520df)),col="red")
legend("bottomright", legend = paste0("P = ",signif(summary(lm(BWG1520df[,1]~Trans(BWG1520df)))$coefficients[2,4],3)),bty = "n")
dev.off()
