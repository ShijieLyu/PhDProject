# SNP additive and dominance effect calculation
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


### The additive effects(AE) calculation for the Top SNP per traits
GetAE <- function(Trait, TopSNP){
  SNPsinfo <- read.table("RawData/SNPsinfo.txt",header=TRUE,sep="\t")
  NHI.Allele <- as.character(SNPsinfo[which(SNPsinfo[,"Markers"]==TopSNP),]["NHI.Allele"][,1])
  WL77.Allele <- as.character(SNPsinfo[which(SNPsinfo[,"Markers"]==TopSNP),]["WL77.Allele"][,1])

  Pheno <- AllCorPheno[[Trait]]
  AddGeno <- as.character(QTLdataAll[which(!is.na(QTLdataAll[,Trait])),TopSNP])
  DomGeno <- as.character(QTLdataAll[which(!is.na(QTLdataAll[,Trait])),TopSNP])

  AddGeno[which(AddGeno == paste0(WL77.Allele,WL77.Allele))] <- 0
  AddGeno[c(which(AddGeno == paste0(NHI.Allele,WL77.Allele)),which(AddGeno == paste0(WL77.Allele,NHI.Allele)))] <- 1
  AddGeno[which(AddGeno == paste0(NHI.Allele,NHI.Allele))] <- 2
  
  DomGeno[which(DomGeno == paste0(WL77.Allele,WL77.Allele))] <- 0
  DomGeno[c(which(DomGeno == paste0(NHI.Allele,WL77.Allele)),which(DomGeno == paste0(WL77.Allele,NHI.Allele)))] <- 1
  DomGeno[which(DomGeno == paste0(NHI.Allele,NHI.Allele))] <- 0

  model <- lm(as.numeric(Pheno) ~ as.numeric(AddGeno) + as.numeric(DomGeno))
  B <- summary(model)$coefficients["as.numeric(Geno)","Estimate"]
  SE <- summary(model)$coefficients["as.numeric(Geno)","Std. Error"]
  Pvalue <- summary(model)$coefficients["as.numeric(Geno)","Pr(>|t|)"]
  result <- cbind(B=B,SE=SE,Pvalue=Pvalue)
  rownames(result) <- paste0("AE:", TopSNP, "-vs-", Trait)
  return(result)
}

(mean(Pheno[which(AddGeno=="2"),]) - mean(Pheno[which(AddGeno=="0"),]))/2
mean(Pheno[which(DomGeno=="1"),])- mean(Pheno[which(DomGeno == "0"),])

### The dominance effects(DE) calculation for the Top SNP per traits

GetDE <- function(Trait, TopSNP){
  SNPsinfo <- read.table("RawData/SNPsinfo.txt",header=TRUE,sep="\t")
  NHI.Allele <- as.character(SNPsinfo[which(SNPsinfo[,"Markers"]==TopSNP),]["NHI.Allele"][,1])
  WL77.Allele <- as.character(SNPsinfo[which(SNPsinfo[,"Markers"]==TopSNP),]["WL77.Allele"][,1])

  Pheno <- AllCorPheno[[Trait]]
  Geno <- as.character(QTLdataAll[which(!is.na(QTLdataAll[,Trait])),TopSNP])
  
  Geno[which(Geno == paste0(WL77.Allele,WL77.Allele))] <- 0
  Geno[c(which(Geno == paste0(NHI.Allele,WL77.Allele)),which(Geno == paste0(WL77.Allele,NHI.Allele)))] <- 1
  Geno[which(Geno == paste0(NHI.Allele,NHI.Allele))] <- 0

  model <- lm(as.numeric(Pheno) ~ as.numeric(Geno))
  B <- summary(model)$coefficients["as.numeric(Geno)","Estimate"]
  SE <- summary(model)$coefficients["as.numeric(Geno)","Std. Error"]
  Pvalue <- summary(model)$coefficients["as.numeric(Geno)","Pr(>|t|)"]
  result <- cbind(B=B,SE=SE,Pvalue=Pvalue)
  rownames(result) <- paste0("DE:",TopSNP, "-vs-", Trait)
  return(result)
}

Trait <- "Gew_10Wo"
TopSNP <- "rs14490774"

GetAE(Trait,TopSNP)
GetDE(Trait,TopSNP)

