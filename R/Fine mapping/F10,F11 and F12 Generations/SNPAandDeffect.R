# SNP additive and dominance effect calculation
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/FineMapping")
CorPheno <- read.table("Analysis/MLM_CorPheno.txt",sep="\t",header=TRUE)

### The additive effects(AE) calculation for the Top SNP per traits
GetAE <- function(Trait, TopSNP){
  SNPsinfo <- read.table("RawData/SNPsinfo.txt",header=TRUE,sep="\t")
  NHI.Allele <- as.character(SNPsinfo[which(SNPsinfo[,"Markers"]==TopSNP),]["NHI.Allele"][,1])
  WL77.Allele <- as.character(SNPsinfo[which(SNPsinfo[,"Markers"]==TopSNP),]["WL77.Allele"][,1])

  Pheno <- CorPheno[,Trait]
  Geno <- as.character(CorPheno[,TopSNP])

  Geno[which(Geno == paste0(WL77.Allele,WL77.Allele))] <- 0
  Geno[c(which(Geno == paste0(NHI.Allele,WL77.Allele)),which(Geno == paste0(WL77.Allele,NHI.Allele)))] <- 1
  Geno[which(Geno == paste0(NHI.Allele,NHI.Allele))] <- 2

  model <- lm(as.numeric(Pheno) ~ as.numeric(Geno))
  B <- summary(model)$coefficients["as.numeric(Geno)","Estimate"]
  SE <- summary(model)$coefficients["as.numeric(Geno)","Std. Error"]
  Pvalue <- summary(model)$coefficients["as.numeric(Geno)","Pr(>|t|)"]
  result <- cbind(B=B,SE=SE,Pvalue=Pvalue)
  rownames(result) <- paste0("AE:", TopSNP, "-vs-", Trait)
  return(result)
}

### The dominance effects(DE) calculation for the Top SNP per traits

GetDE <- function(Trait, TopSNP){
  SNPsinfo <- read.table("RawData/SNPsinfo.txt",header=TRUE,sep="\t")
  NHI.Allele <- as.character(SNPsinfo[which(SNPsinfo[,"Markers"]==TopSNP),]["NHI.Allele"][,1])
  WL77.Allele <- as.character(SNPsinfo[which(SNPsinfo[,"Markers"]==TopSNP),]["WL77.Allele"][,1])

  Pheno <- CorPheno[,Trait]
  Geno <- as.character(CorPheno[,TopSNP])

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


Trait <- "Gew_1d"
TopSNP <- "rs14490774"

GetAE(Trait,TopSNP)
GetDE(Trait,TopSNP)

