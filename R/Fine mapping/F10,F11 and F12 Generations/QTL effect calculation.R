# QTL effect calculation
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/FineMapping")
F10genotypes <- read.table("RawData/F10_Genotypes.txt", na.strings="-", header=TRUE, sep="\t")
F11genotypes <- read.table("RawData/F11_Genotypes.txt", na.strings="-", header=TRUE, sep="\t")
F12genotypes <- read.table("RawData/F12_Genotypes-FILL IN.txt", na.strings="-", header=TRUE, sep="\t")
SlaughterTraits <- read.table("RawData/Slaughter_Data_All_Generations.txt", na.strings=c("NA","?"," "),header=TRUE, sep="\t", colClasses = c("character", rep("factor",9), rep("numeric",111)))
QTLdataAll <- read.table("Analysis/QTLdataAll.txt", sep = "\t",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("character", 9), rep("factor",6)))

F12IDForST <- SlaughterTraits[which(SlaughterTraits[,"Generation"]=="F12"),][,"ID.Nr"]    # Select and order the F12 chicken ID according to the slaughter traits
F12genotypes <- F12genotypes[match(F12IDForST,F12genotypes[,"ID.Nr"]),]

BandP <- QTLdataAll[match(SlaughterTraits[,"ID.Nr"],QTLdataAll[,"ID.Nr"]),][,c("Batch","Parents")] # Add the Batch and Parents information
SlaughterTraits <- cbind(SlaughterTraits, BandP)

SlaughterTraits <- cbind(SlaughterTraits, LegnoSkin = rep(NA, nrow(SlaughterTraits)), LegSkin = rep(NA, nrow(SlaughterTraits)), LegnoSkin_bone= rep(NA, nrow(SlaughterTraits)),
                         BreastnoSkin = rep(NA, nrow(SlaughterTraits)), BreastSkin = rep(NA, nrow(SlaughterTraits)), VisceralFat = rep(NA, nrow(SlaughterTraits)), TotalFat = rep(NA, nrow(SlaughterTraits)),
                         TransNeckFat = rep(NA, nrow(SlaughterTraits)),TransVisceralFat = rep(NA, nrow(SlaughterTraits)),TransTotalFat = rep(NA, nrow(SlaughterTraits)),LegBonemass = rep(NA, nrow(SlaughterTraits))
                         )

SlaughterTraits[,"LegnoSkin"] <- SlaughterTraits[,"Keule.re..ohne.Haut"] + SlaughterTraits[,"Keule.li..ohne.Haut"]
SlaughterTraits[,"LegSkin"] <- SlaughterTraits[,"Keule.Haut.re"] + SlaughterTraits[,"Keule.Haut.li."]
SlaughterTraits[,"LegnoSkin_bone"] <- SlaughterTraits[,"Keule.re..ohne.Haut.und.ohne.Knochen"] + SlaughterTraits[,"Keule.li..ohne.Haut.u.ohne.Knochen"]
SlaughterTraits[,"BreastnoSkin"] <- SlaughterTraits[,"Brust.re..ohne.Haut"] + SlaughterTraits[,"Brust.li..ohne.Haut"]
SlaughterTraits[,"BreastSkin"] <- SlaughterTraits[,"Brust.Haut.re."] + SlaughterTraits[,"Brust.Haut.li."]
SlaughterTraits[,"VisceralFat"] <- SlaughterTraits[,"Fett.Herz."]+SlaughterTraits[,"Fett.Magen."] + SlaughterTraits[,"Fett.Leber."] + SlaughterTraits[,"Fett.Milz."] + SlaughterTraits[,"Abdominalfett"]
SlaughterTraits[,"TotalFat"] <- SlaughterTraits[,"Fett.Herz."]+SlaughterTraits[,"Fett.Magen."] + SlaughterTraits[,"Fett.Leber."] + SlaughterTraits[,"Fett.Milz."] + SlaughterTraits[,"Abdominalfett"] + SlaughterTraits[,"Hals.Fett"]
SlaughterTraits[,"TransNeckFat"] <- (SlaughterTraits[,"Hals.Fett"])^0.5
SlaughterTraits[,"TransVisceralFat"] <- (SlaughterTraits[,"VisceralFat"])^0.5
SlaughterTraits[,"TransTotalFat"] <- (SlaughterTraits[,"TotalFat"])^0.5
SlaughterTraits[,"LegBonemass"] <- SlaughterTraits[,"LegnoSkin"] - SlaughterTraits[,"LegnoSkin_bone"]# Organised the genotypes
SelectColumn <- c("ID.Nr","rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")
genotypes <- rbind(F10genotypes[,SelectColumn],F11genotypes[,SelectColumn],F12genotypes[,SelectColumn])

# Selected slaughter traits
library(lme4)
sst <- c("BW.bratfertig.","Kopf","Hals","Flugel","LegnoSkin","LegnoSkin_bone","BreastnoSkin","Hals.Fett","VisceralFat","TotalFat","LegBonemass")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

AllCorPheno <- vector("list", length(sst)) 
names(AllCorPheno) <- sst

for (Phenotype in sst){
  idx <- which(!is.na(SlaughterTraits[,Phenotype]))
  onlyEnv <- lmer((as.numeric(SlaughterTraits[,Phenotype])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])), REML=FALSE)
  #onlyEnv <- lmer((as.numeric(SlaughterTraits[,Phenotype])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx]))+as.numeric(SlaughterTraits[,"BW.nuchtern"][idx]), REML=FALSE)
  Intercept <- coef(onlyEnv)$`(as.factor(SlaughterTraits[, "Parents"])[idx])`["(Intercept)"]
  corPheno <- NULL
  for (x in 1:length(idx)){
    corPhenoind <- round(resid(onlyEnv)[x] + Intercept[which(rownames(Intercept)==SlaughterTraits[idx[x],"Parents"]),])
    corPheno <- rbind(corPheno,corPhenoind)
  }
  rownames(corPheno) <- SlaughterTraits[idx,"Parents"]
  colnames(corPheno) <- paste0("Corrected",Phenotype)
  AllCorPheno[[Phenotype]] <- corPheno
}

## QTL effect for Leg Bone mass
LBMID <- which(!is.na(SlaughterTraits[,"LegBonemass"]))
LBMGeno <- as.character(genotypes[LBMID,"rs14490774"])
NLBMGeno <- factor(LBMGeno, levels=levels(LBMGeno)<-c("TT","CT","CC")) 
LBM <- AllCorPheno$LegBonemass[LBMID]
LBMdf <- data.frame(as.numeric(LBM),as.character(LBMGeno)) 
colnames(LBMdf) <- c("CorrectedLBM","Geno")
mean(LBMdf[which(LBMdf[,"Geno"]=="CC"),"CorrectedLBM"])-mean(LBMdf[which(LBMdf[,"Geno"]=="TT"),"CorrectedLBM"])

## QTL effect for Leg muscle weight
LNSBID <- which(!is.na(SlaughterTraits[,"LegnoSkin_bone"]))
LNSBGeno <- as.character(genotypes[LNSBID,"rs14490774"])
NLNSBGeno <- factor(LNSBGeno, levels=levels(LNSBGeno)<-c("TT","CT","CC")) 
LNSB <- AllCorPheno$LegnoSkin_bone[LNSBID]
LNSBdf <- data.frame(as.numeric(LNSB),as.character(LNSBGeno)) 
colnames(LNSBdf) <- c("CorrectedLNSB","Geno")
mean(LNSBdf[which(LNSBdf[,"Geno"]=="CC"),"CorrectedLNSB"])-mean(LNSBdf[which(LNSBdf[,"Geno"]=="TT"),"CorrectedLNSB"])






