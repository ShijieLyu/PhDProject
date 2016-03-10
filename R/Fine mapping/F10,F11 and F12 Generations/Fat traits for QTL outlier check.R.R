# Check the fat traits. is there outlier? 
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/FineMapping")
F10genotypes <- read.table("RawData/F10_Genotypes.txt", na.strings="-", header=TRUE, sep="\t")
F11genotypes <- read.table("RawData/F11_Genotypes.txt", na.strings="-", header=TRUE, sep="\t")
F12genotypes <- read.table("RawData/F12_Genotypes-FILL IN.txt", na.strings="-", header=TRUE, sep="\t")
SlaughterTraits <- read.table("RawData/Slaughter_Data_All_Generations-MRI.txt", na.strings=c("NA","?"," "),header=TRUE, sep="\t", colClasses = c("character", rep("factor",9), rep("numeric",86)))
QTLdataAll <- read.table("Analysis/QTLdataAll.txt", sep = "\t",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("character", 9), rep("factor",6)))

F12IDForST <- SlaughterTraits[which(SlaughterTraits[,"Generation"]=="F12"),][,"ID.Nr"]    # Select and order the F12 chicken ID according to the slaughter traits
F12genotypes <- F12genotypes[match(F12IDForST,F12genotypes[,"ID.Nr"]),]

BandP <- QTLdataAll[match(SlaughterTraits[,"ID.Nr"],QTLdataAll[,"ID.Nr"]),][,c("Batch","Parents")] # Add the Batch and Parents information
SlaughterTraits <- cbind(SlaughterTraits, BandP)

SlaughterTraits <- cbind(SlaughterTraits, LegnoSkin = rep(NA, nrow(SlaughterTraits)), LegSkin = rep(NA, nrow(SlaughterTraits)), LegnoSkin_bone= rep(NA, nrow(SlaughterTraits)),
                         BreastnoSkin = rep(NA, nrow(SlaughterTraits)), BreastSkin = rep(NA, nrow(SlaughterTraits)), VisceralFat = rep(NA, nrow(SlaughterTraits)), TotalFat = rep(NA, nrow(SlaughterTraits)),
                         TransNeckFat = rep(NA, nrow(SlaughterTraits)),TransVisceralFat = rep(NA, nrow(SlaughterTraits)),TransTotalFat = rep(NA, nrow(SlaughterTraits)),
                         BreastMRI = rep(NA, nrow(SlaughterTraits)),LegMRI = rep(NA, nrow(SlaughterTraits)),
                         BreastWHC1h = rep(NA, nrow(SlaughterTraits)),BreastWHC24h = rep(NA, nrow(SlaughterTraits)),BreastWHC48h = rep(NA, nrow(SlaughterTraits)),
                         BreastWHC72h = rep(NA, nrow(SlaughterTraits)),BreastWHC6Tage = rep(NA, nrow(SlaughterTraits)),LegWHC1h = rep(NA, nrow(SlaughterTraits)),
                         LegWHC24h = rep(NA, nrow(SlaughterTraits)),LegWHC48h = rep(NA, nrow(SlaughterTraits)),LegWHC72h = rep(NA, nrow(SlaughterTraits)),
                         LegWHC6Tage = rep(NA, nrow(SlaughterTraits)))

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

Neg2NA <- function(dat){dat[dat <= 0] <- NA; return(dat)}
SlaughterTraits[,"BreastMRI"] <- (Neg2NA(SlaughterTraits[,"MRI.Brust.Fett"]) + Neg2NA(SlaughterTraits[,"MRI.Brust.Fett.1"]))/2/Neg2NA((SlaughterTraits[,"Gewicht.TeilstÃ.ck.Brust.fÃ.r.Messung"]))*100
SlaughterTraits[,"LegMRI"] <- (Neg2NA(SlaughterTraits[,"MRI.Keule.ohne.Knochen.Fett"]) + Neg2NA(SlaughterTraits[,"MRI.Keule.ohne.Knochen.Fett.1"]))/2/Neg2NA((SlaughterTraits[,"Gewicht.TeilstÃ.ck.Keule.ohne.Knochen.f..Messung"]))*100

SlaughterTraits[,"BreastWHC1h"] <- SlaughterTraits[,"Gewicht.Brust.li.TeilstÃ.ck.nach.1h"]/SlaughterTraits[,"Gewicht.Brust.li.TeilstÃ.ck"]*100
SlaughterTraits[,"BreastWHC24h"] <- SlaughterTraits[,"Gewicht.Brust.li..24h"]/SlaughterTraits[,"Gewicht.Brust.li.TeilstÃ.ck"]*100
SlaughterTraits[,"BreastWHC48h"] <- SlaughterTraits[,"Gewicht.Brust.li..48h"]/SlaughterTraits[,"Gewicht.Brust.li.TeilstÃ.ck"]*100
SlaughterTraits[,"BreastWHC72h"] <- SlaughterTraits[,"Gewicht.Brust.li..72h"]/SlaughterTraits[,"Gewicht.Brust.li.TeilstÃ.ck"]*100
SlaughterTraits[,"BreastWHC6Tage"] <- SlaughterTraits[,"Gewicht.Brust.li..6Tage"]/SlaughterTraits[,"Gewicht.Brust.li.TeilstÃ.ck"]*100

SlaughterTraits[,"LegWHC1h"] <- SlaughterTraits[,"Gewicht.Keule.li..TeilstÃ.ck.nach.1h"]/SlaughterTraits[,"Gewicht.Keule.li..TeilstÃ.ck"]*100
SlaughterTraits[,"LegWHC24h"] <- SlaughterTraits[,"Gewicht.Keule.li..24h"]/SlaughterTraits[,"Gewicht.Keule.li..TeilstÃ.ck"]*100
SlaughterTraits[,"LegWHC48h"] <- SlaughterTraits[,"Gewicht.Keule.li..48h"]/SlaughterTraits[,"Gewicht.Keule.li..TeilstÃ.ck"]*100
SlaughterTraits[,"LegWHC72h"] <- SlaughterTraits[,"Gewicht.Keule.li..72h"]/SlaughterTraits[,"Gewicht.Keule.li..TeilstÃ.ck"]*100
SlaughterTraits[,"LegWHC6Tage"] <- SlaughterTraits[,"Gewicht.Keule.li..6Tage"]/SlaughterTraits[,"Gewicht.Keule.li..TeilstÃ.ck"]*100

#write.table(SlaughterTraits,"Analysis/SlaughterTraits-All generation.txt",sep="\t",row.names = FALSE,quote = FALSE)

# Organised the genotypes
SelectColumn <- c("ID.Nr","rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")
genotypes <- rbind(F10genotypes[,SelectColumn],F11genotypes[,SelectColumn],F12genotypes[,SelectColumn])

# Get the varance of the fix effect

SNPvarPerc <- function(mF){
  X <- getME(mF,"X")
  TotalFixed <- 0
  for (x in 2:length(fixef(mF))){  
    TotalFixed <- TotalFixed + (fixef(mF)[x] * X[, x])  
  }
  SNPFixed <- 0
  for (x in grep("OneSNP", names(fixef(mF)))){     
    SNPFixed <- SNPFixed + (fixef(mF)[x] * X[, x])
  }
  FixedTotalVar <- var(TotalFixed)
  FixedSNPVar <- var(SNPFixed)
  TotalVar <- FixedTotalVar + sum(data.frame(VarCorr(mF))[,4])
  Perc <- round(FixedSNPVar/TotalVar,3)*100
  return(Perc)
}

#########################################
## correcting the outlier and re-analyze the fat traits with the model no BodyWeight as covariate

library(robustHD)
library(lme4)
FatTraits <- SlaughterTraits[,c("Hals.Fett","VisceralFat","TotalFat")]
ReOutFattraits <- apply(FatTraits,2,winsorize)                                           # Correcting the outliers using winsorize
colnames(ReOutFattraits) <- c("ReOutHals.Fett","ReOutVisceralFat","ReOutTotalFat")
Phenotypes <- cbind(FatTraits,ReOutFattraits)

boxplot(Phenotypes)
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

NoBWLod <- vector("list", length(Phenotypes)) 
names(NoBWLod) <- colnames(Phenotypes)

for (OnePheno in colnames(Phenotypes)){
  EachTraitLod <- NULL                                                     
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(Phenotypes[,OnePheno]))
    model.full <- lmer((as.numeric(Phenotypes[,OnePheno])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])) + (as.character(genotypes[,OneSNP])[idx]), REML=FALSE) # qqnorm(resid(model.full)) 
    model.null <- lmer((as.numeric(Phenotypes[,OnePheno])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])), REML=FALSE)
    res <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[8]]),SNPvar))
    colnames(EachTraitLod) <- c("Residuals","SNP","SNPvar")
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  NoBWLod[[OnePheno]] <- EachTraitLod
}

#########################################
## correcting the outlier and re-analyze the fat traits with the model no BodyWeight as covariate

WithBWLod <- vector("list", length(Phenotypes)) 
names(WithBWLod) <- colnames(Phenotypes)

for (OnePheno in colnames(Phenotypes)){
  EachTraitLod <- NULL                                                     
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(Phenotypes[,OnePheno]))
    model.full <- lmer((as.numeric(Phenotypes[,OnePheno])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])) + (as.numeric(SlaughterTraits[,"BW.nuchtern"])[idx])+(as.character(genotypes[,OneSNP])[idx]), REML=FALSE) # qqnorm(resid(model.full)) 
    model.null <- lmer((as.numeric(Phenotypes[,OnePheno])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])) + (as.numeric(SlaughterTraits[,"BW.nuchtern"])[idx]), REML=FALSE)
    res <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[8]]),SNPvar))
    colnames(EachTraitLod) <- c("Residuals","SNP","SNPvar")
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  WithBWLod[[OnePheno]] <- EachTraitLod
}

#########################################
## correcting the outlier and re-analyze the fat traits with the model no BodyWeight as covariate

FatperBWLod <- vector("list", length(Phenotypes)) 
names(FatperBWLod) <- colnames(Phenotypes)

for (OnePheno in colnames(Phenotypes)){
  EachTraitLod <- NULL                                                     
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(Phenotypes[,OnePheno]))
    model.full <- lmer(((as.numeric(Phenotypes[,OnePheno])[idx])/(as.numeric(SlaughterTraits[,"BW.nuchtern"])[idx])) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])) +(as.character(genotypes[,OneSNP])[idx]), REML=FALSE) # qqnorm(resid(model.full)) 
    model.null <- lmer(((as.numeric(Phenotypes[,OnePheno])[idx])/(as.numeric(SlaughterTraits[,"BW.nuchtern"])[idx])) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])), REML=FALSE)
    res <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[8]]),SNPvar))
    colnames(EachTraitLod) <- c("Residuals","SNP","SNPvar")
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  FatperBWLod[[OnePheno]] <- EachTraitLod
}

#############################################################
### summarize the data and make comparison

ResultList <- vector("list", length(Phenotypes)) 
names(ResultList) <- colnames(Phenotypes)

for (OnePheno in colnames(Phenotypes)){
  AllResults <- cbind(NoBWLod[[OnePheno]],WithBWLod[[OnePheno]],FatperBWLod[[OnePheno]])
  ResultList[[OnePheno]] <- AllResults
}

###############################################################
#### Check the difference between every group






