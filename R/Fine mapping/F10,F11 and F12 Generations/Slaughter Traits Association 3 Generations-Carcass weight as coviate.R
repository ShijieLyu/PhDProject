# Association Analysis using mixed-effects model for the F12 Chicken
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
                         TransNeckFat = rep(NA, nrow(SlaughterTraits)),TransVisceralFat = rep(NA, nrow(SlaughterTraits)),TransTotalFat = rep(NA, nrow(SlaughterTraits)),Bonemass = rep(NA, nrow(SlaughterTraits))
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
SlaughterTraits[,"Bonemass"] <- SlaughterTraits[,"LegnoSkin"] - SlaughterTraits[,"LegnoSkin_bone"]

write.table(SlaughterTraits,"Analysis/SlaughterTraits-All generation.txt",sep="\t",row.names = FALSE,quote = FALSE)

# Organised the genotypes
SelectColumn <- c("ID.Nr","rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")
genotypes <- rbind(F10genotypes[,SelectColumn],F11genotypes[,SelectColumn],F12genotypes[,SelectColumn])

# Association Analysis
AllSTName <- names(SlaughterTraits)[-c(45,54,71,82,96,109)]  # Remove the "ID.Nr" column
AllSTName <- AllSTName[c(11:115,118:128)]                    # Select the phenotype traits

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

############################################################################
library(lme4)
sst <- c("BW.bratfertig.","Kopf","Hals","Flugel","LegnoSkin","LegnoSkin_bone","BreastnoSkin","Hals.Fett","VisceralFat","TotalFat","Bonemass")
Traits  <- sst
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

GTLod <- vector("list", length(Traits)) 
names(GTLod) <- Traits

#options(warn=2)
for (Phenotype in Traits){                                                     
  EachTraitLod <- NULL
  if(length(which(!is.na(SlaughterTraits[,Phenotype]))) == 0){
        EachTraitLod <- rbind(EachTraitLod, NA)
        colnames(EachTraitLod) <- "SNP"
  }else{
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(SlaughterTraits[,Phenotype]))
    model.full <- lmer((as.numeric(SlaughterTraits[,Phenotype])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])) + (as.numeric(SlaughterTraits[,"BW.nuchtern"])[idx]) + (as.character(genotypes[,OneSNP])[idx]), REML=FALSE) # qqnorm(resid(model.full)) 
    model.null <- lmer((as.numeric(SlaughterTraits[,Phenotype])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])) + (as.numeric(SlaughterTraits[,"BW.nuchtern"])[idx]), REML=FALSE)
    res <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[8]]),SNPvar))
    colnames(EachTraitLod) <- c("Residuals","SNP","SNPvar")
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  GTLod[[Phenotype]] <- EachTraitLod
  }
}

#######################################################################################################################################
library(lme4)
sst <- c("BW.bratfertig.","Kopf","Hals","Flugel","LegnoSkin","LegnoSkin_bone","BreastnoSkin","Hals.Fett","VisceralFat","TotalFat")
Traits  <- sst
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

GTLod <- vector("list", length(Traits)) 
names(GTLod) <- Traits

#options(warn=2)
for (Phenotype in Traits){                                                     
  EachTraitLod <- NULL
  if(length(which(!is.na(SlaughterTraits[,Phenotype]))) == 0){
        EachTraitLod <- rbind(EachTraitLod, NA)
        colnames(EachTraitLod) <- "SNP"
  }else{
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(SlaughterTraits[,Phenotype]))
    model.full <- lmer(((as.numeric(SlaughterTraits[,Phenotype])[idx])/(as.numeric(SlaughterTraits[,"BW.nuchtern"])[idx])) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])) + (as.character(genotypes[,OneSNP])[idx]), REML=FALSE) # qqnorm(resid(model.full)) 
    model.null <- lmer(((as.numeric(SlaughterTraits[,Phenotype])[idx])/(as.numeric(SlaughterTraits[,"BW.nuchtern"])[idx])) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])), REML=FALSE) 
    res <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[8]]),SNPvar))
    colnames(EachTraitLod) <- c("Residuals","SNP","SNPvar")
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  GTLod[[Phenotype]] <- EachTraitLod
  }
}