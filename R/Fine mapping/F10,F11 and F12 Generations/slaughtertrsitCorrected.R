# Correct the selected slaughter traits data for the chicken (All 3 Generations)
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
SlaughterTraits[,"LegBonemass"] <- SlaughterTraits[,"LegnoSkin"] - SlaughterTraits[,"LegnoSkin_bone"]

# Organised the genotypes
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
  onlyEnv <- lmer((as.numeric(SlaughterTraits[,Phenotype])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx]))+as.numeric(SlaughterTraits[,"BW.nuchtern"][idx]), REML=FALSE)
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
### QTL effect on Breast muscle 
BNSID <- which(!is.na(SlaughterTraits[,"BreastnoSkin"]))
BNSGeno <- as.character(genotypes[BNSID,"rs14490774"])
BNS <- SlaughterTraits[BNSID,"BreastnoSkin"]
mean(BNS[which(BNSGeno=="CC")])-mean(BNS[which(BNSGeno=="TT")])

### Plot for carcass weight and TotalFat
CWID <- which(!is.na(SlaughterTraits[,"BW.bratfertig."]))
CWGeno <- as.character(genotypes[CWID,"rs14490774"])
NCWGeno <- factor(CWGeno, levels=levels(CWGeno)<-c("TT","CT","CC")) 
CW <- SlaughterTraits[CWID,"BW.bratfertig."]                                   #AllCorPheno[["BW.bratfertig."]]
CWdf <- data.frame(as.numeric(CW),as.character(CWGeno))
lm(CWdf[,1]~as.numeric(CWdf[,2]))
colnames(CWdf) <- c("CorrectedCW","Geno")

TFID <- which(!is.na(SlaughterTraits[,"TotalFat"]))
TFGeno <- as.character(genotypes[TFID,"rs14490774"])
NTFGeno <- factor(TFGeno, levels=levels(TFGeno)<-c("TT","CT","CC")) 
TF <- SlaughterTraits[TFID,"TotalFat"] 
TFdf <- data.frame(as.numeric(TF),as.character(TFGeno))  # wilcox.test, shapiro.test
colnames(TFdf) <- c("CorrectedTF","Geno")

LNSBID <- which(!is.na(SlaughterTraits[,"LegnoSkin_bone"]))
LNSBGeno <- as.character(genotypes[LNSBID,"rs14490774"])
NLNSBGeno <- factor(LNSBGeno, levels=levels(LNSBGeno)<-c("TT","CT","CC")) 
LNSB <- SlaughterTraits[LNSBID,"LegnoSkin_bone"]
LNSBdf <- data.frame(as.numeric(LNSB),as.character(LNSBGeno))  # wilcox.test, shapiro.test
colnames(LNSBdf) <- c("CorrectedLNSB","Geno")

LBMID <- which(!is.na(SlaughterTraits[,"LegBonemass"]))
LBMGeno <- as.character(genotypes[LBMID,"rs14490774"])
NLBMGeno <- factor(LNSBGeno, levels=levels(LBMGeno)<-c("TT","CT","CC")) 
LBM <- SlaughterTraits[LBMID,"LegBonemass"]
LBMdf <- data.frame(as.numeric(LBM),as.character(LBMGeno))  # wilcox.test, shapiro.test
colnames(LBMdf) <- c("CorrectedLBM","Geno")

BNSID <- which(!is.na(SlaughterTraits[,"BreastnoSkin"]))
BNSGeno <- as.character(genotypes[BNSID,"rs14490774"])
NBNSGeno <- factor(BNSGeno, levels=levels(BNSGeno)<-c("TT","CT","CC")) 
BNS <- SlaughterTraits[BNSID,"BreastnoSkin"]
BNSdf <- data.frame(as.numeric(BNS),as.character(BNSGeno))  # wilcox.test, shapiro.test
colnames(BNSdf) <- c("CorrectedBNS","Geno")

Trans <- function(CWdf){
  Geno <- as.character(CWdf[,2])
  Geno[which(Geno == "TT")] <- 1
  Geno[which(Geno == "CT")] <- 2
  Geno[which(Geno == "CC")] <- 3
  return(as.numeric(Geno))
}


tiff("Analysis/4 traits phenotypes.tif", res = 300,width = 2400, height = 2400, compression = "lzw")
par(mfrow=c(2,2))
boxplot(CW~NCWGeno, ylab= "CW (g)", xlab = "Genotypes", sub = "rs14490774", cex.lab=1.3) #text(c(1,2,3),1390, c("n=41","n=13", "n=2"), cex= 0.5) 
abline(lm(CWdf[,1]~Trans(CWdf)),col="red")
legend("bottomright", legend = paste0("P = ",signif(summary(lm(CWdf[,1]~Trans(CWdf)))$coefficients[2,4],3)),bty = "n")

boxplot(TF~NTFGeno, ylab= "WAT (g)", xlab = "Genotypes", sub = "rs14490774", cex.lab=1.3)
abline(lm(TFdf[,1]~Trans(TFdf)),col="red")
legend("bottomright", legend = paste0("P = ",signif(summary(lm(TFdf[,1]~Trans(TFdf)))$coefficients[2,4],3)),bty = "n")

boxplot(LBM~NLBMGeno, ylab= "Leg bone mass (g)", xlab = "Genotypes", sub = "rs14490774", cex.lab=1.3)
abline(lm(LBMdf[,1]~Trans(LBMdf)),col="red")
legend("bottomright", legend = paste0("P = ",signif(summary(lm(LBMdf[,1]~Trans(LBMdf)))$coefficients[2,4],3)),bty = "n")

boxplot(BNS~NBNSGeno, ylab= "Breast muscle weight (g)", xlab = "Genotypes", sub = "rs14490774", cex.lab=1.3)
abline(lm(BNSdf[,1]~Trans(BNSdf)),col="red")
legend("bottomright", legend = paste0("P = ",signif(summary(lm(BNSdf[,1]~Trans(BNSdf)))$coefficients[2,4],3)),bty = "n")
dev.off()




