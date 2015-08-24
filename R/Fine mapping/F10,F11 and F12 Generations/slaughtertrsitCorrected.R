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
                         BreastnoSkin = rep(NA, nrow(SlaughterTraits)), BreastSkin = rep(NA, nrow(SlaughterTraits)), VisceralFat = rep(NA, nrow(SlaughterTraits)), TotalFat = rep(NA, nrow(SlaughterTraits)))

SlaughterTraits[,"LegnoSkin"] <- SlaughterTraits[,"Keule.re..ohne.Haut"] + SlaughterTraits[,"Keule.li..ohne.Haut"]
SlaughterTraits[,"LegSkin"] <- SlaughterTraits[,"Keule.Haut.re"] + SlaughterTraits[,"Keule.Haut.li."]
SlaughterTraits[,"LegnoSkin_bone"] <- SlaughterTraits[,"Keule.re..ohne.Haut.und.ohne.Knochen"] + SlaughterTraits[,"Keule.li..ohne.Haut.u.ohne.Knochen"]
SlaughterTraits[,"BreastnoSkin"] <- SlaughterTraits[,"Brust.re..ohne.Haut"] + SlaughterTraits[,"Brust.li..ohne.Haut"]
SlaughterTraits[,"BreastSkin"] <- SlaughterTraits[,"Brust.Haut.re."] + SlaughterTraits[,"Brust.Haut.li."]
SlaughterTraits[,"VisceralFat"] <- SlaughterTraits[,"Fett.Herz."]+SlaughterTraits[,"Fett.Magen."] + SlaughterTraits[,"Fett.Leber."] + SlaughterTraits[,"Fett.Milz."] + SlaughterTraits[,"Abdominalfett"]
SlaughterTraits[,"TotalFat"] <- SlaughterTraits[,"Fett.Herz."]+SlaughterTraits[,"Fett.Magen."] + SlaughterTraits[,"Fett.Leber."] + SlaughterTraits[,"Fett.Milz."] + SlaughterTraits[,"Abdominalfett"] + SlaughterTraits[,"Hals.Fett"]


# Organised the genotypes
SelectColumn <- c("ID.Nr","rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")
genotypes <- rbind(F10genotypes[,SelectColumn],F11genotypes[,SelectColumn],F12genotypes[,SelectColumn])

# Selected slaughter traits
library(lme4)
sst <- c("BW.bratfertig.","Kopf","Hals","Flugel","LegnoSkin","LegnoSkin_bone","BreastnoSkin","Hals.Fett","VisceralFat","TotalFat")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

AllCorPheno <- vector("list", length(sst)) 
names(AllCorPheno) <- sst

for (Phenotype in sst){
  idx <- which(!is.na(SlaughterTraits[,Phenotype]))
  onlyEnv <- lmer((as.numeric(SlaughterTraits[,Phenotype])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])), REML=FALSE)
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

### Plot for carcass weight and TotalFat
CWID <- which(!is.na(SlaughterTraits[,"BW.bratfertig."]))
CWGeno <- genotypes[CWID,"rs14490774"]
NCWGeno <- factor(CWGeno, levels=levels(CWGeno)<-c("TT","CT","CC")) 
CW <- AllCorPheno[["BW.bratfertig."]]
CWdf <- data.frame(as.numeric(CW),as.character(CWGeno))
colnames(CWdf) <- c("CorrectedCW","Geno")

TFID <- which(!is.na(SlaughterTraits[,"TotalFat"]))
TFGeno <- genotypes[TFID,"rs14490774"]
NTFGeno <- factor(TFGeno, levels=levels(TFGeno)<-c("TT","CT","CC")) 
TF <- AllCorPheno[["TotalFat"]]
TFdf <- data.frame(as.numeric(TF),as.character(TFGeno))  # wilcox.test, shapiro.test
colnames(TFdf) <- c("CorrectedTF","Geno")

LNSBID <- which(!is.na(SlaughterTraits[,"LegnoSkin_bone"]))
LNSBGeno <- genotypes[LNSBID,"rs14490774"]
NLNSBGeno <- factor(LNSBGeno, levels=levels(LNSBGeno)<-c("TT","CT","CC")) 
LNSB <- AllCorPheno[["LegnoSkin_bone"]]
LNSBdf <- data.frame(as.numeric(LNSB),as.character(LNSBGeno))  # wilcox.test, shapiro.test
colnames(LNSBdf) <- c("CorrectedLNSB","Geno")

BNSID <- which(!is.na(SlaughterTraits[,"BreastnoSkin"]))
BNSGeno <- genotypes[BNSID,"rs14490774"]
NBNSGeno <- factor(BNSGeno, levels=levels(BNSGeno)<-c("TT","CT","CC")) 
BNS <- AllCorPheno[["BreastnoSkin"]]
BNSdf <- data.frame(as.numeric(BNS),as.character(BNSGeno))  # wilcox.test, shapiro.test
colnames(BNSdf) <- c("CorrectedBNS","Geno")

tiff("Analysis/4 traits phenotypes.tif", res = 300,width = 2400, height = 2400, compression = "lzw")
par(mfrow=c(2,2))
boxplot(CW~NCWGeno, ylab= "Carcass weight (g)", xlab = "Genotypes", sub = "rs14490774" ) #text(c(1,2,3),1390, c("n=41","n=13", "n=2"), cex= 0.5) 
boxplot(LNSB~NLNSBGeno, ylab= "Leg weight without skin and bone (g)", xlab = "Genotypes", sub = "rs14490774" )
boxplot(BNS~NBNSGeno, ylab= "Breast weight without skin (g)", xlab = "Genotypes", sub = "rs14490774" )
boxplot(TF~NTFGeno, ylab= "Total adipose tissues weight (g)", xlab = "Genotypes", sub = "rs14490774" )
dev.off()




