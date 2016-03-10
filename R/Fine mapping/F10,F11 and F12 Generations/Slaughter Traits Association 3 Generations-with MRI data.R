# Association Analysis using mixed-effects model for the F12 Chicken
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
#########################################################################################################
## No body weight as the covariate
library(lme4)
sst <- c("BW.bratfertig.","Kopf","Hals","Flugel","LegnoSkin","LegnoSkin_bone","BreastnoSkin","Hals.Fett","VisceralFat","TotalFat","BreastMRI","LegMRI",
         "BreastWHC1h","BreastWHC24h","BreastWHC48h","BreastWHC72h","BreastWHC6Tage","LegWHC1h","LegWHC24h","LegWHC48h","LegWHC72h","LegWHC6Tage")
Traits  <- sst
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

GTLod <- vector("list", length(Traits)) 
names(GTLod) <- Traits

for (Phenotype in Traits){                                                     
  EachTraitLod <- NULL
  if(length(which(!is.na(SlaughterTraits[,Phenotype]))) == 0){
        EachTraitLod <- rbind(EachTraitLod, NA)
        colnames(EachTraitLod) <- "SNP"
  }else{
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(SlaughterTraits[,Phenotype]))
    model.full <- lmer((as.numeric(SlaughterTraits[,Phenotype])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])) + (as.character(genotypes[,OneSNP])[idx]), REML=FALSE) # qqnorm(resid(model.full)) 
    model.null <- lmer((as.numeric(SlaughterTraits[,Phenotype])[idx]) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])), REML=FALSE)
    res <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[8]]),SNPvar))
    colnames(EachTraitLod) <- c("Residuals","SNP","SNPvar")
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  GTLod[[Phenotype]] <- EachTraitLod
  }
}

#######################################################################################################
## Body weight as the covariate 

library(lme4)
sst <- c("BW.bratfertig.","Kopf","Hals","Flugel","LegnoSkin","LegnoSkin_bone","BreastnoSkin","Hals.Fett","VisceralFat","TotalFat","BreastMRI","LegMRI",
         "BreastWHC1h","BreastWHC24h","BreastWHC48h","BreastWHC72h","BreastWHC6Tage","LegWHC1h","LegWHC24h","LegWHC48h","LegWHC72h","LegWHC6Tage")
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
    #model.full <- lmer(((as.numeric(SlaughterTraits[,Phenotype])[idx])/(as.numeric(SlaughterTraits[,"BW.nuchtern"])[idx])) ~ (as.factor(SlaughterTraits[,"Batch"])[idx]) + (1|(as.factor(SlaughterTraits[,"Parents"])[idx])) + (as.character(genotypes[,OneSNP])[idx]), REML=FALSE) # qqnorm(resid(model.full)) 
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

####################################################################################################################################################################


if(!file.exists("Analysis/SlaughterTraitsSig-3Generation.txt")){
  for(x in 1:length(GTLod)){
    SigSNP <- which(GTLod[[x]][,"SNP"] > -log10(0.05/9))                                                      # Analyse the resulting profiles, and look at which markers are above the threshold
    if (length(SigSNP) == 0) {
    cat(names(GTLod)[x], "\n", file = "Analysis/SlaughterTraitsSig-3Generation.txt", append=TRUE)
    }else{
    cat(names(GTLod)[x], names(SigSNP), "\n", file = "Analysis/SlaughterTraitsSig-3Generation.txt", append=TRUE)
    }                                                             
  }
}else{
  cat("Loading Slaughted Traits with Significant SNPs from disk\n")
  SigSNP <- read.table("Analysis/SlaughterTraitsSig-3Generation.txt",sep="\t", header=FALSE)
}

SigTraits <- NULL                                                                                      
for (x in 1:nrow(SigSNP)){
  if(length(unlist(strsplit(as.character(SigSNP[x,1]), " "))) > 1) {
  SigTraits <- c(SigTraits, strsplit(as.character(SigSNP[x,1]), " ")[[1]][1])
  }
}

par(mfrow=c(4,5))
  for(x in 101:110){
    plot(x=c(0,10), y=c(0,10), t="n", ylab="LOD Score", xlab="", xaxt="n", main= names(GTLod[x]), cex.main=0.9)
    points(GTLod[[x]][,2],t="l")
    abline(h=-log10(0.05/9), lty=2)    
    axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
  }
title("The association analysis for the Slaughter Traits", outer=TRUE)

#### Selected Slaughter traits
sst <- c("BW.bratfertig.","Kopf","Hals","Flugel","LegnoSkin","LegnoSkin_bone","BreastnoSkin","Hals.Fett","VisceralFat","TotalFat")
GTlodsst <- GTLod[sst]
Threshold <- -log10(0.05/(length(sst)*length(SNPsForAnalysis)))
SigThreshold <- -log10(0.01/(length(sst)*length(SNPsForAnalysis)))

SNPsinfo <- read.table("RawData/SNPsinfo.txt",header=TRUE,sep="\t")

#tiff("Analysis/SlaughterTrait-AllinOne.tif", res = 300,width = 2100, height = 2000, compression = "lzw")
#pdf("Analysis/SlaughterTrait-AllinOne.pdf")
par(mai = c(2.2, 1, 1, 1))
plot(x=c(as.numeric(SNPsinfo[1,3]-100000),as.numeric(SNPsinfo[9,3]+100000)), y=c(0,10), t="n", ylab="LOD Score", xlab="Physical Position (Mb)", xaxt="n")
for(x in 1:10){
  points(x= SNPsinfo[,"Location"], y=GTlodsst[[x]][,2],t="l",col=rainbow(10)[x], lwd=1.8)
}
abline(h = Threshold , lty=2, lwd=1.9)
abline(h = SigThreshold , lty=1, lwd=1.9)
#text(x=SNPsinfo[1,3]+600000,y=Threshold+0.25, labels = "5% threshold")
#text(x=SNPsinfo[1,3]+600000,y=SigThreshold+0.25, labels = "1% threshold")
axis(1, at=seq(69000000,78000000,1000000), c("69","70","71","72","73","74","75","76","77","78"), las=1)

points(x=SNPsinfo[,3], y = rep(-0.3,length(SNPsinfo[,3])), pch=17)
#axis(1, at=SNPsinfo[,3], SNPsinfo[,"Markers"], cex.axis=0.4,las=2)

legend(69250000,-2.8, c("CW","HW","NW","WW","LNS"),lty=1,col=(rainbow(10)[1:5]),horiz = TRUE,xpd = TRUE, bty = "n", lwd=2.2)
legend(68100000,-3.5, c("LNSB","BNS","SubcAT","ViscAT","WAT"),lty=1,col=(rainbow(10)[6:10]),horiz = TRUE,xpd = TRUE, bty = "n", lwd=2.2)
#dev.off()

### Estimate the position of the CI

CIAll <- NULL
TopMarkers <- NULL 
for (x in 1:length(GTLod)){
  TopMarkers <- cbind(TopMarkers ,names(which.max(GTLod[[x]][,2])))
  CIAll <- cbind(CIAll, max(GTLod[[x]][,2])-1.5)
}
colnames(TopMarkers) <- names(GTLod)
colnames(CIAll) <- names(GTLod)

SNPsinfo <- read.table("RawData/SNPsinfo.txt",header=TRUE,sep="\t")
SNPPosition <- SNPsinfo[, c("Markers","Location")]
rownames(SNPPosition) <- SNPPosition[,"Markers"]

CIPosCal <- function(SNP1, SNP2, Trait){
  y1 <- as.numeric(GTLod[[Trait]][SNP1,2]) 
  y2 <- as.numeric(GTLod[[Trait]][SNP2,2]) 
  x1 <- as.numeric(SNPPosition[SNP1, "Location"])
  x2 <- as.numeric(SNPPosition[SNP2, "Location"])
  a <- as.numeric((y1-y2)/(x1-x2))
  b <- y1-a*(x1)
  CIPos <- (CIAll[,Trait] -b)/a
  return(CIPos)
}

CIPosCal("rs318175270", "rs14492508","BreastnoSkin")

par(mfrow=c(3,4))
for(x in 1:10){
  plot(x=c(0,10), y=c(0,10), t="n", ylab="LOD Score", xlab="", xaxt="n", main=names(GTlodsst)[x])
  points(GTlodsst[[x]][,2],t="l")
  abline(h = -log10(0.05/(length(sst)*length(SNPsForAnalysis))) , lty=1)
  abline(h = CIAll[,sst][x], lty=2)  
  axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
}
