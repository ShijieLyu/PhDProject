# Association Analysis using mixed-effects model for the chicken (All 3 Generations)
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/GGR27-fine-mapping") 
QTLdataAll <- read.table("RawData/GGR27-ADL0376-QTLdata.txt", sep = "\t",header=TRUE,na.strings=c("NA","?"," ",""), colClasses=c(rep("character",1), rep("factor",13),rep("character",1), rep("numeric",77)))

QTLdataAll <- cbind(QTLdataAll, LegnoSkin = rep(NA, nrow(QTLdataAll)), LegSkin = rep(NA, nrow(QTLdataAll)), LegnoSkin_bone= rep(NA, nrow(QTLdataAll)),
                    BreastnoSkin = rep(NA, nrow(QTLdataAll)), BreastSkin = rep(NA, nrow(QTLdataAll)), VisceralFat = rep(NA, nrow(QTLdataAll)), TotalFat = rep(NA, nrow(QTLdataAll)),
                    TransNeckFat = rep(NA, nrow(QTLdataAll)),TransVisceralFat = rep(NA, nrow(QTLdataAll)),TransTotalFat = rep(NA, nrow(QTLdataAll)),LegBonemass = rep(NA, nrow(QTLdataAll))
                    )

QTLdataAll[,"LegnoSkin"] <- QTLdataAll[,"Keule.re..ohne.Haut"] + QTLdataAll[,"Keule.li..ohne.Haut"]
QTLdataAll[,"LegSkin"] <- QTLdataAll[,"Keule.Haut.re"] + QTLdataAll[,"Keule.Haut.li."]
QTLdataAll[,"LegnoSkin_bone"] <- QTLdataAll[,"Keule.re..ohne.Haut.und.ohne.Knochen"] + QTLdataAll[,"Keule.li..ohne.Haut.u.ohne.Knochen"]
QTLdataAll[,"BreastnoSkin"] <- QTLdataAll[,"Brust.re..ohne.Haut"] + QTLdataAll[,"Brust.li..ohne.Haut"]
QTLdataAll[,"BreastSkin"] <- QTLdataAll[,"Brust.Haut.re."] + QTLdataAll[,"Brust.Haut.li."]
QTLdataAll[,"VisceralFat"] <- QTLdataAll[,"Fett.Herz."]+QTLdataAll[,"Fett.Magen."] + QTLdataAll[,"Fett.Leber."] + QTLdataAll[,"Fett.Milz."] + QTLdataAll[,"Abdominalfett"]
QTLdataAll[,"TotalFat"] <- QTLdataAll[,"Fett.Herz."]+QTLdataAll[,"Fett.Magen."] + QTLdataAll[,"Fett.Leber."] + QTLdataAll[,"Fett.Milz."] + QTLdataAll[,"Abdominalfett"] + QTLdataAll[,"Hals.Fett"]
QTLdataAll[,"TransNeckFat"] <- (QTLdataAll[,"Hals.Fett"])^0.5
QTLdataAll[,"TransVisceralFat"] <- (QTLdataAll[,"VisceralFat"])^0.5
QTLdataAll[,"TransTotalFat"] <- (QTLdataAll[,"TotalFat"])^0.5
QTLdataAll[,"LegBonemass"] <- QTLdataAll[,"LegnoSkin"] - QTLdataAll[,"LegnoSkin_bone"]


# The function for calculating the 
SNPvarPerc <- function(mF){
  X <- getME(mF,"X")                                               # Get the fixed-effects model matrix
  TotalFixed <- 0
  for (x in 2:length(fixef(mF))){  
    TotalFixed <- TotalFixed + (fixef(mF)[x] * X[, x])  
  }
  SNPFixed <- 0
  for (x in grep("OneSNP", names(fixef(mF)))){     
    SNPFixed <- SNPFixed + (fixef(mF)[x] * X[, x])
  }
  FixedTotalVar <- var(TotalFixed)                                 # Get the total fixed effects variance 
  FixedSNPVar <- var(SNPFixed)                                     # Get the specific fixed effect variance
  TotalVar <- FixedTotalVar + sum(data.frame(VarCorr(mF))[,4])     # Calculate the total variance with fixed and random
  Perc <- round(FixedSNPVar/TotalVar,3)*100
  return(Perc)
}

### Association Analysis using the mixed linear model
library(lme4)
Traits  <- names(QTLdataAll)[16:ncol(QTLdataAll)]
SNPsForAnalysis <- "ADL0376"

GTLod <- vector("list", length(Traits)) 
names(GTLod) <- Traits
for (Phenotype in Traits){                                                     
  EachTraitLod <- NULL
  if(length(which(!is.na(QTLdataAll[,Phenotype]))) == 0){
        EachTraitLod <- rbind(EachTraitLod, NA)
        colnames(EachTraitLod) <- "SNP"
  }else{
  OneSNP <- SNPsForAnalysis
  idx <- which(!is.na(QTLdataAll[,Phenotype]))
  model.full <- lmer(QTLdataAll[,Phenotype][idx] ~ QTLdataAll[,"Batch"][idx] + (1|QTLdataAll[,"Parents"][idx]) + QTLdataAll[,OneSNP][idx], REML=FALSE) # plot(fitted(model.full),residuals(model.full)) qqnorm(resid(model.full)) 
  model.null <- lmer(QTLdataAll[,Phenotype][idx] ~ QTLdataAll[,"Batch"][idx] + (1|QTLdataAll[,"Parents"][idx]), REML=FALSE)
  #model.full <- lmer((as.numeric(QTLdataAll[,Phenotype])[idx]) ~ (as.factor(QTLdataAll[,"Batch"])[idx]) + (1|(as.factor(QTLdataAll[,"Parents"])[idx])) + (as.numeric(QTLdataAll[,"BW.nuchtern"])[idx]) + (as.character(QTLdataAll[,OneSNP])[idx]), REML=FALSE) # qqnorm(resid(model.full)) 
  #model.null <- lmer((as.numeric(QTLdataAll[,Phenotype])[idx]) ~ (as.factor(QTLdataAll[,"Batch"])[idx]) + (1|(as.factor(QTLdataAll[,"Parents"])[idx])) + (as.numeric(QTLdataAll[,"BW.nuchtern"])[idx]), REML=FALSE)
  res <- anova(model.null,model.full)
  SNPvar <- SNPvarPerc(model.full)
  EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[8]]),SNPvar))
  colnames(EachTraitLod) <- c("Residuals","SNP","SNPvar")
  rownames(EachTraitLod) <- SNPsForAnalysis
  GTLod[[Phenotype]] <- EachTraitLod
  }
}
GTLod

### Selecte the traits which above the threshold -log10(0.05)
SigTraits <- NULL
for (x in 1:length(GTLod)){
  if (GTLod[[x]][,"SNP"] > -log10(0.05)){
    rownames(GTLod[[x]]) <- names(GTLod[x])
    SigTraits <- rbind(SigTraits,GTLod[[x]])           #paste(names(GTLod[x]),GTLod[[x]]))
  }
}

png("Analysis/SigTraits.png")
par(mar=c(10,4,4,2))
plot(x=c(0,15), y=c(0,10), t="n", ylab="LOD Score", xlab="", xaxt="n")
xcoor <- 1
for (trait in rownames(SigTraits)){
  points(SigTraits[trait,2],x=xcoor)
  xcoor <- xcoor + 1
}
abline(h = -log10(0.05) , lty=1, col="blue")
abline(h = -log10(0.01) , lty=1, col="red")
TraitNames <- c("BWG1015", "BWG1520", "BW nuchtern","Wings weight","Neck weight","Heart weight","Fat weight on stomach","Spleen weight","left leg muscle pH after 1h","left leg muscle weight after 48h","left shank length","left shank weight","right shank length","right shank weight","Leg bone mass")
axis(1, at=1:15, TraitNames, las=2, cex.axis = 0.77)
dev.off()


ShankTraits <- c("StÃ.nder..links..LÃ.nge.","StÃ.nder..links..Gewicht..","StÃ.nder..rechts..LÃ.nge.","StÃ.nder..rechts..Gewicht..")
png("Analysis/RightShank.png")
ShankMatrix <- QTLdataAll[,c(ShankTraits,"ADL0376")]
par(mfrow=c(2,1))
par(mar=c(2,4,4,2))
#boxplot(ShankMatrix[,"StÃ.nder..links..LÃ.nge."]~ShankMatrix[,"ADL0376"],main="Left Shank Length")
#boxplot(ShankMatrix[,"StÃ.nder..links..Gewicht.."]~ShankMatrix[,"ADL0376"],main="Left Shank Weight")
boxplot(ShankMatrix[,"StÃ.nder..rechts..LÃ.nge."]~ShankMatrix[,"ADL0376"],main="Right Shank Length")
boxplot(ShankMatrix[,"StÃ.nder..rechts..Gewicht.."]~ShankMatrix[,"ADL0376"],main="Right Shank Weight")
dev.off()


GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
sst <- c("BW.bratfertig.","Kopf","Hals","Flugel","LegBonemass","LegnoSkin_bone","BreastnoSkin","Hals.Fett","VisceralFat","TotalFat")

Threshold <- -log10(0.05/((length(GrowthTraits)+length(sst))))
SigThreshold <- -log10(0.01/((length(GrowthTraits)+length(sst))))

png("Analysis/GrowthTraits-2.png")
plot(x=c(0,10), y=c(0,10), t="n", ylab="LOD Score", xlab="", xaxt="n")
xcoor <- 1
for (trait in GrowthTraits){
  points(GTLod[[trait]][,2],x=xcoor)
  xcoor <- xcoor + 1
}
abline(h = -log10(0.05) , lty=1, col="blue")
abline(h = -log10(0.01) , lty=1, col="red")
axis(1, at=1:9, GrowthTraits, las=2, cex.axis = 0.77)
dev.off()

png("Analysis/SlaughterTraits-BW as covariate.png")
plot(x=c(0,10), y=c(0,10), t="n", ylab="LOD Score", xlab="", xaxt="n")
xcoor <- 1
for (trait in sst){
  points(GTLod[[trait]][,2],x=xcoor)
  xcoor <- xcoor + 1
}
abline(h = Threshold , lty=1, col="blue")
abline(h = SigThreshold , lty=1, col="red")
axis(1, at=1:10, sst, las=2, cex.axis = 0.77)
dev.off()