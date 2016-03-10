# Chicken F13 association analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Shijie Lyu
# last modified Feb, 2016
# first written Feb, 2016

setwd("D:/Chicken/Rdata/Feed_Experiment_on_F13")
F13genotypes <- read.table("RawData/F13_genotypes.txt", header= TRUE, sep="\t", check.names=FALSE, colClasses=c("character",rep("factor",4)))   # Load Genotypes
F13phenotypes <- read.table("RawData/F13_phenotypes.txt", header = TRUE, sep="\t",check.names=FALSE, na.strings = "", colClasses= c("character", rep("factor",7), rep("numeric",33)))         # Load Phenotypes

### Association analysis on 13 weeks of age for different diet raisen in single cage
FullData <- cbind(F13phenotypes,F13genotypes[,-1],parents=paste0(F13phenotypes[,"Vater"],F13phenotypes[,"Mutter"]))
F13pheno13w <- FullData[which(FullData[,"Week"]== "13w"),]             # Select the chicken slaughterd at 13 weeks

Normal <- FullData[which(FullData[,"Diet"]=="Standard"),]
Diet0 <- F13pheno13w[which(F13pheno13w[,"Diet"]== 0),]
Diet1 <- F13pheno13w[which(F13pheno13w[,"Diet"]== 1),]
Diet2 <- F13pheno13w[which(F13pheno13w[,"Diet"]== 2),]

GrowthTraits  <- colnames(F13pheno13w)[9:22]
par(mfrow=c(3,5))
for (trait in GrowthTraits){
  boxplot(Normal[,trait],Diet0[,trait],Diet1[,trait],Diet2[,trait], main = trait) 
  axis(1, at=1:4, c("Normal","Diet0","Diet1","Diet2"), las=2, cex.axis = 0.77)
}

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

### Different Diet

library(lme4)
GrowthTraits  <- colnames(F13pheno13w)[9:22]
SNPsForAnalysis <- c("rs16435551", "rs14490774", "rs318175270", "rs14492508")

par(mfrow=c(4,8))
GTLod <- vector("list", length(GrowthTraits)) 
names(GTLod) <- GrowthTraits
for (Phenotype in GrowthTraits){                                                     
  EachTraitLod <- NULL
  if(length(which(!is.na(F13pheno13w[,Phenotype]))) == 0){
        EachTraitLod <- rbind(EachTraitLod, NA)
        colnames(EachTraitLod) <- "SNP"
  }else{
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(F13pheno13w[,Phenotype]))
    model.full <- lmer(F13pheno13w[,Phenotype][idx] ~ F13pheno13w[,"SD"][idx] + F13pheno13w[,"Family"][idx] + F13pheno13w[,"Diet"][idx] + (1|F13pheno13w[,"parents"][idx]) + F13pheno13w[,OneSNP][idx], REML=FALSE) # plot(fitted(model.full),residuals(model.full)) qqnorm(resid(model.full)) 
    plot(fitted(model.full),residuals(model.full))
    model.null <- lmer(F13pheno13w[,Phenotype][idx] ~ F13pheno13w[,"SD"][idx] + F13pheno13w[,"Family"][idx] + F13pheno13w[,"Diet"][idx] + (1|F13pheno13w[,"parents"][idx]), REML=FALSE)
    res <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[8]]),SNPvar))
    colnames(EachTraitLod) <- c("Residuals","SNP","SNPvar")
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  GTLod[[Phenotype]] <- EachTraitLod
  }
}

Threshold <- -log10(0.05/((length(GrowthTraits))*(length(SNPsForAnalysis))))
SigThreshold <- -log10(0.01/((length(GrowthTraits))*(length(SNPsForAnalysis))))

par(mfrow=c(3,5))
for(x in 1:14){
  plot(x=c(0,5), y=c(0,5), t="n", ylab="LOD Score", xlab="", xaxt="n", main=names(GTLod)[x])
  points(GTLod[[x]][,2],t="l")
  abline(h = Threshold , lty=1, col="blue")
  abline(h = SigThreshold , lty=1, col="red")
  axis(1, at=1:4, SNPsForAnalysis, las=2, cex.axis = 0.77)
}

## Corrected the phenotypes using the Slaughtr Data, Family, Diet and parents

library(lme4)
GrowthTraits  <- colnames(F13pheno13w)[9:22]
SNPsForAnalysis <- c("rs16435551", "rs14490774", "rs318175270", "rs14492508")

AllCorPheno <- NULL
AllCorPheno <- vector("list", length(GrowthTraits)) 
names(AllCorPheno) <- GrowthTraits

for (Phenotype in GrowthTraits){
  idx <- which(!is.na(F13pheno13w[,Phenotype]))
  onlyEnv <- lmer(F13pheno13w[,Phenotype][idx] ~ F13pheno13w[,"SD"][idx] + F13pheno13w[,"Family"][idx] + F13pheno13w[,"Diet"][idx] + (1|F13pheno13w[,"parents"][idx]), REML=FALSE) 
  Intercept <- coef(onlyEnv)$`F13pheno13w[, "parents"][idx]`["(Intercept)"]
  corPheno <- NULL
  for (x in 1:length(idx)){
    corPhenoind <- round(resid(onlyEnv)[x] + Intercept[which(rownames(Intercept)==F13pheno13w[idx[x],"parents"]),])
    corPheno <- rbind(corPheno,corPhenoind)
  }
  rownames(corPheno) <- F13pheno13w[idx,"parents"]
  colnames(corPheno) <- paste0("Corrected",Phenotype)
  corPheno <- cbind(ID = F13pheno13w[idx,"ID"], parents = F13pheno13w[idx,"parents"],corPheno)
  AllCorPheno[[Phenotype]] <- corPheno
}

CorPheMa <- matrix(NA, nrow(F13pheno13w), length(GrowthTraits), dimnames = list(F13pheno13w[,"ID"], paste0(GrowthTraits)))
for(x in 1:length(AllCorPheno)){
  colname <- names(AllCorPheno)[x]
  for(y in 1:nrow(AllCorPheno[[x]])){
    CorPheMa[AllCorPheno[[x]][y,1], colname] <- AllCorPheno[[x]][y,3]
  }
}
colnames(CorPheMa) <- paste0("Corr-",GrowthTraits)
write.table(CorPheMa,"Analysis/Corrected_F13_GrowthTraits.txt",row.names = TRUE, col.names=TRUE,sep = "\t")

### Merge the F10, F11 and F12 with the F13 -- Only the 20week with normal diet

QTLdataAll <- read.table("RawData/F10toF13-QTLdataAll.txt", sep = "\t",na.strings = "",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("character", 9), rep("factor",6)))

library(lme4)
GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
SNPsForAnalysis <- c("rs16435551", "rs14490774", "rs318175270", "rs14492508")

GTLod <- vector("list", length(GrowthTraits)) 
names(GTLod) <- GrowthTraits
for (Phenotype in GrowthTraits){                                                     
  EachTraitLod <- NULL
  if(length(which(!is.na(QTLdataAll[,Phenotype]))) == 0){
        EachTraitLod <- rbind(EachTraitLod, NA)
        colnames(EachTraitLod) <- "SNP"
  }else{
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(QTLdataAll[,Phenotype]))
    model.full <- lmer(QTLdataAll[,Phenotype][idx] ~ QTLdataAll[,"Batch"][idx] + (1|QTLdataAll[,"Parents"][idx]) + QTLdataAll[,OneSNP][idx], REML=FALSE) # plot(fitted(model.full),residuals(model.full)) qqnorm(resid(model.full)) 
    model.null <- lmer(QTLdataAll[,Phenotype][idx] ~ QTLdataAll[,"Batch"][idx] + (1|QTLdataAll[,"Parents"][idx]), REML=FALSE)
    res <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[8]]),SNPvar))
    colnames(EachTraitLod) <- c("Residuals","SNP","SNPvar")
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  GTLod[[Phenotype]] <- EachTraitLod
  }
}

Threshold <- -log10(0.05/((length(GrowthTraits))*(length(SNPsForAnalysis))))
SigThreshold <- -log10(0.01/((length(GrowthTraits))*(length(SNPsForAnalysis))))

par(mfrow=c(3,3))
for(x in 1:9){
  plot(x=c(0,5), y=c(0,10), t="n", ylab="LOD Score", xlab="", xaxt="n", main=names(GTLod)[x])
  points(GTLod[[x]][,2],t="l")
  abline(h = Threshold , lty=1, col="blue")
  abline(h = SigThreshold , lty=1, col="red")
  axis(1, at=1:4, SNPsForAnalysis, las=2, cex.axis = 0.77)
}




