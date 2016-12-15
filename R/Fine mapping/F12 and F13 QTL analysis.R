# Association Analysis using mixed-effects model for the F12 + F13 Chicken
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Dec, 2016
# first written Dec, 2016

setwd("D:/Chicken/Rdata/F12_and_F13_Combined_for_growth_traits")## Load F12 data
QTLdata <- read.table("RawData/F12 and F13 Growth traits with genotypiong.txt", header= TRUE, sep="\t", check.names=FALSE, na.strings="NA",
                       colClasses=c(rep("character",1), rep("factor",9), rep("character", 4),rep("numeric",21)))
QTLdata <- cbind(QTLdata,parents=paste0(QTLdata[,"Vater"],QTLdata[,"Mutter"]))

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

library(lme4)
GrowthTraits  <-  c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520") #names(QTLdata)[15:35]
SNPsForAnalysis <- c("rs16435551", "rs14490774", "rs318175270", "rs14492508")

GTLod <- vector("list", length(GrowthTraits)) 
names(GTLod) <- GrowthTraits
for (Phenotype in GrowthTraits){                                                     
  EachTraitLod <- NULL
  if(length(which(!is.na(QTLdata[,Phenotype]))) == 0){
        EachTraitLod <- rbind(EachTraitLod, NA)
        colnames(EachTraitLod) <- "SNP"
  }else{
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(QTLdata[,Phenotype]))
    model.full <- lmer(QTLdata[,Phenotype][idx] ~ QTLdata[,"Batch"][idx] + QTLdata[,"Futter"][idx] +(1|QTLdata[,"parents"][idx]) + QTLdata[,OneSNP][idx], REML=FALSE) # plot(fitted(model.full),residuals(model.full)) qqnorm(resid(model.full)) 
    #model.full <- lmer(QTLdata[,Phenotype][idx] ~ QTLdata[,"Batch"][idx] + (1|QTLdata[,"parents"][idx]) + QTLdata[,OneSNP][idx], REML=FALSE) # plot(fitted(model.full),residuals(model.full)) qqnorm(resid(model.full)) 
    #plot(fitted(model.full),residuals(model.full))
    model.null <- lmer(QTLdata[,Phenotype][idx] ~ QTLdata[,"Batch"][idx] + QTLdata[,"Futter"][idx] +(1|QTLdata[,"parents"][idx]), REML=FALSE)
    #model.null <- lmer(QTLdata[,Phenotype][idx] ~ QTLdata[,"Batch"][idx] + (1|QTLdata[,"parents"][idx]), REML=FALSE)

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

par(mfrow=c(4,5))
for(x in 2:21){
  plot(x=c(0,5), y=c(0,5), t="n", ylab="LOD Score", xlab="", xaxt="n", main=names(GTLod)[x])
  points(GTLod[[x]][,2],t="l")
  abline(h = Threshold , lty=1, col="blue")
  abline(h = SigThreshold , lty=1, col="red")
  axis(1, at=1:4, SNPsForAnalysis, las=2, cex.axis = 0.77)
}





