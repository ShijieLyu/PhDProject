# Association Analysis using mixed-effects model for the chicken (All 3 Generations)
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/FineMapping")
QTLdataAll <- read.table("Analysis/QTLdataAll.txt", sep = "\t",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("character", 9), rep("factor",6)))

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
GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

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
    model.full <- lmer(QTLdataAll[,Phenotype][idx] ~ QTLdataAll[,"Batch"][idx] + (1|QTLdataAll[,"Parents"][idx]) + QTLdataAll[,OneSNP][idx], REML=FALSE) # qqnorm(resid(model.full)) 
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
GTLod$Gew_15Wo

### Threshold under bonferroni correction
Threshold <- -log10(0.05/(length(GrowthTraits)*length(SNPsForAnalysis)))
SigThreshold <- -log10(0.01/(length(GrowthTraits)*length(SNPsForAnalysis)))

### Confidence Interval calculation -- Drop down 1.5-LOD of the Top marker 
CIAll <- NULL
TopMarkers <- NULL 
for (x in 1:length(GTLod)){
  TopMarkers <- cbind(TopMarkers ,names(which.max(GTLod[[x]][,2])))
  CIAll <- cbind(CIAll, max(GTLod[[x]][,2])-1.5)
}
colnames(TopMarkers) <- names(GTLod)
colnames(CIAll) <- names(GTLod)

### Plot the results
pdf("Analysis/MLMPlot.pdf")
par(mfrow=c(3,3))
for(x in 1:9){
  plot(x=c(0,10), y=c(0,10), t="n", ylab="LOD Score", xlab="", xaxt="n", main=names(GTLod)[x])
  points(GTLod[[x]][,2],t="l")
  abline(h = Threshold , lty=1)
  abline(h = CIAll[,x], lty=2)  
  axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
}
dev.off()

### Final Pic.
SNPsinfo <- read.table("RawData/SNPsinfo.txt",header=TRUE,sep="\t")

tiff("Analysis/MLMPlot-AllinOne-growth traits.tif", res = 300,width = 2100, height = 2000, compression = "lzw")
#pdf("Analysis/MLMPlot-AllinOne.pdf")
par(mai = c(2.2, 1, 1, 1))
plot(x=c(as.numeric(SNPsinfo[1,3]-100000),as.numeric(SNPsinfo[9,3]+100000)), y=c(0,10), t="n", ylab="LOD Score", xlab="Physical Position (Mb)", xaxt="n")
for(x in 1:9){
  points(x= SNPsinfo[,"Location"], y=GTLod[[x]][,2],t="l",col=rainbow(9)[x], lwd=1.8)
}
abline(h = Threshold , lty=2, lwd=1.9)
abline(h = SigThreshold , lty=1, lwd=1.9)
#text(x=SNPsinfo[1,3]+600000,y=Threshold+0.25, labels = "5% threshold")
#text(x=SNPsinfo[1,3]+600000,y=SigThreshold+0.25, labels = "1% threshold")
axis(1, at=seq(69000000,78000000,1000000), c("69","70","71","72","73","74","75","76","77","78"), las=1)

points(x=SNPsinfo[,3], y = rep(-0.3,length(SNPsinfo[,3])), pch=17)
#axis(1, at=SNPsinfo[,3], SNPsinfo[,"Markers"], cex.axis=0.4,las=2)

legend(68500000,-2.8, c("BW0","BW5","BW10","BW15","BW20"),lty=1,col=(rainbow(9)[1:5]),horiz = TRUE,xpd = TRUE, bty = "n", lwd=2.2)
legend(68300000,-3.5, c("BWG05","BWG510","BWG1015","BWG1520"),lty=1,col=(rainbow(9)[6:9]),horiz = TRUE,xpd = TRUE, bty = "n", lwd=2.2)
dev.off()


### Estimate the position of the CI

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

CIPosCal("rs16435551","rs14490774","Gew_10Wo")













