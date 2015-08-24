# Association Analysis using mixed-effects model for the F12 Chicken
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/FineMapping")
QTLdataAll <- read.table("Analysis/QTLdataAll.txt", sep = "\t",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("character", 9), rep("factor",6)))
F12QTLdata <- read.table("Analysis/F12_QTLdataWithBWG-Geno NA fill in.txt", na.strings=c("NA","?"),header=TRUE, sep="\t")

F12ST <- F12QTLdata[which(F12QTLdata[,"SchlachtWo"]=="20W" & F12QTLdata[,"Futter"] == "normal"),]     # F12 slaughtered traits measured at 20th week, and the chicken were all with normal diet
F12ST <- cbind(F12ST, Batch = QTLdataAll[which(QTLdataAll[,"ID.Nr"] %in% as.character(F12ST[,"ID.Nr"])),"Batch"], Parents = QTLdataAll[which(QTLdataAll[,"ID.Nr"] %in% as.character(F12ST[,"ID.Nr"])),"Parents"])

GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")

F12STN  <- names(F12ST)[-c(99,104)]                            # Remove "Stander._links_.Durchmesser" and "Stander._rechts_.Durchmesser_", since miss to many chickens                 
F12STN  <- c(GrowthTraits, F12STN[c(37:39,41:143)])            # The names of slaughtered traits which for the association analysis

### Association Analysis using the mixed linear model
library(lme4)
Traits  <- F12STN
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

GTLod <- vector("list", length(Traits)) 
names(GTLod) <- Traits

for (Phenotype in Traits){                                                     
  EachTraitLod <- NULL
  if(length(which(!is.na(F12ST[,Phenotype]))) == 0){
        EachTraitLod <- rbind(EachTraitLod, NA)
        colnames(EachTraitLod) <- "SNP"
  }else{
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(F12ST[,Phenotype]))
    model.full <- lmer((as.numeric(F12ST[,Phenotype])[idx]) ~ (as.factor(F12ST[,"Batch"])[idx]) + (1|(as.factor(F12ST[,"Parents"])[idx])) + (as.character(F12ST[,OneSNP])[idx]), REML=FALSE) # qqnorm(resid(model.full)) 
    model.null <- lmer((as.numeric(F12ST[,Phenotype])[idx]) ~ (as.factor(F12ST[,"Batch"])[idx]) + (1|(as.factor(F12ST[,"Parents"])[idx])), REML=FALSE)
    res <- anova(model.null,model.full)
    EachTraitLod <- rbind(EachTraitLod, -log10(res[[8]]))
    colnames(EachTraitLod) <- c("Residuals","SNP")  
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  GTLod[[Phenotype]] <- EachTraitLod
  }
}

GTLod$Gew_15Wo

#library(marray)
#write.list(GTLod,"Analysis/F12SlaughterTraitsLod.txt")

if(!file.exists("Analysis/F12SlaughterTraitsSignificant.txt")){
  for(x in 1:length(GTLod)){
    SigSNP <- which(GTLod[[x]][,"SNP"] > -log10(0.05/9))                                                      # Analyse the resulting profiles, and look at which markers are above the threshold
    if (length(SigSNP) == 0) {
    cat(names(GTLod)[x], "\n", file = "Analysis/F12SlaughterTraitsSignificant.txt", append=TRUE)
    }else{
    cat(names(GTLod)[x], names(SigSNP), "\n", file = "Analysis/F12SlaughterTraitsSignificant.txt", append=TRUE)
    }                                                             
  }
}else{
  cat("Loading Slaughted Traits with Significant SNPs from disk\n")
  SigSNP <- read.table("Analysis/F12SlaughterTraitsSignificant.txt",sep="\t", header=FALSE)
}

SigTraits <- NULL                                                                                      
for (x in 1:nrow(SigSNP)){
  if(length(unlist(strsplit(as.character(SigSNP[x,1]), " "))) > 1) {
  SigTraits <- c(SigTraits, strsplit(as.character(SigSNP[x,1]), " ")[[1]][1])
  }
}

pdf("Analysis/F12SlaughterTraitsPlot.pdf",paper="USr", width=10, height=8)
  par(mfrow=c(4,5))
  for(SigTrait in SigTraits){
    plot(x=c(0,10), y=c(0,5), t="n", ylab="LOD Score", xlab="", xaxt="n", main= SigTrait, cex.main=0.9)
    points(GTLod[[SigTrait]][,2],t="l")
    abline(h=-log10(0.05/9), lty=2)    
    axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
  }
  title("The association analysis for the F12 Slaughter Traits", outer=TRUE)
dev.off()

### Export all the traits with LOD score-- Generate on total plot every 20 traits

par(mfrow=c(4,5))
  for(x in 1:20){
    plot(x=c(0,10), y=c(0,5), t="n", ylab="LOD Score", xlab="", xaxt="n", main= names(GTLod[x]), cex.main=0.9)
    points(GTLod[[x]][,2],t="l")
    abline(h=-log10(0.05/9), lty=2)    
    axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
  }
title("The association analysis for the F12 Slaughter Traits", outer=TRUE)







