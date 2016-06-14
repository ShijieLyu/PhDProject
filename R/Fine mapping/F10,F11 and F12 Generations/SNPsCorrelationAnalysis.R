# SNP correlation anlaysis for 3 generations
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified June, 2016
# first written June, 2016

setwd("D:/Chicken/Rdata/FineMapping")
QTLdataAll <- read.table("Analysis/QTLdataAll.txt", sep = "\t",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("factor", 9), rep("factor",6)))

Generations <- c("F10","F11","F12")
SNPs <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

library("psych")
GenerationCorr <- vector("list", length(Generations)) 
names(GenerationCorr) <- Generations
for (GR in Generations){
  GRSNP <- QTLdataAll[which(QTLdataAll[,"Generation"] == GR),SNPs]
  SNPmatrix <- NULL
  for (x in 1:ncol(GRSNP)){
    NSNP <- as.numeric(GRSNP[,x])
    SNPmatrix <- cbind(SNPmatrix,NSNP)
  }
  colnames(SNPmatrix) <- colnames(GRSNP)
  GenerationCorr[[GR]] <- corr.test(SNPmatrix,method="pearson",use = "pairwise")[[1]] 
  png(paste("Analysis/SNPsCorrelation/",GR,"-SNPsCorrelation.png",sep=""))
  image(GenerationCorr[[GR]],breaks=c(0, 0.2, 0.4, 0.6, 0.8, 0.99, 1.0), col=c("white", "gold", "goldenrod", "orange", "red","red3"), xaxt='n', yaxt='n', main=paste("LD analysis for generation",GR))
  axis(1, at = seq(0, 1, 1/(length(colnames(SNPcor$r))-1)), colnames(SNPcor$r), las=2, cex.axis=0.7)
  axis(2, at = seq(0, 1, 1/(length(colnames(SNPcor$r))-1)), colnames(SNPcor$r), las=2, cex.axis=0.7)
  abline(h=0.0625 + seq(0, 0.875, 1/(length(colnames(SNPcor$r))-1)), col="white")
  abline(v=0.0625 + seq(0, 0.875, 1/(length(colnames(SNPcor$r))-1)), col="white")
  box()
  dev.off()
}



