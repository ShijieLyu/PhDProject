# Correct the phenotypic data for the chicken (All 3 Generations)
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/FineMapping")
QTLdataAll <- read.table("Analysis/QTLdataAll.txt", sep = "\t",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("character", 9), rep("factor",6)))

library(lme4)
GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

AllCorPheno <- NULL
for (Phenotype in GrowthTraits){
  onlyEnv <- lmer(QTLdataAll[,Phenotype] ~ QTLdataAll[,"Batch"] + (1|QTLdataAll[,"Parents"]), REML=FALSE) 
  Intercept <- coef(onlyEnv)$`QTLdataAll[, "Parents"]`["(Intercept)"]
  corPheno <- NULL
  for (x in 1:nrow(QTLdataAll)){
    corPhenoind <- round(resid(onlyEnv)[x] + Intercept[which(rownames(Intercept)==QTLdataAll[x,"Parents"]),])
    corPheno <- rbind(corPheno,corPhenoind)
  }
  AllCorPheno <- cbind(AllCorPheno, corPheno)
}
colnames(AllCorPheno) <- GrowthTraits
rownames(AllCorPheno) <- NULL
CorQTLdataAll <- cbind(QTLdataAll[,c("ID.Nr","Generation","Family","Schlupf","Parents","Batch")], AllCorPheno ,QTLdataAll[,SNPsForAnalysis])
write.table(CorQTLdataAll,"Analysis/MLM_CorPheno.txt",sep="\t",row.names = FALSE,quote = FALSE)

CorQTLdataAll <- read.table("Analysis/MLM_CorPheno.txt", sep = "\t",header=TRUE)
plot(CorQTLdataAll[,"Gew_15Wo"]~CorQTLdataAll[,"rs318175270"])