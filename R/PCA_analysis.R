# Principle component analysis (PCA)
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/FineMapping")
QTLdataAll <- read.table(file="Analysis/QTLdataAll.txt", sep = "\t",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("character", 9), rep("factor",6)))

factors <- QTLdataAll[,c("Generation","Month","Date","Season")]
factors <- apply(factors,2,function(x){as.numeric(as.factor(x))})
fac.pca <- princomp(factors, cor=TRUE)

summary(fac.pca,loadings=TRUE)
pre <- predict(fac.pca) 

cp1<- pre[,1]
cp2<- pre[,2]
cp3<- pre[,3]

Phenotype <- "Gew_20Wo"
OneSNP <-  "rs14490774"
library(lme4)
model.full <- lmer(QTLdataAll[,Phenotype] ~ cp1 + cp2+ cp3 + (1|QTLdataAll[,"Parents"]) + QTLdataAll[,OneSNP], REML=FALSE)
model.null <- lmer(QTLdataAll[,Phenotype] ~ cp1 + cp2+ cp3 + (1|QTLdataAll[,"Parents"]), REML=FALSE)
anova(model.null,model.full)