# Calculate the variance explaination for weight using F12 and F13
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Dec, 2016
# first written Dec, 2016

setwd("D:/Chicken/Rdata/F12_and_F13_Combined_for_growth_traits")## Load F12 data
PheoAndGeno <- read.table("RawData/F12 and F13 Growth traits with genotypiong.txt", header= TRUE, sep="\t", check.names=FALSE, na.strings="NA",
                       colClasses=c(rep("character",1), rep("factor",9), rep("character", 4),rep("numeric",21)))
PheoAndGeno <- cbind(PheoAndGeno,parents=paste0(PheoAndGeno[,"Vater"],PheoAndGeno[,"Mutter"]))

Generations <- c("F12","F13")
Diets <- c("1","2","normal")
Marker <- "rs14490774"
Traits <- names(PheoAndGeno)[15:35]

VarExpl <- vector("list", length(Diets)) 
names(VarExpl) <- Diets
for (Diet in Diets){
GeVar <- NULL
  for (Generation in Generations){
  Var <- NULL
    for (trait in Traits){
      SelectedInd <- PheoAndGeno[which(PheoAndGeno[,"Generation"]== Generation & PheoAndGeno[,"Futter"]== Diet),]
      model <- anova(lm(SelectedInd[,trait]~SelectedInd[,"parents"] + SelectedInd[,Marker]))
      EveryVar <- model["SelectedInd[, Marker]","Sum Sq"]/sum(model[,"Sum Sq"]) * 100
      Var <- cbind(Var, EveryVar)
    }
  GeVar <- rbind(GeVar, Var)
  }
colnames(GeVar) <- Traits
rownames(GeVar) <- Generations  
VarExpl[[Diet]] <- GeVar
}

## Only normal diet for F13 using linear model
F13Var <- NULL
for (trait in Traits){
  SelectedInd <- PheoAndGeno[which(PheoAndGeno[,"Generation"]== "F13" & PheoAndGeno[,"Futter"]== "normal"),]
  model <- anova(lm(SelectedInd[,trait]~ SelectedInd[,"Batch"] + SelectedInd[,"parents"] + SelectedInd[,Marker]))
  EveryVar <- model["SelectedInd[, Marker]","Sum Sq"]/sum(model[,"Sum Sq"]) * 100
  F13Var <- cbind(F13Var, EveryVar)
}
colnames(F13Var) <- Traits

## Only normal diet for F13 using linear mixed model
SNPvarPerc <- function(mF){
  X <- getME(mF,"X")                                               # Get the fixed-effects model matrix
  TotalFixed <- 0
  for (x in 2:length(fixef(mF))){  
    TotalFixed <- TotalFixed + (fixef(mF)[x] * X[, x])  
  }
  SNPFixed <- 0
  for (x in grep("Marker", names(fixef(mF)))){     
    SNPFixed <- SNPFixed + (fixef(mF)[x] * X[, x])
  }
  FixedTotalVar <- var(TotalFixed)                                 # Get the total fixed effects variance 
  FixedSNPVar <- var(SNPFixed)                                     # Get the specific fixed effect variance
  TotalVar <- FixedTotalVar + sum(data.frame(VarCorr(mF))[,4])     # Calculate the total variance with fixed and random
  Perc <- round(FixedSNPVar/TotalVar,3)*100
  return(Perc)
}

Marker <- c()
library(lme4)
F13Var <- NULL
for (trait in Traits){
  SelectedInd <- PheoAndGeno[which(PheoAndGeno[,"Generation"]== "F13" & PheoAndGeno[,"Futter"]== "normal"),]
  model.full <- lmer(SelectedInd[,trait]~ SelectedInd[,"Batch"] + (1|SelectedInd[,"parents"]) + SelectedInd[,Marker])
  EveryVar <- SNPvarPerc(model.full)
  F13Var <- cbind(F13Var, EveryVar)
}
colnames(F13Var) <- Traits






