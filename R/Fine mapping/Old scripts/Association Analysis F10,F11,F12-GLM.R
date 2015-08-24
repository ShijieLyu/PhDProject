# Association Analysis for the chicken using F10, F11 and F12 Generation
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Apr, 2015
# first written Apr, 2015

setwd("D:/Chicken/Rdata/FineMapping")
F10genotypes <- read.table("RawData/F10_Genotypes.txt", header=TRUE, sep="\t")
F10BodyWeight <- read.table("RawData/F10_BodyWeight_Phenotypes.txt", header=TRUE, sep="\t")
F11genotypes <- read.table("RawData/F11_Genotypes.txt", na.strings="-", header=TRUE, sep="\t")
F11BodyWeight <- read.table("RawData/F11_BodyWeight_Phenotypes.txt", na.strings="",header=TRUE, sep="\t")
F12Traits <- read.table("Analysis/F12_QTLdataWithBWG-Geno NA fill in.txt", na.strings=c("NA","?"),header=TRUE, sep="\t")

SNPnames <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

F11BodyWeight <- cbind(F11BodyWeight, BWG05 = rep(NA, nrow(F11BodyWeight)),BWG510 = rep(NA, nrow(F11BodyWeight)),BWG1015 = rep(NA, nrow(F11BodyWeight)),BWG1520 = rep(NA, nrow(F11BodyWeight)))
for (x in 1:nrow(F11BodyWeight)){
  F11BodyWeight[x,"BWG05"]   <- F11BodyWeight[x,"Gew_5Wo"] - F11BodyWeight[x,"Gew_1d"]
  F11BodyWeight[x,"BWG510"]  <- F11BodyWeight[x,"Gew_10Wo"] - F11BodyWeight[x,"Gew_5Wo"]
  F11BodyWeight[x,"BWG1015"] <- F11BodyWeight[x,"Gew_15Wo"] - F11BodyWeight[x,"Gew_10Wo"]
  F11BodyWeight[x,"BWG1520"] <- F11BodyWeight[x,"Gew_20Wo"] - F11BodyWeight[x,"Gew_15Wo"]
}

F11QTLdata <- cbind(F11BodyWeight,F11genotypes[,SNPnames])

F10BodyWeight <- cbind(F10BodyWeight, BWG05 = rep(NA, nrow(F10BodyWeight)),BWG510 = rep(NA, nrow(F10BodyWeight)),BWG1015 = rep(NA, nrow(F10BodyWeight)),BWG1520 = rep(NA, nrow(F10BodyWeight)))
for (x in 1:nrow(F10BodyWeight)){
  F10BodyWeight[x,"BWG05"]   <- F10BodyWeight[x,"Gew_5Wo"] - F10BodyWeight[x,"Gew_1d"]
  F10BodyWeight[x,"BWG510"]  <- F10BodyWeight[x,"Gew_10Wo"] - F10BodyWeight[x,"Gew_5Wo"]
  F10BodyWeight[x,"BWG1015"] <- F10BodyWeight[x,"Gew_15Wo"] - F10BodyWeight[x,"Gew_10Wo"]
  F10BodyWeight[x,"BWG1520"] <- F10BodyWeight[x,"Gew_20Wo"] - F10BodyWeight[x,"Gew_15Wo"]
}

F10QTLdata <- cbind(F10BodyWeight,F10genotypes[,SNPnames])

### 3-Select the normal diet and individuals with the BodyWeight data of every 5 weeks, should be 46 chickens
F12QTLdata <- F12Traits[,c(names(F11BodyWeight),SNPnames)]
F12QTLdata <- F12QTLdata[which(!is.na(F12QTLdata[,"Gew_20Wo"]) & F12QTLdata[,"Futter"] == "normal"),] 

QTLdataAll <- rbind(F12QTLdata,F11QTLdata,F10QTLdata) 

Parents <- NULL
for (x in 1:nrow(QTLdataAll)){
  parent <- paste0(QTLdataAll[x,"Vater"],QTLdataAll[x,"Mutter"])
  Parents <- c(Parents,parent) 
}
QTLdataAll <- cbind(QTLdataAll,Parents)

getSeason <- function(DATES) {
  mmonths <- as.numeric(unlist(lapply(strsplit(as.character(DATES), ".", fixed=TRUE), "[", 2)))
  ret <- rep(NA, length(mmonths))
  ret[mmonths >= 3 & mmonths <= 5]                  <- "Spring"
  ret[mmonths >= 6 & mmonths <= 8]                  <- "Summer"
  ret[mmonths >= 9 & mmonths <= 11]                 <- "Fall"
  ret[mmonths == 12 | mmonths == 1 | mmonths == 2]  <- "Winter"
  return(ret)
}

getMonth <- function(DATES) {
  mmonths <- as.numeric(unlist(lapply(strsplit(as.character(DATES), ".", fixed=TRUE), "[", 2)))
  return(mmonths)
}

getDate <- function(DATES) {
  Date <- rep(NA, length(DATES))
  for(x in 1:length(DATES)){
  Date[x] <- paste0(strsplit(as.character(DATES[x]), ".", fixed=TRUE)[[1]][c(1,2)], collapse = "/")
  } 
  return(Date)
}

QTLdataAll <- cbind(QTLdataAll, Batch = paste0(QTLdataAll[,"Generation"],QTLdataAll[,"Schlupf"]), GS = paste0(QTLdataAll[,"Generation"],QTLdataAll[,"Schlupf"]), 
                    Season=getSeason(QTLdataAll[,"Schlupf"]), Month=getMonth(QTLdataAll[,"Schlupf"]), Date=getDate(QTLdataAll[,"Schlupf"]))

write.table(QTLdataAll, file="Analysis/QTLdataAll.txt", sep = "\t", row.names=FALSE,quote=FALSE)


# Association Analysis

GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

GTLod <- vector("list", length(GrowthTraits)) 
names(GTLod) <- GrowthTraits

for (Phenotype in GrowthTraits){                                                       # Model after selection
  EachTraitLod <- NULL
  for (OneSNP in SNPsForAnalysis){
    res <- anova(lm(as.numeric(QTLdataAll[,Phenotype]) ~ QTLdataAll[,"Generation"] + as.factor(QTLdataAll[,"Family"]) + QTLdataAll[,"Schlupf"] + QTLdataAll[,OneSNP]))
      Pfactors <- NULL
      for (x in 1:5){
        Pfactors <- c(Pfactors, round((sum(res[x,2])/sum(res[2]))*100, digits = 2))    # percentage(contribution) per factor
      }
    EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[5]]),Pfactors))
    colnames(EachTraitLod) <- c("Generation","Family","Schlupf","SNP", "Residuals","Generation%","Family%","Schlupf%","SNP%", "Residuals%")  
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  GTLod[[Phenotype]] <- EachTraitLod
}
GTLod$Gew_15Wo

## 1000 times Permutation, Got threschlod

GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
Covs  <- c("Generation","Family","Schlupf")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

genotypes    <- QTLdataAll[,SNPsForAnalysis]                                 # Variables
phenotypes   <- QTLdataAll[,GrowthTraits]
covariates   <- QTLdataAll[,Covs]

mapQTLs <- function(genotypes, phenotypes, covariates){                      # Function to map QTLs
  LODS <- matrix(NA, ncol(genotypes), ncol(phenotypes),dimnames=list(colnames(phenotypes), colnames(genotypes)))     # Matrix LOD scores of markers x phenotypes
  for(x in 1:ncol(phenotypes)){  
    for(y in 1:ncol(genotypes)){ 
      res <- anova(lm(as.numeric(phenotypes[,x]) ~ covariates[,"Generation"] + covariates[,"Family"] + covariates[,"Schlupf"] + genotypes[,y]))
      LODS[x, y] <- -log10(res[[5]][length(res[[5]]) - 1])
    }
  }
  return(LODS)
}

LodForAll <- mapQTLs(genotypes, phenotypes, covariates)                       # Map the real data

if(!file.exists("Analysis/Permutation-1000-withF10.txt")){                   # Permutation for getting an overall 95 % confidence interval (corrected for N markers and N phenotypes
  maxes <- NULL
  for(p in 1:1000){
     ngenotypes <- genotypes[sample(nrow(genotypes)), ]
     res <- mapQTLs(ngenotypes, phenotypes, covariates)
     maxes <- c(maxes, max(res))
  }
     write.table(maxes,"Analysis/Permutation-1000-withF10.txt",sep="\t")
  }else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading Permutation information from disk\n")
  maxes <- read.table("Analysis/Permutation-1000-withF10.txt",sep="\t", header=TRUE)
}

Threshold <- sort(maxes[,1])[1000 * 0.95]       # Cut-off for significance 95 % -- Threshold = 3.111283
cat("Threshold =", Threshold, "\n")

## Calculate the 95% Confidence Interval for the top marker of every traits

if(!file.exists("Analysis/CIboots1000-withF10.txt")){
  boots <- NULL
  for(x in 1:1000){
    idx <- sample(nrow(genotypes), nrow(QTLdataAll), replace=TRUE)
    boots <- rbind(boots, mapQTLs(genotypes[idx, ], phenotypes[idx,], covariates[idx,]))
  }
    boots <- cbind(Traits = rownames(boots), boots)
    write.table(boots,"Analysis/CIboots1000-withF10.txt",sep="\t")
  }else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading boots information from disk\n")
  boots <- read.table("Analysis/CIboots1000-withF10.txt",sep="\t", header=TRUE, row.names=NULL)
}

CIAll <- NULL
for (trait in names(phenotypes)){
  TraitBoots <- boots[which(boots[,"Traits"] == trait),]
  topM <- names(which.max(LodForAll[trait,]))
  BootI <- sort(TraitBoots[,topM])[1000 * 0.025]
  CIAll <- cbind(CIAll, BootI)
}
colnames(CIAll) <- names(phenotypes)

pdf("Analysis/QTLplots3Generation.pdf")
par(mfrow=c(3,3))
for(x in 1:9){
  plot(x=c(0,10), y=c(0,10), t="n", ylab="LOD Score", xlab="", xaxt="n", main=names(GTLod)[x])
  points(GTLod[[x]][,4],t="l")
  abline(h = Threshold, lty=1)    
  abline(h = CIAll[x], lty=2)
  axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
}
dev.off()

### Correct the phenotypic data 
if(!file.exists("Analysis/CorQTLdataAll-3Generation.txt")){
  AllCorPheno <- NULL
  for (Phenotype in GrowthTraits){
    onlyEnv <- lm(as.numeric(QTLdataAll[,Phenotype]) ~ QTLdataAll[,"Generation"] + QTLdataAll[,"Family"] + QTLdataAll[,"Schlupf"])
    corPheno <- round(onlyEnv$residuals + onlyEnv$coefficients["(Intercept)"])
    AllCorPheno <- cbind(AllCorPheno, corPheno)
  }
  colnames(AllCorPheno) <- GrowthTraits
  CorQTLdataAll <- cbind(QTLdataAll[,c("ID.Nr","Generation","Family","Schlupf")], AllCorPheno ,QTLdataAll[,SNPsForAnalysis])
  write.table(CorQTLdataAll,"Analysis/CorQTLdataAll-3Generation.txt",sep="\t",row.names = FALSE,quote = FALSE)
}else{                                                                                                        # Or load them from Disk if we already have them
  cat("The corrected phenotypes are already calculated, please find it in the Analysis folder\n")
}

