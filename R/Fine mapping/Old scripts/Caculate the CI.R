# Calculate the Confidence Interval for the Growth Traits of every 5 weeks using F11 and F12
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Apr, 2015
# first written Apr, 2015

setwd("D:/Chicken/Rdata/FineMapping")
F11genotypes <- read.table("RawData/F11_Genotypes.txt", na.strings="-", header=TRUE, sep="\t")
F11BodyWeight <- read.table("RawData/F11_BodyWeight_Phenotypes.txt", na.strings="",header=TRUE, sep="\t")
F12Traits <- read.table("Analysis/F12_QTLdataWithBWG-Geno NA fill in.txt", na.strings=c("NA","?"),header=TRUE, sep="\t")

F11BodyWeight <- cbind(F11BodyWeight, BWG05 = rep(NA, nrow(F11BodyWeight)),BWG510 = rep(NA, nrow(F11BodyWeight)),BWG1015 = rep(NA, nrow(F11BodyWeight)),BWG1520 = rep(NA, nrow(F11BodyWeight)))
for (x in 1:nrow(F11BodyWeight)){
  F11BodyWeight[x,"BWG05"]   <- F11BodyWeight[x,"Gew_5Wo"] - F11BodyWeight[x,"Gew_1d"]
  F11BodyWeight[x,"BWG510"]  <- F11BodyWeight[x,"Gew_10Wo"] - F11BodyWeight[x,"Gew_5Wo"]
  F11BodyWeight[x,"BWG1015"] <- F11BodyWeight[x,"Gew_15Wo"] - F11BodyWeight[x,"Gew_10Wo"]
  F11BodyWeight[x,"BWG1520"] <- F11BodyWeight[x,"Gew_20Wo"] - F11BodyWeight[x,"Gew_15Wo"]
}

F11SNPnames <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")
F11QTLdata <- cbind(F11BodyWeight,F11genotypes[,F11SNPnames])

### 3-Select the normal diet and individuals with the BodyWeight data of every 5 weeks, should be 46 chickens
F12QTLdata <- F12Traits[,c(names(F11BodyWeight),F11SNPnames)]
F12QTLdata <- F12QTLdata[which(!is.na(F12QTLdata[,"Gew_20Wo"]) & F12QTLdata[,"Futter"] == "normal"),] 

### 4-Combine the F11 and F12 data in one file
QTLdata1112 <- rbind(F11QTLdata,F12QTLdata) 

GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
Covs  <- c("Generation","Family","Schlupf")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

# Variables
genotypes    <- QTLdata1112[,SNPsForAnalysis]
phenotypes   <- QTLdata1112[,GrowthTraits]
covariates   <- QTLdata1112[,Covs]

# Function to map QTLs
mapQTLs <- function(genotypes, phenotypes, covariates){
  # Matrix LOD scores of markers x phenotypes
  LODS <- matrix(NA, ncol(phenotypes), ncol(genotypes),dimnames=list(colnames(phenotypes), colnames(genotypes))) # Markers x phenotypes
  for(x in 1:ncol(phenotypes)){  #cat(x,"\n")
    for(y in 1:ncol(genotypes)){ #cat(y,"\n")
      res <- anova(lm(as.numeric(phenotypes[,x]) ~ covariates[,"Generation"] + covariates[,"Family"] + covariates[,"Schlupf"] + genotypes[,y]))
      LODS[x, y] <- -log10(res[[5]][length(res[[5]]) - 1])
    }
  }
  return(LODS)
}

# Map the real data
LodForAll <- mapQTLs(genotypes, phenotypes, covariates)

## Confidence Interval
if(!file.exists("Analysis/CIboots1000-withF10.txt")){
  boots <- NULL
  for(x in 1:1000){
    idx <- sample(nrow(genotypes), 97, replace=TRUE)
    boots <- rbind(boots, mapQTLs(genotypes[idx, ], phenotypes[idx,], covariates[idx,]))
  }
    boots <- cbind(Traits = rownames(boots), boots)
    write.table(boots,"Analysis/CIboots1000-withF10.txt",sep="\t")
  }else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading boots information from disk\n")
  boots <- read.table("Analysis/CIboots1000-withF10.txt",sep="\t", header=TRUE, row.names=FALSE)
}

CIAll <- NULL
for (trait in names(phenotypes)){
  TraitBoots <- boots[which(rownames(boots) == trait),]
  topM <- names(which.max(LodForAll[trait,]))
  BootI <- sort(TraitBoots[,topM])[1000 * 0.025]
  CIAll <- cbind(CIAll, BootI)
}
colnames(CIAll) <- names(phenotypes)

# Calculate for the position of the CI

SNPsinfo <- read.table("RawData/SNPsinfo.txt",header=TRUE,sep="\t")
SNPPosition <- SNPsinfo[, c("Markers","Location")]
rownames(SNPPosition) <- SNPPosition[,"Markers"]


CIPosCal <- function(SNP1, SNP2, Trait){
  y1 <- as.numeric(LodForAll[Trait, SNP1])
  y2 <- as.numeric(LodForAll[Trait, SNP2])
  x1 <- as.numeric(SNPPosition[SNP1, "Location"])
  x2 <- as.numeric(SNPPosition[SNP2, "Location"])
  a <- as.numeric((y1-y2)/(x1-x2))
  b <- y1-a*(x1)
  CIPos <- (CIAll[,Trait] -b)/a
  return(CIPos)
}

CIPosCal("rs318175270","rs14492508","Gew_5Wo")

