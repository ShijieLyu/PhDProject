# Fine Mapping for the chicken using F11 and F12 Generation
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Apr, 2015
# first written Apr, 2015


setwd("D:/Chicken/Rdata/FineMapping")
F10genotypes <- read.table("RawData/F10_Genotypes.txt", header=TRUE, sep="\t")
F10BodyWeight <- read.table("RawData/F10_BodyWeight_Phenotypes.txt", header=TRUE, sep="\t")
F11genotypes <- read.table("RawData/F11_Genotypes.txt", na.strings="-", header=TRUE, sep="\t")
F11BodyWeight <- read.table("RawData/F11_BodyWeight_Phenotypes.txt", na.strings="",header=TRUE, sep="\t")
F12Traits <- read.table("Analysis/F12_QTLdataWithBWG-Geno NA fill in.txt", na.strings=c("NA","?"),header=TRUE, sep="\t")

### 1-linkage disequilibrium format construct for the SHEsis website FOR F11 9 markers

F11SNPnames <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

if(!file.exists("Analysis/F11LDGenotypes.txt")){
  F11LDGenotypes <- F11genotypes[,c("ID.Nr",F11SNPnames)]
  for (mrow in 1:nrow(F11LDGenotypes)){
    onerow <- NULL
    for (mcol in 2:ncol(F11LDGenotypes)){
      alleleA <- strsplit(as.character(F11LDGenotypes[mrow,mcol]),"")[[1]][1]
      alleleB <- strsplit(as.character(F11LDGenotypes[mrow,mcol]),"")[[1]][2]
      onerow <- cbind(onerow, paste(alleleA,alleleB))
    }
  cat(c(as.character(F11LDGenotypes[mrow,"ID.Nr"]),onerow), "\n", file="Analysis/F11LDGenotypes.txt", append = TRUE)
  }
}

F10SNPnames <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

if(!file.exists("Analysis/F10LDGenotypes.txt")){
  F10LDGenotypes <- F10genotypes[,c("ID.Nr",F10SNPnames)]
  for (mrow in 1:nrow(F10LDGenotypes)){
    onerow <- NULL
    for (mcol in 2:ncol(F10LDGenotypes)){
      alleleA <- strsplit(as.character(F10LDGenotypes[mrow,mcol]),"")[[1]][1]
      alleleB <- strsplit(as.character(F10LDGenotypes[mrow,mcol]),"")[[1]][2]
      onerow <- cbind(onerow, paste(alleleA,alleleB))
    }
  cat(c(as.character(F10LDGenotypes[mrow,"ID.Nr"]),onerow), "\n", file="Analysis/F10LDGenotypes.txt", append = TRUE)
  }
}

### 2-Add the bodyweight gain and the genotypes in the F11QTLdata

F11BodyWeight <- cbind(F11BodyWeight, BWG05 = rep(NA, nrow(F11BodyWeight)),BWG510 = rep(NA, nrow(F11BodyWeight)),BWG1015 = rep(NA, nrow(F11BodyWeight)),BWG1520 = rep(NA, nrow(F11BodyWeight)))
for (x in 1:nrow(F11BodyWeight)){
  F11BodyWeight[x,"BWG05"]   <- F11BodyWeight[x,"Gew_5Wo"] - F11BodyWeight[x,"Gew_1d"]
  F11BodyWeight[x,"BWG510"]  <- F11BodyWeight[x,"Gew_10Wo"] - F11BodyWeight[x,"Gew_5Wo"]
  F11BodyWeight[x,"BWG1015"] <- F11BodyWeight[x,"Gew_15Wo"] - F11BodyWeight[x,"Gew_10Wo"]
  F11BodyWeight[x,"BWG1520"] <- F11BodyWeight[x,"Gew_20Wo"] - F11BodyWeight[x,"Gew_15Wo"]
}

F11QTLdata <- cbind(F11BodyWeight,F11genotypes[,F11SNPnames])

F10BodyWeight <- cbind(F10BodyWeight, BWG05 = rep(NA, nrow(F10BodyWeight)),BWG510 = rep(NA, nrow(F10BodyWeight)),BWG1015 = rep(NA, nrow(F10BodyWeight)),BWG1520 = rep(NA, nrow(F10BodyWeight)))
for (x in 1:nrow(F10BodyWeight)){
  F10BodyWeight[x,"BWG05"]   <- F10BodyWeight[x,"Gew_5Wo"] - F10BodyWeight[x,"Gew_1d"]
  F10BodyWeight[x,"BWG510"]  <- F10BodyWeight[x,"Gew_10Wo"] - F10BodyWeight[x,"Gew_5Wo"]
  F10BodyWeight[x,"BWG1015"] <- F10BodyWeight[x,"Gew_15Wo"] - F10BodyWeight[x,"Gew_10Wo"]
  F10BodyWeight[x,"BWG1520"] <- F10BodyWeight[x,"Gew_20Wo"] - F10BodyWeight[x,"Gew_15Wo"]
}

F10QTLdata <- cbind(F10BodyWeight,F10genotypes[,F10SNPnames])

### 3-Select the normal diet and individuals with the BodyWeight data of every 5 weeks, should be 46 chickens
F12QTLdata <- F12Traits[,c(names(F11BodyWeight),F11SNPnames)]
F12QTLdata <- F12QTLdata[which(!is.na(F12QTLdata[,"Gew_20Wo"]) & F12QTLdata[,"Futter"] == "normal"),] 

### 4-Combine the F11 and F12 data in one file
QTLdata1112 <- rbind(F11QTLdata,F12QTLdata) 

### 5-Association analysis for F11 and F12 growth Traits
GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

GTPvalue <- vector("list", length(GrowthTraits)) 
names(GTPvalue) <- GrowthTraits

for (Phenotype in GrowthTraits){                                                       # Model after selection
  EachTraitPvalue <- NULL
  for (OneSNP in SNPsForAnalysis){
    res <- anova(lm(as.numeric(QTLdata1112[,Phenotype]) ~ QTLdata1112[,"Generation"] + QTLdata1112[,"Family"] + QTLdata1112[,"Schlupf"] + QTLdata1112[,OneSNP]))
      Pfactors <- NULL
      for (x in 1:5){
        Pfactors <- c(Pfactors, round((sum(res[x,2])/sum(res[2]))*100, digits = 2))    # percentage(contribution) per factor
      }
    EachTraitPvalue <- rbind(EachTraitPvalue, c(res[[5]],Pfactors))
    colnames(EachTraitPvalue) <- c("Generation","Family","Schlupf","SNP", "Residuals","Generation%","Family%","Schlupf%","SNP%", "Residuals%")  
  }
  rownames(EachTraitPvalue) <- SNPsForAnalysis
  GTPvalue[[Phenotype]] <- EachTraitPvalue
}
GTPvalue$Gew_15Wo
#library(marray)
#write.list(GTPvalue,"Analysis/GrowthTraitsPvalue-every5weeeks-Geno NA fill in.txt")

GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

#pdf("Analysis/AllSnpsPlot.pdf")
#par(mfrow=c(3,3))
#for(x in 1:9){
  #plot(x=c(0,10), y=c(0,10), t="n", ylab="LOD Score", xlab="", xaxt="n", main=names(GTPvalue)[x])
  #points(-log10(GTPvalue[[x]][,4]),t="l")
  #abline(h=-log10(0.05/(9*9)), lty=2)    # Threshold = -log10(0.05/(9*9))
  #axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
#}
#dev.off()

### 6-Correct the phenotypic data 
AllCorPheno <- NULL
for (Phenotype in GrowthTraits){
  onlyEnv <- lm(as.numeric(QTLdata1112[,Phenotype]) ~ QTLdata1112[,"Generation"] + QTLdata1112[,"Family"] + QTLdata1112[,"Schlupf"])
  corPheno <- round(onlyEnv$residuals + onlyEnv$coefficients["(Intercept)"])
  AllCorPheno <- cbind(AllCorPheno, corPheno)
}
colnames(AllCorPheno) <- GrowthTraits
CorQTLdata1112 <- cbind(QTLdata1112[,c("ID.Nr","Generation","Family","Schlupf")], AllCorPheno ,QTLdata1112[,SNPsForAnalysis])
#write.table(CorQTLdata1112,"Analysis/CorQTLdata1112.txt",sep="\t",row.names = FALSE,quote = FALSE)

SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")
pdf("Analysis/rs315966269.pdf")
par(mfrow=c(3,3))
for (Phenotype in GrowthTraits){ 
  plot(CorQTLdata1112[,Phenotype]~ CorQTLdata1112[,"rs315966269"], ylab="Body Weight(g)", xlab= "Genotypes", main= Phenotype)
}
dev.off()

### 7-Check the Allele origin

annotation <- read.table("D:/Chicken/600KSNPchip/RawData/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", na.strings = "---", header=TRUE)                  
SNPsforfinemapping <- read.table("D:/Chicken/600KSNPchip/Analysis/SNPsforfinemapping.txt", sep ="\t", header=TRUE)
chipannotation <- read.table("D:/Chicken/600KSNPchip/RawData/chickenNumber.txt", sep ="\t", header=TRUE)

SNPsForAnalysis <- c("rs10725580", "rs316584759", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

number <- NULL                                                        # the ProbeID for every SNP
for (OneSNP in SNPsForAnalysis){
  number <- c(number, which(annotation[,"dbSNP.RS.ID"]== OneSNP))
}
ProbeIDs <- as.character(annotation[number,"Probe.Set.ID"])

Probenumber <- NULL                                                   # Given the ProbeID, find the genotype for every individual
for (OneProbe in ProbeIDs){
  Probenumber <- c(Probenumber, which(rownames(SNPsforfinemapping) == OneProbe))
}

chipannotation[,"ID_Chip"] <- paste0("X", chipannotation[,"ID_Chip"]) 
WL77 <- chipannotation[which(chipannotation[,"popPk"] == "IB77xx"), "ID_Chip"]
NHI <- chipannotation[which(chipannotation[,"popPk"] == "IBNHxx"), "ID_Chip"]
WL77info <- cbind(annotation[number,c("Probe.Set.ID","dbSNP.RS.ID","Flank","Allele.A","Allele.B")],SNPsforfinemapping[Probenumber,WL77])
NHIinfo <- cbind(annotation[number,c("Probe.Set.ID","dbSNP.RS.ID","Flank","Allele.A","Allele.B")],SNPsforfinemapping[Probenumber,NHI])

### 8- For all traits for the 20week-slaughtered-chickens fed with normal feed

F12AllTraits20Week <- F12Traits[which(F12Traits[,"SchlachtWo"]=="20W" & F12Traits[,"Futter"] == "normal"),]  # F12AT20W (F12AllTraits20Week)

GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

F12AT20W  <- names(F12AllTraits20Week)[-c(99,104)]
F12AT20W  <- c(GrowthTraits, F12AT20W[c(37:39,41:143)])

F12AT20WLod <- vector("list", length(F12AT20W)) 
names(F12AT20WLod) <- F12AT20W

for (Phenotype in F12AT20W){                                                       # Model after selection
  EveryTraitPvalue <- NULL
    if(length(which(!is.na(F12AllTraits20Week[,Phenotype]))) == 0){
        EveryTraitPvalue <- rbind(EveryTraitPvalue, NA)
        colnames(EveryTraitPvalue) <- "SNP"
    }else{
      for (OneSNP in SNPsForAnalysis){
        res <- anova(lm(as.numeric(F12AllTraits20Week[,Phenotype]) ~ F12AllTraits20Week[,"Family"] + F12AllTraits20Week[,"Schlupf"] + F12AllTraits20Week[,OneSNP]))
        Pfactors <- NULL
        for (x in 1:4){
          Pfactors <- c(Pfactors, round((sum(res[x,2])/sum(res[2]))*100, digits = 2))    # percentage(contribution) per factor
        }
      EveryTraitPvalue <- rbind(EveryTraitPvalue, c(-log10(res[[5]]),Pfactors))
      colnames(EveryTraitPvalue) <- c("Family","Schlupf","SNP", "Residuals","Family%","Schlupf%","SNP%", "Residuals%")  
      }
    rownames(EveryTraitPvalue) <- SNPsForAnalysis
    F12AT20WLod[[Phenotype]] <- EveryTraitPvalue
  }
}
F12AT20WLod$Gew_15Wo
threshold <- -log10(0.05/9)
#library(marray)
#write.list(F12AT20WLod,"Analysis/F12AllTraits20WeekLod--Geno NA fill in.txt")

if(!file.exists("Analysis/F12AllTraits20WLinkageAnalysis.txt")){           
  for(x in 1:length(F12AT20WLod)){
    SigSNP <- which(F12AT20WLod[[x]][,"SNP"] > threshold)                                                      # Analyse the resulting profiles, and look at which markers are above the threshold
    cat(names(F12AT20WLod)[x], names(SigSNP),"\n",file = "Analysis/F12AllTraits20WLinkageAnalysis.txt",append=TRUE)                                                                  
  }
  SigSNP <- read.table("Analysis/F12AllTraits20WLinkageAnalysis.txt",sep="\t", header=FALSE)
}else{
  cat("Loading Linkage Analysis information from disk\n")
  SigSNP <- read.table("Analysis/F12AllTraits20WLinkageAnalysis.txt",sep="\t", header=FALSE)   # I remove the Space at the end of each traits manually!!!!!
}

### Make a plot for the Significant traits
SigTraits <- NULL                                                                                      
for (x in 1:nrow(SigSNP)){
  if(length(unlist(strsplit(as.character(SigSNP[x,1]), " "))) > 1) {
  SigTraits <- c(SigTraits, strsplit(as.character(SigSNP[x,1]), " ")[[1]][1])
  }
}

pdf("Analysis/SigBodycomposition.pdf")
par(mfrow=c(4,4))
for(SigTrait in SigTraits){
  plot(x=c(0,10), y=c(0,5), t="n", ylab="LOD Score", xlab="", xaxt="n", main= SigTrait)
  points(F12AT20WLod[[SigTrait]][,3],t="l")
  abline(h=-log10(0.05/9), lty=2)    
  axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
}
dev.off()

### 9- Correlation for all traits
F12TraitsName <- names(F12Traits)[20:149]
F12TraitsOnly <- F12Traits[,F12TraitsName]

for (EachTrait in F12TraitsName){
  if(length(which(!is.na(F12TraitsOnly[,EachTrait]))) == 0){
    F12TraitsOnly[,EachTrait] <- NULL
  }
}
pdf("Analysis/CorForAllF12Traits.pdf")
F12TraitsCor <- cor(data.matrix(F12TraitsOnly),method="spearman",use = "pairwise")
#image(F12TraitsCor)
plot(hclust(dist(t(F12TraitsCor))), cex=0.3)
dev.off()

cor(F12TraitsOnly[,"Abdominalfett"],F12TraitsOnly[, "BW.nuchtern"], method="spearman", use = "pairwise")

### 10- Catch the genes in the target Region

library("biomaRt")                                                                                              # How many genes in different regions for different traits
Mart <- useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
GenesForGew_5Wo_BWG05 <- getBM(attributes = c("chromosome_name","ensembl_gene_id","start_position","end_position","external_gene_name","description"), filters = "chromosomal_region", values = "4:69352655:73629630", mart = Mart)
write.table(GenesForGew_5Wo_BWG05, row.names = FALSE, "Analysis/GenesForGew_5Wo_BWG05.txt", sep="\t") 

GenesForGew_10Wo__BWG510 <- getBM(attributes = c("chromosome_name","ensembl_gene_id","start_position","end_position","external_gene_name","description"), filters = "chromosomal_region", values = "4:73629630:77032597", mart = Mart)
write.table(GenesForGew_10Wo__BWG510, row.names = FALSE, "Analysis/GenesForGew_10Wo__BWG510.txt", sep="\t") 

GenesForGew_15Wo_Gew_20Wo <- getBM(attributes = c("chromosome_name","ensembl_gene_id","start_position","end_position","external_gene_name","description"), filters = "chromosomal_region", values = "4:69352655:77032597", mart = Mart)
write.table(GenesForGew_15Wo_Gew_20Wo, row.names = FALSE, "Analysis/GenesForGew_15Wo_Gew_20Wo.txt", sep="\t") 

GenesForBWG1015 <- getBM(attributes = c("chromosome_name","ensembl_gene_id","start_position","end_position","external_gene_name","description"), filters = "chromosomal_region", values = "4:72301597:74439992", mart = Mart)
write.table(GenesForBWG1015, row.names = FALSE, "Analysis/GenesForBWG1015.txt", sep="\t") 


### 11- 1d-20W slaughtered chickens with the bodyweight for every week --should be 41 chickens
F12allweeks <- F12Traits[which(!is.na(F12Traits[,"Gew_14Wo"])),] 
F12allweeks <- F12allweeks[,c(1:40,146:150)]

GrowthTraitsAllweek  <- names(F12allweeks[20:44])
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

GTPvalueAllWeeks <- vector("list", length(GrowthTraitsAllweek)) 
names(GTPvalueAllWeeks) <- GrowthTraitsAllweek

for (Phenotype in GrowthTraitsAllweek){
  EachTraitPvalue <- NULL
  if (Phenotype != "Gew_1d"){# && Phenotype != "Gew_1Wo" && Phenotype != "Gew_2Wo" && Phenotype != "Gew_3Wo" && Phenotype != "Gew_4Wo" && Phenotype != "Gew_5Wo"){
    for (OneSNP in SNPsForAnalysis){
      res <- anova(lm(as.numeric(F12allweeks[,Phenotype]) ~ F12allweeks[,"Family"] + F12allweeks[,"Futter"] + F12allweeks[,"Schlupf"] + F12allweeks[,OneSNP]))
      EachTraitPvalue <- rbind(EachTraitPvalue, res[[5]])
      #colnames(EachTraitPvalue) <- c("Family", "Futter", "Schlupf", "SNP", "Residuals")  
    }
  }else{
    for (OneSNP in SNPsForAnalysis){
      res <- anova(lm(as.numeric(F12allweeks[,Phenotype]) ~ F12allweeks[,"Family"] + F12allweeks[,"Schlupf"] + F12allweeks[,OneSNP]))
      EachTraitPvalue <- rbind(EachTraitPvalue, res[[5]])
      #colnames(EachTraitPvalue) <- c("Family","Schlupf", "SNP", "Residuals")  
    }
  }
  rownames(EachTraitPvalue) <- SNPsForAnalysis
  GTPvalueAllWeeks[[Phenotype]] <- EachTraitPvalue
}


### 12- FDR
GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

GTPvalue <- vector("list", length(GrowthTraits)) 
names(GTPvalue) <- GrowthTraits

for (Phenotype in GrowthTraits){                                                       # Model after selection
  EachTraitPvalue <- NULL
  for (OneSNP in SNPsForAnalysis){
    res <- anova(lm(as.numeric(QTLdata1112[,Phenotype]) ~ QTLdata1112[,"Generation"] + QTLdata1112[,"Family"] + QTLdata1112[,"Schlupf"] + QTLdata1112[,OneSNP]))
      Pfactors <- NULL
      for (x in 1:5){
        Pfactors <- c(Pfactors, round((sum(res[x,2])/sum(res[2]))*100, digits = 2))    # percentage(contribution) per factor
      }
    EachTraitPvalue <- rbind(EachTraitPvalue, c(p.adjust(res[[5]], "BH",81),Pfactors))
    colnames(EachTraitPvalue) <- c("Generation","Family","Schlupf","SNP", "Residuals","Generation%","Family%","Schlupf%","SNP%", "Residuals%")  
  }
  rownames(EachTraitPvalue) <- SNPsForAnalysis
  GTPvalue[[Phenotype]] <- EachTraitPvalue
}
GTPvalue$Gew_15Wo

GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

pdf("Analysis/AllSnpsPlot-FDR.pdf")
par(mfrow=c(3,3))
for(x in 1:9){
  plot(x=c(0,10), y=c(0,6), t="n", ylab="LOD Score", xlab="", xaxt="n", main=names(GTPvalue)[x])
  points(-log10(GTPvalue[[x]][,4]),t="l")
  abline(h=-log10(0.05), lty=2)    
  axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
}
dev.off()

### 13- Permutation and the CI

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
  LODS <- matrix(NA, ncol(genotypes), ncol(phenotypes),dimnames=list(colnames(phenotypes), colnames(genotypes))) # Markers x phenotypes
  for(x in 1:ncol(phenotypes)){  #cat(x,"\n")
    for(y in 1:ncol(genotypes)){ #cat(y,"\n")
      res <- anova(lm(as.numeric(phenotypes[,x]) ~ covariates[,"Generation"] + covariates[,"Family"] + covariates[,"Schlupf"] + genotypes[,y]))
      LODS[x, y] <- -log10(res[[5]][length(res[[5]]) - 1])
    }
  }
  return(LODS)
}

# Map the real data
realData <- mapQTLs(genotypes, phenotypes, covariates)


# Permutation for getting an overall 95 % confidence interval (corrected for N markers and N phenotypes
if(!file.exists("Analysis/Permutation-1000-withF10.txt")){
  maxes <- NULL
  for(p in 1:1000){
     ngenotypes <- genotypes[sample(nrow(genotypes)), ]
     res <- mapQTLs(ngenotypes, phenotypes, covariates)
     maxes <- c(maxes, max(res))
     #cat("Finished permutation", p,"\n")
  }
     write.table(maxes,"Analysis/Permutation-1000-withF10.txt",sep="\t")
  }else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading Permutation information from disk\n")
  maxes <- read.table("Analysis/Permutation-1000-withF10.txt",sep="\t", header=TRUE)
}

# Null-distribution
hist(maxes[,1])

# Cut-off for significance 95 %
Threshold <- sort(maxes[,1])[1000 * 0.95]       #<- Threshold = 3.179255

pdf("Analysis/AllSnpsPlot-1000perputation-withF10.pdf")
par(mfrow=c(3,3))
for(x in 1:9){
  plot(x=c(0,10), y=c(0,6), t="n", ylab="LOD Score", xlab="", xaxt="n", main=names(GTPvalue)[x])
  points((-log10(GTPvalue[[x]][,4]),t="l")
  abline(h=Threshold, lty=2)    
  axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
}
dev.off()

## Confidence Interval
if(!file.exists("Analysis/CIboots.txt")){
  boots <- NULL
  for(x in 1:1000){
    idx <- sample(nrow(genotypes), 97, replace=TRUE)
    boots <- rbind(boots, mapQTLs(genotypes[idx, ], phenotypes[idx,], covariates[idx,]))
  }
    write.table(boots,"Analysis/CIboots.txt",sep="\t")
  }else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading boots information from disk\n")
  boots <- read.table("Analysis/CIboots.txt",sep="\t", header=TRUE, row.names=NULL)
}

# Get the 10 week data
G10wB <- boots[which(boots[,"row.names"] == "Gew_10Wo"),]

# What is the top marker
topM <- names(which.max(realData["Gew_10Wo",]))

# What is the LOD score below we are out of the confidence interval
BootI <- sort(G10wB[,topM])[1000 * 0.025]

plot(c(0,9), c(0,10), t = 'n')
points(realData["Gew_10Wo",], t = 'l')
abline(h=BootI)


Allbot <- NULL
for (trait in GrowthTraits){
  Eachbot <- apply(boots[which(boots[,"row.names"] == trait),], 2, function(x){sort(x)[1000 * 0.025]})
  Allbot <- rbind(Allbot ,Eachbot) 
}
rownames(Allbot) <- GrowthTraits
Allbot <- Allbot[,-(1:2)]

#top <- apply(boots, 2, function(x){sort(x)[1000 * 0.975]})
#bot <- apply(boots, 2, function(x){sort(x)[1000 * 0.025]})

pdf("Analysis/AllTraitwithCI.pdf")
par(mfrow=c(3,3))
for(x in 1:9){
  plot(x=c(0,10), y=c(0,6), t="n", ylab="LOD Score", xlab="", xaxt="n", main=names(GTPvalue)[x])
  points(-log10(GTPvalue[[x]][,4]),t="l")
  abline(h = 3.179255, lty=1)    # Threshold = sort(maxes[,1])[1000 * 0.95] = 3.179255
  abline(h = Allbot[x,which.max(-log10(GTPvalue[[x]][,4]))], lty=2)
  axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
}
dev.off()



