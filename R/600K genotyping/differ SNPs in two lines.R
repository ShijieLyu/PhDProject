# selecting SNPs which are different in two lines
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Feb, 2016
# first written Feb, 2016

setwd("D:/Chicken/600KSNPchip") 

genotypes <- read.table("RawData/600kSNPgenotypes.txt", sep="\t", na.strings = "-1", header=TRUE, check.names = FALSE)                                 # Load the genotypes
annotation <- read.table("RawData/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", na.strings = "---", header=TRUE)                  # Load the annotation
chipannotation <- read.table("RawData/chickenNumber.txt", sep ="\t", header=TRUE, check.names = FALSE)

WLpop <- as.character(chipannotation[which(chipannotation[,"popPk"] == "IB77xx"), "ID_Chip"])
NHpop <- as.character(chipannotation[which(chipannotation[,"popPk"] == "IBNHxx"), "ID_Chip"])

differSNPs <- NULL
for(x in 1:nrow(genotypes)){                                  # Decide if genotypes[,WLpop] != genotypes[,NHpop]
  if(!any(genotypes[x,] == 1, na.rm=TRUE)){                   # 0 - No heterozygous animals
    tblPop1 <- table(as.numeric(genotypes[x, WLpop]))
    if(length(tblPop1) == 1){                                 # 1 - Check if genotypes[,WLpop] are consistent
      tblPop2 <- table(as.numeric(genotypes[x, NHpop]))
      if(length(tblPop2) == 1){                               # 2 - Check if genotypes[,NHpop] are consistent
        if(names(tblPop1) != names(tblPop2)){                 # 3 - Check if genotypes[,WLpop] != genotypes[,NHpop]
          differSNPs <- rbind(differSNPs, genotypes[x,])
        }
      }
    }
  }
}

dim(differSNPs)                                    
rownames(differSNPs) <- differSNPs[,"ID"]

write.table(differSNPs[,c(WLpop,NHpop)], "Analysis/differSNPs.txt", sep="\t")

AnnoSNPs <- annotation[which(annotation[,"Probe.Set.ID"] %in% differSNPs[,"ID"]),]   # The differ SNPs with annotation