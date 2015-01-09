# 600K SNP-chip data analysis, selecting SNPs for fine mapping of the QTL found by M
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Dec, 2014
# first written Dec, 2014

setwd("C:/Codes_R_analysis_for_Phd/600KSNPchip")

genotypes <- read.table("RawData/600kSNPgenotypes.txt", sep="\t", na.strings = "-1", header=TRUE)                                 # Load the genotypes
annotation <- read.table("RawData/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", na.strings = "---", header=TRUE)                  # Load the annotation
chipannotation <- read.table("RawData/chickenNumber.txt", sep ="\t", header=TRUE)

chipannotation[,"ID_Chip"] <- paste0("X", chipannotation[,"ID_Chip"])            # Names in genotypes start with an X, change the column also

startOn4 <- 59553432    # Start of our region
endOn4   <- 84774762    # End of our region

inChr4Region    <- which(annotation[,"Chromosome"] == 4 &                                                 # On chromosome 4
                         as.numeric(as.character(annotation[,"Physical.Position"])) > startOn4 &          # Larger then our start
                         as.numeric(as.character(annotation[,"Physical.Position"])) < endOn4)             # Smaller then our end`
annotation      <- annotation[inChr4Region, ]                                                             # Only use the annotation in our region

genotypesOn4    <- which(genotypes[,"ID"] %in% annotation[,"Probe.Set.ID"])                               # Which genotypes do we have ?
genotypes <- genotypes[genotypesOn4, ]                                                                    # ONly use the ones which are in our region

dim(annotation)       # Dimensions of the annotation in our chromosome 4 region
dim(genotypes)        # Dimensions of our new genotype matrix  "a total 12232 SNPs in this region"

pop1 <- chipannotation[which(chipannotation[,"popPk"] == "IB77xx"), "ID_Chip"]
pop2 <- chipannotation[which(chipannotation[,"popPk"] == "IBNHxx"), "ID_Chip"]

snpsForFinemapping <- NULL
for(x in 1:nrow(genotypes)){                                  # Decide if genotypes[,pop1] != genotypes[,pop2]
  if(!any(genotypes[x,] == 1, na.rm=TRUE)){                   # 0 - No heterozygous animals
    tblPop1 <- table(as.numeric(genotypes[x, pop1]))
    if(length(tblPop1) == 1){                                 # 1 - Check if genotypes[,pop1] are consistent
      tblPop2 <- table(as.numeric(genotypes[x, pop2]))
      if(length(tblPop2) == 1){                               # 2 - Check if genotypes[,pop2] are consistent
        if(names(tblPop1) != names(tblPop2)){                 # 3 - Check if genotypes[,pop1] != genotypes[,pop2]
          snpsForFinemapping <- rbind(snpsForFinemapping, genotypes[x,])
        }
      }
    }
  }
  cat("Done",x,"/", nrow(genotypes),"\n")
}

dim(snpsForFinemapping)                                      # "a total of 3855 SNPs"
rownames(snpsForFinemapping) <- snpsForFinemapping[,"ID"]

write.table(snpsForFinemapping[,c(pop1,pop2)], "Analysis/SNPsforfinemapping.txt", sep="\t")