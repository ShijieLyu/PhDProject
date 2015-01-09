# Calculate the proportion of the missing data of every SNP in the two different lines
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Dec, 2014
# first written Dec, 2014

setwd("C:/Codes_R_analysis_for_Phd/600KSNPchip")

SNPsforfinemapping <- read.table("Analysis/SNPsforfinemapping.txt", sep="\t", header=TRUE)
chipannotation     <- read.table("RawData/chickenNumber.txt", sep ="\t", header=TRUE)

chipannotation[,"ID_Chip"] <- paste0("X", chipannotation[,"ID_Chip"])            # Names in genotypes start with an X, change the column also

pop1 <- chipannotation[which(chipannotation[,"popPk"] == "IB77xx"), "ID_Chip"]
pop2 <- chipannotation[which(chipannotation[,"popPk"] == "IBNHxx"), "ID_Chip"]

missingdata <- matrix(NA, nrow(SNPsforfinemapping), 2)
for(x in 1:nrow(SNPsforfinemapping)){
  perc1 <- (sum(is.na(SNPsforfinemapping[x, pop1])) / length(pop1)) * 100
  perc2 <- (sum(is.na(SNPsforfinemapping[x, pop2])) / length(pop2)) * 100
  missingdata[x,] <- c(perc1, perc2)
}

rownames(missingdata) <- rownames(SNPsforfinemapping)
colnames(missingdata) <- c("IB77xx", "IBNHxx")

write.table(missingdata, "Analysis/MissingDataMatrix.txt", sep = "\t")


