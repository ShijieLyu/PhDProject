# Annotated SNPs which will be used for fine mapping
#
# copyright (c) 2014-2020 - Shijie Lyu
# last modified Dec, 2014
# first written Dec, 2014

setwd("C:/Data/600KSNPchip")

snps <- read.table("Analysis/SNPsforfinemapping.txt", sep = "\t", header=TRUE)                                     # load the all the SNPs for fine mapping
annotation <- read.table("RawData/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", na.strings = "---", header=TRUE)   # load the annotation file

targetSNPsNr <- which(annotation[,"Probe.Set.ID"] %in% rownames(snps))                                             # Give the Fine Mapping SNPS the annotation
finemappingSNPswithannotation <- annotation[targetSNPsNr,]
write.table(finemappingSNPswithannotation, "Analysis/finemappingSNPswithannotation.txt", sep="\t")
