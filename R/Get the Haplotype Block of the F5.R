# F5 generation haplotype analysis
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Jun, 2015
# first written Jun, 2015

setwd("D:/Chicken/60KSNPchip")
genotypes <- read.table("RawData/Full _Data_Table.txt", header=TRUE, sep="\t",na.strings = "--",check.names=FALSE, colClasses="character")
annotation <- read.table("RawData/Chicken_60K_allResults.txt", header = TRUE, sep="\t",colClasses="character")
F5chicken <- read.table("RawData/F5phenotypeDataforQTL.txt", header = TRUE, sep="\t",na.strings = ".")

annotation <- annotation[,c("SNP.ID","flanking.seq","chr..galGal4.","position..galGal4.")]
annotation[,"position..galGal4."] <- as.numeric(gsub(",","", annotation[,"position..galGal4."]))
annotation <- annotation[-which(is.na(annotation[,"position..galGal4."])),]

annotation <- annotation[which(annotation[,"SNP.ID"] %in% genotypes[,"Name"]),]
genotypes <- genotypes[which(genotypes[,"Name"] %in% annotation[,"SNP.ID"]),]

Reorder <- match(annotation[,"SNP.ID"], genotypes[,"Name"])
genotypes <- genotypes[Reorder,]

GenoUpdataPos <- cbind(Name = genotypes[,"Name"], Position = annotation[,"position..galGal4."], genotypes[,7:ncol(genotypes)])

