# Run of homozygosity analysis in the parents and Identical by descent analysis
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Shijie Lyu
# last modified Feb, 2016
# first written Feb, 2016

setwd("D:/Chicken/Rdata/F5_QTL_Analysis")
annotation <- read.table("RawData/Chicken_60K_annotations.txt", header = TRUE, sep="\t",colClasses="character")   # chicken 60k-SNP chip annotations
genotypes <- read.table("RawData/F5_60K_genotyping_results.txt", header= TRUE, sep="\t",na.strings = "--",check.names=FALSE, colClasses="character")   # Load Genotypes
phenotypes <- read.table("RawData/F5_organised_phenotypes.txt", header = TRUE, sep="\t",na.strings = "", colClasses= c(rep("character",3), rep("factor",5), rep("numeric",66)))         # Load Phenotypes

### Select the probes which have annotation  
annotation <- annotation[,c("SNP.ID","flanking.seq","chr..galGal4.","position..galGal4.")]
annotation[,"position..galGal4."] <- as.numeric(gsub(",","", annotation[,"position..galGal4."]))
annotation <- annotation[-which(is.na(annotation[,"position..galGal4."])),]
annotation <- annotation[which(annotation[,"SNP.ID"] %in% genotypes[,"Name"]),]
genotypes <- genotypes[which(genotypes[,"Name"] %in% annotation[,"SNP.ID"]),]

### Select the F5 chicken in the genotyping file
Ind_Nr <- NULL
for(x in 1:nrow(phenotypes)){
 PerInd_Nr <- paste0(phenotypes[x,"ID.Nr"], ".Top Alleles")    # Get the individual name in the genotyping file
 Ind_Nr <- c(Ind_Nr, PerInd_Nr)
}
genotypes <- t(genotypes[,c("Name",Ind_Nr)])
colnames(genotypes) <- genotypes["Name",]
genotypes <- cbind(ID.Nr=as.character(phenotypes[,"ID.Nr"]),genotypes[-1,])    # attention: remove the first row of the regotypes first which is the same whih the rowname
