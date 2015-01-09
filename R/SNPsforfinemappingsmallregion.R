# The location for the SNPs which will be used for fine mapping(the small region and the big region)
#
# copyright (c) 2014-2020 - Shijie Lyu
# last modified Dec.4, 2014
# first written Dec.4, 2014

setwd("C:/Data/600KSNPchip")

bigregionSNPs <- read.table("Analysis/finemappingSNPswithannotation.txt", sep = "\t", header=TRUE)    # load the SNPs for fine mapping which are in the big region(3855 SNPs)

smallregionstart <- 69583407
smallregionend   <- 78715886

smallregionind <- which(as.numeric(as.character(bigregionSNPs[,"Physical.Position"])) > smallregionstart &          # Larger than the start
                     as.numeric(as.character(bigregionSNPs[,"Physical.Position"])) < smallregionend)                # Smaller then the end

smallregionSNPs <- bigregionSNPs[smallregionind, ]           # SNPs in the small region(1246 SNPs)
write.table(smallregionSNPs, "Analysis/finemappingSNPinsmallregion.txt", sep="\t")