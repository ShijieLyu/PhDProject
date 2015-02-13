# Analyzing what the GAPs are in the big region for finemapping which is observed through the 600k SNP chip
# copyright (c) 2014-2020 - Danny Arends & Shijie Lyu
# last modified Feb, 2015
# first written Feb, 2015

setwd("C:/Data/600KSNPchip")

## Find the Gaps in the target region, the GAPs means a region with no the parental snps
snpsbigregion <- read.table("Analysis/finemappingSNPswithannotation.txt", sep = "\t", header=TRUE)  # load the full data
snpslocation <- snpsbigregion[,"Physical.Position"]                                                 # load the col-name which means location 
 
windowsize <- 20000                                                                                 # create the size of the window
windowstarts <- seq(min(snpslocation), max(snpslocation)-windowsize, by = windowsize)                     # calculate the start for exery window
windownumber <- length(windowstarts)                                                                      # how many windows will use 
SNPsNrinwindows <- NULL
for(x in 1:windownumber){                                                                           # creat the sliding windows with the 20kb size
  SNPsNrperwindow <- length(which(snpslocation <= windowstarts[x] + (windowsize-1) & snpslocation >= windowstarts[x]))
  SNPsNrinwindows <- c(SNPsNrinwindows, SNPsNrperwindow)
} 
noSNPswindow <- which(SNPsNrinwindows == 0)                                                         # the Region with no SNPs which is parental SNPs

## Distinguish and annotate all the SNPs in the GAPs

genotypes <- read.table("RawData/600kSNPgenotypes.txt", sep="\t", na.strings = "-1", header=TRUE)                                 # Load the genotypes
annotation <- read.table("RawData/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", na.strings = "---", header=TRUE)                  # Load the annotation
chipannotation <- read.table("RawData/chickenNumber.txt", sep ="\t", header=TRUE) 
chipannotation[,"ID_Chip"] <- paste0("X", chipannotation[,"ID_Chip"])            # Names in genotypes start with an X, change the column also

SNPsinGap <- NULL
for (x in 1:length(noSNPswindow)) {
  StartPerGap <- windowstarts[noSNPswindow[x]]
  EndPerGap <- windowstarts[noSNPswindow[x]] + 19999
  SNPinOneGap <- which(annotation[,"Chromosome"] == 4 &                                                 
                         as.numeric(as.character(annotation[,"Physical.Position"])) > StartPerGap &        
                         as.numeric(as.character(annotation[,"Physical.Position"])) < EndPerGap) 
  SNPsinGap <- c(SNPsinGap, SNPinOneGap)
}

annotation <- annotation[SNPsinGap, ]
ArraysinGap <- genotypes[which(genotypes[,"ID"] %in% annotation[,"Probe.Set.ID"]),]

pop1 <- chipannotation[which(chipannotation[,"popPk"] == "IB77xx"), "ID_Chip"]
pop2 <- chipannotation[which(chipannotation[,"popPk"] == "IBNHxx"), "ID_Chip"]

write.table(ArraysinGap[,c("ID",pop1,pop2)], "Analysis/ArraysinGap.txt", sep="\t")
