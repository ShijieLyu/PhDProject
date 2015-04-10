# 5 Target SNPs sequence for Primer design +/- 100bp  
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Feb, 2015
# first written Feb, 2015

setwd("C:/Data/TargetSNPsSequence-Primer") 
SNPs <- read.table("5-More-SNP-ID.txt", header=TRUE)  # load the 5 Target SNPs info. (5 SNPs were selected based on the SNP selection results-- a total of 37, and the knowledge, literature)

# source("http://bioconductor.org/biocLite.R")
# biocLite("BSgenome.Ggallus.UCSC.galGal4")      # download the gallus4 genomics

suppressMessages(library(BSgenome.Ggallus.UCSC.galGal4))
SNPsSequence <- NULL
for (x in 1:nrow(SNPs)){
  chromosome <- SNPs[x,"Chr"]                    # the Chr Number
  startPos <- SNPs[x,"Location"] - 100           # Get +/- 100bp of the SNP for primer design
  endPos <- SNPs[x,"Location"] + 100
  msequence <- toString(subseq(Ggallus[[chromosome]], startPos, endPos))  # Get the sequence from the Genomics
  looseLetters <- strsplit(msequence,"")
  ref <- looseLetters[[1]][101]                  # Reference mutation
  alt <- as.character(SNPs[x,"Alt"])             # Altered mutation
  nsequence <- paste0(paste0(looseLetters[[1]][1:100], collapse=""),"[",ref,"/",alt,"]", paste0(looseLetters[[1]][102:201], collapse=""))  # Mark the SNP 
  SequenceperSNP <- c(as.character(SNPs[x,"SNP.ID"]),SNPs[x,"Location"], paste0("4:",startPos,":",endPos), nsequence)
  SNPsSequence <- rbind(SNPsSequence, SequenceperSNP)
}
colnames(SNPsSequence) <- c("SNP.ID", "SNP_Location", "SequenceLocation", "Sequence")
write.table(SNPsSequence,"5_More.Target_SNPs_Sequence.txt", row.names =FALSE, sep="\t")
