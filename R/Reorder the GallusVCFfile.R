# Re-order the Gallus_gallusVCF.file into a new order for DNA-desequencing data analysis 
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Jan, 2015
# first written Jan, 2015


setwd("D:/Chicken")
gallusVCF <- read.table("Gallus_gallus.vcf", header=FALSE, skip = 13, sep="\t")
chr <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 32, "Z", "W", "MT")

neworder <- NULL
for(x in chr){
    neworder <- rbind(neworder, gallusVCF[which(gallusVCF[,1] == x),])
}
header <- readLines("Gallus_gallus.vcf",n=13)

cat(header, file="Gallus_gallus.reordered.vcf", sep="\n")
write.table(neworder, "Gallus_gallus.reordered.vcf", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)


