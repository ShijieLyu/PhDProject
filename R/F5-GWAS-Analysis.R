# Chicken F5 GWAS analysis using 60K SNP-Chip
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

### Remove the marker which has no segregating
OnlyOneGeno <- which(apply(genotypes,2,function(x){length(unique(x[which(!is.na(x))]))})==1)
OnlyNa <- which(apply(genotypes,2,function(x){length(unique(x[which(!is.na(x))]))})==0)
Ngenotypes <- genotypes[,-c(OnlyOneGeno,OnlyNa)]                                                  # Genotypes after filtering

### Association analysis
phenotypes <- cbind(phenotypes, parents=paste0(phenotypes[,"V"],phenotypes[,"M"]))
traits <- colnames(phenotypes)[9:74]    # colnume 9 to 74 are the measured phenotypes
markers <- colnames(Ngenotypes)[-1]     # the first colnume is the "ID", so remove it from the marker names.

if(!file.exists("Analysis/pvaluesGWASforF5chickens.txt")){
  pvalues <- matrix(NA, length(traits), length(markers), dimnames = list(traits, markers)) #options(warn=2)
  for (phe in traits){
    cat("Computing GWAS results for:", phe, "\n")
    pvalues[phe,] <- apply(Ngenotypes[,-1], 2, function(marker){
      idx <- which(!is.na(phenotypes[,phe]))
      tryCatch(res <- anova(lm(as.numeric(phenotypes[,phe][idx]) ~ as.factor(marker[idx])))[[5]][1], error = function(e){res <<- NA})
      return(res)
    })
  }
  write.table(pvalues, "Analysis/pvaluesGWASforF5chickens.txt", sep="\t")
}else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading GWAS pvalues from disk\n")
  pvalues <- read.table("Analysis/pvaluesGWASforF5chickens.txt",sep="\t", header=TRUE)
}
### Plot the QTL

MarkerInfo <- annotation[which(annotation[,"SNP.ID"] %in% markers),]
rownames(MarkerInfo) <- MarkerInfo[,1]
MarkerInfo <- MarkerInfo[markers,]
pvalues <- pvalues[,MarkerInfo[,"SNP.ID"]]

phe <- "BW20W_g"
chrcols <- 1+ (as.numeric(as.factor(MarkerInfo[,"chr..galGal4."])) %% 2)
plot(x = c(1, ncol(pvalues)), y = c(0,10), t='n', main=phe, xlab="LOD Score", ylab = "Markers")
points(x=1:ncol(pvalues),y=-log10(pvalues[phe,]), pch = 19, cex = 0.5, col=chrcols)
abline(h = -log10(0.1/ncol(pvalues)), col="gray")

### Generate a Q-Q plot

phe <- "BW20W_g"

observed <- sort(pvalues[phe,])
lobs <- -(log10(observed))

expected <- c(1:length(observed)) 
lexp <- -(log10(expected / (length(expected))))


#pdf("qqplot.pdf", width=6, height=6)
plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,7), ylim=c(0,7), las=1, xaxs="i", yaxs="i", bty="l")
points(lexp, lobs, pch=23, cex=.4, bg="black") 
#dev.off()







### Select the SNPs which have alternative allele in the target region on Chr4
SNPsof60K <- read.table("D:/Chicken/Rdata/SNPsinFMRegion/Analysis/SNPsof60K.txt",sep="\t", header=TRUE, check.names=FALSE)
OrderedSNP <- SNPsof60K[order(SNPsof60K[,"Position"]),]          # Order the SNP accoring to the position
AllMarkers <- genotypes[,as.character(OrderedSNP[,"Name"])]      # All the SNP markers with genotypes infoemation for the 31 F5 individuals
write.table(AllMarkers, "Analysis/AllMarkersinQTLregion.txt", sep="\t")








