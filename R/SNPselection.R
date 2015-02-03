# SNPs for fine mapping selection and filter 
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

setwd("C:/Data/SNPselection") 
genes <- read.table("RawData/genename.txt", sep="\t", header=TRUE)                                 # load the gene which is selected for the further analysis
SNPsinbigregion <- read.table("RawData/finemappingSNPswithannotation.txt", sep="\t", header=TRUE)  # load the SNPs which was found using 600K SNP-Chip in the big region 
SNPsinBigRegionDNAseq <- read.table("RawData/SNPsinBigRegionDNAseq.txt", sep="\t", header=TRUE)    # load the SNPs which was found using the DNA-resequencing Data
SNPs_Called_DNAseq <- read.csv("RawData/output.raw.snps.indels.3samples.vcf",skip=59,header=TRUE,sep="\t")

### Filter the SNPs, the aim is to select the target SNPs which no SNP is near to within 60bp to design the Allele specific Primer for genotyping
## 1st, get all the dbSNPs in big region in the database 

if(!file.exists("Analysis/dbSNPsinbigregion.txt")){                                                # Obtain all the dbSNPs in the big region using biomart
  library("biomaRt")
    Mart <- useMart(biomart="snp", dataset="ggallus_snp")
    dbSNPs <- getBM(attributes = c("refsnp_id","chrom_start","chrom_end","chrom_strand","ensembl_gene_stable_id","sift_score","sift_prediction"), filters = "chromosomal_region", values = "4:59553432:84774762", mart = Mart)
    write.table(dbSNPs, row.names = FALSE, "Analysis/dbSNPsinbigregion.txt", sep="\t")
}else{                                                                                             # Or load them from Disk if we already have them
  cat("Loading all dbSNPs from disk\n")
  dbSNPs <- read.table("Analysis/dbSNPsinbigregion.txt",sep="\t", header=TRUE)
}

## 2nd, get all the SNPs using the DNA-resequencing data which is the parents(6929, 6954 & 8425)

SNPs_Called_DNAseq <- cbind(SNPs_Called_DNAseq, GT6929 = rep(NA, nrow(SNPs_Called_DNAseq)),GT6954 = rep(NA, nrow(SNPs_Called_DNAseq)),GT8425 = rep(NA, nrow(SNPs_Called_DNAseq)))

Split6929 <- strsplit(as.character(SNPs_Called_DNAseq[,"X6929"]), ":")
Split6954 <- strsplit(as.character(SNPs_Called_DNAseq[,"X6954"]), ":")
Split8425 <- strsplit(as.character(SNPs_Called_DNAseq[,"X8425."]), ":")

for(x in 1:length(Split6929)){                                 # Split6954 has the same length as Split6929
  SNPs_Called_DNAseq[x,"GT6929"] <- substr(Split6929[x],4,6)   # Get the genotype in the strings, 4 and 6 are the start location and Stop location of the genotype in the strings
  SNPs_Called_DNAseq[x,"GT6954"] <- substr(Split6954[x],4,6)
  SNPs_Called_DNAseq[x,"GT8425"] <- substr(Split8425[x],4,6)  
}

SNPs_Called_DNAseq <- SNPs_Called_DNAseq[-which(SNPs_Called_DNAseq[,"FILTER"] == "LowQual"),c("X.CHROM","POS","ID","REF","ALT","GT6929","GT6954","GT8425")] # Remove the SNPs with low quality
SNPsinChr4 <- SNPs_Called_DNAseq[which(SNPs_Called_DNAseq[,"X.CHROM"] == "4"),]                               # Get the SNPs on Chr4

startOn4 <- 59553432                                                              # Start of our region
endOn4   <- 84774762                                                              # End of our region

SNPsinBigRegionDNAseq <- SNPsinChr4[which(SNPsinChr4[,"POS"] > startOn4 & SNPsinChr4[,"POS"] < endOn4),]  # The SNPs in our big region
SNPsinBigRegionDNAseq <- SNPsinChr4[-which(SNPsinChr4[,"GT6929"] == "0/0" & SNPsinChr4[,"GT6954"] == "0/0" & SNPsinChr4[,"GT8425"] == "0/0"),]

## Based on the dbSNPs and DNA-seq-SNPs in the big region, filter and select the target SNPs

filteredSNPs <- NULL                                                                               # Test whether there is SNP near the target SNP within 60bp distance
for (x in 1:nrow(SNPsinbigregion)){
  beforeSNPstart <- SNPsinbigregion[x,"Physical.Position"]-1
  beforeSNPend <- SNPsinbigregion[x,"Physical.Position"]-60
  afterSNPstart <- SNPsinbigregion[x,"Physical.Position"]+1
  afterSNPend <- SNPsinbigregion[x,"Physical.Position"]+60
  if(!any(c(beforeSNPstart:beforeSNPend, afterSNPstart:afterSNPend) %in% c(dbSNPs[,"chrom_end"],SNPsinBigRegionDNAseq[,"POS"]))){
    filteredSNPs <- rbind(filteredSNPs, SNPsinbigregion[x,c("dbSNP.RS.ID","Physical.Position")])
  }
}
filteredSNPswithinfo <-dbSNPs[which(dbSNPs[,"refsnp_id"] %in% filteredSNPs[,"dbSNP.RS.ID"]),]
write.table(filteredSNPswithinfo, row.names = FALSE, "Analysis/filteredSNPswithinfoDNAseqfiltered.txt", sep="\t")  

### Select the SNPs which are in or near the target Genes, 5 interesting genes which selected according the microarry and gene description
SNPsnearGenes <- NULL 
for (x in 1:nrow(filteredSNPswithinfo)){
  if(any(filteredSNPswithinfo[x,"ensembl_gene_stable_id"] %in% genes[,"Ensembl_gene_id"])){
    SNPsnearGenes <- rbind(SNPsnearGenes, filteredSNPswithinfo[x,])
  }
}
write.table(SNPsnearGenes, row.names = FALSE, "Analysis/SNPsnearGenesDNAseqFiltered.txt", sep="\t")
