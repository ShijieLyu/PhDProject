# SNP selection and filter
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Jan, 2015
# first written Jan, 2015

setwd("C:/Data/SNPselection") 
genes <- read.table("RawData/genename.txt", sep="\t", header=TRUE)                                 # load the gene which is selected for the further analysis
SNPsinbigregion <- read.table("RawData/finemappingSNPswithannotation.txt", sep="\t", header=TRUE)  # load the SNPs which was found using 600K SNP-Chip in the big region 

### Filter the SNPs, the aim is to select the target SNPs which no SNP is near to within 50bp.

if(!file.exists("Analysis/dbSNPsinbigregion.txt")){                                                # Obtain all the dbSNPs in the big region using biomart
  library("biomaRt")
    Mart <- useMart(biomart="snp", dataset="ggallus_snp")
    dbSNPs <- getBM(attributes = c("refsnp_id","chrom_start","chrom_end","chrom_strand","ensembl_gene_stable_id","sift_score","sift_prediction"), filters = "chromosomal_region", values = "4:59553432:84774762", mart = Mart)
    write.table(dbSNPs, row.names = FALSE, "Analysis/dbSNPsinbigregion.txt", sep="\t")
}else{                                                                                             # Or load them from Disk if we already have them
  cat("Loading all dbSNPs from disk\n")
  dbSNPs <- read.table("Analysis/dbSNPsinbigregion.txt",sep="\t", header=TRUE)
}

filteredSNPs <- NULL                                                                               # Test whether there is SNP near the target SNP within 50bp distance
for (x in 1:nrow(SNPsinbigregion)){
  beforeSNPstart <- SNPsinbigregion[x,"Physical.Position"]-1
  beforeSNPend <- SNPsinbigregion[x,"Physical.Position"]-60
  afterSNPstart <- SNPsinbigregion[x,"Physical.Position"]+1
  afterSNPend <- SNPsinbigregion[x,"Physical.Position"]+60
  if(!any(c(beforeSNPstart:beforeSNPend, afterSNPstart:afterSNPend) %in% dbSNPs[,"chrom_end"])){
    filteredSNPs <- rbind(filteredSNPs, SNPsinbigregion[x,c("dbSNP.RS.ID","Physical.Position")])
  }
}
filteredSNPswithinfo <-dbSNPs[which(dbSNPs[,"refsnp_id"] %in% filteredSNPs[,"dbSNP.RS.ID"]),]
write.table(filteredSNPswithinfo, row.names = FALSE, "Analysis/filteredSNPswithinfo.txt", sep="\t")  

### Select the SNPs which are in or near the target Genes
SNPsnearGenes <- NULL 
for (x in 1:nrow(filteredSNPswithinfo)){
  if(any(filteredSNPswithinfo[x,"ensembl_gene_stable_id"] %in% genes[,"Ensembl_gene_id"])){
    SNPsnearGenes <- rbind(SNPsnearGenes, filteredSNPswithinfo[x,])
  }
}
write.table(SNPsnearGenes, row.names = FALSE, "Analysis/SNPsnearGenes.txt", sep="\t")
