# Using the DNA-reseq , 600K and 60k to select the interesting SNPs
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Jun, 2015
# first written Jun, 2015

setwd("D:/Chicken/Rdata/SNPsinFMRegion")                              

###### The SNPs of 600K-chip in new 95% CI

if(!file.exists("Analysis/SNPsof600K.txt")){

  genotypes <- read.table("RawData/600kSNPgenotypes.txt", sep="\t", na.strings = "-1", header=TRUE, check.names=FALSE)              # Load the genotypes
  annotation <- read.table("RawData/Axiom_GW_GT_Chicken.na35.annot.csv", sep=",", na.strings = "---", header=TRUE)                    # Load the annotation
  chipannotation <- read.table("RawData/chickenNumber.txt", sep ="\t", header=TRUE)

  startOn4 <- 73629630    # Start of our region
  endOn4   <- 77032597    # End of our region

  inChr4Region    <- which(annotation[,"Chromosome"] == 4 &                                                 # On chromosome 4
                         as.numeric(as.character(annotation[,"Physical.Position"])) > startOn4 &            # Larger then our start
                         as.numeric(as.character(annotation[,"Physical.Position"])) < endOn4)               # Smaller then our end`
  annotation      <- annotation[inChr4Region, ]                                                             # Only use the annotation in our region

  genotypesOn4    <- which(genotypes[,"ID"] %in% annotation[,"Probe.Set.ID"])                               # Which genotypes do we have ?
  genotypes <- genotypes[genotypesOn4, ]                                                                    # ONly use the ones which are in our region

  pop1 <- as.character(chipannotation[which(chipannotation[,"popPk"] == "IB77xx"), "ID_Chip"])  # pop1 is the WL77
  pop2 <- as.character(chipannotation[which(chipannotation[,"popPk"] == "IBNHxx"), "ID_Chip"])  # pop2 is the NHI

  snpsForFinemapping <- NULL
  for(x in 1:nrow(genotypes)){                                  # Decide if genotypes[,pop1] != genotypes[,pop2]
    if(!any(genotypes[x,] == 1, na.rm=TRUE)){                   # 0 - No heterozygous animals
      tblPop1 <- table(as.numeric(genotypes[x, pop1]))
      if(length(tblPop1) == 1){                                 # 1 - Check if genotypes[,pop1] are consistent
        tblPop2 <- table(as.numeric(genotypes[x, pop2]))
        if(length(tblPop2) == 1){                               # 2 - Check if genotypes[,pop2] are consistent
          if(names(tblPop1) != names(tblPop2)){                 # 3 - Check if genotypes[,pop1] != genotypes[,pop2]
            snpsForFinemapping <- rbind(snpsForFinemapping, genotypes[x,])
          }
        }
      }
    }
  }
  rownames(snpsForFinemapping) <- snpsForFinemapping[,"ID"]
  snpsForFinemapping <- snpsForFinemapping[,c(pop1,pop2)]
  targetSNPsNr <- which(annotation[,"Probe.Set.ID"] %in% rownames(snpsForFinemapping))                                             # Give the Fine Mapping SNPS the annotation
  finemappingSNPswithannotation <- annotation[targetSNPsNr,]
  SNPsof600K <- cbind(finemappingSNPswithannotation, snpsForFinemapping)
  write.table(SNPsof600K, "Analysis/SNPsof600K.txt", sep="\t")
}else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading SNPs of 600K-Chip from disk\n")
  SNPsof600K <- read.table("Analysis/SNPsof600K.txt",sep="\t", header=TRUE, check.names=FALSE)
}

###### The SNPs of 60K-chip in new 95% CI

if(!file.exists("Analysis/SNPsof60K.txt")){
  genotypes <- read.table("RawData/FullDataTable_3parents.txt", sep="\t", na.strings = "--", header=TRUE, check.names=FALSE)
  annotation <- read.table("RawData/Chicken_60K_allResults.txt", header = TRUE, sep="\t",colClasses="character")
  
  annotation[,"position..galGal4."] <- as.numeric(gsub(",","", annotation[,"position..galGal4."]))
  annotation <- annotation[-which(is.na(annotation[,"position..galGal4."])),]

  annotation <- annotation[which(annotation[,"SNP.ID"] %in% genotypes[,"Name"]),]
  genotypes <- genotypes[which(genotypes[,"Name"] %in% annotation[,"SNP.ID"]),]

  Reorder <- match(annotation[,"SNP.ID"], genotypes[,"Name"])
  genotypes <- genotypes[Reorder,]

  Ngenotypes <- cbind(Name = genotypes[,"Name"], genotypes[,7:ncol(genotypes)], Chr = annotation[,"chr..galGal4."], Position = annotation[,"position..galGal4."], flanking.seq = annotation[,"flanking.seq"])  # give genotypes correct location
  
  startCIon4 <- 73629630    # Start of our region
  endCIon4   <- 77032597    # End of our region
    
  Ngenotypes <- Ngenotypes[which(Ngenotypes[,"Chr"] == 4 & as.numeric(as.character(Ngenotypes[,"Position"])) > startCIon4 &            # Larger then our start
                                 as.numeric(as.character(Ngenotypes[,"Position"])) < endCIon4 &    
                                 Ngenotypes[,"6929.Top Alleles"] == Ngenotypes[,"6954.Top Alleles"] & 
                                 Ngenotypes[,"8425.Top Alleles"] != Ngenotypes[,"6929.Top Alleles"]),]

  for (x in 1:nrow(Ngenotypes)){                                                                                                       # Label heterozygous as NA
    if (length(unique(strsplit(as.character(Ngenotypes[x,"6929.Top Alleles"]),"")[[1]])) == 2) {Ngenotypes[x,"6929.Top Alleles"] <- NA}
  
  }
  for (x in 1:nrow(Ngenotypes)){
    if (length(unique(strsplit(as.character(Ngenotypes[x,"6954.Top Alleles"]),"")[[1]])) == 2) {Ngenotypes[x,"6954.Top Alleles"] <- NA}
  
  }
  for (x in 1:nrow(Ngenotypes)){
    if (length(unique(strsplit(as.character(Ngenotypes[x,"8425.Top Alleles"]),"")[[1]])) == 2) {Ngenotypes[x,"8425.Top Alleles"] <- NA}
  
  }
  Ngenotypes <- Ngenotypes[which(!is.na(Ngenotypes[,"6929.Top Alleles"])),]
  Ngenotypes <- Ngenotypes[which(!is.na(Ngenotypes[,"6954.Top Alleles"])),]
  Ngenotypes <- Ngenotypes[which(!is.na(Ngenotypes[,"8425.Top Alleles"])),]
  SNPsof60K <- Ngenotypes
  write.table(SNPsof60K, "Analysis/SNPsof60K.txt", sep="\t")
  
}else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading SNPs of 60K-Chip from disk\n")
  SNPsof60K <- read.table("Analysis/SNPsof60K.txt",sep="\t", header=TRUE, check.names=FALSE)
}

## The SNPs of the DNA-resequencing data in new 95% CI

if(!file.exists("Analysis/SNPsofDNAseq.txt")){
  SNPs_Called_DNAseq <- read.csv("RawData/output.raw.snps.indels.3samples.vcf",skip=59,header=TRUE,sep="\t")
  SNPs_Called_DNAseq <- cbind(SNPs_Called_DNAseq, DP = rep(NA, nrow(SNPs_Called_DNAseq)), GT6929 = rep(NA, nrow(SNPs_Called_DNAseq)),GT6954 = rep(NA, nrow(SNPs_Called_DNAseq)),GT8425 = rep(NA, nrow(SNPs_Called_DNAseq)))

  Split6929 <- strsplit(as.character(SNPs_Called_DNAseq[,"X6929"]), ":")
  Split6954 <- strsplit(as.character(SNPs_Called_DNAseq[,"X6954"]), ":")
  Split8425 <- strsplit(as.character(SNPs_Called_DNAseq[,"X8425."]), ":")

  for(x in 1:length(Split6929)){                                 # Split6954 has the same length as Split6929
    SNPs_Called_DNAseq[x,"GT6929"] <- substr(Split6929[x],4,6)   # Get the genotype in the strings, 4 and 6 are the start location and Stop location of the genotype in the strings
    SNPs_Called_DNAseq[x,"GT6954"] <- substr(Split6954[x],4,6)
    SNPs_Called_DNAseq[x,"GT8425"] <- substr(Split8425[x],4,6)
    SNPs_Called_DNAseq[x,"DP"] <- sub(".*?DP=(.*?);.*", "\\1", as.character(SNPs_Called_DNAseq[x,"INFO"]))
  }
    
  SNPs_Called_DNAseq <- SNPs_Called_DNAseq[-which(SNPs_Called_DNAseq[,"FILTER"] == "LowQual"),c("X.CHROM","POS","ID","REF","ALT","QUAL","DP","GT6929","GT6954","GT8425")] # Remove the SNPs with low quality
  SNPsinChr4 <- SNPs_Called_DNAseq[which(SNPs_Called_DNAseq[,"X.CHROM"] == "4"),]                               # Get the SNPs on Chr4

  startOn4 <- 73629630                                                              # Start of our region
  endOn4   <- 77032597                                                              # End of our region

  SNPsinBigRegion <- SNPsinChr4[which(SNPsinChr4[,"POS"] > startOn4 & SNPsinChr4[,"POS"] < endOn4),]
  NSNPsinBigRegion <- SNPsinBigRegion[which(SNPsinBigRegion[,"GT6929"] == SNPsinBigRegion[,"GT6954"] & 
                                            SNPsinBigRegion[,"GT6929"] != SNPsinBigRegion[,"GT8425"]),]
  NSNPsinBigRegion <- NSNPsinBigRegion[c(which(NSNPsinBigRegion[,"GT6929"]=="0/0"),which(NSNPsinBigRegion[,"GT6929"]=="1/1"),which(NSNPsinBigRegion[,"GT6929"]=="2/2")),]
  NSNPsinBigRegion <- NSNPsinBigRegion[which(NSNPsinBigRegion[,"GT6954"]==NSNPsinBigRegion[,"GT6929"]),]
  NSNPsinBigRegion <- NSNPsinBigRegion[c(which(NSNPsinBigRegion[,"GT8425"]=="0/0"),which(NSNPsinBigRegion[,"GT8425"]=="1/1"),which(NSNPsinBigRegion[,"GT8425"]=="2/2")),]
  SNPsofDNAseq <- NSNPsinBigRegion
  write.table(SNPsofDNAseq, "Analysis/SNPsofDNAseq.txt", sep="\t",row.names=FALSE)
  
}else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading SNPs of DNA-seq from disk\n")
  SNPsofDNAseq <- read.table("Analysis/SNPsofDNAseq.txt",sep="\t", header=TRUE)
}

# Merge the three types data with the same format
SNPname60k <- as.character(SNPsof60K[,"Name"])
Nname <- NULL
for (x in 1: length(SNPname60k)){
  Nname <- c(Nname, substr(SNPname60k[x],5,length(strsplit(SNPname60k[x],"")[[1]])))
}
SNPsof60K <- cbind(SNPsof60K,Nname)
SNPsof60K <- SNPsof60K[-which(SNPsof60K[,"Nname"] %in% SNPsof600K[,"dbSNP.RS.ID"]),]
BriefSNPsof60K <- SNPsof60K[,c("Nname","Position","6929.Top Alleles","8425.Top Alleles")]
BriefSNPsof60K <- cbind(BriefSNPsof60K, WL77_Allele = rep(NA, nrow(BriefSNPsof60K)), NHI_Allele = rep(NA, nrow(BriefSNPsof60K)))
for (x in 1:nrow(BriefSNPsof60K)){
  BriefSNPsof60K[x,"WL77_Allele"] <- unique(strsplit(as.character(BriefSNPsof60K[x,"6929.Top Alleles"]),"")[[1]])
  BriefSNPsof60K[x,"NHI_Allele"] <- unique(strsplit(as.character(BriefSNPsof60K[x,"8425.Top Alleles"]),"")[[1]])
}
BriefSNPsof60K <- BriefSNPsof60K[,c("Nname","Position","WL77_Allele","NHI_Allele")]
colnames(BriefSNPsof60K) <- c("SNP-ID","Position","WL77_Allele","NHI_Allele")

SNPsof600K <- cbind(SNPsof600K, WL77_Allele = rep(NA, nrow(SNPsof600K)), NHI_Allele = rep(NA, nrow(SNPsof600K)))
for (x in 1:nrow(SNPsof600K)){
  if (SNPsof600K[x,"5149"] == 0) {SNPsof600K[x,"WL77_Allele"] = as.character(SNPsof600K[x,"Allele.A"]); SNPsof600K[x,"NHI_Allele"] = as.character(SNPsof600K[x,"Allele.B"])}
  if (SNPsof600K[x,"5149"] == 2) {SNPsof600K[x,"WL77_Allele"] = as.character(SNPsof600K[x,"Allele.B"]); SNPsof600K[x,"NHI_Allele"] = as.character(SNPsof600K[x,"Allele.A"])}
}
BriefSNPsof600K <- SNPsof600K[,c("dbSNP.RS.ID","Physical.Position","WL77_Allele","NHI_Allele")]
colnames(BriefSNPsof600K) <- c("SNP-ID","Position","WL77_Allele","NHI_Allele")

SNPsofDNAseq <- cbind(SNPsofDNAseq, WL77_Allele = rep(NA, nrow(SNPsofDNAseq)), NHI_Allele = rep(NA, nrow(SNPsofDNAseq)))
for (x in 1:nrow(SNPsofDNAseq)){
  if (SNPsofDNAseq[x,"GT6929"] == "0/0") {SNPsofDNAseq[x,"WL77_Allele"] = as.character(SNPsofDNAseq[x,"REF"]); SNPsofDNAseq[x,"NHI_Allele"] = as.character(SNPsofDNAseq[x,"ALT"])}
  if (SNPsofDNAseq[x,"GT6929"] == "1/1") {SNPsofDNAseq[x,"WL77_Allele"] = as.character(SNPsofDNAseq[x,"ALT"]); SNPsofDNAseq[x,"NHI_Allele"] = as.character(SNPsofDNAseq[x,"REF"])}
  if (SNPsofDNAseq[x,"GT6929"] == "2/2" && SNPsofDNAseq[x,"GT8425"] == "1/1") {SNPsofDNAseq[x,"WL77_Allele"] = strsplit(as.character(SNPsofDNAseq[x,"ALT"]),",")[[1]][2]; SNPsofDNAseq[x,"NHI_Allele"] = strsplit(as.character(SNPsofDNAseq[x,"ALT"]),",")[[1]][1]}
  if (SNPsofDNAseq[x,"GT6929"] == "2/2" && SNPsofDNAseq[x,"GT8425"] == "0/0") {SNPsofDNAseq[x,"WL77_Allele"] = strsplit(as.character(SNPsofDNAseq[x,"ALT"]),",")[[1]][2]; SNPsofDNAseq[x,"NHI_Allele"] = as.character(SNPsofDNAseq[x,"REF"])}
}
BriefSNPsofDNAseq <- SNPsofDNAseq[,c("ID","POS","WL77_Allele","NHI_Allele")]
colnames(BriefSNPsofDNAseq) <- c("SNP-ID","Position","WL77_Allele","NHI_Allele")

MergeSNPs <- rbind(BriefSNPsofDNAseq, BriefSNPsof60K, BriefSNPsof600K)

# Add the SNP ID for the SNP with no ID
if(!file.exists("Analysis/MergeSNPswithRefreshID.txt")){
  library("biomaRt")                                                                                              # How many genes in different regions for different traits
  mart <- useMart(biomart="snp", dataset="ggallus_snp")
  pos <- MergeSNPs[,"Position"]
  UpID <- NULL
  for (x in 1:length(pos)){
    ID <-getBM(c("refsnp_id"), filters = c('chr_name','start','end'), values = list("4",pos[x]-1,pos[x]+1), mart=mart)
    UpID <- c(UpID,as.character(ID))
  }
  NMergeSNPs <- cbind(MergeSNPs,RefreshID = UpID)
  write.table(NMergeSNPs, "Analysis/MergeSNPswithRefreshID.txt", sep="\t",row.names=FALSE)
  }else{                                                                                                        # Or load them from Disk if we already have them
  cat("See the results in analysis folder\n")
}





