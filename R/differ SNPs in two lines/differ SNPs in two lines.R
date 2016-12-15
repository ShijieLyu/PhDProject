# selecting SNPs which are different in two lines
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Feb, 2016
# first written Feb, 2016

setwd("D:/Chicken/DifferSNPsinTwoLines") 

### Differ SNPs between NHI and WL77 in 600K data
genotypes <- read.table("D:/Chicken/600KSNPchip/RawData/600kSNPgenotypes.txt", sep="\t", na.strings = "-1", header=TRUE, check.names = FALSE)                                 # Load the genotypes
annotation <- read.table("D:/Chicken/600KSNPchip/RawData/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", na.strings = "---", header=TRUE)                  # Load the annotation
chipannotation <- read.table("D:/Chicken/600KSNPchip/RawData/chickenNumber.txt", sep ="\t", header=TRUE, check.names = FALSE)

WLpop <- as.character(chipannotation[which(chipannotation[,"popPk"] == "IB77xx"), "ID_Chip"])
NHpop <- as.character(chipannotation[which(chipannotation[,"popPk"] == "IBNHxx"), "ID_Chip"])

differSNPs <- NULL
for(x in 1:nrow(genotypes)){                                  # Decide if genotypes[,WLpop] != genotypes[,NHpop]
  if(!any(genotypes[x,] == 1, na.rm=TRUE)){                   # 0 - No heterozygous animals
    tblPop1 <- table(as.numeric(genotypes[x, WLpop]))
    if(length(tblPop1) == 1){                                 # 1 - Check if genotypes[,WLpop] are consistent
      tblPop2 <- table(as.numeric(genotypes[x, NHpop]))
      if(length(tblPop2) == 1){                               # 2 - Check if genotypes[,NHpop] are consistent
        if(names(tblPop1) != names(tblPop2)){                 # 3 - Check if genotypes[,WLpop] != genotypes[,NHpop]
          differSNPs <- rbind(differSNPs, genotypes[x,])
        }
      }
    }
  }
}

dim(differSNPs)                                    
rownames(differSNPs) <- differSNPs[,"ID"]

write.table(differSNPs[,c(WLpop,NHpop)], "Analysis/differSNPs600K.txt", sep="\t")

### Annotate the differ SNPs

differSNPs <- read.table("Analysis/differSNPs600K.txt", sep="\t",header=TRUE,check.names = FALSE)
chipannotation <- read.table("D:/Chicken/600KSNPchip/RawData/chickenNumber.txt", sep ="\t", header=TRUE, check.names = FALSE)

WLpop <- paste0("WL",as.character(chipannotation[which(chipannotation[,"popPk"] == "IB77xx"), "ID_Chip"]))
NHpop <- paste0("NH",as.character(chipannotation[which(chipannotation[,"popPk"] == "IBNHxx"), "ID_Chip"]))

colnames(differSNPs) <- c(WLpop,NHpop)                                                                          # Add "WL" and "NH" sign in front of the ID number

annotation <- read.table("D:/Chicken/600KSNPchip/RawData/Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", na.strings = "---", header=TRUE)# Load the annotation
AnnoSNPs <- annotation[which(annotation[,"Probe.Set.ID"] %in% rownames(differSNPs)),]                                # The differ SNPs with annotation

SNPEnsemblInput <- cbind(AnnoSNPs[,c("Chromosome","Physical.Position","Position.End")],Allele = NA,Strand=AnnoSNPs[,"Strand"])

TakeChange <- function(Sequence){
  s <- grep("[[]",strsplit(as.character(Sequence),"")[[1]])
  e <- grep("[]]",strsplit(as.character(Sequence),"")[[1]])
  paste0(strsplit(as.character(Sequence),"")[[1]][(s+1):(e-1)],collapse = "")
}

for(x in 1:nrow(AnnoSNPs)){ 
  SNPEnsemblInput[x,"Allele"] <- TakeChange(AnnoSNPs[x,"Flank"])
}

write.table(SNPEnsemblInput,"Analysis/DifferSNP-VEP-EnsemblInput600K.txt", sep="\t",row.names = FALSE, quote=FALSE)

## Differ SNPs between NHI and WL77 in 60K Data
genotypes <- read.table("D:/Chicken/Rdata/SNPsinFMRegion/RawData/FullDataTable_3parents.txt", sep="\t", na.strings = "--", header=TRUE, check.names=FALSE)
annotation <- read.table("D:/Chicken/Rdata/SNPsinFMRegion/RawData/Chicken_60K_allResults.txt", header = TRUE, sep="\t",colClasses="character")
  
annotation[,"position..galGal4."] <- as.numeric(gsub(",","", annotation[,"position..galGal4."]))
annotation <- annotation[-which(is.na(annotation[,"position..galGal4."])),]

annotation <- annotation[which(annotation[,"SNP.ID"] %in% genotypes[,"Name"]),]
genotypes <- genotypes[which(genotypes[,"Name"] %in% annotation[,"SNP.ID"]),]

Reorder <- match(annotation[,"SNP.ID"], genotypes[,"Name"])
genotypes <- genotypes[Reorder,]

Ngenotypes <- cbind(Name = genotypes[,"Name"], genotypes[,7:ncol(genotypes)], Chr = annotation[,"chr..galGal4."], Position = annotation[,"position..galGal4."], flanking.seq = annotation[,"flanking.seq"])  # give genotypes correct location
    
Ngenotypes <- Ngenotypes[which(Ngenotypes[,"6929.Top Alleles"] == Ngenotypes[,"6954.Top Alleles"] & 
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
write.table(SNPsof60K, "Analysis/DifferSNPsof60K.txt", sep="\t")
  
### The format for input of the Biomart website
#for (x in 1:nrow(SNPEnsemblInput60K)){
  #cat(paste0(SNPEnsemblInput60K[x,"Chr"],":",SNPEnsemblInput60K[x,"Position"],":",SNPEnsemblInput60K[x,"Position"],":","1"), "\n", file = "Analysis/60KSNPinputBiomart.txt",append=TRUE)
#}

### Get the SNP name based on the location of the Chip info.
SNPEnsemblInput60K <- cbind(SNPsof60K[,c("Name","Chr","Position","Position")],SNPID = NA)

library("biomaRt")
mart <- useMart("ENSEMBL_MART_SNP","ggallus_snp")
for (x in 1:nrow(SNPEnsemblInput60K)){
  ID <-getBM(c("refsnp_id"), filters = c('chr_name','start','end'), values = list(SNPEnsemblInput60K[x,"Chr"],SNPEnsemblInput60K[x,"Position"],SNPEnsemblInput60K[x,"Position"]), mart=mart)
  cat(paste0(SNPEnsemblInput60K[x,"Name"]," ",as.character(ID)), "\n", file = "Analysis/60KSNPIDBiomart.txt",append=TRUE)
}

write.table(SNPEnsemblInput60K,"Analysis/DifferSNP-VEP-EnsemblInput60K-SNPIDs.txt", sep="\t",row.names = FALSE, quote=FALSE)

### Differ Mutations in NHI and WL77 of DNA-seq Data

SNPs_Called_DNAseq <- read.csv("D:/Chicken/Rdata/SNPsinFMRegion/RawData/output.raw.snps.indels.3samples.vcf",skip=59,header=TRUE,sep="\t")
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

NSNPs_Called_DNAseq <- SNPs_Called_DNAseq[which(SNPs_Called_DNAseq[,"GT6929"] == SNPs_Called_DNAseq[,"GT6954"] & 
                                            SNPs_Called_DNAseq[,"GT6929"] != SNPs_Called_DNAseq[,"GT8425"]),]
NSNPs_Called_DNAseq <- NSNPs_Called_DNAseq[c(which(NSNPs_Called_DNAseq[,"GT6929"]=="0/0"),which(NSNPs_Called_DNAseq[,"GT6929"]=="1/1"),which(NSNPs_Called_DNAseq[,"GT6929"]=="2/2")),]
NSNPs_Called_DNAseq <- NSNPs_Called_DNAseq[which(NSNPs_Called_DNAseq[,"GT6954"]==NSNPs_Called_DNAseq[,"GT6929"]),]
NSNPs_Called_DNAseq <- NSNPs_Called_DNAseq[c(which(NSNPs_Called_DNAseq[,"GT8425"]=="0/0"),which(NSNPs_Called_DNAseq[,"GT8425"]=="1/1"),which(NSNPs_Called_DNAseq[,"GT8425"]=="2/2")),]
SNPsofDNAseq <- NSNPs_Called_DNAseq

write.table(SNPsofDNAseq, "Analysis/SNPsofDNAseq.txt", sep="\t",row.names=FALSE)


SNPEnsemblInput6DNAseq <- cbind(SNPsofDNAseq[,c("X.CHROM","POS")],NoName1=".",SNPsofDNAseq[,c("REF","ALT")],NoName2=".",NoName3=".",NoName4=".")
write.table(SNPEnsemblInput6DNAseq,"Analysis/DifferSNP-VEP-EnsemblInputDNAseq.VCF", sep="\t",row.names = FALSE, quote=FALSE)

















