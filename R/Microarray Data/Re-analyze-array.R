# Re-analysis of the micro array data from high fat and low fat chickens
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Dec, 2014
# first written Dec, 2014

setwd("C:/Data/Array")
arrays <- read.table("RawData/arrays.txt", sep="\t", header=TRUE, colClasses = "character")

lowfat  <- arrays[which(arrays[,"fat"] == "low"),"filename"]
highfat <- arrays[which(arrays[,"fat"] == "high"),"filename"]

group1 <- paste0("lowfat",1:4)
group2 <- paste0("highfat",1:4)

if(!file.exists("Analysis/fatdata.txt")){
  fatdata <- NULL
  for(filename in c(lowfat,highfat)){
    mdata <- read.csv(paste0("RawData/", filename), sep = "\t", skip = 9, header = TRUE)
    fatdata <- cbind(fatdata, mdata$gProcessedSignal)
  }
  colnames(fatdata) <- c(group1, group2)

### QC
  cor(fatdata[,c(group1,group2)], method="spearman")
  plot(fatdata[,group1][,c(1,4)])
  plot(hclust(dist(t(fatdata))))

### Log2 and Quantile normalisation of the data
  fatdataLOG <- apply(fatdata, 2, log2)
  library(preprocessCore)
  fatdata <- normalize.quantiles(fatdataLOG)
  colnames(fatdata) <- colnames(fatdataLOG)    
  rownames(fatdata) <- rownames(fatdataLOG)

  fatdata <- cbind(mdata[,c("ProbeName", "Sequence")], fatdata)
  colnames(fatdata) <- c("ProbeName", "Sequence", group1, group2)
  write.table(fatdata, "Analysis/fatdata.txt", sep="\t")
}else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading fatdata from disk\n")
  fatdata <- read.table("Analysis/fatdata.txt",sep="\t", header=TRUE)
}
### Quality control !!!
plot(hclust(dist(t(fatdata[,c(group1,group2)]))))

### Create a fasta file, and blast to the genome, re annotate the data that we have
if(!file.exists("Analysis/probes.fasta")){
  cat("", file = "Analysis/probes.fasta") 
  for(x in 1:nrow(fatdata)){
    if(as.character(fatdata[x,"Sequence"]) != ""){
      cat(">", as.character(fatdata[x,"ProbeName"]), "\n", sep = "", file="Analysis/probes.fasta", append = TRUE)
      cat(as.character(fatdata[x,"Sequence"]), "\n", sep = "", file="Analysis/probes.fasta", append = TRUE)
    }
  }
}

### Blast all the probes, using the reference made by Danny, and this step was ran in Danny's Desktop.
# blastn -task blastn -query e:\Chicken\RNA\Arrays_for_fat-2011-01\Analysis\probes.fasta -db Gallus_gallus.Galgal4.74.dna.db -perc_identity 95 -outfmt 6 -evalue 0.1 -num_alignments 5 -out E:\Chicken\RNA\Arrays_for_fat-2011-01\Analysis\ProbeLocations.txt
# Here we will get the "ProbeLocations.txt"

### Call the Uniqueprobes (menas which has unique location after blast)

probes <- read.table("Analysis/ProbeLocations.txt", sep="\t", header=FALSE)
colnames(probes) <- c("ProbeName", "Chr", "Ident", "Length", "U1", "U2", "U3", "Match", "Start", "Stop", "evalue", "Score")

if(!file.exists("Analysis/Uniqueprobes.txt")){
  notunique <- which(probes[,"ProbeName"] %in% unique(probes[duplicated(probes[,"ProbeName"]),"ProbeName"]))  # Remove the probes which are not unique
  probes <- probes[-notunique,]
  write.table(probes, row.names = FALSE, "Analysis/Uniqueprobes.txt", sep="\t")
}else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading Unique probes information from disk\n")
  probes <- read.table("Analysis/Uniqueprobes.txt",sep="\t", header=TRUE)
}

### How many unique probes and genes in the QTL region

if(!file.exists("Analysis/UniqueprobesinQTLregion.txt")){
  startOn4 <- 59553432                                                              # Start of our region
  endOn4   <- 84774762                                                              # End of our region
  starts <- apply(probes[,c("Start","Stop")],1,min)
  stops <- apply(probes[,c("Start","Stop")],1,max)
  onChr4 <- which(probes[,"Chr"] == 4 & starts > startOn4 & stops < endOn4)
  probes <- probes[onChr4,]
  write.table(probes, row.names = FALSE, "Analysis/UniqueprobesinQTLregion.txt", sep="\t")                      # In this case, 545 unique probes were found
}else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading Unique probes in QTL region information from disk\n")
  probes <- read.table("Analysis/UniqueprobesinQTLregion.txt",sep="\t", header=TRUE)
}

if(!file.exists("Analysis/GenesinQTLregion.txt")){
  library("biomaRt")                                                                            # How many genes in the QTL region
  Mart <- useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
  genes <- getBM(attributes = c("chromosome_name","ensembl_gene_id","start_position","end_position","external_gene_name","description"), filters = "chromosomal_region", values = "4:59553432:84774762", mart = Mart)
  write.table(genes, row.names = FALSE, "Analysis/GenesinQTLregion.txt", sep="\t")                           # In this case, 276 genes in the QTL region were found 
}else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading Genes in QTL region information from disk\n")
  genes <- read.table("Analysis/GenesinQTLregion.txt",sep="\t", header=TRUE)
}

### Match the ProbesInfo to probes, and get the probeinfo (Sequence and gProcessedSignal)

if(!file.exists("Analysis/ProbesinGenes.txt")){
  ProbesInfo <- read.table("Analysis/fatdata.txt", sep="\t", header=TRUE) 
  ProbesNr <- match(probes[,"ProbeName"], ProbesInfo[,"ProbeName"])        # Match the Unique probes in QTL region to the original data to get the Nr. of the probes
  ProbeInformation <- ProbesInfo[ProbesNr, c("Sequence", group1, group2)]
  probes <- cbind(probes, ProbeInformation)                                # Give the Unique Probes the Sequence info.

### Get the info. between the probes and the genes
  starts <- apply(probes[,c("Start","Stop")],1,min)  # The start location of a probe
  stops <- apply(probes[,c("Start","Stop")],1,max)   # The end location of a probe

  ProbesinGenes <- NULL
  for (y in 1:nrow(genes)){                          # Every gene
    for (x in 1:nrow(probes)){                       # every probe
    if (starts[x] > genes[y,"start_position"] & stops[x] < genes[y,"end_position"]) # Every time, let all probes compare to the gene from 1st to last
      ProbesinGenes <- rbind(ProbesinGenes, cbind(probes[x,],genes[y,]))
    }
  }
  ProbesinGenes <- ProbesinGenes[,c("ProbeName", "Chr", "Sequence","Start", "Stop", "ensembl_gene_id", "start_position", "end_position","description","external_gene_name",group1,group2)]
  write.table(ProbesinGenes, row.names = FALSE, "Analysis/ProbesinGenes.txt", sep="\t") 
}else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading Probes in Genes from disk\n")
  ProbesinGenes <- read.table("Analysis/ProbesinGenes.txt",sep="\t", header=TRUE)
}

### difference comparison between the low-fat group and high-fat group
ProbesinGenes <- cbind(ProbesinGenes, P.value = rep(NA, nrow(ProbesinGenes)), Mean = rep(NA, nrow(ProbesinGenes)), Var = rep(NA, nrow(ProbesinGenes)), MeanxVar= rep(NA, nrow(ProbesinGenes)))
for(x in 1:nrow(ProbesinGenes)){
  ProbesinGenes[x,"P.value"] = t.test(ProbesinGenes[x,group1], ProbesinGenes[x,group2])$p.value
  ProbesinGenes[x,"Mean"] = mean(as.numeric(ProbesinGenes[x,c(group1,group2)]))
  ProbesinGenes[x,"Var"] = var(as.numeric(ProbesinGenes[x,c(group1,group2)]))
  ProbesinGenes[x,"MeanxVar"] = ProbesinGenes[x,"Mean"] * ProbesinGenes[x,"Var"]
}
write.table(ProbesinGenes, row.names = FALSE, "Analysis/ProbesinGenesTtest&Mean.txt", sep="\t") 

### Get the TOP25 Probes with high Means and make a plot  
### Since correlation between individuals is too high, and no significant genes are found, So, just select the TOP25 using the mean value

snpsbigregion <- read.table("C:/Codes_R_analysis_for_Phd/600KSNPchip/Analysis/finemappingSNPswithannotation.txt", sep = "\t", header=TRUE)  # load the full data
snpsmallregion <- read.table("C:/Codes_R_analysis_for_Phd/600KSNPchip/Analysis/finemappingSNPinsmallregion.txt", sep = "\t", header=TRUE)  # load the full data
 
SortMean <- sort(ProbesinGenes[,"Mean"],decreasing = TRUE, index.return = TRUE) 
SortProbesusingMeans <- ProbesinGenes[SortMean$ix[1:353],]                                        # here can change the Number  TOP 30 or 40 or...100
write.table(SortProbesusingMeans, row.names = FALSE, "Analysis/SortProbesusingMeans.txt", sep="\t") 

Bigregionstart <- min(snpsbigregion[,"Physical.Position"])
Bigregionend <- max(snpsbigregion[,"Physical.Position"])
Smallregionstart <- min(snpsmallregion[,"Physical.Position"])
Smallregionend <- max(snpsmallregion[,"Physical.Position"])

ablinestart <- seq(Bigregionstart, Bigregionend, by= 1000000)

regions <- cbind(seq(Bigregionstart,Bigregionend,1000000), seq(Bigregionstart+1000000,Bigregionend+1000000,1000000))
colnames(regions) <- c("Start","End")
regions <- cbind(regions, Gene = NA)
regions <- cbind(regions, Middle = NA)
regions <- cbind(regions, Description = NA)

pdf("Analysis/TheGenesinBigregion.pdf")   
plot(x=c(Bigregionstart,Bigregionend), y = c(-1,1.6), t = "n", xlab = "The Big Region for the QTL", yaxt='n')
abline(v=ablinestart, col="gray", lwd=0.5)

nG <- nrow(SortProbesusingMeans)
genesDraw <- NULL
for(x in 1:nG){
  middle <- (SortProbesusingMeans[x,"start_position"]+SortProbesusingMeans[x,"end_position"]) / 2
  inRegion <- which(regions[,"Start"] < middle & regions[,"End"] > middle)
  if(!(SortProbesusingMeans[x,"ensembl_gene_id"] %in% genesDraw)){
    if(is.na(regions[inRegion,"Gene"])){
      lines(x = c(SortProbesusingMeans[x,"start_position"], SortProbesusingMeans[x,"end_position"]), y = c(0.4, 0.4) - x/(nG*0.6),lwd=5)
      text(SortProbesusingMeans[x,"start_position"]+400000,0.36 - x/(nG*0.6),SortProbesusingMeans[x,"external_gene_name"], cex=0.7)
      regions[inRegion,"Gene"] <- SortProbesusingMeans[x,"external_gene_name"]
      regions[inRegion,"Middle"] <- middle
      regions[inRegion,"Description"] <- SortProbesusingMeans[x,"description"]
      genesDraw <- c(genesDraw, SortProbesusingMeans[x,"ensembl_gene_id"])
    }else{
      lines(x = c(SortProbesusingMeans[x,"start_position"], SortProbesusingMeans[x,"end_position"]), y = c(0.6, 0.6) + x/(nG*1.1),lwd=1)
      text(SortProbesusingMeans[x,"start_position"]+200000,0.61 + x/(nG*1.1),SortProbesusingMeans[x,"external_gene_name"], cex=0.2)
      genesDraw <- c(genesDraw, SortProbesusingMeans[x,"ensembl_gene_id"])
    }
  }
}

points(x= snpsbigregion[,"Physical.Position"] , y= rep(0.5, length(snpsbigregion[,"Physical.Position"])), pch="|")
points(x=c(Smallregionstart,Smallregionend), y=rep(0.5,2), pch = "|", col="red", cex=2)

dev.off()

mean(diff(as.numeric(regions[,"Middle"])),na.rm=T)

### Caculate the threshold and get the Differential Expression Probes
threshold <- 0.05 / length(unique(ProbesinGenes[,"ProbeName"]))
ProbesInGenesDE <- ProbesinGenes[which(ProbesinGenes[,"P.value"] < threshold),]
cat("Finally", nrow(ProbesInGenesDE),"DE Probes were found")

