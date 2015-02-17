# Analysis of hetrozygousness in the 600k SNP chips
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Feb, 2015
# first written Feb, 2015

setwd("D:/Chicken/600KSNPchip/RawData")

arrayAnnotation <- read.table("Axiom_GW_GT_Chicken.na34.annot.csv", sep=",", header=TRUE, colClasses="character", na.string="---")
arrayAnnotation <- arrayAnnotation[-which(is.na(arrayAnnotation[,"Physical.Position"])),c("Probe.Set.ID","Chromosome", "Physical.Position")]

calldata <- read.table("600kSNPgenotypes.txt", header= TRUE, na.strings="-1")
calldata <- calldata[which(calldata[,"ID"] %in% arrayAnnotation[,"Probe.Set.ID"]),]

chipannotation <- read.table("chickenNumber.txt", sep ="\t", header=TRUE)
chipannotation[,"ID_Chip"] <- paste0("X", chipannotation[,"ID_Chip"]) 
pop1 <- chipannotation[which(chipannotation[,"popPk"] == "IB77xx"), "ID_Chip"]
pop2 <- chipannotation[which(chipannotation[,"popPk"] == "IBNHxx"), "ID_Chip"]

calldatapop1 <- calldata[,c("ID",pop1)]
calldatapop2 <- calldata[,c("ID",pop2)]

testHetrozygousity <- function(snpdata, annotation, regionwidth = 100000){
  output <- NULL
  chromosomes <- as.character(unique(annotation[,"Chromosome"]))
  for(chromosome in chromosomes){
    chrAnnotation <- annotation[annotation[,"Chromosome"] == chromosome,]
	ChrSnps <- snpdata[which(snpdata[,"ID"] %in% chrAnnotation[,"Probe.Set.ID"]),]
    maxlength <- max(chrAnnotation[,"Physical.Position"])
    regions <- cbind(start = seq(1,as.numeric(maxlength),regionwidth), end = seq(regionwidth,as.numeric(maxlength)+regionwidth,regionwidth))
    for(x in 1:nrow(regions)){
      snpInRegion <- chrAnnotation[which(as.numeric(chrAnnotation[,"Physical.Position"]) > regions[x,"start"] & as.numeric(chrAnnotation[,"Physical.Position"]) < regions[x,"end"]),"Probe.Set.ID"]
      if(length(snpInRegion) > 0){
        isHe <- sum(ChrSnps[which(ChrSnps[,"ID"] %in% snpInRegion),] == 1,na.rm=TRUE)
        isHo <- sum(ChrSnps[which(ChrSnps[,"ID"] %in% snpInRegion),] != 1,na.rm=TRUE)
      }else{
        isHe <- 0 ; isHo <- 1
      }
      output <- rbind(output, c(chromosome, regions[x,"start"], regions[x,"end"], round((isHe/isHo) * 100, d=2)))
    }
    cat("Done chromosome",chromosome,"\n")
  }
  colnames(output) <- c("chromosome","start","end","score")
  write.table(output, "D:/Chicken/600KSNPchip/Analysis/HetroRatio.txt", sep="\t")
  return(output)
}

permuteHetrozygousity <- function(snpdata, annotation, regionwidth = 100000, nperm = 100){
  chromosomes <- as.character(unique(annotation[,"Chromosome"]))
  thresholdperChr <- NULL
  for(chromosome in chromosomes){
    permutations <- NULL
    for(y in 1:nperm){
	  annotation[,"Probe.Set.ID"] <- annotation[sample(nrow(annotation)), "Probe.Set.ID"]    # Shuffle the ProbeSetIDs to new locations
      chrAnnotation <- annotation[annotation[,"Chromosome"] == chromosome,]
	  ChrSnps <- snpdata[which(snpdata[,"ID"] %in% chrAnnotation[,"Probe.Set.ID"]),]
      maxlength <- max(chrAnnotation[,"Physical.Position"])
      regions <- cbind(start = seq(1,as.numeric(maxlength),regionwidth), end = seq(regionwidth,as.numeric(maxlength)+regionwidth,regionwidth))
	  output <- NULL
	  Max <- NULL
      for(x in 1:nrow(regions)){
	    snpInRegion <- chrAnnotation[which(as.numeric(chrAnnotation[,"Physical.Position"]) > regions[x,"start"] & as.numeric(chrAnnotation[,"Physical.Position"]) < regions[x,"end"]),"Probe.Set.ID"]
        if(length(snpInRegion) > 0){
          isHe <- sum(ChrSnps[which(ChrSnps[,"ID"] %in% snpInRegion),] == 1,na.rm=TRUE)
          isHo <- sum(ChrSnps[which(ChrSnps[,"ID"] %in% snpInRegion),] != 1,na.rm=TRUE)
        }else{
          isHe <- 0 ; isHo <- 1
        }
        output <- rbind(output, c(chromosome, regions[x,"start"], regions[x,"end"], round((isHe/isHo) * 100, d=2)))
		Max <- max(output[,4])
	  }
	  permutations <- c(permutations, Max)
    }
  permutations <- permutations[sort(permutations,index.return=TRUE)$ix]
  threshold <- permutations[round(length(permutations) * .95, d = 0)]
  thresholdperChr <- c(thresholdperChr,c(chromosome,threshold))
  cat("Done Threshold caculation chromosome",chromosome,"\n")
  }
  write.table(thresholdperChr, "D:/Chicken/600KSNPchip/Analysis/thresholdperChr.txt", sep="\t")
  return(thresholdperChr)
}

regionwidth = 250000
nperm       = 1000

realresults <- testHetrozygousity(calldatapop1, arrayAnnotation, regionwidth)
permutations <- permuteHetrozygousity(calldatapop1, arrayAnnotation, regionwidth, nperm)