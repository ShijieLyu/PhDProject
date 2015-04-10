# F5 generation hetrozygousness analysis
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Feb, 2015
# first written Feb, 2015

setwd("D:/Chicken/60KSNPchip")
genotypes <- read.table("RawData/Full _Data_Table.txt", header=TRUE, sep="\t",na.strings = "--",check.names=FALSE, colClasses="character")
annotation <- read.table("RawData/Chicken_60K_allResults.txt", header = TRUE, sep="\t",colClasses="character")
F5chicken <- read.table("RawData/F5phenotypeDataforQTL.txt", header = TRUE, sep="\t",na.strings = ".")

annotation <- annotation[,c("SNP.ID","flanking.seq","chr..galGal4.","position..galGal4.")]
annotation[,"position..galGal4."] <- as.numeric(gsub(",","", annotation[,"position..galGal4."]))
annotation <- annotation[-which(is.na(annotation[,"position..galGal4."])),]

annotation <- annotation[which(annotation[,"SNP.ID"] %in% genotypes[,"Name"]),]
genotypes <- genotypes[which(genotypes[,"Name"] %in% annotation[,"SNP.ID"]),]

#Reorder <- match(annotation[,"SNP.ID"], genotypes[,"Name"])
#genotypes <- genotypes[Reorder,]


Ind_Nr <- NULL
for(x in 1:nrow(F5chicken)){
 PerInd_Nr <- paste0(F5chicken[x,"Tier.Nrn"], ".Top Alleles")    # Get the individual name in the genotyping file
 Ind_Nr <- c(Ind_Nr, PerInd_Nr)
}
genotypes <- genotypes[,c("Name",Ind_Nr)]

if(!file.exists("Analysis/Genotypes60KSNP.txt")) {
  for (x in 1:nrow(genotypes)) {
    for (y in 2:ncol(genotypes)) {
      if (length(unique(strsplit(genotypes[x,y],"")[[1]])) == 2) {
         genotypes[x,y] <- 1
      }else{
         if(is.na(genotypes[x,y])){genotypes[x,y] <- NA}
         else{genotypes[x,y] <- 0}
      }
    }
  }
write.table(genotypes, "D:/Chicken/60KSNPchip/Analysis/Genotypes60KSNP.txt", sep="\t")
}else{
  cat("Loading the genotypes file from disk\n")
  genotypes <- read.table("Analysis/Genotypes60KSNP.txt", sep="\t",header= TRUE)
}


testHetrozygousity <- function(snpdata, annotation, regionwidth = 100000){
  output <- NULL
  chromosomes <- as.character(unique(annotation[,"chr..galGal4."]))
  chromosomes <- chromosomes[c(1:29,258:261)]
  for(chromosome in chromosomes){
    chrAnnotation <- annotation[annotation[,"chr..galGal4."] == chromosome,]
    ChrSnps <- snpdata[which(snpdata[,"Name"] %in% chrAnnotation[,"SNP.ID"]),]
    maxlength <- max(chrAnnotation[,"position..galGal4."])
    regions <- cbind(start = seq(1,as.numeric(maxlength),regionwidth), end = seq(regionwidth,as.numeric(maxlength)+regionwidth,regionwidth))
    for(x in 1:nrow(regions)){
      snpInRegion <- chrAnnotation[which(as.numeric(chrAnnotation[,"position..galGal4."]) > regions[x,"start"] & as.numeric(chrAnnotation[,"position..galGal4."]) < regions[x,"end"]),"SNP.ID"]
      if(length(snpInRegion) > 0){   
        isHe <- sum(ChrSnps[which(ChrSnps[,"Name"] %in% snpInRegion),] == 1,na.rm=TRUE)
        isHo <- sum(ChrSnps[which(ChrSnps[,"Name"] %in% snpInRegion),] != 1,na.rm=TRUE)
      }else{
        isHe <- 0 ; isHo <- 1
      }
      output <- rbind(output, c(chromosome, regions[x,"start"], regions[x,"end"], round((isHe/isHo) * 100, d=2)))
    }
    cat("Done chromosome",chromosome,"\n")
  }
  colnames(output) <- c("chromosome","start","end","score")
  write.table(output, "D:/Chicken/60KSNPchip/Analysis/HetroRatio.txt", sep="\t")
  return(output)
}

permuteHetrozygousity <- function(snpdata, annotation, regionwidth = 100000, nperm = 100){
  chromosomes <- as.character(unique(annotation[,"chr..galGal4."]))
  chromosomes <- chromosomes[c(1:29,258:261)]
  thresholdperChr <- NULL
  for(chromosome in chromosomes){
    permutations <- NULL
    for(y in 1:nperm){
	  chrAnnotation <- annotation[annotation[,"chr..galGal4."] == chromosome,]
	  chrAnnotation[,"SNP.ID"] <- chrAnnotation[sample(nrow(chrAnnotation)), "SNP.ID"]
	  ChrSnps <- snpdata[which(snpdata[,"Name"] %in% chrAnnotation[,"SNP.ID"]),]
      maxlength <- max(chrAnnotation[,"position..galGal4."])
      regions <- cbind(start = seq(1,as.numeric(maxlength),regionwidth), end = seq(regionwidth,as.numeric(maxlength)+regionwidth,regionwidth))
	  output <- NULL
	  Max <- NULL
      for(x in 1:nrow(regions)){
	    snpInRegion <- chrAnnotation[which(as.numeric(chrAnnotation[,"position..galGal4."]) > regions[x,"start"] & as.numeric(chrAnnotation[,"position..galGal4."]) < regions[x,"end"]),"SNP.ID"]
        if(length(snpInRegion) > 0){
          isHe <- sum(ChrSnps[which(ChrSnps[,"Name"] %in% snpInRegion),] == 1,na.rm=TRUE)
          isHo <- sum(ChrSnps[which(ChrSnps[,"Name"] %in% snpInRegion),] != 1,na.rm=TRUE)
        }else{
          isHe <- 0 ; isHo <- 1
        }
        output <- rbind(output, c(chromosome, regions[x,"start"], regions[x,"end"], round((isHe/isHo) * 100, d=2)))
		Max <- max(as.numeric(output[,4]))
	  }
	  permutations <- c(permutations, Max)
    }
  permutations <- permutations[sort(permutations,index.return=TRUE)$ix]
  threshold <- permutations[round(length(permutations) * .95, d = 0)]
  thresholdperChr <- c(thresholdperChr,c(chromosome,threshold))
  cat("Done Threshold caculation chromosome",chromosome,"\n")
  }
  write.table(thresholdperChr, "D:/Chicken/60KSNPchip/Analysis/thresholdperChr.txt", sep="\t")
  return(thresholdperChr)
}

regionwidth = 250000
nperm       = 1000

realresults <- testHetrozygousity(genotypes, annotation, regionwidth)

permutations <- permuteHetrozygousity(genotypes, annotation, regionwidth, nperm)



