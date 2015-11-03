# Analysis of heterozygous genotypes
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Mar, 2015
# first written Feb, 2015

scan.heterozygous <- function(genotypes, map, chrinfo, window.size = 750000, step.size = window.size/2, nperms = 10, verbose = FALSE){
  he <- apply(genotypes, 1, function(x){ sum(x == "H", na.rm=TRUE) })                                                                 # Calculate the number of heterozygous alleles
  ho <- apply(genotypes, 1, function(x){ sum(x != "H", na.rm=TRUE) })                                                                 # Calculate the number of homozygous alleles

  bins <- NULL
  for(chr in chrinfo[,"Chr"]){
    chr.length <- chrinfo[chrinfo[,"Chr"] == chr, "Length"]                                                                           # Length of the chromosome
    chr.bins <- cbind(Start = seq(1, chr.length, step.size), Stop = seq(1, chr.length, step.size) + window.size)                      # Create our bins
    smap <- map[map[,"Chr"] == chr, ]                                                                                                 # Create a copy of the map (only this chromosome)

    stats <- t(apply(chr.bins, 1, function(x){
      chr.markers <- rownames(smap)[as.numeric(smap[,"Position"]) > x["Start"] & as.numeric(smap[,"Position"]) < x["Stop"]]           # Which markers are in my current bin
      return(c(length(chr.markers), sum(he[chr.markers]) / sum(ho[chr.markers])))                                                     # Calculate HE / HO ratio
    }))
    chr.bins <- cbind(chr.bins, nSNPs = stats[, 1], Score = stats[, 2])                                                               # Remember the scores

    for(perm in 1:nperms){                                                                                                            # permute this chromosome
      permstats <- t(apply(chr.bins, 1, function(x){
        if(as.numeric(x["nSNPs"]) != 0){                                                                                              # If there are SNPs in the bin
          chr.markers <- rownames(smap)[sample(length(rownames(smap)), as.numeric(x["nSNPs"]))]                                       # Draw N SNPs at random from this chromosome
          return(c(length(chr.markers), sum(he[chr.markers]) / sum(ho[chr.markers])))                                                 # Calculate HE / HO ratio
        }else{
          return(c(0, NaN))                                                                                                           # No SNPs, just return no score
        }
      }))
      chr.bins <- cbind(chr.bins, permstats[, 2])                                                                                     # Add the scores observed during permutation
    }
    if(verbose) cat(paste0("Done with chromosome ", chr, "/", length(chrinfo[,"Chr"]),"\n"))
    chr.bins <- cbind(Chr = chr, chr.bins)
    bins <- rbind(bins, chr.bins)
  }
  colnames(bins)[which(colnames(bins) == "")] <- paste0("P" , 1:length(which(colnames(bins) == "")))

  invisible(return(list(results = bins[,-grep("^P", colnames(bins))], permutations = apply(bins[,grep("^P", colnames(bins))], 2, as.numeric))))
}

test.scan.heterozygous <- function(){
  setwd("D:/Github/heterozygous")
  chrinfo     <- read.table("chromosomeinfo.txt", colClasses=c("character", "numeric"), header=TRUE)                                # Chromosome information
  genotypes   <- read.table("genotypes.txt", sep="\t", check.names=FALSE, header=TRUE)
  map         <- read.table("map.txt", sep="\t", colClasses=c("character"), header=TRUE)                                            # Chromosome information

  results     <- scan.heterozygous(genotypes, map, chrinfo, nperms = 5)
}