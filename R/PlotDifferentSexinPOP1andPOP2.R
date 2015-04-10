# Get the significant region with high Hetro/Homo ratio and plot it in different sex of POP1,POP2
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Feb, 2015
# first written Feb, 2015


setwd("C:/Data/600KSNPchip")

### Load the POP1() and POP2() Analysis Data -- Ratio and Threshold

POP1 <- read.table("Analysis/HetroRatioPOP1-100.txt", header=TRUE)    # load the Population 1 Hetro Ratio POP1=WL77
POP2 <- read.table("Analysis/HetroRatioPOP2-100.txt", header=TRUE)    # load the Population 2 Hetro Ratio POP2=NHI
POP1Thr <- read.table("Analysis/thresholdperChr-POP1-26.02.15.txt", header=TRUE)   # Load the threshold for Population 1 under 1000 times permutation
POP2Thr <- read.table("Analysis/thresholdperChr-POP2-26.02.15.txt", header=TRUE)   # Load the threshold for Population 2 under 1000 times permutation
pop1sex0 <- read.table("Analysis/SexHetro/HetroRatio-pop1sex0.txt",header = TRUE)
pop1sex1 <- read.table("Analysis/SexHetro/HetroRatio-pop1sex1.txt",header = TRUE)
pop2sex0 <- read.table("Analysis/SexHetro/HetroRatio-pop2sex0.txt",header = TRUE)
pop2sex1 <- read.table("Analysis/SexHetro/HetroRatio-pop2sex1.txt",header = TRUE)
pop1sex0Thr <- read.table("Analysis/SexHetro/thresholdperChr-pop1sex0-100.txt",header = TRUE)
pop1sex1Thr <- read.table("Analysis/SexHetro/thresholdperChr-pop1sex1-100.txt",header = TRUE)
pop2sex0Thr <- read.table("Analysis/SexHetro/thresholdperChr-pop2sex0-100.txt",header = TRUE)
pop2sex1Thr <- read.table("Analysis/SexHetro/thresholdperChr-pop2sex1-100.txt",header = TRUE)

### Re-organised the Threshold format

ThresholdReformat <- function(POP1Thr){
  POP1matrix <- as.matrix(1:32)                                         # The format of the POP1Thr is nor correct, here I re-organised it
  POP1ThrChr <- as.character(POP1Thr[seq(1,64,2),])
  POP1ThrScore <- as.numeric(as.character(POP1Thr[seq(2,64,2),]))
  rownames(POP1matrix) <- POP1ThrChr
  POP1matrix[,1] <- POP1ThrScore
  POP1Thr <- POP1matrix
  colnames(POP1Thr) <- c("threshold")
  return(POP1Thr)
}

POP1Thr <-ThresholdReformat(POP1Thr)
POP2Thr <-ThresholdReformat(POP2Thr)
pop1sex0Thr <-ThresholdReformat(pop1sex0Thr)
pop1sex1Thr <-ThresholdReformat(pop1sex1Thr)
pop2sex0Thr <-ThresholdReformat(pop2sex0Thr)
pop2sex1Thr <-ThresholdReformat(pop2sex1Thr)

###Select the ratio which higher than the threshold for each Population

TargetRegion <- function(POP1,POP1Thr) {                             
  chromosomes <- as.character(unique(POP1[,"chromosome"]))
  Region <- NULL
  for (chr in chromosomes){ 
    ScorePerChr <- POP1[which(POP1[,"chromosome"] == chr),]
    thresholdPerChr <- as.numeric(POP1Thr[which(rownames(POP1Thr) == chr),"threshold"])
    TargetRegionPerchr <- ScorePerChr[which(as.numeric(ScorePerChr[,"score"]) > thresholdPerChr),]
    Region <- rbind(Region, TargetRegionPerchr)
  }
  return(Region)
}

TargetRegionPOP1 <- TargetRegion(POP1,POP1Thr)
TargetRegionPOP2 <- TargetRegion(POP2,POP2Thr)
TargetRegionpop1sex0 <- TargetRegion(pop1sex0,pop1sex0Thr)
TargetRegionpop1sex1 <- TargetRegion(pop1sex1,pop1sex1Thr)
TargetRegionpop2sex0 <- TargetRegion(pop2sex0,pop2sex0Thr)
TargetRegionpop2sex1 <- TargetRegion(pop2sex0,pop2sex1Thr)

### Plot chromosomes and add the region information to it---VERITAL

## 1. Plot the POP1 and different Sex in POP1

chrinfo <- read.table("C:/Data/600KSNPchip/RawData/chromosomeinfo.txt", header=TRUE)     # Load the Chromosome information

chromosomes <- chrinfo[,"Chromosome"]
mlength <- max(chrinfo[,"Length"]) * 1.10              

pdf("Analysis/POP1DifferentSexHetro.pdf")                 
plot(y=c(0,  mlength), x=c(0, nrow(chrinfo)) + 0.5, t='n', yaxt="n", xlab="Chromosome", ylab="Location (Mbp)", xaxt="n")     # Make a frame

invisible(apply(chrinfo, 1, function(mrow){                          # Plot the Chromosomes
  xloc <- which(chromosomes == mrow["Chromosome"])
  lines(c(xloc, xloc),y=c(0, mrow["Length"]), lwd=4)
}))

invisible(apply(TargetRegionPOP1, 1, function(mrow){                 # Plot the Significant region for population 1  
  xloc <- which(chromosomes == mrow["chromosome"])
  if(length(xloc) > 0) lines(c(xloc, xloc), y=c(mrow["start"], mrow["end"]), t="l", lwd=4, col="green",lty=1)
  #return(as.numeric(mrow["start"]) + abs(as.numeric(mrow["start"]) - as.numeric(mrow["end"])) / 2)
}))

invisible(apply(TargetRegionpop1sex0, 1, function(mrow){                 # Plot the Significant region for pop1sex0
  xloc <-which(chromosomes == mrow["chromosome"])
  if(length(xloc) > 0) lines(c(xloc, xloc)-0.25, y=c(mrow["start"], mrow["end"]), t="l", lwd=4, col="red",lty=1)
}))

invisible(apply(TargetRegionpop1sex1, 1, function(mrow){                 # Plot the Significant region for pop1sex1
  xloc <- which(chromosomes == mrow["chromosome"])
  if(length(xloc) > 0) lines(c(xloc, xloc)+0.25, y=c(mrow["start"], mrow["end"]), t="l", lwd=4, col="blue",lty=1)
  #return(as.numeric(mrow["start"]) + abs(as.numeric(mrow["start"]) - as.numeric(mrow["end"])) / 2)
}))

axis(1, at=1:nrow(chrinfo), chromosomes, las=3)
axis(2, at=seq(0, mlength, 10000000), seq(0,mlength,10000000)/ 1000000, las=2)
legend("topright", c("POP1","POP1Sex0","POP1Sex1"),lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c("green","red","blue")) 
dev.off()

## 2. Plot the POP2 and different Sex in POP2

chrinfo <- read.table("C:/Data/600KSNPchip/RawData/chromosomeinfo.txt", header=TRUE)     # Load the Chromosome information

chromosomes <- chrinfo[,"Chromosome"]
mlength <- max(chrinfo[,"Length"]) * 1.10              

pdf("Analysis/POP2DifferentSexHetro.pdf")                 
plot(y=c(0,  mlength), x=c(0, nrow(chrinfo)) + 0.5, t='n', yaxt="n", xlab="Chromosome", ylab="Location (Mbp)", xaxt="n")     # Make a frame

invisible(apply(chrinfo, 1, function(mrow){                          # Plot the Chromosomes
  xloc <- which(chromosomes == mrow["Chromosome"])
  lines(c(xloc, xloc),y=c(0, mrow["Length"]), lwd=4)
}))

invisible(apply(TargetRegionPOP2, 1, function(mrow){                 # Plot the Significant region for population 2  
  xloc <- which(chromosomes == mrow["chromosome"])
  if(length(xloc) > 0) lines(c(xloc, xloc), y=c(mrow["start"], mrow["end"]), t="l", lwd=4, col="green",lty=1)
  #return(as.numeric(mrow["start"]) + abs(as.numeric(mrow["start"]) - as.numeric(mrow["end"])) / 2)
}))

invisible(apply(TargetRegionpop2sex0, 1, function(mrow){                 # Plot the Significant region for pop2sex0
  xloc <-which(chromosomes == mrow["chromosome"])
  if(length(xloc) > 0) lines(c(xloc, xloc)-0.25, y=c(mrow["start"], mrow["end"]), t="l", lwd=4, col="red",lty=1)
}))

invisible(apply(TargetRegionpop2sex1, 1, function(mrow){                 # Plot the Significant region for pop2sex1
  xloc <- which(chromosomes == mrow["chromosome"])
  if(length(xloc) > 0) lines(c(xloc, xloc)+0.25, y=c(mrow["start"], mrow["end"]), t="l", lwd=4, col="blue",lty=1)
  #return(as.numeric(mrow["start"]) + abs(as.numeric(mrow["start"]) - as.numeric(mrow["end"])) / 2)
}))

axis(1, at=1:nrow(chrinfo), chromosomes, las=3)
axis(2, at=seq(0, mlength, 10000000), seq(0,mlength,10000000)/ 1000000, las=2)
legend("topright", c("POP2","POP2Sex0","POP2Sex1"),lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c("green","red","blue")) 
dev.off()