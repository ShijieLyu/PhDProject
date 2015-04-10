# Get the significant region with high Hetro/Homo ratio and plot it in F5, POP1 and POP2
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Feb, 2015
# first written Feb, 2015

setwd("C:/Data/60KSNPchip")

### Load the F5 generation Analysis Data Ratio and Threshold

F5Generation <- read.table("Analysis/HetroRatio-F5.txt", header=TRUE)           # load the Population 1 Hetro Ratio
F5GenerationThr <- read.table("Analysis/thresholdperChr-F5.txt", header=TRUE)   # Load the threshold for Population 1 under 1000 times permutation

F5Generationmatrix <- as.matrix(1:33)                                           # The format of the F5GenerationThr is nor correct, here I re-organised it
F5GenerationThrChr <- as.character(F5GenerationThr[seq(1,66,2),])
F5GenerationThrScore <- as.numeric(as.character(F5GenerationThr[seq(2,66,2),]))
rownames(F5Generationmatrix) <- F5GenerationThrChr
F5Generationmatrix[,1] <- F5GenerationThrScore
F5GenerationThr <- F5Generationmatrix
colnames(F5GenerationThr) <- c("threshold")

chromosomes <- as.character(unique(F5Generation[,"chromosome"]))
TargetRegion <- function(F5Generation,F5GenerationThr) {                             # Select the ratio which higher than the threshold for each Population
  Region <- NULL
  for (chr in chromosomes){
    ScorePerChr <- F5Generation[which(F5Generation[,"chromosome"] == chr),]
    thresholdPerChr <- as.numeric(F5GenerationThr[which(rownames(F5GenerationThr) == chr),"threshold"])
    TargetRegionPerchr <- ScorePerChr[which(as.numeric(ScorePerChr[,"score"]) > thresholdPerChr),]
    Region <- rbind(Region, TargetRegionPerchr)
  }
  return(Region)
}
TargetRegionF5 <- TargetRegion(F5Generation,F5GenerationThr)

### Load the POP1() and POP2() Analysis Data -- Ratio and Threshold

POP1 <- read.table("C:/Data/600KSNPchip/Analysis/HetroRatioPOP1-100.txt", header=TRUE)    # load the Population 1 Hetro Ratio POP1=WL77
POP2 <- read.table("C:/Data/600KSNPchip/Analysis/HetroRatioPOP2-100.txt", header=TRUE)    # load the Population 2 Hetro Ratio POP2=NHI
POP1Thr <- read.table("C:/Data/600KSNPchip/Analysis/thresholdperChr-POP1-26.02.15.txt", header=TRUE)   # Load the threshold for Population 1 under 1000 times permutation
POP2Thr <- read.table("C:/Data/600KSNPchip/Analysis/thresholdperChr-POP2-26.02.15.txt", header=TRUE)   # Load the threshold for Population 2 under 1000 times permutation

POP1matrix <- as.matrix(1:32)                                         # The format of the POP1Thr is nor correct, here I re-organised it
POP1ThrChr <- as.character(POP1Thr[seq(1,64,2),])
POP1ThrScore <- as.numeric(as.character(POP1Thr[seq(2,64,2),]))
rownames(POP1matrix) <- POP1ThrChr
POP1matrix[,1] <- POP1ThrScore
POP1Thr <- POP1matrix

POP2matrix <- as.matrix(1:32)
POP2ThrChr <- as.character(POP2Thr[seq(1,64,2),])
POP2ThrScore <- as.numeric(as.character(POP2Thr[seq(2,64,2),]))
rownames(POP2matrix) <- POP2ThrChr
POP2matrix[,1] <- POP2ThrScore
POP2Thr <- POP2matrix

colnames(POP1Thr) <- c("threshold")
colnames(POP2Thr) <- c("threshold")

chromosomes <- as.character(unique(POP1[,"chromosome"]))

TargetRegion <- function(POP1,POP1Thr) {                             # Select the ratio which higher than the threshold for each Population
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

### Plot chromosomes and add the region information to it---VERITAL

chrinfo <- read.table("C:/Data/600KSNPchip/RawData/chromosomeinfo.txt", header=TRUE)     # Load the Chromosome information

chromosomes <- chrinfo[,"Chromosome"]
mlength <- max(chrinfo[,"Length"]) * 1.10              

pdf("Analysis/TwoBreedsAndandF5HetroComparison_Verital.pdf")                 
plot(y=c(0,  mlength), x=c(0, nrow(chrinfo)) + 0.5, t='n', yaxt="n", xlab="Chromosome", ylab="Location (Mbp)", xaxt="n")     # Make a frame

invisible(apply(chrinfo, 1, function(mrow){                          # Plot the Chromosomes
  xloc <- which(chromosomes == mrow["Chromosome"])
  lines(c(xloc, xloc),y=c(0, mrow["Length"]), lwd=4)
}))

invisible(apply(TargetRegionPOP1, 1, function(mrow){                 # Plot the Significant region for population 1  
  xloc <- which(chromosomes == mrow["chromosome"])
  if(length(xloc) > 0) lines(c(xloc, xloc)-0.15, y=c(mrow["start"], mrow["end"]), t="l", lwd=4, col="red",lty=1)
  #return(as.numeric(mrow["start"]) + abs(as.numeric(mrow["start"]) - as.numeric(mrow["end"])) / 2)
}))

invisible(apply(TargetRegionPOP2, 1, function(mrow){                 # Plot the Significant region for population 2
  xloc <-which(chromosomes == mrow["chromosome"])
  if(length(xloc) > 0) lines(c(xloc, xloc)+0.15, y=c(mrow["start"], mrow["end"]), t="l", lwd=4, col="blue",lty=1)
}))

invisible(apply(TargetRegionF5, 1, function(mrow){                 # Plot the Significant region for F5
  xloc <- which(chromosomes == mrow["chromosome"])
  if(length(xloc) > 0) lines(c(xloc, xloc), y=c(mrow["start"], mrow["end"]), t="l", lwd=4, col="green",lty=1)
  #return(as.numeric(mrow["start"]) + abs(as.numeric(mrow["start"]) - as.numeric(mrow["end"])) / 2)
}))

axis(1, at=1:nrow(chrinfo), chromosomes, las=3)
axis(2, at=seq(0, mlength, 10000000), seq(0,mlength,10000000)/ 1000000, las=2)
legend("topright", c("POP1","POP2","F5"),lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c("red","blue","green")) 
dev.off()
