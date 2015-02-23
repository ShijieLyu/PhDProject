# Get the significant region with high Hetro/Homo ratio and plot it
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Feb, 2015
# first written Feb, 2015

setwd("C:/Data/600KSNPchip")

POP1 <- read.table("Analysis/HetroRatioPOP1-100.txt", header=TRUE)    # load the Population 1 Hetro Ratio
POP2 <- read.table("Analysis/HetroRatioPOP2-100.txt", header=TRUE)    # load the Population 2 Hetro Ratio
POP1Thr <- read.table("Analysis/thresholdperChrPOP1-1000.txt", header=TRUE)   # Load the threshold for Population 1 under 1000 times permutation
POP2Thr <- read.table("Analysis/thresholdperChrPOP2-1000.txt", header=TRUE)   # Load the threshold for Population 2 under 1000 times permutation

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
    thresholdPerChr <- POP1Thr[which(rownames(POP1Thr) == chr),"threshold"]
    TargetRegionPerchr <- ScorePerChr[which(ScorePerChr[,"score"] > thresholdPerChr),]
    Region <- rbind(Region, TargetRegionPerchr)
  }
  return(Region)
}

TargetRegionPOP1 <- TargetRegion(POP1,POP1Thr)
TargetRegionPOP2 <- TargetRegion(POP2,POP2Thr)


### Plot chromosomes and add the region information to it---VERITAL

chrinfo <- read.table("RawData/chromosomeinfo.txt", header=TRUE)     # Load the Chromosome information

chromosomes <- chrinfo[,"Chromosome"]
mlength <- max(chrinfo[,"Length"]) * 1.10              

pdf("Analysis/TwoBreedsHetroComparison_Verital_1000.pdf")                 
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

axis(1, at=1:nrow(chrinfo), chromosomes)
axis(2, at=seq(0, mlength, 10000000), seq(0,mlength,10000000)/ 1000000, las=2)
dev.off()

### ### Plot chromosomes and add the region information to it---HORIZONTAL

chrinfo <- read.table("RawData/chromosomeinfo.txt", header=TRUE)
chrinfo <- chrinfo[rev(seq_len(nrow(chrinfo))),]           # Reverse the Chr order

chromosomes <- chrinfo[,"Chromosome"]
mlength <- max(chrinfo[,"Length"]) * 1.10

pdf("Analysis/TwoBreedsHetroComparison_Horizontal_1000.pdf")
plot(x=c(0,mlength), y=c(0,nrow(chrinfo)), type= "n", xaxt="n", xaxs='i', yaxt="n", xlab="Location(Mbp)", ylab="Chromosome")

yStep <- 1
invisible(apply(chrinfo,1,function(mrow){
  lines(x=c(0, mrow["Length"]), y = c(yStep,yStep), type="l", lwd=2.5)
  yStep <<- yStep + 1
}))

invisible(apply(TargetRegionPOP1,1,function(mrow){
  yStep <- which(chromosomes == mrow["chromosome"])
  lines(x=c(mrow["start"], mrow["end"]), y = c(yStep,yStep)+0.15, t="l", lwd=2.5, col="red",lty=1)
}))

invisible(apply(TargetRegionPOP2,1,function(mrow){
  yStep <- which(chromosomes == mrow["chromosome"])
  lines(x=c(mrow["start"], mrow["end"]), y = c(yStep,yStep)-0.15, t="l", lwd=2.5, col="blue",lty=1)
}))

axis(1, at=seq(0, mlength, 10000000), seq(0,mlength,10000000)/ 1000000)
axis(2, at=1:nrow(chrinfo), chromosomes, las=2)
abline(v=seq(0, mlength, 10000000), col = "lightgray", lty = "dotted")
dev.off()

