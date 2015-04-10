# Feed efficiency QTL Map
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Feb, 2015
# first written Feb, 2015

setwd("C:/Data/Feed_Efficiency_Map")

FCR <- read.table("FeedConversionRatio.txt", header= TRUE, sep="\t")
RFI <- read.table("ResidualFeedIntake.txt", header= TRUE, sep="\t")
FE <- read.table("FeedEfficiency.txt", header= TRUE, sep="\t")
chrinfo <- read.table("C:/Data/600KSNPchip/RawData/chromosomeinfo.txt", header=TRUE)     # Load the Chromosome information

chromosomes <- as.character(chrinfo[,"Chromosome"])
mlength <- max(chrinfo[,"Length"]) * 1.10   

## FCR Picture           
pdf("FeedConversionRatio.pdf")       
plot(x=c(0, nrow(chrinfo)) + 0.5, y=c(0,  mlength), t='n', yaxt="n", xlab="Chromosome", ylab="Location (Mbp)", main= "The Feed Conversion Ratio QTL Map", xaxt="n")     # Make a frame

invisible(apply(chrinfo, 1, function(mrow){                          # Plot the Chromosomes
  xloc <- which(chromosomes == mrow["Chromosome"])
  lines(c(xloc, xloc),y=c(0, mrow["Length"]), lwd=4)
}))

invisible(apply(FCR[1:7,],1,function(mrow){
  xloc <- which(chromosomes == as.numeric(as.character(mrow["Chr"])))
  lines(x=c(xloc, xloc)+0.15, y=c(mrow["Start"], mrow["End"]), t="l", lwd=5, col="red",lty=1)
}))

invisible(apply(FCR[7,],1,function(mrow){
  xloc <- which(chromosomes == mrow["Chr"])
  if(length(xloc) > 0) lines(x=c(xloc, xloc)+0.15, y=c(mrow["Start"], mrow["End"]), t="l", lwd=1.8, col="red",lty=1)
}))

invisible(apply(FCR[8,],1,function(mrow){
  xloc <- which(chromosomes == mrow["Chr"])
  if(length(xloc) > 0) lines(x=c(xloc, xloc)+0.2, y=c(mrow["Start"], mrow["End"]), t="l", lwd=1.8, col="blue",lty=1)
  #; text(xloc, y= mrow["Start"],"Ref1",cex=0.5)
}))

invisible(apply(FCR[9:10,],1,function(mrow){
  xloc <- which(chromosomes == mrow["Chr"])
  if(length(xloc) > 0) lines(x=c(xloc, xloc)+0.25, y=c(mrow["Start"], mrow["End"]), t="l", lwd=1.8, col="green",lty=1)
}))

axis(1, at=1:nrow(chrinfo), chromosomes, las=3)
axis(2, at=seq(0, mlength, 10000000), seq(0,mlength,10000000)/ 1000000, las=2)
legend("topright", c("De Koning,D.,et al.(2004)","Nones,K.,et al.(2006)","Sheng,Z.,et al.(2013)"),lty=c(1,1,1), lwd=c(2.5,2.5,2.5),col=c("red","blue","green")) 

dev.off()

## RFI Picture
pdf("ResidualFeedIntake.pdf")
plot(x=c(0, nrow(chrinfo)) + 0.5, y=c(0,  mlength), t='n', yaxt="n", xlab="Chromosome", ylab="Location (Mbp)",main= "The Residual Feed Intake QTL Map", xaxt="n")     # Make a frame

invisible(apply(chrinfo, 1, function(mrow){                          # Plot the Chromosomes
  xloc <- which(chromosomes == mrow["Chromosome"])
  lines(c(xloc, xloc),y=c(0, mrow["Length"]), lwd=4)
}))

invisible(apply(RFI[1:12,],1,function(mrow){
  xloc <- which(chromosomes == mrow["Chr"])
  if(length(xloc) > 0) lines(x=c(xloc, xloc)+0.15, y=c(mrow["Start"], mrow["End"]), t="l", lwd=1.8, col="red",lty=1)
}))

invisible(apply(RFI[13:15,],1,function(mrow){
  xloc <- which(chromosomes == mrow["Chr"])
  if(length(xloc) > 0) lines(x=c(xloc, xloc)+0.2, y=c(mrow["Start"], mrow["End"]), t="l", lwd=1.8, col="blue",lty=1)
}))

axis(1, at=1:nrow(chrinfo), chromosomes,las=3)
axis(2, at=seq(0, mlength, 10000000), seq(0,mlength,10000000)/ 1000000, las=2)
legend("topright", c("De Koning,D.,et al.(2004)","Wolc,A.,et al.(2013)"),lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue")) 

dev.off()

## FE Picture
pdf("FeedEfficiency.pdf")
plot(x=c(0, nrow(chrinfo)) + 0.5, y=c(0,  mlength), t='n', yaxt="n", xlab="Chromosome", ylab="Location (Mbp)", main= "The Feed Efficiency QTL Map", xaxt="n")     # Make a frame

invisible(apply(chrinfo, 1, function(mrow){                          # Plot the Chromosomes
  xloc <- which(chromosomes == mrow["Chromosome"])
  lines(c(xloc, xloc),y=c(0, mrow["Length"]), lwd=4)
}))

invisible(apply(FE[1,],1,function(mrow){
  xloc <- which(chromosomes == mrow["Chr"])
  if(length(xloc) > 0) lines(x=c(xloc, xloc)+0.15, y=c(mrow["Start"], mrow["End"]), t="l", lwd=1.8, col="red",lty=1)
}))

invisible(apply(FE[2,],1,function(mrow){
  xloc <- which(chromosomes == mrow["Chr"])
  if(length(xloc) > 0) lines(x=c(xloc, xloc)+0.2, y=c(mrow["Start"], mrow["End"]), t="l", lwd=1.8, col="blue",lty=1)
}))

axis(1, at=1:nrow(chrinfo), chromosomes, las=3)
axis(2, at=seq(0, mlength, 10000000), seq(0,mlength,10000000)/ 1000000, las=2)
legend("topright", c("Hansen,C.,et al.(2005)","Rosario,M.,et al.(2014) "),lty=c(1,1), lwd=c(2.5,2.5),col=c("red","blue")) 

dev.off()