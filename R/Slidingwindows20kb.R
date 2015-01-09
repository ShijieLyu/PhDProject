# show the location of the SNPs and construct a sliding window of 20kb
#
# copyright (c) 2014-2020 - Shijie Lyu
# last modified Dec, 2014
# first written Dec, 2014

setwd("C:/Data/600KSNPchip")    

snpsbigregion <- read.table("Analysis/finemappingSNPswithannotation.txt", sep = "\t", header=TRUE)  # load the full data
snpslocation <- snpsbigregion[,"Physical.Position"]                                                 # load the col-name which means location 
plot(x= snpslocation, y= rep(1, length(snpslocation)), pch="|")                                     # plot the location
 
windowsize <- 20000                                                                                 # create the size of the window
starts <- seq(min(snpslocation), max(snpslocation)-windowsize, by = windowsize)                     # calculate the start for exery window
windownumber <- length(starts)                                                                      # how many windows will use 
SNPsNrinwindows <- NULL
for(x in 1:windownumber){                                                                           # creat the sliding windows with the 20kb size
  SNPsNrperwindow <- length(which(snpslocation <= starts[x] + (windowsize-1) & snpslocation >= starts[x]))
  SNPsNrinwindows <- c(SNPsNrinwindows, SNPsNrperwindow)
} 
noSNPswindow <- which(SNPsNrinwindows == 0)  
noSNPswindow

pdf("Analysis/rplot.pdf")                                                                                    # output as pdf
plot(starts,SNPsNrinwindows,type="l",xlab="SNPs Location", ylab="SNPs in windows", main= "The density of the SNPs for fine mapping(20kb sliding window)")
points(x= snpslocation, y= rep(0, length(snpslocation)), pch="|")
dev.off()





