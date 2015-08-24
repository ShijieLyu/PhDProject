# correlation and clustering for F12 slaughter traits
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/FineMapping")
QTLdataAll <- read.table("Analysis/QTLdataAll.txt", sep = "\t",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("character", 9), rep("factor",6)))
F12QTLdata <- read.table("Analysis/F12_QTLdataWithBWG-Geno NA fill in.txt", na.strings=c("NA","?"),header=TRUE, sep="\t")

F12ST <- F12QTLdata[which(F12QTLdata[,"SchlachtWo"]=="20W" & F12QTLdata[,"Futter"] == "normal"),]       # F12 slaughtered traits measured at 20th week, and the chicken were all with normal diet
F12ST <- cbind(F12ST, Batch = QTLdataAll[which(QTLdataAll[,"ID.Nr"] %in% as.character(F12ST[,"ID.Nr"])),"Batch"], Parents = QTLdataAll[which(QTLdataAll[,"ID.Nr"] %in% as.character(F12ST[,"ID.Nr"])),"Parents"])

GrowthTraits  <- c("Gew_1d","Gew_5Wo","Gew_10Wo","Gew_15Wo","Gew_20Wo","BWG05","BWG510","BWG1015","BWG1520")

F12STN <- names(F12ST)[-c(99,104)]                              # Remove "Stander._links_.Durchmesser" and "Stander._rechts_.Durchmesser_", since miss to many chickens                 
F12STN <- c(GrowthTraits, F12STN[c(37:39,41:143)])              # The names of slaughtered traits which for the association analysis

library("psych")

# For all traits
F12AllSTonly <- F12QTLdata[,F12STN]

for (EachTrait in F12STN){
  if(length(which(!is.na(F12AllSTonly[,EachTrait]))) == 0){
    F12AllSTonly[,EachTrait] <- NULL
  }
}

CorAllST <- corr.test(data.matrix(F12AllSTonly),method="pearson",use = "pairwise") 
#write.table(CorAllST[[1]], "Analysis/CorForAllSlaughterTraits.txt",sep="\t", row.names = TRUE, col.names = TRUE)

#mat <- data.matrix(F12AllSTonly)
#for (x in 1: ncol(mat)){
  #for (y in (x+1):ncol(mat)){
   # r <- round(cor.test(mat[,x], mat[,y],method="pearson", use = "pairwise")$estimate,2)
   # pvalue <- cor.test(mat[,x], mat[,y],method="pearson", use = "pairwise")$p.value
   # png(paste("Analysis/CorrelationPlots/", colnames(mat)[x],"-and-",colnames(mat)[y],".png",sep=""))
   # plot(mat[,x]~mat[,y],ylab=colnames(mat)[x],xlab=colnames(mat)[y],main=paste("Correlation between",colnames(mat)[x],"and",colnames(mat)[y]),cex.main=0.9)
   # abline(lm(mat[,x]~mat[,y]),col="red",cex=1.2)
   # legend("bottomright", c(paste("r=",r),paste("p=",pvalue)),cex=0.9)
   # dev.off()
  #}
#}

# For seleted slaughter(carcass) traits

F12STS <- F12STN[c(13,14,15,16,17,18,20,21,23,24,26,27,29,31,33,34,35,36,37,38,39,40,41,42,43)]                    # Select the traits which is necessary to analyse
CorSST <- corr.test(data.matrix(F12AllSTonly[,F12STS]),method="pearson",use = "pairwise")                         # Correlation analysis for the selected slaughter traits

#write.table(CorSST[[1]], "Analysis/CorForSlaughterTraits.txt",sep="\t", row.names = TRUE, col.names = TRUE)
#write.table(CorSST[[4]], "Analysis/CorForSlaughterTraits-Pvalue.txt",sep="\t", row.names = TRUE, col.names = TRUE)

plot(hclust(dist(CorSST[[1]])),cex=0.7,h=-1, main= "Clustering of Slaughter Traits", xlab= "Traits")               # Generate the dendrogram




