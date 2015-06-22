# Calculate the Confidence Interval for the Growth Traits of every 5 weeks using F11 and F12
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Apr, 2015
# first written Apr, 2015

setwd("D:/Chicken/Rdata/FineMapping")
F12Traits <- read.table("Analysis/F12_QTLdataWithBWG-Geno NA fill in.txt", na.strings=c("NA","?"),header=TRUE, sep="\t")

F12AllTraits20Week <- F12Traits[which(F12Traits[,"SchlachtWo"]=="20W" & F12Traits[,"Futter"] == "normal"),]  # F12AT20W (F12AllTraits20Week)

SigSNP <- read.table("Analysis/F12AllTraits20WLinkageAnalysis.txt",sep="\t", header=FALSE)

SigTraits <- NULL                                                                                      
for (x in 1:nrow(SigSNP)){
  if(length(unlist(strsplit(as.character(SigSNP[x,1]), " "))) > 1) {
  SigTraits <- c(SigTraits, strsplit(as.character(SigSNP[x,1]), " ")[[1]][1])
  }
}

Covs  <- c("Family","Schlupf")
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

# Variables
genotypes    <- F12AllTraits20Week[,SNPsForAnalysis]
phenotypes   <- F12AllTraits20Week[,SigTraits]
covariates   <- F12AllTraits20Week[,Covs]

mapQTLs <- function(genotypes, phenotypes, covariates){
  # Matrix LOD scores of markers x phenotypes
  LODS <- matrix(NA, ncol(phenotypes), ncol(genotypes),dimnames=list(colnames(phenotypes), colnames(genotypes))) # Markers x phenotypes
  for(x in 1:ncol(phenotypes)){  #cat(x,"\n")
    for(y in 1:ncol(genotypes)){ #cat(y,"\n")
      res <- anova(lm(as.numeric(phenotypes[,x]) ~ covariates[,"Family"] + covariates[,"Schlupf"] + genotypes[,y]))
      LODS[x, y] <- -log10(res[[5]][length(res[[5]]) - 1])
    }
  }
  return(LODS)
}

# Map the real data
LodForAll <- mapQTLs(genotypes, phenotypes, covariates)

## Confidence Interval
if(!file.exists("Analysis/CIboots1000-20WeekNormal.txt")){
  boots <- NULL
  for(x in 1:1000){
    idx <- sample(nrow(genotypes), 46, replace=TRUE)
    boots <- rbind(boots, mapQTLs(genotypes[idx, ], phenotypes[idx,], covariates[idx,]))
  }
    write.table(boots,"Analysis/CIboots1000-20WeekNormal.txt",sep="\t")
  }else{                                                                                                        # Or load them from Disk if we already have them
  cat("Loading boots information from disk\n")
  boots <- read.table("Analysis/CIboots1000-20WeekNormal.txt",sep="\t", header=TRUE, row.names=NULL)
}

CIAll <- NULL
for (trait in names(phenotypes)){
  TraitBoots <- boots[which(boots[,"row.names"] == trait),]
  topM <- names(which.max(LodForAll[trait,]))
  BootI <- sort(TraitBoots[,topM])[1000 * 0.025]
  CIAll <- cbind(CIAll, BootI)
}
colnames(CIAll) <- names(phenotypes)

pdf("Analysis/SigBodycomposition.pdf")

par(mfrow=c(4,4))
for(Trait in rownames(LodForAll)){
  plot(x=c(0,10), y=c(0,5), t="n", ylab="LOD Score", xlab="", xaxt="n", main= Trait)
  points(LodForAll[Trait, ],t="l")
  abline(h=-log10(0.05/9), lty=1) 
  abline(h=CIAll[,Trait], lty=2)  
  axis(1, at=1:9, SNPsForAnalysis, las=2, cex.axis = 0.77)
}

dev.off()



