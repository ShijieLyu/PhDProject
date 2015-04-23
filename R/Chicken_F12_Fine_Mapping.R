# Fine Mapping for the chicken F12 Generation
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Apr, 2015
# first written Apr, 2015

### Reorder the slaughter data as the body weight rank and uniform the format of two phenotype data 

setwd("D:/Chicken/Rdata/FineMapping")
bodyweight <- read.table("RawData/F12_BodyWeight_Phenotypes.txt", na.strings="",header=TRUE, sep="\t")
slaughtdata <- read.table("RawData/F12_Slaught_Phenotypes.txt", na.strings="",header=TRUE, sep="\t")
genotypes <- read.table("RawData/F12_Genotypes.txt", na.strings="-", header=TRUE, sep="\t")

genotypes <- genotypes[match(bodyweight[,"ID.Nr"], genotypes[,"ID.Nr"]),]                  # Reorder the genotypes as bodyweight
slaughtdata <- slaughtdata[match(bodyweight[,"ID.Nr"], slaughtdata[,"ID.Nr"]),]            # Reorder the slaughtdata as bodyweight
slaughtdata <- cbind(slaughtdata, Family = bodyweight[,"Family"])                          # Add Family information for every individual
slaughtdata <- cbind(slaughtdata[,1:10], Family = slaughtdata[,124],slaughtdata[,12:123])  # Reorder the column
slaughtdata <- slaughtdata[,6:123]
colnames(slaughtdata)[3:4] <- c("Vatter","Mutter")  # Give the bodyweight and the slaughterdata same colnames

bodyweight[3]<-NULL  # Remove the SD column
bodyweight <- bodyweight[,c(1,2,4,5,3,6:28)]        # Give the bodyweight and the slaughterdata same colnames

bodyweight[which(bodyweight[,"ID.Nr"] == 430),]     # check a samples
slaughtdata[which(slaughtdata[,"ID.Nr"] == 430),1:15]

# write.table(slaughtdata,"Analysis/F12_Slaught_Phenotypes_withfamily.txt",sep="\t",row.names = FALSE,quote = FALSE)

### Check the info. of body weight and slaughter data is same or not
MessUprow <- NULL
for (mrow in 1:nrow(bodyweight)){
  if (any((bodyweight[mrow,1:7] == slaughtdata[mrow,1:7]) == FALSE)) {
    MessUprow <- c(MessUprow, mrow)
  }
}

bodyweight[MessUprow,1:7]
slaughtdata[MessUprow,1:7]  # ID 320 430 448 468-1 483-1 the father or mother is messed up

# write.table(c(bodyweight[MessUprow,1:7],slaughtdata[MessUprow,1:7]),"MessUP.txt",sep="\t",row.names = FALSE,quote = FALSE)

### Combine the Slaughter Data, Bodyweight Data and the Genotypes to one matrix

QTLdata <- cbind(genotypes[,-(1:2)],bodyweight,slaughtdata[,-c(1:7,42,51,68,79,93,106)])
QTLdata <- QTLdata[,c(6:12,1:5,13:ncol(QTLdata))]

### Association Analyisis for the Markers and the Phenotypes

LinkageAnalysis <- function(QTLdata, Phennotype = "Gew_1d" ){
  SNPsforQTL <- colnames(QTLdata[,8:12])                         # Load the Markers and genotypes
  PhenotypesforQTL <- colnames(QTLdata[,13:ncol(QTLdata)])       # Load the Phenotypes Data
  LODscores <- NULL
  for (OneSNP in SNPsforQTL){
    if(length(which(!is.na(QTLdata[,Phennotype]))) == 0){
      LODscores <- rbind(LODscores, NA)
      colnames(LODscores) <- "SNP_Name"
    }else if(length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"Futter"])) == 1 && length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"SchlachtWo"])) == 1 && length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"Family"])) == 1){
      res <- anova(lm(as.numeric(QTLdata[,Phennotype]) ~ QTLdata[,OneSNP]))
      LODscores <- rbind(LODscores, -log10(res[[5]]))
      colnames(LODscores) <- c("SNP_Name", "Residuals")    
    }else if(length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"Futter"])) == 1 && length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"SchlachtWo"])) == 1){
      res <- anova(lm(as.numeric(QTLdata[,Phennotype]) ~ QTLdata[,"Family"] + QTLdata[,OneSNP]))
      LODscores <- rbind(LODscores, -log10(res[[5]]))
      colnames(LODscores) <- c("Family", "SNP_Name", "Residuals")  
    }else if(length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"Futter"])) == 1 && length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"Family"])) == 1){
      res <- anova(lm(as.numeric(QTLdata[,Phennotype]) ~ QTLdata[,"SchlachtWo"] + QTLdata[,OneSNP]))
      LODscores <- rbind(LODscores, -log10(res[[5]]))
      colnames(LODscores) <- c("SchlachtWo", "SNP_Name", "Residuals")  
    }else if(length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"SchlachtWo"])) == 1 && length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"Family"])) == 1){
      res <- anova(lm(as.numeric(QTLdata[,Phennotype]) ~ QTLdata[,"Futter"] + QTLdata[,OneSNP]))
      LODscores <- rbind(LODscores, -log10(res[[5]]))
      colnames(LODscores) <- c("Futter", "SNP_Name", "Residuals")  
    }else if(length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"SchlachtWo"])) == 1){
      res <- anova(lm(as.numeric(QTLdata[,Phennotype]) ~ QTLdata[,"Futter"] + QTLdata[,"Family"] + QTLdata[,OneSNP]))
      LODscores <- rbind(LODscores, -log10(res[[5]]))
      colnames(LODscores) <- c("Futter", "Family", "SNP_Name", "Residuals")
    }else if(length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"Family"])) == 1){
      res <- anova(lm(as.numeric(QTLdata[,Phennotype]) ~ QTLdata[,"Futter"] + QTLdata[,"SchlachtWo"] + QTLdata[,OneSNP]))
      LODscores <- rbind(LODscores, -log10(res[[5]]))
      colnames(LODscores) <- c("Futter", "SchlachtWo", "SNP_Name", "Residuals")
    }else if(length(unique(QTLdata[which(!is.na(QTLdata[,Phennotype])),"Futter"])) == 1){
      res <- anova(lm(as.numeric(QTLdata[,Phennotype]) ~ QTLdata[,"Family"] + QTLdata[,"SchlachtWo"] + QTLdata[,OneSNP]))
      LODscores <- rbind(LODscores, -log10(res[[5]]))
      colnames(LODscores) <- c("Family", "SchlachtWo", "SNP_Name", "Residuals")
    }else{
      res <- anova(lm(as.numeric(QTLdata[,Phennotype]) ~ QTLdata[,"Futter"] + QTLdata[,"Family"] + QTLdata[,"SchlachtWo"] + QTLdata[,OneSNP]))
      LODscores <- rbind(LODscores, -log10(res[[5]]))
      colnames(LODscores) <- c("Futter", "Family", "SchlachtWo", "SNP_Name", "Residuals")  
    }    
  }
  rownames(LODscores) <- SNPsforQTL
  return(LODscores)
}

LinkageLOD <- vector("list", length(PhenotypesforQTL)) 
names(LinkageLOD) <- PhenotypesforQTL
for (Phennotype in PhenotypesforQTL) {
  LinkageLOD[[Phennotype]] <- LinkageAnalysis(QTLdata, Phennotype)
}

# threshold <- -log10(0.05/ (length(PhenotypesforQTL) * length(SNPsforQTL)))

threshold <- -log10(0.05/ length(SNPsforQTL))
for(x in 1:length(LinkageLOD)){
  SigSNP <- which(LinkageLOD[[x]][,"SNP_Name"] > threshold)                                                      # Analyse the resulting profiles, and look at which markers are above the threshold
  cat(names(LinkageLOD)[x], names(SigSNP),"\n")                                                                  
}



##############################################################################################################################################
### Organise the Data to the r/QTL package required

for (mcol in 3:7){
  for (mrow in 1:nrow(genotypes)){
    if unique(genotypes[mrow,mcol]
  }

}

genotypes











