# correlation analysis foe the fat content using the wet and dry methods
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Shijie Lyu
# last modified Sep, 2015
# first written Sep, 2015


setwd("D:/Chicken/Rdata/Chicken Fat measurement")
F12SD <- read.table("D:/Chicken/Rdata/FineMapping/RawData/F12_Slaught_Phenotypes.txt", na.strings=c("","?"),header=TRUE, sep="\t")   # F12 Slaughter data
FDM <- read.table("RawData/Fettmessung.txt", na.strings="", header=TRUE, sep="\t")                                                   # Fat measurement with Dry method
FDM <- cbind(FDM, BMFO = rep(NA, nrow(FDM)), BTFO = rep(NA, nrow(FDM)), BMAFO = rep(NA, nrow(FDM)), BMBFO = rep(NA, nrow(FDM)), 
                  LUpFO = rep(NA, nrow(FDM)), LLoFO = rep(NA, nrow(FDM)), LMAFO = rep(NA, nrow(FDM)), LMBFO = rep(NA, nrow(FDM)),
                  BMF = rep(NA, nrow(FDM)), BTF = rep(NA, nrow(FDM)), BMAF = rep(NA, nrow(FDM)), BMBF = rep(NA, nrow(FDM)), 
                  LUpF = rep(NA, nrow(FDM)), LLoF = rep(NA, nrow(FDM)), LMAF = rep(NA, nrow(FDM)), LMBF = rep(NA, nrow(FDM)))  
                  
toN <- function(x){as.numeric(paste0(strsplit(as.character(x[1]),"")[[1]][1:length(strsplit(as.character(x[1]),"")[[1]])-1],collapse = ""))}   # Remove the %, transform to numericfor (x in 1:nrow(FDM)){
  
for (x in 1:nrow(FDM)){
  FDM[x,"BMFO"]  <- mean(c(toN(FDM[x,"Brust.Mitte.1.1"]),toN(FDM[x,"Brust.Mitte.2.1"]),toN(FDM[x,"Brust.Mitte.3.1"])),na.rm = TRUE)*0.01      # BMFO = Original breast middle fat content, which means to use the percentage in dry muscle directly.
  FDM[x,"BTFO"]  <- mean(c(toN(FDM[x,"Brust.Spitze.1.1"]),toN(FDM[x,"Brust.Spitze.2.1"]),toN(FDM[x,"Brust.Spitze.3.1"])),na.rm = TRUE)*0.01   # BTFO = Original breast tip fat
  FDM[x,"BMAFO"] <- mean(c(toN(FDM[x,"Brust.MusA.1.1"]),toN(FDM[x,"Brust.MusA.2.1"]),toN(FDM[x,"Brust.MusA.3.1"])),na.rm = TRUE)*0.01         # BMAFO = Original breast muscle A fat  
  FDM[x,"BMBFO"] <- mean(c(toN(FDM[x,"Brust.MusB.1.1"]),toN(FDM[x,"Brust.MusB.2.1"]),toN(FDM[x,"Brust.MusB.3.1"])),na.rm = TRUE)*0.01         # BMBFO = Original breast muscle B fat 
  FDM[x,"LUpFO"] <- mean(c(toN(FDM[x,"Keule.oben.1.1"]),toN(FDM[x,"Keule.oben.2.1"]),toN(FDM[x,"Keule.oben.3.1"])),na.rm = TRUE)*0.01         # LUpFO = Original leg upper part fat 
  FDM[x,"LLoFO"] <- mean(c(toN(FDM[x,"Keule.unten.1.1"]),toN(FDM[x,"Keule.unten.2.1"]),toN(FDM[x,"Keule.unten.3.1"])),na.rm = TRUE)*0.01      # LLoFO = Original leg lower part fat
  FDM[x,"LMAFO"] <- mean(c(toN(FDM[x,"Keule.MusA.1.1"]),toN(FDM[x,"Keule.MusA.2.1"]),toN(FDM[x,"Keule.MusA.3.1"])),na.rm = TRUE)*0.01         # LMAFO = Original leg muscle A fat  
  FDM[x,"LMBFO"] <- mean(c(toN(FDM[x,"Keule.MusB.1.1"]),toN(FDM[x,"Keule.MusB.2.1"]),toN(FDM[x,"Keule.MusB.3.1"])),na.rm = TRUE)*0.01         # LMBFO = Original leg muscle B fat     
  FDM[x,"BMF"]   <- mean(c((toN(FDM[x,"Brust.Mitte.1.1"])*0.01*FDM[x,"Brust.Mitte.1"])/FDM[x,"Brust.Mitte"],(toN(FDM[x,"Brust.Mitte.2.1"])*0.01*FDM[x,"Brust.Mitte.2"])/FDM[x,"Brust.Mitte"],(toN(FDM[x,"Brust.Mitte.3.1"])*0.01*FDM[x,"Brust.Mitte.3"])/FDM[x,"Brust.Mitte"]),na.rm = TRUE)           # BMF = breast middle fat in wet muscle
  FDM[x,"BTF"]   <- mean(c((toN(FDM[x,"Brust.Spitze.1.1"])*0.01*FDM[x,"Brust.Spitze.1"])/FDM[x,"Brust.Spitze"],(toN(FDM[x,"Brust.Spitze.2.1"])*0.01*FDM[x,"Brust.Spitze.2"])/FDM[x,"Brust.Spitze"],(toN(FDM[x,"Brust.Spitze.3.1"])*0.01*FDM[x,"Brust.Spitze.3"])/FDM[x,"Brust.Spitze"]),na.rm = TRUE)  # BTF = breast tip fat in wet muscle
  FDM[x,"BMAF"]  <- mean(c((toN(FDM[x,"Brust.MusA.1.1"])*0.01*FDM[x,"Brust.MusA.1"])/FDM[x,"Brust.MusA"],(toN(FDM[x,"Brust.MusA.2.1"])*0.01*FDM[x,"Brust.MusA.2"])/FDM[x,"Brust.MusA"],(toN(FDM[x,"Brust.MusA.3.1"])*0.01*FDM[x,"Brust.MusA.3"])/FDM[x,"Brust.MusA"]),na.rm = TRUE)                    # BMAF = breast muscle A fat in wet muscle
  FDM[x,"BMBF"]  <- mean(c((toN(FDM[x,"Brust.MusB.1.1"])*0.01*FDM[x,"Brust.MusB.1"])/FDM[x,"Brust.MusB"],(toN(FDM[x,"Brust.MusB.2.1"])*0.01*FDM[x,"Brust.MusB.2"])/FDM[x,"Brust.MusB"],(toN(FDM[x,"Brust.MusB.3.1"])*0.01*FDM[x,"Brust.MusB.3"])/FDM[x,"Brust.MusB"]),na.rm = TRUE)                    # BMBF = breast muscle B fat in wet muscle
  FDM[x,"LUpF"]  <- mean(c((toN(FDM[x,"Keule.oben.1.1"])*0.01*FDM[x,"Keule.oben.1"])/FDM[x,"Keule.oben"],(toN(FDM[x,"Keule.oben.2.1"])*0.01*FDM[x,"Keule.oben.2"])/FDM[x,"Keule.oben"],(toN(FDM[x,"Keule.oben.3.1"])*0.01*FDM[x,"Keule.oben.3"])/FDM[x,"Keule.oben"]),na.rm = TRUE)                    # LUpF = leg upper part fat in wet muscle
  FDM[x,"LLoF"]  <- mean(c((toN(FDM[x,"Keule.unten.1.1"])*0.01*FDM[x,"Keule.unten.1"])/FDM[x,"Keule.unten"],(toN(FDM[x,"Keule.unten.2.1"])*0.01*FDM[x,"Keule.unten.2"])/FDM[x,"Keule.unten"],(toN(FDM[x,"Keule.unten.3.1"])*0.01*FDM[x,"Keule.unten.3"])/FDM[x,"Keule.unten"]),na.rm = TRUE)           # LLoF = leg lower part fat in wet muscle
  FDM[x,"LMAF"]  <- mean(c((toN(FDM[x,"Keule.MusA.1.1"])*0.01*FDM[x,"Keule.MusA.1"])/FDM[x,"Keule.MusA"],(toN(FDM[x,"Keule.MusA.2.1"])*0.01*FDM[x,"Keule.MusA.2"])/FDM[x,"Keule.MusA"],(toN(FDM[x,"Keule.MusA.3.1"])*0.01*FDM[x,"Keule.MusA.3"])/FDM[x,"Keule.MusA"]),na.rm = TRUE)                    # LMAF = Original leg lower part fat
  FDM[x,"LMBF"]  <- mean(c((toN(FDM[x,"Keule.MusB.1.1"])*0.01*FDM[x,"Keule.MusB.1"])/FDM[x,"Keule.MusB"],(toN(FDM[x,"Keule.MusB.2.1"])*0.01*FDM[x,"Keule.MusB.2"])/FDM[x,"Keule.MusB"],(toN(FDM[x,"Keule.MusB.3.1"])*0.01*FDM[x,"Keule.MusB.3"])/FDM[x,"Keule.MusB"]),na.rm = TRUE)        # LMBF = leg lower part fat     
}

Mean <- apply(FDM[,64:ncol(FDM)],2,mean,na.rm=TRUE)
plot(Mean,type = "h",lwd = 15, lend = 2,xaxt="n",ylim=c(0,0.2))
SD <-apply(FDM[,64:ncol(FDM)],2,sd,na.rm=TRUE)
for (x in 1:ncol(FDM[,64:ncol(FDM)])){
  lines(c(x,x),c(Mean[x],Mean[x]+SD[x]))
}
axis(1, at=1:ncol(FDM[,64:ncol(FDM)]), names(Mean), las=3)


F12SD <- cbind(F12SD, BF = rep(NA, nrow(F12SD)), LF = rep(NA, nrow(F12SD)))
for (x in 1:nrow(F12SD)){
  F12SD[x,"BF"] <- mean(c(F12SD[x,"MRI.Brust.Fett"]/F12SD[x,"Gewicht.Teilstuck_Brust_fur_Messung"], F12SD[x,"MRI.Brust.Fett.1"]/F12SD[x,"Gewicht.Teilstuck_Brust_fur_Messung"]),na.rm=TRUE)   # BF=Breast Fat 
  F12SD[x,"LF"] <- mean(c(F12SD[x,"MRI.Keule_ohne_Knochen.Fett"]/F12SD[x,"Gewicht.Teilstuck.Keule_ohne_Knochen_f._Messung"], F12SD[x,"MRI.Keule_ohne_Knochen.Fett.1"]/F12SD[x,"Gewicht.Teilstuck.Keule_ohne_Knochen_f._Messung"]),na.rm=TRUE)   # LF= Leg Fat  
}

FDMF12 <- FDM[which(FDM[,"Gene.ration"] == "F12"),]                     # Select the F12 generation
F12SD <- F12SD[match(FDMF12[,"Tier.Nr."],F12SD[,"ID.Nr"]),]             # filter the same individuals
F12SDN <- F12SD[which(!is.na(F12SD[,"ID.Nr"])),]                        # remove the NA, which can not match in the two groups
FDMF12N <- FDMF12[which(FDMF12[,"Tier.Nr."] %in% F12SDN[,"ID.Nr"]),]     # make the two group have the same individuals

# Correlation analysis

TraitsMatrix <- cbind(FDMF12N[,62:ncol(FDMF12N)],F12SDN[,124:125])
Cormatrix <- cor(TraitsMatrix,method="pearson",use = "pairwise")

PAT <- NULL                                                             # Pvalue all traits
for (x in 1:ncol(TraitsMatrix)){
  PPT <- NULL                                                           # Pvalue per trait compare with the other traits
  for (n in 0:(ncol(TraitsMatrix)-1)){
    y <- n + 1
    PPT <- c(PPT,format.pval(cor.test(TraitsMatrix[,x],TraitsMatrix[,y],method = "pearson")$p.value))
  }
  PAT <- rbind(PAT,PPT)
}
rownames(PAT) <- names(TraitsMatrix)
colnames(PAT) <- names(TraitsMatrix)

write.table(Cormatrix, file = "Analysis/correlation coefficient.txt", sep = "\t")
write.table(PAT, file = "Analysis/correlation pvalue.txt", sep = "\t")