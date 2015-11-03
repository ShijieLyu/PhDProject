# correlation analysis foe the fat content using the wet and dry methods
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Shijie Lyu
# last modified Sep, 2015
# first written Sep, 2015

setwd("D:/Chicken/Rdata/Chicken Fat measurement")
F12SD <- read.table("D:/Chicken/Rdata/FineMapping/RawData/F12_Slaught_Phenotypes.txt", na.strings=c("","?"),header=TRUE, sep="\t")   # F12 Slaughter data
FDM <- read.table("RawData/Fettmessung-withTransformed.txt", na.strings="", header=TRUE, sep="\t")#, colClasses=c(rep("character",4), rep("numeric",8),rep("character",1),rep("numeric",72)))                                                   # Fat measurement with Dry method
FDM <- cbind(FDM, BMFOT = rep(NA, nrow(FDM)), BTFOT = rep(NA, nrow(FDM)), BMAFOT = rep(NA, nrow(FDM)), BMBFOT = rep(NA, nrow(FDM)), 
                  LUpFOT = rep(NA, nrow(FDM)), LLoFOT = rep(NA, nrow(FDM)), LMAFOT = rep(NA, nrow(FDM)), LMBFOT = rep(NA, nrow(FDM)),
                  BMFT = rep(NA, nrow(FDM)), BTFT = rep(NA, nrow(FDM)), BMAFT = rep(NA, nrow(FDM)), BMBFT = rep(NA, nrow(FDM)), 
                  LUpFT = rep(NA, nrow(FDM)), LLoFT = rep(NA, nrow(FDM)), LMAFT = rep(NA, nrow(FDM)), LMBFT = rep(NA, nrow(FDM)))  
                  
toN <- function(x){as.numeric(x)}  
  
for (x in 1:nrow(FDM)){
  FDM[x,"BMFOT"]  <- mean(c(toN(FDM[x,"Breast.Mitte.Trans1.fat."]),toN(FDM[x,"Breast.Mitte.Trans2.fat."]),toN(FDM[x,"Breast.Mitte.Trans3.fat."])),na.rm = TRUE)*0.01      # BMFO = Original breast middle fat content, which means to use the percentage in dry muscle directly.
  FDM[x,"BTFOT"]  <- mean(c(toN(FDM[x,"Breast.Spitze.Trans1.fat."]),toN(FDM[x,"Breast.Spitze.Trans2.fat."]),toN(FDM[x,"Breast.Spitze.Trans3.fat."])),na.rm = TRUE)*0.01   # BTFO = Original breast tip fat
  FDM[x,"BMAFOT"] <- mean(c(toN(FDM[x,"Breast.HomoA.Trans1.fat."]),toN(FDM[x,"Breast.HomoA.Trans2.fat."]),toN(FDM[x,"Breast.HomoA.Trans3.fat."])),na.rm = TRUE)*0.01         # BMAFO = Original breast muscle A fat  
  FDM[x,"BMBFOT"] <- mean(c(toN(FDM[x,"Breast.HomoB.Trans1.fat."]),toN(FDM[x,"Breast.HomoB.Trans2.fat."]),toN(FDM[x,"Breast.HomoB.Trans3.fat."])),na.rm = TRUE)*0.01         # BMBFO = Original breast muscle B fat 
  FDM[x,"LUpFOT"] <- mean(c(toN(FDM[x,"Leg.Up.Trans1.fat."]),toN(FDM[x,"Leg.Up.Trans2.fat."]),toN(FDM[x,"Leg.Up.Trans3.fat."])),na.rm = TRUE)*0.01         # LUpFO = Original leg upper part fat 
  FDM[x,"LLoFOT"] <- mean(c(toN(FDM[x,"Leg.Lower.Trans1.fat."]),toN(FDM[x,"Leg.Lower.Trans2.fat."]),toN(FDM[x,"Leg.Lower.Trans3.fat."])),na.rm = TRUE)*0.01      # LLoFO = Original leg lower part fat
  FDM[x,"LMAFOT"] <- mean(c(toN(FDM[x,"Leg.HomoA.Trans1.fat."]),toN(FDM[x,"Leg.HomoA.Trans2.fat."]),toN(FDM[x,"Leg.HomoA.Trans3.fat."])),na.rm = TRUE)*0.01         # LMAFO = Original leg muscle A fat  
  FDM[x,"LMBFOT"] <- mean(c(toN(FDM[x,"Leg.HomoB.Trans1.fat."]),toN(FDM[x,"Leg.HomoB.Trans2.fat."]),toN(FDM[x,"Leg.HomoB.Trans3.fat."])),na.rm = TRUE)*0.01         # LMBFO = Original leg muscle B fat     
  FDM[x,"BMFT"]   <- mean(c((toN(FDM[x,"Breast.Mitte.Trans1.fat."])*0.01*FDM[x,"Brust.Mitte.1"])/FDM[x,"Brust.Mitte"],(toN(FDM[x,"Breast.Mitte.Trans2.fat."])*0.01*FDM[x,"Brust.Mitte.2"])/FDM[x,"Brust.Mitte"],(toN(FDM[x,"Breast.Mitte.Trans3.fat."])*0.01*FDM[x,"Brust.Mitte.3"])/FDM[x,"Brust.Mitte"]),na.rm = TRUE)           # BMF = breast middle fat in wet muscle
  FDM[x,"BTFT"]   <- mean(c((toN(FDM[x,"Breast.Spitze.Trans1.fat."])*0.01*FDM[x,"Brust.Spitze.1"])/FDM[x,"Brust.Spitze"],(toN(FDM[x,"Breast.Spitze.Trans2.fat."])*0.01*FDM[x,"Brust.Spitze.2"])/FDM[x,"Brust.Spitze"],(toN(FDM[x,"Breast.Spitze.Trans3.fat."])*0.01*FDM[x,"Brust.Spitze.3"])/FDM[x,"Brust.Spitze"]),na.rm = TRUE)  # BTF = breast tip fat in wet muscle
  FDM[x,"BMAFT"]  <- mean(c((toN(FDM[x,"Breast.HomoA.Trans1.fat."])*0.01*FDM[x,"Brust.MusA.1"])/FDM[x,"Brust.MusA"],(toN(FDM[x,"Breast.HomoA.Trans2.fat."])*0.01*FDM[x,"Brust.MusA.2"])/FDM[x,"Brust.MusA"],(toN(FDM[x,"Breast.HomoA.Trans3.fat."])*0.01*FDM[x,"Brust.MusA.3"])/FDM[x,"Brust.MusA"]),na.rm = TRUE)                    # BMAF = breast muscle A fat in wet muscle
  FDM[x,"BMBFT"]  <- mean(c((toN(FDM[x,"Breast.HomoB.Trans1.fat."])*0.01*FDM[x,"Brust.MusB.1"])/FDM[x,"Brust.MusB"],(toN(FDM[x,"Breast.HomoB.Trans2.fat."])*0.01*FDM[x,"Brust.MusB.2"])/FDM[x,"Brust.MusB"],(toN(FDM[x,"Breast.HomoB.Trans3.fat."])*0.01*FDM[x,"Brust.MusB.3"])/FDM[x,"Brust.MusB"]),na.rm = TRUE)                    # BMBF = breast muscle B fat in wet muscle
  FDM[x,"LUpFT"]  <- mean(c((toN(FDM[x,"Leg.Up.Trans1.fat."])*0.01*FDM[x,"Keule.oben.1"])/FDM[x,"Keule.oben"],(toN(FDM[x,"Leg.Up.Trans2.fat."])*0.01*FDM[x,"Keule.oben.2"])/FDM[x,"Keule.oben"],(toN(FDM[x,"Leg.Up.Trans3.fat."])*0.01*FDM[x,"Keule.oben.3"])/FDM[x,"Keule.oben"]),na.rm = TRUE)                    # LUpF = leg upper part fat in wet muscle
  FDM[x,"LLoFT"]  <- mean(c((toN(FDM[x,"Leg.Lower.Trans1.fat."])*0.01*FDM[x,"Keule.unten.1"])/FDM[x,"Keule.unten"],(toN(FDM[x,"Leg.Lower.Trans2.fat."])*0.01*FDM[x,"Keule.unten.2"])/FDM[x,"Keule.unten"],(toN(FDM[x,"Leg.Lower.Trans3.fat."])*0.01*FDM[x,"Keule.unten.3"])/FDM[x,"Keule.unten"]),na.rm = TRUE)           # LLoF = leg lower part fat in wet muscle
  FDM[x,"LMAFT"]  <- mean(c((toN(FDM[x,"Leg.HomoA.Trans1.fat."])*0.01*FDM[x,"Keule.MusA.1"])/FDM[x,"Keule.MusA"],(toN(FDM[x,"Leg.HomoA.Trans2.fat."])*0.01*FDM[x,"Keule.MusA.2"])/FDM[x,"Keule.MusA"],(toN(FDM[x,"Leg.HomoA.Trans3.fat."])*0.01*FDM[x,"Keule.MusA.3"])/FDM[x,"Keule.MusA"]),na.rm = TRUE)                    # LMAF = Original leg lower part fat
  FDM[x,"LMBFT"]  <- mean(c((toN(FDM[x,"Leg.HomoB.Trans1.fat."])*0.01*FDM[x,"Keule.MusB.1"])/FDM[x,"Keule.MusB"],(toN(FDM[x,"Leg.HomoB.Trans2.fat."])*0.01*FDM[x,"Keule.MusB.2"])/FDM[x,"Keule.MusB"],(toN(FDM[x,"Leg.HomoB.Trans3.fat."])*0.01*FDM[x,"Keule.MusB.3"])/FDM[x,"Keule.MusB"]),na.rm = TRUE)        # LMBF = leg lower part fat     
}

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

TraitsMatrix <- cbind(FDMF12N[,86:ncol(FDMF12N)],F12SDN[,124:125])
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

write.table(Cormatrix, file = "Analysis/correlation coefficient-transformed.txt", sep = "\t")
write.table(PAT, file = "Analysis/correlation pvalue-transformed.txt", sep = "\t")