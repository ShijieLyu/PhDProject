# correlation analysis foe the fat content using the wet and dry methods
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Shijie Lyu
# last modified Sep, 2015
# first written Sep, 2015


setwd("D:/Chicken/Rdata/FineMapping")
FDM <- read.table("RawData/Fettmessung-QTl.txt", na.strings="", header=TRUE, sep="\t")
QTLdataAll <- read.table("Analysis/QTLdataAll.txt", sep = "\t",header=TRUE, colClasses=c(rep("character",1), rep("factor",9),rep("numeric",9), rep("character", 9), rep("factor",6)))                                                   # Fat measurement with Dry method
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

FDM10 <- FDM[which(FDM[,"Gene.ration"] == "F10"),]
F10NID <- paste0("C",FDM10[,"Tier.Nr."])
FDM10 <- cbind(ID.Nr = F10NID, FDM10)

FDM11 <- FDM[which(FDM[,"Gene.ration"] == "F11"),]
F11NID <- paste0("B",FDM11[,"Tier.Nr."])
FDM11 <- cbind(ID.Nr = F11NID, FDM11)

FDM12 <- FDM[which(FDM[,"Gene.ration"] == "F12"),]
FDM12 <- cbind(ID.Nr = FDM12[,"Tier.Nr."], FDM12)

FDMAll <- rbind(FDM12,FDM11,FDM10)
FDMAll <- FDMAll[match(QTLdataAll[,"ID.Nr"],FDMAll[,"ID.Nr"]),]
FDMAll <- FDMAll[which(!is.na(FDMAll[,"ID.Nr"])),]

QTLdataAll <- QTLdataAll[which(QTLdataAll[,"ID.Nr"] %in% FDMAll[,"ID.Nr"]),]
QTLdataAll <- cbind(QTLdataAll,FDMAll[,64:ncol(FDMAll)])
QTLdataAll <- QTLdataAll[which(!is.na(QTLdataAll[,"Gew_20Wo"])),]

# Association analysis
SNPvarPerc <- function(mF){
  X <- getME(mF,"X")                                               # Get the fixed-effects model matrix
  TotalFixed <- 0
  for (x in 2:length(fixef(mF))){  
    TotalFixed <- TotalFixed + (fixef(mF)[x] * X[, x])  
  }
  SNPFixed <- 0
  for (x in grep("OneSNP", names(fixef(mF)))){     
    SNPFixed <- SNPFixed + (fixef(mF)[x] * X[, x])
  }
  FixedTotalVar <- var(TotalFixed)                                 # Get the total fixed effects variance 
  FixedSNPVar <- var(SNPFixed)                                     # Get the specific fixed effect variance
  TotalVar <- FixedTotalVar + sum(data.frame(VarCorr(mF))[,4])     # Calculate the total variance with fixed and random
  Perc <- round(FixedSNPVar/TotalVar,3)*100
  return(Perc)
}

library(lme4)
FatContentTraits  <- names(FDMAll)[64:ncol(FDMAll)]
SNPsForAnalysis <- c("rs10725580", "rs315966269", "rs313283321", "rs16435551", "rs14490774", "rs314961352", "rs318175270", "rs14492508","rs312839183")

FCLod <- vector("list", length(FatContentTraits))                      # FC = Fat content
names(FCLod) <- FatContentTraits
for (Phenotype in FatContentTraits){                                                     
  EachTraitLod <- NULL
  if(length(which(!is.na(QTLdataAll[,Phenotype]))) == 0){
        EachTraitLod <- rbind(EachTraitLod, NA)
        colnames(EachTraitLod) <- "SNP"
  }else{
  for (OneSNP in SNPsForAnalysis){
    idx <- which(!is.na(QTLdataAll[,Phenotype]))
    model.full <- lmer(QTLdataAll[,Phenotype][idx] ~ QTLdataAll[,"Batch"][idx] + (1|QTLdataAll[,"Parents"][idx]) + (as.numeric(QTLdataAll[,"Gew_20Wo"])[idx]) + QTLdataAll[,OneSNP][idx], REML=FALSE) # qqnorm(resid(model.full)) 
    model.null <- lmer(QTLdataAll[,Phenotype][idx] ~ QTLdataAll[,"Batch"][idx] + (1|QTLdataAll[,"Parents"][idx]) + (as.numeric(QTLdataAll[,"Gew_20Wo"])[idx]), REML=FALSE)
    res <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    EachTraitLod <- rbind(EachTraitLod, c(-log10(res[[8]]),SNPvar))
    colnames(EachTraitLod) <- c("Residuals","SNP","SNPvar")
  }
  rownames(EachTraitLod) <- SNPsForAnalysis
  FCLod[[Phenotype]] <- EachTraitLod
  }
}


Threshold <- -log10(0.05/(length(FCLod)*length(SNPsForAnalysis)))
SigThreshold <- -log10(0.01/(length(FCLod)*length(SNPsForAnalysis)))

SNPsinfo <- read.table("RawData/SNPsinfo.txt",header=TRUE,sep="\t")

#tiff("Analysis/SlaughterTrait-AllinOne.tif", res = 300,width = 2100, height = 2000, compression = "lzw")
#pdf("Analysis/SlaughterTrait-AllinOne.pdf")
par(mai = c(2.2, 1, 1, 1))
plot(x=c(as.numeric(SNPsinfo[1,3]-100000),as.numeric(SNPsinfo[9,3]+100000)), y=c(0,10), t="n", ylab="LOD Score", xlab="Physical Position (Mb)", xaxt="n")
for(x in 1:10){
  points(x= SNPsinfo[,"Location"], y=FCLod[[x]][,2],t="l",col=rainbow(10)[x], lwd=1.8)
}
abline(h = Threshold , lty=2, lwd=1.9)
abline(h = SigThreshold , lty=1, lwd=1.9)
#text(x=SNPsinfo[1,3]+600000,y=Threshold+0.25, labels = "5% threshold")
#text(x=SNPsinfo[1,3]+600000,y=SigThreshold+0.25, labels = "1% threshold")
axis(1, at=seq(69000000,78000000,1000000), c("69","70","71","72","73","74","75","76","77","78"), las=1)

points(x=SNPsinfo[,3], y = rep(-0.3,length(SNPsinfo[,3])), pch=17)
#axis(1, at=SNPsinfo[,3], SNPsinfo[,"Markers"], cex.axis=0.4,las=2)

legend(69250000,-2.8, c("CW","HW","NW","WW","LNS"),lty=1,col=(rainbow(10)[1:5]),horiz = TRUE,xpd = TRUE, bty = "n", lwd=2.2)
legend(68100000,-3.5, c("LNSB","BNS","SubcAT","ViscAT","WAT"),lty=1,col=(rainbow(10)[6:10]),horiz = TRUE,xpd = TRUE, bty = "n", lwd=2.2)
#dev.off()


