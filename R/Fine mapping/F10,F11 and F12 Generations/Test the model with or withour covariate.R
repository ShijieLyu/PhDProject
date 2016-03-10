phenotype <- rnorm(188, 63.56, 30)     #SlaughterTraits[,"TotalFat"] # 
  
LodOri <- NULL  
for (OneSNP in SNPsForAnalysis){
    model.full <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"]))+ (as.character(genotypes[,OneSNP])), REML=FALSE)
    model.null <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])), REML=FALSE)
    res <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    LodOri <- rbind(LodOri, c(-log10(res[[8]]),SNPvar))
    colnames(LodOri) <- c("Residuals","SNP","SNPvar")
  }
  rownames(LodOri) <- SNPsForAnalysis

Lodcov <- NULL  
for (OneSNP in SNPsForAnalysis){
    model.full <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])) + (as.numeric(SlaughterTraits[,"BW.nuchtern"])) + (as.character(genotypes[,OneSNP])), REML=FALSE)
    model.null <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])) + (as.numeric(SlaughterTraits[,"BW.nuchtern"])), REML=FALSE)
    resnew <- anova(model.null,model.full)
    SNPvar <- SNPvarPerc(model.full)
    Lodcov <- rbind(Lodcov, c(-log10(res[[8]]),SNPvar))
    colnames(Lodcov) <- c("Residuals","SNP","SNPvar")
  }
  rownames(Lodcov) <- SNPsForAnalysis
  
OneSNP <- "rs14490774"

phenotype <- rnorm(188, 63.5, 29.2)  #SlaughterTraits[,"TotalFat"], SlaughterTraits[,"LegnoSkin"],rnorm(188, 10, 0.1) 
genotype <- as.factor(genotypes[,OneSNP])

anova(lm(as.numeric(SlaughterTraits[,"BW.nuchtern"]) ~ genotype))
anova(lm(phenotype ~ genotype))
anova(lm(phenotype ~ as.numeric(SlaughterTraits[,"BW.nuchtern"])))
anova(lm(phenotype ~ as.numeric(SlaughterTraits[,"BW.nuchtern"]) + genotype))
anova(lm(phenotype ~ genotype + as.numeric(SlaughterTraits[,"BW.nuchtern"])))

anova(lm(lm(phenotype ~ as.numeric(SlaughterTraits[,"BW.nuchtern"]))$residuals ~ as.factor(genotypes[,OneSNP])))


phenotype <- as.numeric(SlaughterTraits[,"LegnoSkin"])
covariate <- as.numeric(SlaughterTraits[,"BW.nuchtern"])
model.full <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])) + covariate + (as.factor(genotypes[,OneSNP])), REML=FALSE)
model.null <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])) + covariate , REML=FALSE)
anova(model.null,model.full)

phenotype <- SlaughterTraits[,"BW.nuchtern"]
model.full <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])) + (as.factor(genotypes[,OneSNP])), REML=FALSE)
model.null <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])), REML=FALSE)
anova(model.null,model.full)

phenotype <- SlaughterTraits[,"BW.nuchtern"]
covariate <- as.numeric(SlaughterTraits[,"TotalFat"])
model.full <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])) + covariate + (as.factor(genotypes[,OneSNP])), REML=FALSE)
model.null <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])) + covariate , REML=FALSE)
anova(model.null,model.full)

phenotype <- SlaughterTraits[,"TotalFat"]
covariate <- as.numeric(SlaughterTraits[,"BW.nuchtern"])
model.full <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])) + covariate + (as.factor(genotypes[,OneSNP])), REML=FALSE)
model.null <- lmer(phenotype ~ as.factor(SlaughterTraits[,"Batch"]) + (1|as.factor(SlaughterTraits[,"Parents"])) + covariate , REML=FALSE)
anova(model.null,model.full)