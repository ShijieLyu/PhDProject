# Why is the additive genetic relationship of full brothers and full sisters (full sibs) 0.5?
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends
# last modified Apr, 2016
# first written Apr, 2016

createGenotype <- function(markers = 1000){
  genotype <- NULL
  for(x in 1:markers){
    if(runif(1) < 0.5){
      genotype <- c(genotype, "F1")
    }else{
      genotype <- c(genotype, "F2")
    }
    if(runif(1) < 0.5){
      genotype <- c(genotype, "M1")
    }else{
      genotype <- c(genotype, "M2")
    }
  }
  return(genotype)
}

equality <- NULL
for(x in 1:1000){
  c1 <- createGenotype()
  c2 <- createGenotype()
  equality <- c(equality, length(which(c1 == c2)) / length(c1))
}
hist(equality)

