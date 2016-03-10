# Winsorising outliers in your data
#
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Shijie Lyu
# last modified Feb, 2016
# first written Feb, 2016

library(robustHD)
Data <- cbind(ColA=c(rnorm(10),4),ColB=c(rnorm(10),200))
Data
Data <- apply(Data,2,winsorize)  
Data