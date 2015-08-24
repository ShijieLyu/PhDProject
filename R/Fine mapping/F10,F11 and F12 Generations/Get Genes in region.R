# Get all genes between rs16435551 and rs14492508(flanking markers of the block region)
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified July, 2015
# first written July, 2015

setwd("D:/Chicken/Rdata/FineMapping")
library("biomaRt")                                                                                              # How many genes in different regions for different traits
Mart <- useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
GenesInRegion <- getBM(attributes = c("chromosome_name","ensembl_gene_id","start_position","end_position","external_gene_name","description"), filters = "chromosomal_region", values = "4:73629630:77032597", mart = Mart)
write.table(GenesInRegion, row.names = FALSE, "Analysis/Genes_in_region.txt", sep="\t") 

