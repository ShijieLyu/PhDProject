# HeatMap construction
# copyright (c) 2014-2020 - Brockmann group - HU Berlin, Danny Arends & Shijie Lyu
# last modified Feb, 2015
# first written Feb, 2015

setwd("C:/Data") 
genes <- read.table("frequence.txt", sep="\t", header=TRUE)   
rnames <- genes[,1]
frequence_matrix <- data.matrix(genes[,2:ncol(genes)])
rownames(frequence_matrix) <- rnames
frequence_matrix <- frequence_matrix[rev(seq_len(nrow(frequence_matrix))),]
col <- heat.colors(256)
col <- col[order(col,decreasing=T)]
pdf("frequence_heatmap.pdf")   
frequence_heatmap <- heatmap(frequence_matrix, Rowv=NA, Colv=NA, col = col, scale="column", margins=c(1,20),cexCol = 0.95)
dev.off()