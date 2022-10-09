#!/usr/bin/Rscript

library(NMF)
data = read.table("hla-calls.txt",sep="\t",header=TRUE,row.names=1)
# drop RPEM24 
data = data[!grepl("RPEM24",rownames(data)),]
samples = read.table("sample-order-for-cnv.txt",sep="\t",header=TRUE,row.names=1)
data2 = data[match(samples$normal,rownames(data)),]
myhla = c(rep("A",13),rep("B",22),rep("C",16),rep("DP",18),rep("DQ",22),rep("DR",21))
pdf("hla-gene-usage-heatmap.pdf",height=7,width=20)
aheatmap(data2,Colv=NA,Rowv=NA,annRow = samples$group,col="PuBu:100",annCol=myhla)
dev.off()