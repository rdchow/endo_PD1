#!/usr/bin/Rscript
# plot gene expression violins at cell neighborhood-level

library(miloR)
library(ggbeeswarm)
library(ggplot2)
library(dplyr)
library(data.table)

# the expression data is from the milo-exec.preTx.cellIsolatep2.R script.
# expression files too large to include in this repo
data=data.frame(fread("tdat.milo_mutR-NR_preTx.selectGenes.expr.txt",sep="\t",header=TRUE),row.names=1)
info = read.table("mutR-NR.preTx.milo-results.txt",sep="\t",header=TRUE,row.names=1)

#data=data.frame(fread("tdat.milo_mutR-NR_postTx.selectGenes.expr.txt",sep="\t",header=TRUE),row.names=1)
#info = read.table("mutR-NR.postTx.milo-results.txt",sep="\t",header=TRUE,row.names=1)

info = info[info$finalIdent == "CD8 Activated" ,]
rownames(info) = paste("X",rownames(info),sep="")
data = data[,colnames(data) %in% rownames(info)]
info = info[match(colnames(data),rownames(info)),]
all(rownames(info) == colnames(data))

mygenes = c("KLRG1","EOMES","LAG3","PDCD1","CTLA4","KLRC1") #pre-tx
#mygenes = c("GZMB","GNLY","TNF","HAVCR2","PDCD1","LAG3") #post-tx
dataf = data[mygenes,]
info$sig = "no"
info[info$SpatialFDR < 0.1 & info$logFC > 0,]$sig = "yes"

library(reshape2)
meltdata = melt(t(dataf))
colnames(meltdata) = c("nhood","gene","expr")
meltdata$sig = info[meltdata$nhood,"sig"]

library(ggpubr)
library(ggplot2)

pdf("mutR-vs-NR.preTx.selectGenes.violin.pdf",height=5,width=6,useDingbats=FALSE)
#pdf("mutR-vs-NR.postTx.selectGenes.violin.pdf",height=5,width=6,useDingbats=FALSE)
p = ggviolin(meltdata,x="sig",y="expr",color="sig",palette=c("gray66","orangered3"),draw_quantiles=c(0.5),trim=TRUE,add="mean_se")
facet(p,facet.by="gene",scales="free_y",ncol=3)
dev.off()


