#!/usr/bin/Rscript
# boxplots of cell frequencies in pre vs post

library(magrittr)
library(dplyr)
library(Seurat)
library(reshape2)
meta = read.table("filtdata-all-coarseT_reann.final.metadata.txt",sep="\t",header=TRUE,row.names=1)

#remove PEM23C1, PEM20C5, and PEM21C1 because low cell numbers (< 1000 cells)
removesamples = c("PEM23C1","PEM20C5","PEM21C1")
meta = meta[! meta$orig.ident %in% removesamples,]

#binarize to pre-tx vs post-tx
meta$timeBin = ""
meta[meta$timepoint == 1,]$timeBin = "pre"
meta[meta$timepoint == 3 | meta$timepoint == 5,]$timeBin = "post"
meta[meta$timepoint == "h",]$timeBin = "h"

meta$condition = paste(meta$timeBin,meta$clinical,sep="_")
meta$condition = factor(meta$condition,levels=c("h_Healthy","pre_NR","post_NR","pre_epiR","post_epiR","pre_mutR","post_mutR"))

meta$sample = paste(meta$patient,meta$timepoint,meta$clinical,sep="_")

etable = table(meta$sample,meta$finalIdent)
etable = etable[,colSums(etable) > 0]
etablePct = prop.table(etable, margin = 1)
etablePct = as.data.frame.matrix(etablePct)
class = strsplit(rownames(etablePct),"_") %>% sapply(extract2,3)

#binarize timepoints
timepoints = strsplit(rownames(etablePct),"_") %>% sapply(extract2,2)
timepoints[timepoints %in% c(3,5)] = "post"
timepoints[timepoints %in% c(1)] = "pre"

etablePct$condition = paste(timepoints,class,sep="_")

emelt = melt(etablePct)
colnames(emelt)=c("Condition","Cluster","Frequency")

# filter out healthy samples
emelt = emelt[emelt$Condition != "h_Healthy",]
emelt$Condition=factor(emelt$Condition,levels=c("pre_NR","post_NR","pre_epiR","post_epiR","pre_mutR","post_mutR"))

#coarse T cells
emelt$Cluster = factor(emelt$Cluster,levels=c('CD4 Naive','CD4 Activated','Treg','CD8 Naive','CD8 Activated','CD8 Exhausted','Cycling','GDT','MAIT','Naive B Cells','Activated B Cells','Memory B Cells','Plasma Cells','CD16 NK','CD56 NK','CD14 Monocytes','Intermediate Monocytes','CD16 Monocytes','cDC','pDC','Hematopoetic Progenitors'))

library(ggpubr)
pdf("cell-frequency.boxplots.ttest.coarseT.pdf",height=10,width=20,useDingbats=FALSE)
ggboxplot(emelt,x="Condition",y="Frequency",fill="Condition",facet.by="Cluster",ncol=7,scales="free_y",palette=c("#A3D0D6","#6EB1DA","#FFD479","#C8760E","#F37C79","#B6312F"),xlab=FALSE,add="jitter",shape=21,size=0.1)+stat_compare_means(method="t.test",comparisons=list(c("pre_NR","post_NR"),c("pre_epiR","post_epiR"),c("pre_NR","pre_epiR"),c("pre_mutR","post_mutR"),c("pre_NR","pre_mutR")))
dev.off()
