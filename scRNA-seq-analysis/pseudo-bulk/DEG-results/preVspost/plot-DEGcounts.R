#!/usr/bin/Rscript
#compared DEGcounts in LR vs NR :: NLR vs NR
# NLR = epiR
# LR = mutR
# NR = NR

library(NMF)
library(reshape2)
library(ggplot2)
library(ggpubr)
data=read.table("preVspost.DEGcounts.order1.txt",sep="\t",header=TRUE,row.names=1)
data$cell = rownames(data)
meltdata = melt(data)

meltdata$cell = factor(meltdata$cell,levels=c('CD4-Naive','CD4-Activated','Treg','CD8-Naive','CD8-Activated','CD8-Exhausted','Cycling','GDT','MAIT','Naive-B-Cells','Activated-B-Cells','Memory-B-Cells','Plasma-Cells','CD16-NK','CD56-NK','CD14-Monocytes','Intermediate-Monocytes','CD16-Monocytes','cDC','pDC','Hematopoetic-Progenitors'))


wantcells = c('CD4-Naive','CD4-Activated','Treg','CD8-Naive','CD8-Activated','CD8-Exhausted','CD16-NK','CD56-NK')

meltdataF = meltdata[meltdata$cell %in% wantcells,]

pdf("preVspost.all.DEGcts.pdf",height=8,width=5,useDingbats=FALSE)
ggbarplot(meltdataF,x="variable",y="value",fill="variable",facet.by="cell",ncol=2,scales="free_y",palette=c("#5B0D28","#BD4230","#BB652B","#1F3D78","#4186C7","#83A9C3"),xlab=FALSE)+ scale_y_continuous(limits = c(0,NA))
dev.off()
