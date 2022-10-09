#!/usr/bin/Rscript
#compared DEGcounts in mutR vs NR :: epiR vs NR

library(NMF)
library(reshape2)
library(ggplot2)
library(ggpubr)
data=read.table("merged-DEGcounts-vs-NR.order1.txt",sep="\t",header=TRUE,row.names=1)
data$cell = rownames(data)
predata = data[,c(1:4,9)]
postdata = data[,5:9]

premelt = melt(predata)
postmelt = melt(postdata)

premelt$cell = factor(premelt$cell,levels=c('CD4-Naive','CD4-Activated','Treg','CD8-Naive','CD8-Activated','CD8-Exhausted','Cycling','GDT','MAIT','Naive-B-Cells','Activated-B-Cells','Memory-B-Cells','Plasma-Cells','CD16-NK','CD56-NK','CD14-Monocytes','Intermediate-Monocytes','CD16-Monocytes','cDC','pDC'))

postmelt$cell = factor(postmelt$cell,levels=c('CD4-Naive','CD4-Activated','Treg','CD8-Naive','CD8-Activated','CD8-Exhausted','Cycling','GDT','MAIT','Naive-B-Cells','Activated-B-Cells','Memory-B-Cells','Plasma-Cells','CD16-NK','CD56-NK','CD14-Monocytes','Intermediate-Monocytes','CD16-Monocytes','cDC','pDC'))#'Hematopoetic-Progenitors'))#'Hematopoetic-Progenitors'))


wantcells = c('CD4-Naive','CD4-Activated','Treg','CD8-Naive','CD8-Activated','CD8-Exhausted','CD16-NK','CD56-NK')

premeltF =premelt[premelt$cell %in% wantcells,]
postmeltF = postmelt[postmelt$cell %in% wantcells,]

pdf("pre.vsNR.DEGcts.pdf",height=4,width=7,useDingbats=FALSE)
ggbarplot(premeltF,x="variable",y="value",fill="variable",facet.by="cell",ncol=4,scales="free_y",palette=c("#BD4230","#BB652B","#4186C7","#83A9C3"),xlab=FALSE)+ scale_y_continuous(limits = c(0,NA))
dev.off()

pdf("post.vsNR.DEGcts.pdf",height=4,width=7,useDingbats=FALSE)
ggbarplot(postmeltF,x="variable",y="value",fill="variable",facet.by="cell",ncol=4,scales="free_y",palette=c("#BD4230","#BB652B","#4186C7","#83A9C3"),xlab=FALSE)+ scale_y_continuous(limits = c(0,NA))
dev.off()