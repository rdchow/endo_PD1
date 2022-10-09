#!/usr/bin/Rscript
# boxplots of cell frequencies in pre vs post
# using fine_t cell clusters

library(magrittr)
library(dplyr)
library(reshape2)
meta = read.table("filtdata-all-coarseT_reann.final.metadata.txt",sep="\t",header=TRUE,row.names=1)
info = read.table("scrna-sample-info.txt",sep="\t",header=TRUE,row.names=1)

cellcts = as.data.frame(table(meta$orig.ident))
cellcts = cellcts[!grepl("eric",cellcts$Var1) & !grepl("HD",cellcts$Var1),]
cellcts = cellcts[order(cellcts$Freq),]
cellcts$sample = info[as.vector(cellcts$Var1),"simpleCategory"]

timepoints = strsplit(as.character(cellcts$Var1),"C") %>% sapply(extract2,2)
timepoints[timepoints %in% c(3,5)] = "post"
timepoints[timepoints %in% c(1)] = "pre"
cellcts$timepoint = timepoints

cellcts$category = paste(cellcts$sample,cellcts$timepoint,sep="_")
cellcts$category = factor(cellcts$category,levels=c("NR_pre","NR_post","epiR_pre","epiR_post","mutR_pre","mutR_post"))
mycols = c("#A3D0D6","#6EB1DA","#FFD479","#C8760E","#F37C79","#B6312F")

library(ggpubr)
library(ggplot2)
pdf("barplot-cell-counts.pdf",width=9,height=5,useDingbats=FALSE)
ggbarplot(cellcts,x="Var1",y="Freq",fill="category",color=NA,palette=mycols,order=cellcts$Var1,ylab="Cell count")+ theme(axis.text.x=element_text(angle=90,hjust=1))+geom_hline(yintercept = 1000,lty = 2,color="gray44")
dev.off()

