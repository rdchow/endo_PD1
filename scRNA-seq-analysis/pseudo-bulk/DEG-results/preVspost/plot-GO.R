#!/usr/bin/Rscript
# NLR = epiR
# LR = mutR
# NR = NR

#plots KEGG enrichment
library(ggplot2)
library(reshape2)
library(NMF)

data=read.table("preVspost-select-KEGG.results.rf.txt",sep="\t",header=TRUE)
data$cellGroup = paste(data$cellType,data$group,sep=":")
dataf = data[data$Pathway != "Ribosome",]
dataf$cellGroup = factor(dataf$cellGroup,levels=c(
    "CD4-Naive:NR","CD4-Naive:NLR","CD4-Naive:LR",
    "CD4-Activated:NR","CD4-Activated:NLR","CD4-Activated:LR",
    "CD8-Naive:NR","CD8-Naive:NLR","CD8-Naive:LR",
    "CD8-Activated:NR","CD8-Activated:NLR","CD8-Activated:LR",
    "CD16-NK:NR","CD16-NK:NLR","CD16-NK:LR",
    "CD56-NK:NR","CD56-NK:NLR","CD56-NK:LR"))

#only keep significant results
dataf = dataf[abs(dataf$nlp) > log10(0.05)*-1,]

#compress dynamic range to max abs(q) = 3
dataf[dataf$nlp < -3,]$nlp = -3
dataf[dataf$nlp > 3,]$nlp = 3
dataf$nnlp = abs(dataf$nlp)

# cluster the pathways first to get the desired pathway order
unmelt = dcast(dataf,Pathway ~ cellGroup,value.var="nlp",fill=0)
rownames(unmelt) = unmelt$Pathway
unmelt[,1] = NULL
a=aheatmap(unmelt,Colv=NA,distfun="manhattan")
pathwayOrder = rownames(unmelt)[a$rowInd]
dataf$Pathway = factor(dataf$Pathway,levels=pathwayOrder)

#drop stuff
#dataf = dataf[dataf$Pathway != "Th1 and Th2 cell differentiation",]

pdf("preVspost.KEGG.dot.pdf",height=7,width=10,useDingbats=FALSE)
ggplot(dataf,aes(x=cellGroup,y=Pathway,size=nnlp,color=nlp))+geom_point()+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))+scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE)+scale_color_distiller(palette="BrBG",limits=c(-3,3))
dev.off()
