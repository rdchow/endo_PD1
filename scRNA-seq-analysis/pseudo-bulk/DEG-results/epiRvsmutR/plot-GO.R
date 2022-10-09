#!/usr/bin/Rscript
# NLR = epiR
# LR = mutR
# NR = NR

#plots KEGG enrichment
library(ggplot2)
library(reshape2)
library(NMF)

data=read.table("epiR-vs-mutR-select-KEGG.results.rf.txt",sep="\t",header=TRUE)
data$cellGroup = paste(data$cellType,data$Time,sep=":")
dataf = data[data$Pathway != "Ribosome",]
dataf$cellGroup = factor(dataf$cellGroup,levels=c(
    "CD4-Naive:pre","CD4-Naive:post","CD4-Activated:pre","CD4-Activated:post","CD8-Naive:pre","CD8-Naive:post","CD8-Activated:pre","CD8-Activated:post","CD16-NK:pre","CD16-NK:post","CD56-NK:pre","CD56-NK:post"))

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
removepathways = c("Antigen processing and presentation","Apoptosis","Th1 and Th2 cell differentiation","Mitophagy - animal","Cell adhesion molecules","Phosphatidylinositol signaling system")
dataf2 = dataf[!dataf$Pathway %in% removepathways,]

pdf("epiRvsmutR-KEGG.dot.pdf",height=6,width=9,useDingbats=FALSE)
ggplot(dataf2,aes(x=cellGroup,y=Pathway,size=nnlp,color=nlp))+geom_point()+theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))+scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=FALSE)+scale_color_distiller(palette="BrBG",limits=c(-3,3))
dev.off()
