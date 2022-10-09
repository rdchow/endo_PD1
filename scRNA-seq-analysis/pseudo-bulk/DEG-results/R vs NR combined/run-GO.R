#!/usr/bin/Rscript
#runs enrichment analysis on gene lists: idea is to compare NLR and LR differences vs NR

# NLR = epiR
# LR = mutR
# NR = NR

library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)

ref = read.table("kegg-categories-to-plot.txt",sep="\t",header=TRUE)
#mycell = "CD56-NK"
wantcells = c("CD4-Naive","CD4-Activated","CD8-Naive","CD8-Activated","CD16-NK","CD56-NK")
time = "post"
myref = ref[ref$time == time,]

#mylist = vector(mode="list", length=length(wantcells)*4)
mylist = list()
lappend <- function (lst, ...){
    lst <- c(lst, list(...))
    return(lst)
}  

counter = 1
for (mycell in wantcells){
    # epiR vs NR: 1
    data1 = read.table(paste("../DEG-results/epiRvsNR/",time,"/",time,"_NLR-vs-NR.",mycell,"_deg.txt",sep=""),sep="\t",header=TRUE,row.names=1)
    data1 = data1[data1$padj < 0.05,]
    updata1 = data1[data1$log2FoldChange > 0,]
    downdata1 = data1[data1$log2FoldChange < 0,]
    uplist1 = rownames(updata1)
    downlist1 = rownames(downdata1)
    uplist1 = bitr(uplist1,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
    downlist1 = bitr(downlist1,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

    # mutR vs NR: 2
    data2 = read.table(paste("../DEG-results/mutvsNR/",time,"/",time,"_LR-vs-NR.",mycell,"_deg.txt",sep=""),sep="\t",header=TRUE,row.names=1)
    data2 = data2[data2$padj < 0.05,]
    updata2 = data2[data2$log2FoldChange > 0,]
    downdata2 = data2[data2$log2FoldChange < 0,]
    uplist2 = rownames(updata2)
    downlist2 = rownames(downdata2)
    uplist2 = bitr(uplist2,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
    downlist2 = bitr(downlist2,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

    mylist = lappend(mylist,uplist1$ENTREZID, uplist2$ENTREZID, downlist1$ENTREZID,downlist2$ENTREZID)
}

names(mylist) = c(
    "CD4-Naive:NLR_up","CD4-Naive:LR_up","CD4-Naive:NLR_down","CD4-Naive:LR_down",
    "CD4-Activated:NLR_up","CD4-Activated:LR_up","CD4-Activated:NLR_down","CD4-Activated:LR_down",
    "CD8-Naive:NLR_up","CD8-Naive:LR_up","CD8-Naive:NLR_down","CD8-Naive:LR_down",
    "CD8-Activated:NLR_up","CD8-Activated:LR_up","CD8-Activated:NLR_down","CD8-Activated:LR_down",
    "CD16-NK:NLR_up","CD16-NK:LR_up","CD16-NK:NLR_down","CD16-NK:LR_down",
    "CD56-NK:NLR_up","CD56-NK:LR_up","CD56-NK:NLR_down","CD56-NK:LR_down")

ck = compareCluster(mylist,fun="enrichKEGG",pvalueCutoff=1)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
ckdf = as.data.frame(ck)
ckF = ckdf[ckdf$Description %in% myref$category,]
ckF$nlp = log10(ckF$qvalue)*-1
ckF = ckF[order(ckF$nlp,decreasing=FALSE),]

write.table(ckF,paste(time,"-select-KEGG.results.txt",sep=""),sep="\t",row.names=FALSE)
