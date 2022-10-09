#!/usr/bin/Rscript
#runs enrichment analysis on gene lists: idea is to compare pre vs post, in each pt category
# NLR = epiR
# LR = mutR
# NR = NR

library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)

ref = read.table("kegg-categories-to-plot.txt",sep="\t",header=TRUE)
#mycell = "CD56-NK"
wantcells = c("CD4-Naive","CD4-Activated","CD8-Naive","CD8-Activated","CD16-NK","CD56-NK")

mylist = list()
lappend <- function (lst, ...){
    lst <- c(lst, list(...))
    return(lst)
}  

for (mycell in wantcells){

    # NR: 1
    data1 = read.table(paste("./NR/NR_post-vs-pre.",mycell,"_deg.txt",sep=""),sep="\t",header=TRUE,row.names=1)
    data1 = data1[data1$padj < 0.05,]
    updata1 = data1[data1$log2FoldChange > 0,]
    downdata1 = data1[data1$log2FoldChange < 0,]
    uplist1 = rownames(updata1)
    downlist1 = rownames(downdata1)
    uplist1 = bitr(uplist1,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
    downlist1 = bitr(downlist1,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

    # NLR: 2
    data2 = read.table(paste("./epiR/NLR_post-vs-pre.",mycell,"_deg.txt",sep=""),sep="\t",header=TRUE,row.names=1)
    data2 = data2[data2$padj < 0.05,]
    updata2 = data2[data2$log2FoldChange > 0,]
    downdata2 = data2[data2$log2FoldChange < 0,]
    uplist2 = rownames(updata2)
    downlist2 = rownames(downdata2)
    uplist2 = bitr(uplist2,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
    downlist2 = bitr(downlist2,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

    # LR: 3
    data3 = read.table(paste("./mutR/LR_post-vs-pre.",mycell,"_deg.txt",sep=""),sep="\t",header=TRUE,row.names=1)
    data3 = data3[data3$padj < 0.05,]
    updata3 = data3[data3$log2FoldChange > 0,]
    downdata3 = data3[data3$log2FoldChange < 0,]
    uplist3 = rownames(updata3)
    downlist3 = rownames(downdata3)
    uplist3 = bitr(uplist3,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
    downlist3 = bitr(downlist3,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

    mylist = lappend(mylist,uplist1$ENTREZID, uplist2$ENTREZID, uplist3$ENTREZID, downlist1$ENTREZID, downlist2$ENTREZID,downlist3$ENTREZID)
}

names(mylist) = c(
    "CD4-Naive:NR_up","CD4-Naive:NLR_up","CD4-Naive:LR_up","CD4-Naive:NR_down","CD4-Naive:NLR_down","CD4-Naive:LR_down",
    "CD4-Activated:NR_up","CD4-Activated:NLR_up","CD4-Activated:LR_up","CD4-Activated:NR_down","CD4-Activated:NLR_down","CD4-Activated:LR_down",
    "CD8-Naive:NR_up","CD8-Naive:NLR_up","CD8-Naive:LR_up","CD8-Naive:NR_down","CD8-Naive:NLR_down","CD8-Naive:LR_down",
    "CD8-Activated:NR_up","CD8-Activated:NLR_up","CD8-Activated:LR_up","CD8-Activated:NR_down","CD8-Activated:NLR_down","CD8-Activated:LR_down",
    "CD16-NK:NR_up","CD16-NK:NLR_up","CD16-NK:LR_up","CD16-NK:NR_down","CD16-NK:NLR_down","CD16-NK:LR_down",
    "CD56-NK:NR_up","CD56-NK:NLR_up","CD56-NK:LR_up","CD56-NK:NR_down","CD56-NK:NLR_down","CD56-NK:LR_down")

ck = compareCluster(mylist,fun="enrichKEGG",pvalueCutoff=1)
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
ckdf = as.data.frame(ck)
ckF = ckdf[ckdf$Description %in% ref$Category,]
ckF$nlp = log10(ckF$qvalue)*-1
ckF = ckF[order(ckF$nlp,decreasing=FALSE),]

write.table(ckF,"preVspost-select-KEGG.results.txt",sep="\t",row.names=FALSE)
