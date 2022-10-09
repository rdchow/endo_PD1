#!/usr/bin/Rscript

# NLR = epiR
# LR = mutR
# NR = NR
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)

myref = read.table("kegg-categories-to-plot.signif.txt",sep="\t",header=TRUE)
wantcells = c("CD4-Naive","CD4-Activated","CD8-Naive","CD8-Activated","CD16-NK","CD56-NK")

mylist = list()
lappend <- function (lst, ...){
    lst <- c(lst, list(...))
    return(lst)
}  

for (mycell in wantcells){
    # NLR vs LR: pre [1]
    data1 = read.table(paste("pre/pre_NLR-vs-LR.",mycell,"_deg.txt",sep=""),sep="\t",header=TRUE,row.names=1)
    data1 = data1[data1$padj < 0.05,]
    updata1 = data1[data1$log2FoldChange > 0,]
    downdata1 = data1[data1$log2FoldChange < 0,]
    uplist1 = rownames(updata1)
    downlist1 = rownames(downdata1)
    uplist1 = bitr(uplist1,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
    downlist1 = bitr(downlist1,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")

    # NLR vs LR: post [2]
    data2 = read.table(paste("post/post_NLR-vs-LR.",mycell,"_deg.txt",sep=""),sep="\t",header=TRUE,row.names=1)
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
    "CD4-Naive:NLR_pre_up","CD4-Naive:NLR_post_up","CD4-Naive:NLR_pre_down","CD4-Naive:NLR_post_down",
    "CD4-Activated:NLR_pre_up","CD4-Activated:NLR_post_up","CD4-Activated:NLR_pre_down","CD4-Activated:NLR_post_down",
    "CD8-Naive:NLR_pre_up","CD8-Naive:NLR_post_up","CD8-Naive:NLR_pre_down","CD8-Naive:NLR_post_down",
    "CD8-Activated:NLR_pre_up","CD8-Activated:NLR_post_up","CD8-Activated:NLR_pre_down","CD8-Activated:NLR_post_down",
    "CD16-NK:NLR_pre_up","CD16-NK:NLR_post_up","CD16-NK:NLR_pre_down","CD16-NK:NLR_post_down",
    "CD56-NK:NLR_pre_up","CD56-NK:NLR_post_up","CD56-NK:NLR_pre_down","CD56-NK:NLR_post_down")

ck = compareCluster(mylist,fun="enrichKEGG",pvalueCutoff=1,organism="hsa")
ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
#dotplot(ck)
ckdf = as.data.frame(ck)
    ckdf$nlp = log10(ckdf$qvalue)*-1
    ckdf = ckdf[order(ckdf$nlp,decreasing=TRUE),]
    write.table(ckdf,"epiR-vs-mutR-all-KEGG.results.txt",sep="\t",row.names=FALSE)

ckF = ckdf[ckdf$Description %in% myref$category,]
ckF$nlp = log10(ckF$qvalue)*-1
ckF = ckF[order(ckF$nlp,decreasing=FALSE),]

write.table(ckF,"epiR-vs-mutR-select-KEGG.results.txt",sep="\t",row.names=FALSE)
