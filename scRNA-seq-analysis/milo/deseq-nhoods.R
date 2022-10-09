#!/usr/bin/Rscript

library(ggplot2)
library(dplyr)
library(data.table)

#pre-tx
#data=data.frame(fread("tdat.milo_mutR-NR_preTx.selectGenes.expr.txt",sep="\t",header=TRUE),row.names=1)
#info = read.table("mutR-NR.preTx.milo-results.txt",sep="\t",header=TRUE,row.names=1)

#post-tx
data=data.frame(fread("tdat.milo_mutR-NR_postTx.selectGenes.expr.txt",sep="\t",header=TRUE),row.names=1)
info = read.table("mutR-NR.postTx.milo-results.txt",sep="\t",header=TRUE,row.names=1)

#subset to CD8 activated T cells
info = info[info$finalIdent == "CD8 Activated" ,]
rownames(info) = paste("X",rownames(info),sep="")
data = data[,colnames(data) %in% rownames(info)]
info = info[match(colnames(data),rownames(info)),]
all(rownames(info) == colnames(data))

# Binary variable for whether the neighborhood was significantly enriched in mutR (vs NR)
info$sig = "no"
info[info$SpatialFDR < 0.1 & info$logFC > 0,]$sig = "yes"

# run wilcox test on each gene, comparing nonsig vs significant cell neighborhoods
mydf = matrix(nrow=nrow(data),ncol=3)
for (i in 1:nrow(data)){
    mygene = rownames(data)[i]
    mydata = as.numeric(data[i,])
    wt = wilcox.test(mydata~info$sig)
    pval = wt$p.value
    tt = t.test(mydata~info$sig)
    fc = tt$estimate[2] / tt$estimate[1]
    lfc = log2(fc)

    mydf[i,1] = mygene
    mydf[i,2] = lfc
    mydf[i,3] = pval
}
colnames(mydf)=c("gene","log2FC","pval")
mydf=as.data.frame(mydf)
mydf$adjp = p.adjust(mydf$pval,method="BH")
mydf = mydf[order(mydf$adjp),]

#filter out TRAV, TRBV, TRGV
library(stringr)
mydf2 = mydf[str_detect(mydf$gene,"^TRB",negate=TRUE) & str_detect(mydf$gene,"^TRA",negate=TRUE) & str_detect(mydf$gene, "^TRG",negate=TRUE) & str_detect(mydf$gene, "^IGK",negate=TRUE) & str_detect(mydf$gene, "^IGL",negate=TRUE) & str_detect(mydf$gene, "^IGH",negate=TRUE),]

#write.table(mydf2,"tdat.milo_mutR-NR_preTx.SigCD8T.DEG.wilcox.txt",sep="\t",row.names=FALSE)
write.table(mydf2,"tdat.milo_mutR-NR_postTx.SigCD8T.DEG.wilcox.txt",sep="\t",row.names=FALSE)