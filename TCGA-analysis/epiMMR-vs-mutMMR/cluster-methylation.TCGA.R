#!/usr/bin/Rscript

data=read.table("mlh1-methylation-data.TCGA.txt",sep="\t",header=TRUE,row.names=1)

subdata = data[c(16:44),]
subdata = subdata[,colSums(is.na(subdata)) < nrow(subdata)] # remove all NA samples
mymeans = colMeans(subdata,na.rm=TRUE)

mydf=as.data.frame(sort(mymeans))
mydf$index = seq.int(nrow(mydf))
colnames(mydf)[1] = "MLH1_beta"
library(ggpubr)
library(ggplot2)

expr = read.table("mlh1-expression.TCGA.txt",sep="\t",header=TRUE)
#expr = expr[expr$sample %in% rownames(mydf),]
#expr = expr[match(rownames(mydf),expr$sample),]
rownames(expr)=expr$sample

cutoff = 0.5
mydf$expr = expr[rownames(mydf),"MLH1_expr"]
mydf$mlhbin = "no"
mydf[mydf$MLH1_beta >= cutoff & grepl('_01$',rownames(mydf)),]$mlhbin = "yes"
mydf[grepl('_11$',rownames(mydf)),]$mlhbin = "normal"
write.table(mydf,"avg-mlh1-methylscore.txt",sep="\t",row.names=TRUE)

mydf$mlhbin = factor(mydf$mlhbin,levels=c("normal","no","yes"))

pdf("UCEC-MLH.methylation.pdf",height=4,width=4,useDingbats=FALSE)
#ggscatter(mydf,x="index",y="MLH1_beta",color="MLH1_beta",facet.by="mlhbin",shape=21,fill=NA)+geom_hline(yintercept=cutoff)+scale_color_distiller(palette="Spectral")
ggboxplot(mydf,x="mlhbin",y="MLH1_beta",color="mlhbin",shape=21,fill=NA,add="jitter",palette=c("gray44","dodgerblue3","red3"))+geom_hline(yintercept=cutoff,lty=2,color="gray77")+stat_compare_means(comparisons=list(c("normal","no"),c("no","yes"),c("normal","yes")))
dev.off()

pdf("UCEC-MLH.methylation.expr.boxplot.pdf",height=4,width=4,useDingbats=FALSE)
ggboxplot(mydf,x="mlhbin",y="expr",color="mlhbin",add="jitter",palette=c("gray44","dodgerblue3","red3"))+stat_compare_means(comparisons=list(c("normal","no"),c("no","yes"),c("normal","yes")))
dev.off()

subdf = mydf[mydf$mlhbin %in% c("yes","normal"),]
wilcox.test(subdf$expr~subdf$mlhbin)$p.value