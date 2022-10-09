#!/usr/bin/Rscript
data=read.table("tumor-sample-info.ann.mmr.TCGA.txt",sep="\t",header=TRUE,row.names=1)
library(ggpubr)

data$condition = "non-MSI"
data[data$MSI == "yes",]$condition = "other-MSI"
data[data$MMR_mut=="yes"& data$MSI == "yes",]$condition = "mut-MMRd"
data[data$MLH_meth == "yes" & data$MSI == "yes",]$condition = "epi-MMRd"
# note that methylation status is prioritized over mutation status

table(data$MSI,data$condition)

data$mutRate = data$NonsilentMutRate+1
#data$pMHC = data$pMHC+1
data$condition = factor(data$condition,levels=c("non-MSI","other-MSI","epi-MMRd","mut-MMRd"))
write.table(data,"tumor-sample-info.ann.mmr.4grp.txt",sep="\t",row.names=TRUE)

mycols = c("#787777","#dbcf5e","#C8760E","#B6312F")

g1=ggboxplot(data,y="mutRate",x="condition",group="condition",fill="condition",add="jitter",palette=mycols)+stat_compare_means(comparisons=list(c("non-MSI","epi-MMRd"),c("epi-MMRd","mut-MMRd"),c("non-MSI","mut-MMRd")))+scale_y_log10()

g2=ggboxplot(data,y="pMHC",x="condition",group="condition",fill="condition",add="jitter",palette=mycols)+stat_compare_means(comparisons=list(c("non-MSI","epi-MMRd"),c("epi-MMRd","mut-MMRd"),c("non-MSI","mut-MMRd")))+scale_y_log10()

pdf("TCGA-MMR-types.boxplots.pdf",height=5,width=10,useDingbats=FALSE)
ggarrange(g1,g2,ncol=2)
dev.off()

# MutRate
subinfo = data[data$condition %in% c("non-MSI","epi-MMRd"),]
wilcox.test(subinfo$mutRate~subinfo$condition)$p.value
subinfo = data[data$condition %in% c("non-MSI","mut-MMRd"),]
wilcox.test(subinfo$mutRate~subinfo$condition)$p.value
subinfo = data[data$condition %in% c("mut-MMRd","epi-MMRd"),]
wilcox.test(subinfo$mutRate~subinfo$condition)$p.value

# pMHC
subinfo = data[data$condition %in% c("non-MSI","epi-MMRd"),]
wilcox.test(subinfo$pMHC~subinfo$condition)$p.value
subinfo = data[data$condition %in% c("non-MSI","mut-MMRd"),]
wilcox.test(subinfo$pMHC~subinfo$condition)$p.value


# determine stage association with categories
stageTab= as.data.frame.matrix(table(data$Stage,data$condition))
write.table(stageTab,"UCEC-TCGA-stage-by-MMRcategory.txt",sep="\t",row.names=TRUE)