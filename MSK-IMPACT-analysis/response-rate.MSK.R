#!/usr/bin/Rscript

library(ggpubr)
data = read.table("valero-nc-data.ann.txt",sep="\t",header=TRUE)
data = data[data$Cancer.type == "Endometrial",]
as.data.frame.matrix(table(data$Response.to.ICI,data$MMR_canon))
as.data.frame.matrix(table(data$Response.to.ICI,data$ARID1A_mut))
table(data$Response.to.ICI,data$ARID1A_mut,data$MSI.type)

ggboxplot(data,x="MMR_canon",y="TMB..Mutations.Mb.",fill="MMR_canon")+scale_y_log10()

data$category = paste(data$MMR_canon,data$Response.to.ICI,sep="_")
data$category = factor(data$category,levels=c("no_No","no_Yes","other-MMRd_No","other-MMRd_Yes","mut-MMRd_No","mut-MMRd_Yes"))
ggboxplot(data,x="category",y="TMB..Mutations.Mb.",add="jitter")+scale_y_log10()

data$logTMB = log10(data$TMB..Mutations.Mb.+1)
pdf("barplot-TMB-by-responseMSIcat.pdf",height=4,width=4.5,useDingbats=FALSE)
ggbarplot(data,x="category",y="logTMB",add=c("mean_se","jitter"),fill="category")+stat_compare_means(comparisons = list(c("no_No","no_Yes"),c("other-MMRd_No","other-MMRd_Yes"),c("mut-MMRd_No","mut-MMRd_Yes")))
dev.off()

ggbarplot(data,x="MMR_canon",y="logTMB",add=c("mean_se","jitter"),fill="MMR_canon")+stat_compare_means(comparisons=list(c("no","other-MMRd"),c("other-MMRd","mut-MMRd"),c("no","mut-MMRd")))


m1 = lm(logTMB ~ MMR_canon+Response.to.ICI,data=data)
data$responseBool = as.numeric(data$Response.to.ICI == "Yes")
m1 = glm(responseBool ~ MMR_canon+logTMB,data=data,family=binomial(link='logit'))


# all cancer types, ARID1A
library(ggpubr)
data = read.table("valero-nc-data.ann.txt",sep="\t",header=TRUE)
data = data[data$MSI.type == "Unstable" & data$Cancer.type == "Endometrial",]
as.data.frame.matrix(table(data$Response.to.ICI,data$ARID1A_mut))

data = read.table("valero-nc-data.ann.txt",sep="\t",header=TRUE)
data = data[data$MSI.type != "Unstable" & data$Cancer.type == "Endometrial",]
as.data.frame.matrix(table(data$Response.to.ICI,data$ARID1A_mut))