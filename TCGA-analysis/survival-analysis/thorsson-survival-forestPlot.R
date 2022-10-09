#!/usr/bin/Rscript
# find tumor/immune features associated with survival
library(survival)
library(survminer)
library(glmnet)
library(sjPlot)
library(sjlabelled)
library(sjmisc)

library(data.table)
data = read.table("ucec-thorsson-immune.simple.MSIsensor.txt",sep="\t",header=TRUE,row.names=1) # data from Thorsson et al, with MSIsensor scores from cBioPortal
data = data[complete.cases(data),]
data$MSI = 0
data[data$TCGA.Subtype %in% c("UCEC.MSI","UCEC.POLE"),]$MSI = 1 #MSI categorization

data$OS.Time = data$OS.Time/30.417 #convert days to months
survdat = data[,c("OS.Time","OS")]
mysurv = Surv(survdat$OS.Time,survdat$OS)

dropcols = c("PFI","PFI.Time","OS.Time","OS")

dataf = data[,!colnames(data) %in% dropcols]
msi = dataf$MSI
dataf$MSI = NULL
dataf$TCGA.Subtype = NULL

dataf = as.data.frame(scale(dataf))
dataf$MSI = msi

m1 <- coxph(mysurv ~ Homologous.Recombination.Defects+TGF.beta.Response + indel_num + IFN.gamma.Response + MSIsensor + T.Cells.CD8 + NK.Cells.Activated + Nonsilent.Mutation.Rate,dataf)
ggforest(m1)
summary(m1)

#forest plot
theme_set(theme_sjplot())
pdf("survival-forestplot.pdf",height=4,width=6,useDingbats=FALSE)
ggforest(m1)
dev.off()
