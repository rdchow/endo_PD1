#!/usr/bin/Rscript
#uses CD16 NK gene signatures in epiR patients to stratify survival in TCGA UCEC cohort

library(data.table)
library(matrixStats)

info = read.table("tumor-sample-info.ann.txt",sep="\t",header=TRUE,row.names=1)
data=data.frame(fread("UCEC.merged.tpm.tumor.txt",sep="\t"),row.names=1) #obtain from GDAC; not uploaded to code repo
data = data[,colnames(data) %in% rownames(info)]
data = data[,match(rownames(info),colnames(data))]
all(rownames(info) == colnames(data))

#subset to CD16 NK signature genes
genes = read.table("CD16-NK-signature.epiR.txt",sep="\t",header=TRUE)
htseq = read.table("htseq-gene-names-translated.txt",sep="\t",header=TRUE)
wantgenes = c(intersect(genes$Pre_up_shared,genes$Post_up_shared),intersect(genes$Pre_down_shared,genes$Post_down_shared))

wantgenes = sort(wantgenes[wantgenes != ""])
htseq = htseq[htseq$Name %in% wantgenes,]
dataf = data[rownames(data) %in% htseq$ENSG,]
htseqF = htseq[match(rownames(dataf),htseq$ENSG),]
all(htseqF$ensembl == rownames(dataf))
rownames(dataf) = make.unique(htseqF$Name)
dataf = dataf[rowSums(dataf) >0,]

#######################################
# survival analysis
library(survival)
library(survminer)
library(glmnet)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(data.table)

# run lasso regression
hdata = as.data.frame(t(dataf))
sdata = as.data.frame(scale(hdata))
mysurv = Surv(info$OS_MONTHS,info$OS_STATUS)
x <- model.matrix(~ ., data=sdata)
set.seed(35)
fit <- glmnet(x, mysurv, family="cox")
cv.fit <- cv.glmnet(x, mysurv, family="cox",nfolds=10)
plot(cv.fit)

myres=data.frame(as.matrix(coef(cv.fit, s = "lambda.min")))
myres$coefs = myres$X1

#filter out genes where the directionality of association is inverted
#basically, we keep genes with negative coefs if they are UP in NK
# and keep genes with positive coefs if they are DOWN in NK
upgenes = unique(c(genes$Pre_up_shared,genes$Post_up_shared))
downgenes = unique(c(genes$Pre_down_shared,genes$Post_down_shared))

myres2 = myres[myres$coefs != 0,]
myres2$direction = "no"
myres2$gene = rownames(myres2)
myres2[myres2$coefs >0 & myres2$gene %in% downgenes,]$direction = "yes"
myres2[myres2$coefs <0 & myres2$gene %in% upgenes,]$direction = "yes"

myres3 = myres2[myres2$direction == "yes",]
myvars = rownames(myres3)

subd = as.data.frame(sdata[,colnames(sdata) %in% myvars])
m1 <- coxph(mysurv ~ ., data=subd)
pdf("survival-forestplot.NK16degs.pdf",height=4,width=6,useDingbats=FALSE)
ggforest(m1)
dev.off()
summary(m1)

#################################################
# generate a "meta-gene" signature score
# CD63, PPIB, CEBPB, LDOC1
myvars=c("CD63","PPIB","CEBPB","LDOC1")
hdata = as.data.frame(t(dataf))
sdata = as.data.frame(scale(hdata))
subd = as.data.frame(sdata[,colnames(sdata) %in% myvars])
mysurv = Surv(info$OS_MONTHS,info$OS_STATUS)

m1 <- coxph(mysurv ~ ., data=subd)
summary(m1)

#PCA analysis to generate a composite score of all 4 genes
loadings = prcomp(subd,scale.=FALSE)$rotation
loadings = as.data.frame(loadings)
loadings = loadings[order(loadings$PC1,decreasing=TRUE),]
loadings$gene = rownames(loadings)
write.table(loadings,"epiR-nk4-PC1-loadings.txt",sep="\t",row.names=TRUE)

pdf("epiR-nk4-PC1-loadings.pdf",height=6,width=4)
ggbarplot(loadings,x="gene",y="PC1",fill="PC1")+geom_hline(yintercept=0,lty=2,color="gray44")+scale_fill_distiller(palette="RdBu")
dev.off()

pcres = as.data.frame(prcomp(subd,scale.=FALSE)$x)
info$PC1 = pcres$PC1
info$PC2 = pcres$PC2
info$PC3 = pcres$PC3

#Binarize the principal components
library(gtools)
info$PC1bin = quantcut(info$PC1, q = 2)
info$PC2bin = quantcut(info$PC2, q = 2)
info$PC3bin = quantcut(info$PC3, q = 2)
mysurv = Surv(info$OS_MONTHS,info$OS_STATUS)

# CoxPH survival, grouped by PC1 (epiR-NK4 score)
m1 <- coxph(mysurv ~ PC1bin,data=info)
ggforest(m1)
summary(m1)

m2 <- survfit(mysurv ~ PC1bin,info)
msd = survdiff(mysurv ~ PC1bin,info)
mypval = broom::glance(msd)$p.value

#CoxPH survival, considering PCs1-3
m1 <- coxph(mysurv ~ PC1bin+PC2bin+PC3bin,data=info)
ggforest(m1)
summary(m1)

#PCA Biplot for visualization
library(ggfortify)
pdf("PCA-biplot.pdf",height=5.5,width=7.5,useDingbats=FALSE)
autoplot(prcomp(subd,scale.=FALSE),loadings=TRUE,loadings.label=TRUE,colour="PC1bin",data=info)+theme_bw()
dev.off()


###########################
#stratify TCGA UCEC cohort by MSI status
all(rownames(info) == rownames(subd))

highdata = subd[info$MSI == "yes",] #MSI-H
highinfo = info[info$MSI == "yes",] 
lowdata = subd[info$MSI == "no",] #non MSI
lowinfo = info[info$MSI == "no",] 

# MSI-H samples
highfit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ PC1bin,highinfo)
highmsd = survdiff(Surv(OS_MONTHS,OS_STATUS) ~ PC1bin,highinfo)
highpval = broom::glance(highmsd)$p.value

pdf("4gene-CD16NK.marker.survival.MSIhi.pdf",height=4,width=5,useDingbats=FALSE)
ggsurvplot(highfit,pval=highpval)
dev.off()

###########################
#stratify TCGA UCEC cohort by NK score
all(rownames(info) == rownames(subd))
info$NKbin = quantcut(info$NKscore, q = 2)

highdata = subd[info$NKbin == "(0.0591,0.247]",] #nk high
highinfo = info[info$NKbin == "(0.0591,0.247]",] 
lowdata = subd[info$NKbin == "[0,0.0591]",] #nk low
lowinfo = info[info$NKbin == "[0,0.0591]",] 

# NK high first
highfit <- survfit(Surv(OS_MONTHS,OS_STATUS) ~ PC1bin,highinfo)
highmsd = survdiff(Surv(OS_MONTHS,OS_STATUS) ~ PC1bin,highinfo)
highpval = broom::glance(highmsd)$p.value

pdf("4gene-CD16NK.marker.survival.NKhi.pdf",height=4,width=5)
ggsurvplot(highfit,pval=highpval)
dev.off()