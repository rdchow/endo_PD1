#!/usr/bin/Rscript
# CCFs from Pyclone or by manual calculation (Tarabichi, Nat Methods 2021)
# scaled VAF = VAFs divided by maximum VAF detected in each sample
data = read.table("JAK1-mutation-analysis-CCFs.new.txt",sep="\t",header=TRUE)
data$myvar = paste(data$PatientID,data$Start_Position,data$Reference_Allele,data$Tumor_Seq_Allele2,sep=":")

dataf = data[,c("PatientID","timepoint","VAF","ScaleVAF","Manual_CCF","PyClone_CCF","myvar")]
# add new data for the non-mutant matched ones
dataf[nrow(dataf)+1,] = c("PEM07","pre",0,0,0,0,"PEM07:65330630:T:-")
dataf[nrow(dataf)+1,] = c("PEM10","pre",0,0,0,0,"PEM10:65306997:T:-")
dataf[nrow(dataf)+1,] = c("PEM10","pre",0,0,0,0,"PEM10:65325832:-:G")
dataf[nrow(dataf)+1,] = c("PEM12","pre",0,0,0,0,"PEM12:65325833:G:-")
dataf[nrow(dataf)+1,] = c("PEM19","post",0,0,0,0,"PEM19:65332591:-:GT")
dataf$PatientID=factor(dataf$PatientID,levels=c("PEM01","PEM07","PEM09","PEM10","PEM16","PEM11","PEM12","PEM13","PEM19","PEM20","PEM06","PEM23","PEM25"))

library(ggplot2)
library(ggpubr)
dataf$ScaleVAF = as.numeric(dataf$ScaleVAF)
dataf$Manual_CCF = as.numeric(dataf$Manual_CCF)
dataf$PyClone_CCF = as.numeric(dataf$PyClone_CCF)

#PyClone CCFs
pdf("JAK1-mutant-samples-CCFs.pyclone.dotplot.pdf",height=3,width=5.4,useDingbats=FALSE)
ggplot(dataf,aes(x=PatientID,y=PyClone_CCF,color=timepoint))+geom_point()+geom_line(aes(group=myvar))+theme_bw()+scale_color_manual(values=c("#1F8F00","#9026CD"))
dev.off()

#Manually calculated CCFs
pdf("JAK1-mutant-samples-CCFs.manual.dotplot.pdf",height=3,width=5.4,useDingbats=FALSE)
ggplot(dataf,aes(x=PatientID,y=Manual_CCF,color=timepoint))+geom_point()+geom_line(aes(group=myvar))+theme_bw()+scale_color_manual(values=c("#1F8F00","#9026CD"))
dev.off()