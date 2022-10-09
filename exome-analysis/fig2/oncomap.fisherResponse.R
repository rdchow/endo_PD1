#!/usr/bin/Rscript

library(maftools)
mydata.maf = annovarToMaf(annovar = "merged-all-pathogenic.vars.tumors.txt", refBuild = 'hg19',tsbCol = 'Sample', table = 'refGene')
write.table(mydata.maf,"merged-all-pathogenic.vars.tumors.maf",sep="\t",row.names=FALSE,quote=FALSE)

# filter to pre-treatment tumors
info = read.table("sample-info.txt",sep="\t",header=TRUE)
info$Tumor_Sample_Barcode = info$sampleID
info2 = info[info$keep == "y" & info$timepoint == "pre",]
mygene = read.table("mut-gene-list.txt",sep="\t",header=TRUE)
mygenes = mygene$gene

dataf = mydata.maf[mydata.maf$Tumor_Sample_Barcode %in% info2$sampleID,]
write.table(dataf,"merged-all-pathogenic.vars.tumors.preTx.maf",sep="\t",row.names=FALSE,quote=FALSE) #maf with variants in all mutated genes

dataf2 = dataf[dataf$Hugo_Symbol %in% mygenes,]
write.table(dataf2,"merged-top-pathogenic.vars.tumors.preTx.maf",sep="\t",row.names=FALSE,quote=FALSE) #maf with variants in top mutated genes

#  create oncoplot
mymaf = read.maf("merged-all-pathogenic.vars.tumors.preTx.maf")
info = read.table("sample-info.txt",sep="\t",header=TRUE)
info$Tumor_Sample_Barcode = info$sampleID
info2 = info[info$keep == "y" & info$timepoint == "pre",]

mygene = read.table("mut-gene-list.txt",sep="\t",header=TRUE)
mygenes = mygene$gene
library(MetBrewer)
library(RColorBrewer)
library(ggsci)
library(scales)

mycols = pal_aaas()(10)
show_col(mycols)
mycols = mycols[c(1,8,3,4,5,6)]
mycols[2] = "orangered2"
mycols[5] = "maroon3"
mycols[6]="gray44"

names(mycols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Frame_Shift_Ins',
  'Nonsense_Mutation',
  'Splice_Site',
  'Multi_Hit'
)
anncolors = c("#A3D0D6","#FED479","#F37B79")
names(anncolors) = c("Nonresponder","epiR","mutR")
anncolors2 = list(group = anncolors)

pdf("exome-oncoplot.pretxOnly.pdf",height=12,width=8)
oncoplot(mymaf,top=30,genes=mygenes,colors=mycols,annotationDat = info2,clinicalFeatures="group",annotationColor=anncolors2,showTumorSampleBarcodes=TRUE,sortByAnnotation=TRUE)
dev.off()

###########################
# Statistical analysis of mutations and ICB response
#get binary matrix of mutations
mymatrix = mutCountMatrix(mymaf)
mymatrix[mymatrix >1] = 1
mygene = read.table("mut-gene-list.txt",sep="\t",header=TRUE)
mygenes = mygene$gene

#only genes shown in oncoplot
mymatrixF = mymatrix[rownames(mymatrix) %in% mygenes,]

#calculate whether any gene mutations are enriched in responders
info = read.table("sample-info.txt",sep="\t",header=TRUE)
info$Tumor_Sample_Barcode = info$sampleID
info2 = info[info$keep == "y" & info$timepoint == "pre",]
info2$response = "R";
info2[info2$group == "Nonresponder",]$response = "NR";
info2 = info2[match(colnames(mymatrixF),info2$sampleID),]
#write.table(info2,"sample-info.reorder.txt",sep="\t",row.names=FALSE)

all(colnames(mymatrixF) == info2$sampleID)
colnames(mymatrixF) = info2$patient.ID
write.table(mymatrixF,"binary-mut-matrix.topgenes.txt",sep="\t",row.names=TRUE)

#fisher tests
mydf = matrix(nrow=nrow(mymatrixF),ncol=3)
for (i in 1:nrow(mymatrixF)){
  gene = rownames(mymatrixF)[i]
  mydat = as.numeric(mymatrixF[gene,])
  mylabs = info2$response
  ft = fisher.test(table(mylabs,mydat))
  pval = ft$p.value
  or = ft$estimate
  mydf[i,1] = gene
  mydf[i,2] = or
  mydf[i,3] = pval
}
colnames(mydf) = c("Gene","estimate","fisherPval")
write.table(mydf,"fisher-test.responseMutations.topgenes.txt",sep="\t",row.names=FALSE,quote=FALSE)


# Visualize results of Fisher tests in volcano plot
mydf = read.table("fisher-test.responseMutations.topgenes.txt",sep="\t",header=TRUE)
mydf$lOR = log(mydf$estimate)
mydf$nlp = log(mydf$fisherPval)*-1

library(ggpubr)
library(ggplot2)
pdf("fisher-test-responseMutations.scatter.pdf",height=6,width=5.5,useDingbats=FALSE)
ggscatter(mydf,x="lOR",y="nlp",fill="lOR",shape=21,size=6,label="Gene",repel=TRUE,font.label=c(3,"plain"))+scale_fill_distiller(palette="RdBu")+geom_hline(yintercept=log(0.05)*-1,lty=2,color="gray44")+geom_hline(yintercept=log(0.1)*-1,lty=2,color="gray44")
dev.off()
