#!/usr/bin/Rscript

library(ggplot2)
library(ggpubr)

nrpre = read.table("NR-JAK1mutvsWT-pre.milo-results.txt",sep="\t",header=TRUE)
nrpost = read.table("NR-JAK1mutvsWT-post.milo-results.txt",sep="\t",header=TRUE)

epipre = read.table("epiR-JAK1mutvsWT-pre.milo-results.txt",sep="\t",header=TRUE)
epipost = read.table("epiR-JAK1mutvsWT-post.milo-results.txt",sep="\t",header=TRUE)

mutpre = read.table("mutR-JAK1mutvsWT-pre.milo-results.txt",sep="\t",header=TRUE)
mutpost = read.table("mutR-JAK1mutvsWT-post.milo-results.txt",sep="\t",header=TRUE)

# add nlp column
nrpre$nlp = log10(nrpre$FDR)*-1
nrpre$cat = "nr_pre"
nrpost$nlp = log10(nrpost$FDR)*-1
nrpost$cat = "nr_post"

epipre$nlp = log10(epipre$FDR)*-1
epipre$cat = "epi_pre"
epipost$nlp = log10(epipost$FDR)*-1
epipost$cat = "epi_post"

mutpre$nlp = log10(mutpre$FDR)*-1
mutpre$cat = "mut_pre"
mutpost$nlp = log10(mutpost$FDR)*-1
mutpost$cat = "mut_post"

mydf=as.data.frame(rbind(nrpre,nrpost,epipre,epipost,mutpre,mutpost))

cutoff = log10(0.1)*-1

mydf$cat = factor(mydf$cat,levels=c("nr_pre","epi_pre","mut_pre","nr_post","epi_post","mut_post"))

pdf("JAK1-miloR-volcano.pdf",height=7,width=8.5,useDingbats=FALSE)
ggscatter(data=mydf,x="logFC",y="nlp",facet.by="cat",scales="free_x",shape=21,fill=NA)+geom_hline(yintercept=1,lty=2)
dev.off()

