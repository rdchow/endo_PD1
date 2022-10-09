library(DESeq2)
library(data.table)

data = readRDS("pseudo-bulked.5rep.rds")
info = read.table("pseudo-bulked.sample.info.annot.txt",sep="\t",header=TRUE,row.names=1)
all(rownames(info) == colnames(data))

celltypes = unique(info$cellType)
groups = unique(info$Group)
groups = groups[-1] # drop healthy

#compare within each group of patients, compare pre vs post
#generate subset of data for each cell type
for (cell in celltypes){
	for (group in groups){ # for each group of patients
	
	#drop healthy data
	subinfo = info[info$cellType == cell & info$Group == group,]
	subdata = data[,colnames(data) %in% rownames(subinfo)]
	subinfo$timepoint = factor(subinfo$timepoint,levels=c("pre","post"))

	#make sure there are >=5 metacells for each timepoint, for this celltype, for this group
	samplects = as.numeric(table(subinfo$timepoint))

	if (all(samplects>=5)){
	#keep genes with >0 expression in 50%+ of metacells
	subdataf = subdata[rowSums(subdata>0) >= 0.5*ncol(subdata),]

	# first, find DEGs between pre vs post, within each group.
	dds = DESeqDataSetFromMatrix(countData = subdataf,colData = subinfo,design=~timepoint)
	dds = DESeq(dds)

	res <- results(dds,contrast = c("timepoint","post","pre"), alpha = 0.05)
	#res = res[order(res$padj),]
	res <- lfcShrink(dds,coef = "timepoint_post_vs_pre",res=res)

	resf = res[!is.na(res$padj),]
	resf = resf[order(resf$padj),]
	celln = gsub("\\s","-",cell)

	fwrite(as.data.frame(resf),paste("./pseudobulk-DEGs/preVspost/",group,"/",group,"_post-vs-pre.",celln,"_deg.txt",sep=""),sep="\t",row.names=TRUE)
	}
	}
}


