library(DESeq2)
library(data.table)

data = readRDS("pseudo-bulked.5rep.rds")
info = read.table("pseudo-bulked.sample.info.annot.txt",sep="\t",header=TRUE,row.names=1)
all(rownames(info) == colnames(data))

celltypes = unique(info$cellType)
timepoints = unique(info$timepoint)
timepoints = timepoints[-1] # drop healthy

#compare within pre-tx or post-tx samples, comparing between groups
#generate subset of data for each cell type
for (cell in celltypes){
	for (time in timepoints){
	#subset to one of the timepoints 
	subinfo = info[info$cellType == cell & info$timepoint == time,]
	subdata = data[,colnames(data) %in% rownames(subinfo)]
	subinfo$Group = factor(subinfo$Group,levels=c("mutR","NR","epiR"))
	
	#check that each sample group has at least 5 samples for this cell type
	samplects = as.numeric(table(subinfo$Group))

	if (all(samplects >= 5)){ #only analyze cell types with >=5 pseudocells

	#keep genes with >0 expression in 50%+ of metacells
	subdataf = subdata[rowSums(subdata>0) >= 0.5*ncol(subdata),]
	all(rownames(subinfo)==colnames(subdataf))
	
	dds = DESeqDataSetFromMatrix(countData = subdataf,colData = subinfo,design=~Group)
	dds = DESeq(dds)
	
	# mutR vs NR ###
	#res <- results(dds,contrast = c("Group","mutR","NR"), alpha = 0.05)
	#res = res[order(res$padj),]
	#res <- lfcShrink(dds,coef = "Group_mutR_vs_NR",res=res)

	#resf = res[!is.na(res$padj),]
	#resf = resf[order(resf$padj),]
	#celln = gsub("\\s","-",cell)
	#fwrite(as.data.frame(resf),paste("./pseudobulk-DEGs/mutRvsNR/",time,"/",time,"_mutR-vs-NR.",celln,"_deg.txt",sep=""),sep="\t",row.names=TRUE)


	# epiR vs NR ###

	#res <- results(dds,contrast = c("Group","epiR","NR"), alpha = 0.05)
        #res = res[order(res$padj),]
        #res <- lfcShrink(dds,coef = "Group_epiR_vs_NR",res=res)

        #resf = res[!is.na(res$padj),]
        #resf = resf[order(resf$padj),]
        #celln = gsub("\\s","-",cell)
        #fwrite(as.data.frame(resf),paste("./pseudobulk-DEGs/epiRvsNR/",time,"/",time,"_epiR-vs-NR.",celln,"_deg.txt",sep=""),sep="\t",row.names=TRUE)	
	
	# epiR vs mutR ###
        res <- results(dds,contrast = c("Group","epiR","mutR"), alpha = 0.05)
        #res = res[order(res$padj),]
        res <- lfcShrink(dds,coef = "Group_epiR_vs_mutR",res=res)

        resf = res[!is.na(res$padj),]
        resf = resf[order(resf$padj),]
        celln = gsub("\\s","-",cell)
        fwrite(as.data.frame(resf),paste("./pseudobulk-DEGs/epiRvsmutR/",time,"/",time,"_epiR-vs-mutR.",celln,"_deg.txt",sep=""),sep="\t",row.names=TRUE)
	
	}
	}
}	


