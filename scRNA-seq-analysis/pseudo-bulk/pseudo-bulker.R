library(Seurat)
#generate pseudobulk profiles for each cell type in each sample (get 5 pseudobulk profiles per cell type, per sample)

# scRNA-seq data is on GEO:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212217

tdat = readRDS("filtdata-all-coarseT_reann.final.rds")
DefaultAssay(tdat) = "RNA"
Idents(tdat) = "finalIdent"

#filter out samples with low cell numbers (< 1000)
cellnums = table(tdat$orig.ident)
keepSamples = names(cellnums[cellnums>=1000])
tdat = subset(tdat,subset = orig.ident %in% keepSamples)

set.seed(1)
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

mymeta = tdat@meta.data
myindexes = data.frame()

#create new variable with patient:cell:randomSample# 
samples = unique(tdat$orig.ident)
for (sample in samples){
	submeta = mymeta[mymeta$orig.ident == sample,]
	celltypes = unique(submeta$finalIdent) # list of cell types in this subset

	for (cell in celltypes){
		mycells = rownames(submeta[submeta$finalIdent == cell,])  # cell IDs with the given cell type
		
		if (length(mycells) >= 10){ # only analyze if at least 10 cells of this type in this sample 
			s = sample(mycells)
			rs = chunk2(s,5) # list of 5 random samples for this cell type, in this sample
			mydf = as.data.frame(c(rs[[1]],rs[[2]],rs[[3]],rs[[4]],rs[[5]]))
			mydf$sample = c(rep(1,length(rs[[1]])), rep(2,length(rs[[2]])), rep(3,length(rs[[3]])), rep(4,length(rs[[4]])), rep(5,length(rs[[5]])))
			mydf$cell = cell
			colnames(mydf)[1] = "cellID"
			mydf$cellSample = paste(sample,mydf$cell,mydf$sample,sep="_")	
			myindexes = rbind(myindexes,mydf)
		}
	}
}

#aggregate counts based on new grouping variable
library(Matrix.utils)
dataf = subset(tdat,cells = myindexes$cellID) # since we dropped cell types in a sample if < 10 cells
myindexes2 = myindexes[match(colnames(dataf),myindexes$cellID),] # reorder the grouping dataframe to match the matrix

origcts = GetAssayData(dataf,assay="RNA",slot="data")
all(colnames(origcts) == myindexes2$cellID)

aggmat <- aggregate.Matrix(t(origcts),groupings = myindexes2$cellSample, fun = "sum")
aggmatT = as.data.frame(t(aggmat))

saveRDS(aggmatT,"pseudo-bulked.5rep.rds")
write.table(colnames(aggmatT),"pseudo-bulked.sample.info.txt",sep="\t")

	



		

