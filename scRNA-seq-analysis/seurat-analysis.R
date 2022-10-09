#!/usr/bin/Rscript
library(Seurat)

# Processed RDS object is available on GEO:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212217

# ============================
# plot UMAPs

# All cells
data = readRDS("filtdata-all-coarseT_reann.final.rds") #from GEO
# Full Metadata file:
metadata = read.table("filtdata-all-coarseT_reann.final.metadata.txt",sep="\t",header=TRUE,row.names=1)
data@meta.data=metadata

Idents(data) = "cluster_fine"
options(bitmapType='cairo') 	
png('all-umap.png',width=9,height=6,units='in',res=600,bg='transparent')
DimPlot(data,group.by="cluster_fine",cols=c('#845BA6','#FBC74C','#DFDE4B','#EF4875','#F389A7','#F48268','#8EA33B','#D0A836','#524CA8','#A380BA','#AC1F64','#452B59','#4F7EB7'),raster=F)
dev.off()

# T cells
tdat = subset(data,subset=cluster_fine == "T cells")
Idents(tdat) = "finalIdent"
options(bitmapType='cairo')
png('tcells-umap.png',width=9,height=6,units='in',res=600,bg='transparent')
DimPlot(tdat,group.by="finalIdent",cols = c('#395DA7','#72A5D9','#359E48','#125C2F','#91CD88','#A01F6A','#6B8739','#7F7F7F','#9172B3'),raster=F)
dev.off()

# ============================
# plot Violins with key marker genes
library(Seurat)
data = readRDS("filtdata-all-coarseT_reann.final.rds")
# Full Metadata file:
metadata = read.table("filtdata-all-coarseT_reann.final.metadata.txt",sep="\t",header=TRUE,row.names=1)
data@meta.data=metadata

# All cells
data$cluster_fine=factor(data$cluster_fine,levels=rev(c("T Cells","Naive B Cells","Activated B Cells","Memory B Cells","Plasma Cells","CD56 NK","CD16 NK","CD16 Monocytes","Intermediate Monocytes","CD14 Monocytes","cDC","pDC","Hematopoetic Progenitors")))
Idents(data) = "cluster_fine"
DefaultAssay(data) = "RNA"
mygenes = c("CD3D","CD19","MS4A1","CCR7","IGHD","IGKC","IGHG1","MZB1","GNLY","NCAM1","KLRF1","PRF1","FCGR3A","S100A9","CD14","FCER1A","LILRA4","KIT","CD34")
pdf("all-cells-violin.pdf",height=30,width=35,useDingbats=FALSE)
VlnPlot(data,features=mygenes,pt.size=0,log=TRUE,stack=TRUE)
dev.off()

# T cells
tdat = subset(data,subset=cluster_fine == "T cells")

tdat$finalIdent=factor(tdat$finalIdent,levels=rev(c("CD8 Naive","CD8 Activated","CD8 Exhausted","CD4 Naive","CD4 Activated","Treg","Cycling","GDT","MAIT")))
Idents(tdat)="finalIdent"
DefaultAssay(tdat) = "RNA"
mygenes = c("CD3D","CD8A","CD4","CCR7","TCF7","GZMA","TOX","HAVCR2","CTLA4","MKI67","TRDC","TRAV1-2")
pdf("tCoarse-cells-violin.pdf",height=30,width=35,useDingbats=FALSE)
VlnPlot(tdat,features=mygenes,pt.size=0,log=TRUE,stack=TRUE)
dev.off()