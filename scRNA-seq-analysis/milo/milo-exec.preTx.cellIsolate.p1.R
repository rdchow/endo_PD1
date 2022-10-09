library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(Seurat)

#load data
tdat = readRDS("filtdata-all-coarseT_reann.final.rds")

tdat$condition = "pre"
tdat@meta.data[tdat$timepoint %in% c(3,5),]$condition="post"
tdat$clinRN = ""
tdat@meta.data[tdat$clinical == "epiR",]$clinRN = "epiR"
tdat@meta.data[tdat$clinical == "mutR",]$clinRN = "mutR"
tdat@meta.data[tdat$clinical == "NR",]$clinRN = "NR"
tdat@meta.data[tdat$clinical == "Healthy",]$clinRN = "healthy"

#filter out samples with low cell numbers (< 1000), and patients that do not have matched pre vs post-tx
cellnums = table(tdat$orig.ident)
keepSamples = names(cellnums[cellnums>=1000])

#filter to pre-tx samples only
#filter to two pt grps at a time
myptgrp1 = "mutR"
#myptgrp1 = "epiR"
myptgrp2 = "NR"
print(paste("Analyzing pre-tx data for",myptgrp1,"and",myptgrp2))

subdat = subset(tdat,subset = ((condition == "pre") & (orig.ident %in% keepSamples) & (clinRN %in% c(myptgrp1,myptgrp2))))
Idents(subdat) = "finalIdent"
DefaultAssay(subdat) = "RNA"
subdat = NormalizeData(subdat)

tdat_sce = as.SingleCellExperiment(subdat)
umaps = Embeddings(subdat,reduction="umap")
pcas = Embeddings(subdat,reduction="pca")
reducedDim(tdat_sce,type="UMAP") = umaps
reducedDim(tdat_sce,type="PCA") = pcas

tdat.milo = Milo(tdat_sce)
tdat.milo = buildGraph(tdat.milo,k=200,d=50,reduced.dim="PCA")
tdat.milo = makeNhoods(tdat.milo,prop=0.1,refined=TRUE,k=200,d=50,reduced_dims="PCA")
plotNhoodSizeHist(tdat.milo)

tdat.milo <- countCells(tdat.milo, meta.data = as.data.frame(colData(tdat.milo)), sample="orig.ident")
tdat.design <- data.frame(colData(tdat.milo))[,c("orig.ident","clinRN")]

tdat.design <- distinct(tdat.design)
rownames(tdat.design) <- tdat.design$orig.ident

tdat.milo <- calcNhoodDistance(tdat.milo, d=50, reduced.dim = "PCA")
saveRDS(tdat.milo,paste("tdat.milo_",myptgrp1,"-",myptgrp2,"_preTx.rds"))
# stopping point
# continue with part 2 in milo-exec.preTx.cellIsolatep2.R

