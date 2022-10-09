# R4.0.3
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(Seurat)

tdat.milo = readRDS("tdat.milo_mutR-NR_preTx.rds") # load the .rds created by part 1 script
myptgrp1 = "mutR"
myptgrp2 = "NR"
clin1 = paste("clinRN",myptgrp1,sep="")
clin2 = paste("clinRN",myptgrp2,sep="")
mycontrast = paste(clin1,"-",clin2,sep="")

tdat.design <- data.frame(colData(tdat.milo))[,c("orig.ident","clinRN")]
tdat.design <- distinct(tdat.design)
rownames(tdat.design) <- tdat.design$orig.ident

#export the nhoodexpression object
# used later for plotting violins of gene expression in neighborhoods
library(scran)
dec = modelGeneVar(tdat.milo)
hvgs = getTopHVGs(dec, n=5000)
hvgs = c(hvgs,"IL12A","IL10")
head(hvgs)
testgenes = c("IFNG","TNF","IL10","IL13","IL12A","IL2","IL4","CXCL13")

tdat.milo = calcNhoodExpression(tdat.milo,subset.row=hvgs)
bbbt = nhoodExpression(tdat.milo)
colnames(bbbt)=seq(1:ncol(bbbt))
write.table(bbbt,"tdat.milo_mutR-NR_preTx.selectGenes.expr.txt",sep="\t",row.names=TRUE)

# continue running Milo
da_results <- testNhoods(tdat.milo, design = ~ 0 + clinRN, design.df = tdat.design,model.contrasts = mycontrast)
head(da_results[order(da_results$PValue),])

tdat.milo <- buildNhoodGraph(tdat.milo)
#umap_pl <- plotReducedDim(tdat.milo, dimred = "UMAP", colour_by="clinRN", text_by = "finalIdent", text_size = 3, point_size=0.5) +  guides(fill="none")
#nh_graph_pl <- plotNhoodGraphDA(tdat.milo, da_results, layout="UMAP",alpha=0.1)
#umap_pl + nh_graph_pl + plot_layout(guides="collect")

da_results <- annotateNhoods(tdat.milo, da_results, coldata_col = "finalIdent")
head(da_results[order(da_results$PValue),])

da_results$finalIdent <- ifelse(da_results$finalIdent_fraction < 0.7, "Mixed", da_results$finalIdent)
da_results$finalIdent = factor(da_results$finalIdent)
write.table(da_results[order(da_results$PValue),],paste(myptgrp1,"-",myptgrp2,".preTx.milo-results.txt",sep=""),sep="\t",row.names=TRUE)

# Will plot these results with beeswarm-plots.R
