library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(Seurat)

#load data
#tdat = readRDS("filtdata-Tcell.rds")
tdat = readRDS("filtdata-all-coarseT_reann.final.rds")
info = read.table("scrna-jak1-annot.txt",sep="\t",header=TRUE,row.names=1)

tdat$JAK1 = info[tdat$orig.ident,"JAK1"]
tdat$condition = "pre"
tdat@meta.data[tdat$timepoint %in% c(3,5),]$condition="post"
tdat$clinRN = ""
tdat@meta.data[tdat$clinical == "Nonlynch-like responder",]$clinRN = "NLR"
tdat@meta.data[tdat$clinical == "Lynch-like responder",]$clinRN = "LR"
tdat@meta.data[tdat$clinical == "Nonresponder",]$clinRN = "NR"
tdat@meta.data[tdat$clinical == "Healthy",]$clinRN = "healthy"

#filter out samples with low cell numbers (< 1000), and patients that do not have matched pre vs post-tx
cellnums = table(tdat$orig.ident)
keepSamples = names(cellnums[cellnums>=1000])

#filter to pre-tx samples only
#filter to one patient group
myptgrp = "LR"
mytime = "post"

pdf(paste(myptgrp,"-JAK1mutvsWT-",mytime,"-milo-analysis.allCells.pdf",sep=""))
subdat = subset(tdat,subset = ((condition == mytime) & (orig.ident %in% keepSamples) & (clinRN %in% myptgrp)))
Idents(subdat) = "finalIdent"

tdat_sce = as.SingleCellExperiment(subdat)
umaps = Embeddings(subdat,reduction="umap")
pcas = Embeddings(subdat,reduction="pca")

tdat.milo = Milo(tdat_sce)
tdat.milo = buildGraph(tdat.milo,k=200,d=50,reduced.dim="PCA")
tdat.milo = makeNhoods(tdat.milo,prop=0.1,refined=TRUE,k=200,d=50,reduced_dims="PCA")
plotNhoodSizeHist(tdat.milo)

tdat.milo <- countCells(tdat.milo, meta.data = as.data.frame(colData(tdat.milo)), sample="orig.ident")
tdat.design <- data.frame(colData(tdat.milo))[,c("orig.ident","JAK1")]
#tdat.design$patient = as.factor(tdat.design$patient)
tdat.design <- distinct(tdat.design)
rownames(tdat.design) <- tdat.design$orig.ident

tdat.milo <- calcNhoodDistance(tdat.milo, d=50, reduced.dim = "PCA")
clin1 = "JAK1mut"
clin2 = "JAK1wt"
mycontrast = paste(clin1,"-",clin2,sep="")

da_results <- testNhoods(tdat.milo, design = ~ 0 + JAK1, design.df = tdat.design,model.contrasts = mycontrast)
head(da_results[order(da_results$PValue),])

tdat.milo <- buildNhoodGraph(tdat.milo)
umap_pl <- plotReducedDim(tdat.milo, dimred = "UMAP", colour_by="JAK1", text_by = "finalIdent", text_size = 3, point_size=0.5) +  guides(fill="none")
nh_graph_pl <- plotNhoodGraphDA(tdat.milo, da_results, layout="UMAP",alpha=0.1)

umap_pl + nh_graph_pl + plot_layout(guides="collect")
da_results <- annotateNhoods(tdat.milo, da_results, coldata_col = "finalIdent")
head(da_results[order(da_results$PValue),])

#ggplot(da_results, aes(finalIdent_fraction)) + geom_histogram(bins=50)
#ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + geom_point() + geom_hline(yintercept = 1) 

da_results$finalIdent <- ifelse(da_results$finalIdent_fraction < 0.7, "Mixed", da_results$finalIdent)
da_results$finalIdent = factor(da_results$finalIdent)
write.table(da_results[order(da_results$PValue),],paste(myptgrp,"-JAK1mutvsWT-",mytime,".milo-results.txt",sep=""),sep="\t",row.names=TRUE)

#custom plotDAbeeswarm function
plotDAbeeswarm <- function(da.res, group.by=NULL, alpha=0.1, subset.nhoods=NULL){
  if (!is.null(group.by)) {
    if (!group.by %in% colnames(da.res)) {
      stop(group.by, " is not a column in da.res. Have you forgot to run annotateNhoods(x, da.res, ", group.by,")?")
    }
    if (is.numeric(da.res[,group.by])) {
      stop(group.by, " is a numeric variable. Please bin to use for grouping.")
    }
    da.res <- mutate(da.res, group_by = da.res[,group.by])
  } else {
    da.res <- mutate(da.res, group_by = "g1")
  }

  if (!is.factor(da.res[,"group_by"])) {
    message("Converting group.by to factor...")
    da.res <- mutate(da.res, factor(group_by, levels=unique(group_by)))
    # anno_vec <- factor(anno_vec, levels=unique(anno_vec))
  }

  if (!is.null(subset.nhoods)) {
    da.res <- da.res[subset.nhoods,]
  }

  da.res %>%
    mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
    mutate(logFC_color = ifelse(is_signif==1, logFC, NA)) %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
    ggplot(aes(group_by, logFC, color=logFC_color)) +
    scale_color_gradient2(high="#B31A2C",mid="white",low="#2E6AB5") +
    guides(color="none") +
    xlab(group.by) + ylab("Log Fold Change") +
    geom_quasirandom(alpha=1) +
    coord_flip() +
    theme_bw(base_size=22) +
    theme(strip.text.y =  element_text(angle=0))

}
library(ggbeeswarm)
plotDAbeeswarm(da_results, group.by = "finalIdent",alpha=0.1)+theme_bw()

dev.off()

