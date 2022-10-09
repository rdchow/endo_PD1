#!/usr/bin/Rscript

library(miloR)
library(ggbeeswarm)
library(ggplot2)
library(dplyr)
# plots beeswarms for milo results

plotDAbeeswarmNoSig <- function(da.res, group.by=NULL, alpha=0.1, subset.nhoods=NULL){
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
    mutate(logFC_color = ifelse(is_signif==1, "logFC", NA)) %>%
    arrange(group_by) %>%
    mutate(Nhood=factor(Nhood, levels=unique(Nhood))) %>%
    #ggplot(aes(group_by, logFC, color=logFC_color)) +
    ggplot(aes(group_by, logFC)) +
    scale_color_manual("gray44") + 
    #scale_color_gradient2(high="#B31A2C",mid="white",low="#2E6AB5") +
    guides(color="none") +
    xlab(group.by) + ylab("Log Fold Change") +
    geom_quasirandom(alpha=1,color="#7F7F7F") +
    coord_flip() +
    theme_bw(base_size=22) +
    theme(strip.text.y =  element_text(angle=0))

}

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

data=read.table("mutR-NR.postTx.milo-results.txt",sep="\t",header=TRUE,row.names=1)
pdf("mutR-vs-NR.post.beeswarm.milo.pdf",height=8,width=7,useDingbats=FALSE)
plotDAbeeswarm(data, group.by = "finalIdent",alpha=0.1)+theme_bw()
dev.off()

data=read.table("epiR-NR.postTx.milo-results.txt",sep="\t",header=TRUE,row.names=1)
pdf("epiR-vs-NR.post.beeswarm.milo.pdf",height=8,width=7,useDingbats=FALSE)
plotDAbeeswarmNoSig(data, group.by = "finalIdent",alpha=0.1)+theme_bw()
dev.off()

#same thing with pre-Tx files 
