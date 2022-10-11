library(dplyr)
library(purrr)
library(Matrix)
library(patchwork)
library(Seurat)
library(future)
options(future.globals.maxSize= 5*1024^4, future.rng.onMisuse="ignore")
setwd('~/scratch60/pem/')

ldat <- list(
  readRDS('main.rds'),
  readRDS('eric.rds'),
  readRDS('hd.rds')
)

ggplot(dat@meta.data)+geom_density(aes(nCount_RNA))+scale_x_log10()+geom_vline(xintercept = 400)
ggplot(dat@meta.data)+geom_density(aes(nFeature_RNA))+scale_x_log10()+geom_vline(xintercept = 150)
ggplot(hdat@meta.data)+geom_density(aes(propMt))+geom_vline(xintercept = 0.2)

c916 <- readRDS('barcodes_9_16.rds') # preidentified high mt-reads cells that disrupted clustering
ldat <- map(1:3,function(i){
  dat <- ldat[[i]]
  dat$propMt <- PercentageFeatureSet(dat,'^MT-')/100
  dat$temp2 <- if_else(rownames(dat@meta.data) %in% c916[[i]],'n','y')
  dat <- subset(dat, subset=((nCount_RNA>400) & (nFeature_RNA>150) & (propMt<0.2) & (temp2=='y')))
  dat <- SCTransform(dat,vars.to.regress=c('nCount_RNA','nFeature_RNA','propMt'))
  dat
})
saveRDS(ldat,'ldat1.rds')
