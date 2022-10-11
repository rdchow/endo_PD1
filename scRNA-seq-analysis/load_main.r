library(dplyr)
library(purrr)
library(Matrix)
library(Seurat)
library(future)

options(future.globals.maxSize= 5*1024^4, future.rng.onMisuse="ignore") 

# registering URL/Wherever they are
setwd('~/scratch60/pem/data/raw_counts/santin')
file_list = list.dirs(full.names=FALSE, recursive = FALSE)

tmpList <- map(file_list, function(tmp2){
  tmp <- Read10X(data.dir = tmp2)
  tmp <- CreateSeuratObject(count=tmp, project=tmp2)
  tmp[['file_name']] <- tmp2
  return(tmp)
})

tmpList <- merge(tmpList[[1]], y=tmpList[2:length(tmpList)] )
saveRDS(tmpList,'~/scratch60/pem/main.rds')
