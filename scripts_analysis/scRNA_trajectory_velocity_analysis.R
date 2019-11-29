##########################################################################
##########################################################################
# Project: Aleks single cell RNA-seq analysis
# Script purpose: predict the connection and directionality for clusters from seurat
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Nov 28 13:16:06 2019
##########################################################################
##########################################################################

#Function for soursing functions
source.my.script <- function(name.of.function){
  tryCatch(path <- rstudioapi::getSourceEditorContext()$path,
           error = function(e){
             install.packages("rstudioapi")
             path <-  rstudioapi::getSourceEditorContext()$path})
  source.path <- sub(basename(path), "", path)
  source(paste0(source.path,name.of.function))
}

## set up the paths for the data and results
tryCatch(path <- rstudioapi::getSourceEditorContext()$path, 
         error = function(e){
           install.packages("rstudioapi")
           path <-  rstudioapi::getSourceEditorContext()$path})
source.path <- sub(basename(path), "", path)


user <- "results_jiwang/"
setwd(paste0("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/",user))

version.DATA = 'all_batches'
version.analysis =  paste0(version.DATA, '_20191115')
dataDir = paste0("../data/gene_counts/")
resDir = paste0("results/", version.analysis)
tabDir = paste0("results/", version.analysis, "/tables/")

RdataDir = paste0("results/", version.analysis, "/Rdata/")
if(!dir.exists("results/")){dir.create("results/")}
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}


########################################################
########################################################
# Section : Quick clustering 
# 
########################################################
########################################################
load(file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_timingEst_Normed_bc_Seurat.Rdata'))

library(Seurat)
library(SeuratWrappers)
library(cowplot)
library(ggplot2)

ElbowPlot(ms)
ms <- FindNeighbors(object = ms, reduction = 'pca', dims = 1:20)
ms <- FindClusters(object = ms, resolution = 6, algorithm = 4)
#DimPlot(ms, reduction = "umap") + ggtitle('Leiden')
ms1 = FindClusters(object = ms, resolution = 6, algorithm = 1)
#ms2 = FindClusters(object = ms, resolution = 6, algorithm = 2)
#ms <- RunUMAP(object = ms, dims = 1:20)
p0 = DimPlot(ms, reduction = "umap") + ggtitle('Leiden')
p1 = DimPlot(ms1, reduction = "umap") + ggtitle('louvain')
#p2 = DimPlot(ms2, reduction = "umap") + ggtitle('slm')
plot_grid(p0, p1, ncol = 2)


library(loomR)

## save seurat as loom file
## there are some errors trigged here (https://github.com/mojaveazure/loomR/issues/40)
## when I have NA values in the Seurat object metadata. 
## Weirdly it seems to happen when I have NAs present in a column of type character and not of type numeric.
## so here we filter many columns of metadata to avoid the error
md = ms@meta.data
md = md[, c(1:6, 70:74, 121:130)]
ms@meta.data = md
ms@graphs <- list()

ms.loom = as.loom(ms, assay = 'SCT', filename = paste0(resDir, "/Seurat_tmp.loom"), overwrite = TRUE, verbose = FALSE)

library("tidyverse")
library("reticulate")
#Sys.which("python")

# condaenv is not working
#use_condaenv(condaenv = 'dca', conda = "auto", required = FALSE)
use_python("/Users/jiwang/anaconda3/envs/dca/python")
#Sys.which("python")

sc <- import("scanpy", convert = FALSE)
loomfile = paste0(resDir, "/Seurat_tmp.loom")
adata = sc$read_loom(filename = loomfile, sparse=True, cleanup=False, X_name='SCT', obs_names='CellID', var_names='Gene')


