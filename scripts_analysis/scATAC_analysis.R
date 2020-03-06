##########################################################################
##########################################################################
# Project: Aleks' single-cell MS project 
# Script purpose: analyze scATAC-seq with cisTopic
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Mar  2 14:43:04 2020
##########################################################################
##########################################################################
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

version.DATA = 'scATAC_earlyEmbryo'
version.analysis =  paste0(version.DATA, '_20200302')
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
# Section : Input data 
# 
########################################################
########################################################
library(data.table)
library(Matrix)
library(tictoc)
DownSample.mtx = FALSE

#loadRDS(filter.out, file = paste0(output_dir, '/EmptyDrop_obj.rds'))
input_mtx_dir = '/Volumes/groups/cochella/jiwang/Projects/Aleks/R8898_scATAC/output/filtered_matrix'
#metrics <- paste0(pathTo10X, 'atac_v1_pbmc_5k_singlecell.csv')

mat = readMM(paste0(input_mtx_dir, "/matrix.mtx"))
features = fread(paste0(input_mtx_dir, '/features.txt'), header = F)
barcodes = fread(paste0(input_mtx_dir, '/barcodes.txt'), header = F)
rownames(mat) = features$V1
colnames(mat) = barcodes$V1

if(DownSample.mtx){
  ## downsample
  set.seed(1)
  input = mat[, sample(c(1:ncol(mat)), 2000, replace = FALSE)]
}else{
  input = mat
}

peaknames = rownames(input)
peaknames = sapply(peaknames, function(x) {xx = unlist(strsplit(x, '-')); return(paste0(xx[1], ':', xx[2], "-", xx[3]))} )
rownames(input) = peaknames

##########################################
# cisTopic
##########################################
library(cisTopic)

cisTopicObject <- createcisTopicObject(count.matrix = input,  project.name='earlyEmbro')

tic('run model for cisTopic')
#toc()
nb.topcis = c(seq(10, 100, 10))
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=nb.topcis, seed=987, nCores=6, iterations = 200, addModels=FALSE)
toc()

save(cisTopicObject, file = paste0(RdataDir, 'cisTopicOject_runWarpLDAmodels_localrunning.Rdata'))

#logLikelihoodByIter(cisTopicObject, select=nb.topcis)

##########################################
# model selection and interpretation
##########################################
load(file = 'results/scATAC_earlyEmbryo_20200302/Rdata/cisTopicOject_runWarpLDAmodels.Rdata')

par(mfrow=c(3,3))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
cisTopicObject <- selectModel(cisTopicObject, type='derivative')

cisTopicObject <- runPCA(cisTopicObject, target='cell', seed=123, method='Probability')
#cisTopicObject <- runtSNE(cisTopicObject, target='cell', seed=123, pca = TRUE, method='Probability')
cisTopicObject <- runDM(cisTopicObject, target='cell', seed=123, pca=TRUE, method='Probability')

nb.pcs = 20; n.neighbors = 20; min.dist = 0.2;
#ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, pca = TRUE, method='Probability', n.neighbors = n.neighbors, min.dist = min.dist)


cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
dim(cellassign)
cellassign[1:5,1:5]

#sum(colnames(cellassign) == rownames(metadata))
### make quick clustering using tSNE and densityClust
set.seed(123)
library(Rtsne)
DR <- Rtsne(t(cellassign), pca=F)
DRdist <- dist(DR$Y)
library(densityClust)
dclust <- densityClust(DRdist,gaussian=T)

cutoff.rho = 100
cutoff.delta = 5
dclust <- findClusters(dclust, rho = cutoff.rho, delta = cutoff.delta)
par(mfrow=c(1,1))
options(repr.plot.width=6, repr.plot.height=6)
plot(dclust$rho,dclust$delta, pch=20,cex=0.6,xlab='rho', ylab='delta')
points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20,cex=0.8)
text(dclust$rho[dclust$peaks]-0.2,dclust$delta[dclust$peaks]+0.2,labels=dclust$clusters[dclust$peaks])
abline(v=cutoff.rho)
abline(h=cutoff.delta)

densityClust <- dclust$clusters
densityClust <- as.data.frame(densityClust)
rownames(densityClust) <- cisTopicObject@cell.names
colnames(densityClust) <- 'densityClust'
densityClust[,1] <- as.factor(densityClust[,1])
cisTopicObject <- addCellMetadata(cisTopicObject, densityClust)


par(mfrow=c(1,1))
plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, 
             colorBy=c('densityClust'), 
             cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, 
             col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)

plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, 
             colorBy=c('densityClust'), 
             cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, 
             col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)

plotFeatures(cisTopicObject, method='DM', target='cell', topic_contr=NULL, 
             colorBy=c('densityClust'), 
             cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, 
             col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)


########################################################
########################################################
# Section : test Seurat for scATAC-seq data
# 
########################################################
########################################################
library(Signac)
library(Seurat)
#library(GenomeInfoDb)
#library(EnsDb.Hsapiens.v75)
library(ggplot2)
set.seed(1234)

pbmc <- CreateSeuratObject(
  counts = mat,
  assay = 'peaks',
  project = 'scATAC_earlyEmbro',
  min.cells = 10000,
  meta.data = NULL
)

pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q75')
pbmc <- RunSVD(object = pbmc, n = 100,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi',
  irlba.work=500
)

pbmc = ScaleData(pbmc, features = VariableFeatures(pbmc))
#ms <- RunPCA(object = ms, features = VariableFeatures(ms), verbose = FALSE)
pbmc = RunPCA(pbmc, npcs = 100, assay = 'peaks', features = VariableFeatures(pbmc), verbose = FALSE)

pbmc <- FindNeighbors(object = pbmc, reduction = 'pca', dims = 2:nb.pcs)
pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3, resolution = 1.0)

nb.pcs = 60; n.neighbors = 20; min.dist = 0.1;
pbmc <- RunUMAP(object = pbmc, reduction = 'pca', dims = 2:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(object = pbmc, label = TRUE, reduction = 'umap') + NoLegend()

#pbmc = RunTSNE(pbmc, reduction = 'lsi', dims = 2:nb.pcs)
#DimPlot(object = pbmc, label = TRUE, reduction = 'tsne') + NoLegend()
save(pbmc, file = paste0(RdataDir, 'pbmc_Seurat_TFIDtransform_pca_umap.Rdata'))

