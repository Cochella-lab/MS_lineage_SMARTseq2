##########################################################################
##########################################################################
# Project: Aleks' singlce cell MS lineage 
# Script purpose: functions related to RNA velocity
# Usage example:  
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Nov 15 12:31:15 2019
##########################################################################
##########################################################################
########################################################
########################################################
# Section : test the normalization SCTransform in Seurat
# for outputs from feature counts (ours) and velocity.py output
########################################################
########################################################
library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(ggplot2)

Test.Velocity.py.output = TRUE
if(Test.Velocity.py.output){
  # download an example for test
  #curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', destfile= '~/Downloads/SCG71.loom')
  #ldat <- ReadVelocity(file = "~/Downloads/SCG71.loom")
  ldat <- ReadVelocity(file = "../data/raw_ngs_data/merged.loom")
  
  bm <- as.Seurat(x = ldat)
  bm <- SCTransform(object = bm, assay = "spliced")
  
  bm <- RunPCA(object = bm, verbose = FALSE)
  bm <- FindNeighbors(object = bm, dims = 1:20)
  bm <- FindClusters(object = bm)
  bm <- RunUMAP(object = bm, dims = 1:20)
  DimPlot(bm, reduction = "umap")
  
}

ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:20, n.neighbors = 30)
DimPlot(ms, reduction = "umap", group.by = 'request')

ms <- RunUMAP(object = ms, reduction = 'MNN', dims = 1:20, n.neighbors = 30)
DimPlot(ms, reduction = "umap")

ms.logtransform <- NormalizeData(ms0, assay = "RNA" )
ms.logtransform <- FindVariableFeatures(ms.logtransform, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(ms.logtransform)
ms.logtransform <- ScaleData(ms.logtransform, features = all.genes)

ms.logtransform <- RunPCA(object = ms.logtransform, verbose = FALSE)

ms.logtransform <- FindNeighbors(object = ms.logtransform, dims = 1:20)
ms.logtransform <- FindClusters(object = ms.logtransform)

ms.logtransform <- RunUMAP(object = ms.logtransform, reduction = 'pca', dims = 1:20, n.neighbors = 30)
DimPlot(ms.logtransform, reduction = "umap", group.by = 'request')


########################################################
########################################################
# Section : add veclocity
# 
########################################################
########################################################

ldat <- ReadVelocity(file = "/Volumes/groups/cochella/Aleks/bioinformatics/raw_ngs_data/global.loom")
ms <- as.Seurat(x = ldat)
ms <- SCTransform(object = ms, assay = "spliced")

ms <- RunVelocity(object = ms, deltaT = 1, kCells = 25, fit.quantile = 0.02)

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = ms)))
names(x = ident.colors) <- levels(x = ms)
cell.colors <- ident.colors[Idents(object = ms)]
names(x = cell.colors) <- colnames(x = ms)
show.velocity.on.embedding.cor(emb = Embeddings(object = ms, reduction = "umap"), vel = Tool(object = ms, 
                                                                                                  slot = "RunVelocity"), n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)


FeaturePlot(ms, features = c("nhr-67", "pha-4", "hnd-1", "unc-120", "let-381", "sfrp-1", "cwn-2", "ceh-13", "F28B4.4", "end-1", "tbx-37", "tbx-35"))
FeaturePlot(ms, features = c("unc-30", "let-381", "sfrp-1", "ceh-27", "ceh-32", "ceh-34"))

FeaturePlot(ms, features = c("cutl-2", "D1005.2", "K08B4.2", "noah-2", "let-4"))
FeaturePlot(ms, features = c("pha-4", "hnd-1"))
FeaturePlot(ms, features = c("B0310.2"))

FeaturePlot(ms, features = c("fbxa-81", "fbxa-137"))


DimPlot(ms)
DoHeatmap(ms, features = VariableFeatures(ms)[1:200], size = 4, angle = 90)  
NoLegend()

ms[["percent.mt"]] <- PercentageFeatureSet(ms, pattern = "^MT-")
VlnPlot(ms, features = c( "percent.mt"), ncol = 1)
head(ms@meta.data, 5)

ms.copy <- ms

ms.copy <- RunUMAP(ms.copy, dims = 1:15, n.components = 3)

x.ms <- ms.copy@reductions$umap
y.ms <- ms.copy@reductions$umap[,2]
z.ms <- ms.copy@reductions$umap[,3]

x.ms$umap$UMAP_1

marker.genes.ms <- row.names(ms)[grep("^fbx.", row.names(ms))]


FeaturePlot(ms, features = c(marker.genes.ms))
