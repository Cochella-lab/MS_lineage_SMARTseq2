##########################################################################
##########################################################################
# Project:
# Script purpose: test Seurat package: normalization and visulization
# Usage example:  
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Nov 15 12:31:15 2019
##########################################################################
##########################################################################
path2AleksFolder = '/Volumes/groups/cochella/Aleks/bioinformatics/GitHub/scRNAseq_MS_lineage'
version.DATA = 'scRNA_8613_full'
version.analysis =  paste0(version.DATA, '_20191029')

dataDir = paste0("../data/")
resDir = paste0("../results/", version.analysis)
tabDir = paste0("../results/", version.analysis, "/tables/")
RdataDir = paste0("../results/", version.analysis, "/Rdata/")
RdataDirfromAleks = paste0(path2AleksFolder, "/results/", version.analysis, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

library(Seurat)
library(velocyto.R)
library(SeuratWrappers)
library(ggplot2)
########################################################
########################################################
# Section : test the normalization SCTransform in Seurat
# for outputs from feature counts (ours) and velocity.py output
########################################################
########################################################
load(file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_bcMNN.Rdata'))

Test.Velocity.py.output = FALSE
if(Test.Velocity.py.output){
  # download an example for test
  #curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', destfile= '~/Downloads/SCG71.loom')
  #ldat <- ReadVelocity(file = "~/Downloads/SCG71.loom")
  ldat <- ReadVelocity(file = "/Volumes/groups/cochella/Aleks/bioinformatics/raw_ngs_data/global.loom")
  ms <- as.Seurat(x = ldat)
  ms <- SCTransform(object = ms, assay = "spliced")
}

ms0 = as.Seurat(sce, counts = 'counts', data = 'logcounts', assay = "RNA")

ms <- SCTransform(object = ms) # new normalization from Seurat

ms <- RunPCA(object = ms, verbose = FALSE)
ms <- FindNeighbors(object = ms, dims = 1:20)
ms <- FindClusters(object = ms)

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
# Section : test batch correction in Seurat
# 1) MNN calling from Seurat 
# 2) CCA from Seurat
########################################################
########################################################
library(scater)
library(SingleCellExperiment)
library(scran)

sce <- as.SingleCellExperiment(ms)
batches = sce$request 
bc.uniq = unique(batches)
sce$batches <- batches

Norm.Vars.per.batch = TRUE # HVGs for each batch or not 
Rescale.Batches = FALSE # scaling data in each batch or not 
k.mnn = 20
cos.norm = TRUE
nb.pcs = 50

batch.sequence.to.merge = c('R7130', 'R8612', 'R8526', 'R7926', # 3 plates for each request
                            'R6875','R7116','R8613','R8348') # 1 plate for each request
order2correct = match(batch.sequence.to.merge, bc.uniq) 

source("scRNAseq_functions.R")
#HVGs = find.HVGs(sce, Norm.Vars.per.batch = Norm.Vars.per.batch, method = "scran", ntop = 2000)
HVGs = VariableFeatures(ms)
gene.chosen = match(HVGs, rownames(sce))


cat("nb of HGV : ", length(gene.chosen), "\n")

original = list()
fscs = c()
#original0 = list()
for(n in 1:length(bc.uniq)){
  #xx = nout[[n]];
  if(Rescale.Batches){
    original[[n]] = logcounts((nout[[n]][gene.chosen, ]))
  }else{
    original[[n]] = logcounts((sce[gene.chosen, which(sce$batches == bc.uniq[n])])) 
  }
  fscs = c(fscs, sce$FSC_log2[which(sce$batches == bc.uniq[n])])
}

set.seed(1001)
mnn.out <- do.call(fastMNN, c(original, list(k=k.mnn, cos.norm = cos.norm, d=nb.pcs, auto.order=order2correct,
                                             approximate=TRUE)))

reducedDim(sce, "sct_MNN") <- mnn.out$corrected
sce$mnn_Batch <- as.character(mnn.out$batch)
sce

pbmc = as.Seurat(sce)


ms <- RunUMAP(object = pbmc, reduction = 'PCA', dims = 1:20, n.neighbors = 30)
DimPlot(ms, reduction = "umap", group.by = 'request')

ms <- RunUMAP(object = pbmc, reduction = 'MNN', dims = 1:20, n.neighbors = 30)
DimPlot(ms, reduction = "umap", group.by = 'request')

ms <- RunUMAP(object = pbmc, reduction = 'sct_MNN', dims = 1:20, n.neighbors = 30)
DimPlot(ms, reduction = "umap", group.by = 'request')


pbmc.list <- SplitObject(ms, split.by = "request")
pbmcsca <- RunFastMNN(object.list = pbmc.list, features = 2000, reduction.name = 'mnn_sct', 
                      k=20, cos.norm=TRUE, ndist=3, d=50, approximate=FALSE, auto.order = order2correct)

pbmcsca <- RunUMAP(pbmcsca, reduction = "mnn", dims = 1:30)
pbmcsca <- FindNeighbors(pbmcsca, reduction = "mnn", dims = 1:30)
pbmcsca <- FindClusters(pbmcsca)
DimPlot(pbmcsca, group.by = c("Method", "ident", "CellType"), ncol = 3)


pbmc.list <- SplitObject(ms, split.by = "request")
for (i in names(pbmc.list)) {
  pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = FALSE)
}

pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features)
k.filter <- min(sapply(pbmc.list, ncol)) 
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", 
                                       anchor.features = pbmc.features, k.filter = k.filter)

pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT")

pbmc.integrated <- RunPCA(object = pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunUMAP(object = pbmc.integrated, dims = 1:30)

DimPlot(pbmc.integrated, reduction = "umap", group.by = 'request')

plots <- DimPlot(pbmc.integrated, group.by = c("Method", "CellType"), combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
                                                                                                              byrow = TRUE, override.aes = list(size = 2.5))))
CombinePlots(plots)


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
ms
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


pdf(paste0("~/Desktop/linage_identification/","ms_fbxx.pdf"), width=8, height = 8)

for(n in 1:length(marker.genes.ms)) {
  
  FeaturePlot(ms, features = c(marker.genes.ms[n]))
  
}
dev.off()

FeaturePlot(ms, features = c(marker.genes.ms))
