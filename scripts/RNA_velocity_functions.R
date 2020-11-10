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
# Section : test velocyto.R with velocity.py output being input
# sctransform HGV from seurat and functio all from seurat too
########################################################
########################################################
import.loomFiles.velocytoOut = function(seurat.obj) # import the loom file and process the data according to seurat.obj
{
  # download an example for test
  #curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', destfile= '~/Downloads/SCG71.loom')
  #ldat <- ReadVelocity(file = "~/Downloads/SCG71.loom")
  ldat <- ReadVelocity(file = "../data/velocyto_all_merged.loom")
  bm <- as.Seurat(x = ldat)
  
  DefaultAssay(bm) = 'unspliced'
  samples = colnames(bm)
  cells = basename(samples)
  cells = gsub('Aligned.sortedByCoord.out.bam', '', cells)
  cells = sapply(cells, function(x) unlist(strsplit(as.character(x), ':'))[2])
  cells = gsub('[#]', '_', cells)
  #cells = sapply(cells, function(x) unlist(strsplit(as.character(x), '[.]'))[2])
  bm$cells = cells
  
  Idents(bm) = bm$cells
  bm$old.ident = bm$orig.ident
  bm$orig.ident = bm$cells
  
  # extract metadata and unspliced matrix
  #unspliced = as.matrix(bm@assays$unspliced@counts)
  #metadata = as.data.frame(bm@meta.data)
  #RenameCells(bm, new.names = cells)
  
  mm = match(seurat.obj$samples, cells)
  jj = match(rownames(seurat.obj), rownames(bm))
  
  length(which(is.na(mm)))
  length(which(is.na(jj)))
  
  #cell.sel = cells[mm];
  #cell.sel = cells[which(!is.na(cell.sel))]
  pm = bm[jj[which(!is.na(jj))], mm[which(!is.na(mm))]]
  
  mm = match(pm$cells, seurat.obj$samples)
  pm$annoted.ids = seurat.obj$manual.annot.ids[mm]
  pm$BWM.cells = seurat.obj$BWM.cells[mm]
  pm$seurat_clusters = seurat.obj$seurat_clusters[mm]
  
  saveRDS(pm, file = paste0(RdataDir, 'unsplicedData_processed_by_velocyto.rds'))
  
}

process.unspliced.mRNA.check.cell.identities = function(pm)
{
  pm = readRDS(file = paste0(RdataDir, 'unsplicedData_processed_by_velocyto.rds'))
  DefaultAssay(pm) = 'unspliced'
  
  ## use sran normalization
  library(SingleCellExperiment)
  library(scater)
  sce = as.SingleCellExperiment(pm, assay = 'unspliced')
  colnames(sce) = sce$cells
  
  sce$library.size = apply(counts(sce), 2, sum)
  qclust <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters = qclust)
  sce <- logNormCounts(sce, log = TRUE, pseudo_count = 1)
  #sce <- logNormCounts(sce, log = FALSE, size_factors = NULL)
  par(mfrow = c(1, 1))
  plot(sce$library.size/1e6, sizeFactors(sce), log="xy", xlab="Library size (millions)", 
       ylab="Size factor", cex = 0.5)
  
  pm = as.Seurat(sce, counts = 'counts', data = 'logcounts', assay = "unspliced") # scran normalized data were kept in Seurat
  pm = pm[, which(!is.na(pm$BWM.cells))]
  
  #pm <- SCTransform(object = pm, assay = "unspliced", variable.features.n = 3000)
  
  #FeatureScatter(pm, feature1 = "nCount_unspliced", feature2 = "nFeature_unspliced")
  #pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  #pm <- NormalizeData(pm, normalization.method = "LogNormalize", scale.factor = median(pm$nCount_unspliced))
  
  #ids.bwm = names(table(pm$annoted.ids[!is.na(pm$BWM.cells)], useNA = 'ifany'))
  #cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.bwm))])
  #sub.obj = subset(seurat.obj, cells = cells.sels)
  #VlnPlot(pm, features = c("nFeature_unspliced", "nCount_unspliced"), ncol = 2)
  
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  require(tictoc)
  tic()
  test.umap.params.for.BWM.cells(pm, 
                                 pdfname = 'UMAP_premRNA_param_TEST_BWM_all_scranNormalization.pdf',
                                 group.by = 'annoted.ids', with_legend = FALSE, weight.by.var = TRUE,
                                 nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30, 50),
                                 n.neighbors.sampling = c(10, 30, 50), 
                                 min.dist.sampling = c(0.01, 0.1)
  )
  
  toc()
  
  
  pm <- FindVariableFeatures(pm, selection.method = "vst", nfeatures = 3000)
  VariableFeaturePlot(pm)
  pm <- ScaleData(pm, features = rownames(pm))
  
  pm = RunPCA(pm, features = VariableFeatures(object = pm), verbose = FALSE, weight.by.var = TRUE)
  
  
  nb.pcs = 30 # nb of pcs depends on the considered clusters or ids 
  n.neighbors = 50;
  min.dist = 0.1
  pm <- RunUMAP(object = pm, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs), 
                     spread = 1, n.neighbors = n.neighbors,
                     min.dist = min.dist, verbose = TRUE)
  
  DimPlot(pm, group.by = 'annoted.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) + 
    NoLegend()
  
  saveRDS(pm, file = paste0(RdataDir, 'unsplicedData_processed_by_velocyto_scranNormalized_BWM.rds'))
  
}

run.RNAvelocity.with.velocyto.in.R = function(pm)
{
  library(Seurat)
  library(velocyto.R)
  library(SeuratWrappers)
  library(ggplot2)
  
  pm = readRDS(file = paste0(RdataDir, 'unsplicedData_processed_by_velocyto.rds'))
  pm = RenameCells(pm, new.names = pm$cells)
  #colnames(pm) = pm$cells
  #pm = subset(pm, cells = WhichCells(pm, ))
  
  pm = pm[c(1:nrow(pm)), which(!is.na(pm$BWM.cells))]
  
  DefaultAssay(pm) = 'spliced'
  Idents(pm) = pm$annoted.ids
  
  #FeatureScatter(pm, feature1 = "nCount_unspliced", feature2 = "nFeature_unspliced")
  #pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  pm <- NormalizeData(pm, normalization.method = "LogNormalize", scale.factor = median(pm$nCount_spliced))
  
  pm <- FindVariableFeatures(pm, selection.method = "vst", nfeatures = 3000)
  VariableFeaturePlot(pm)
  pm <- ScaleData(pm, features = rownames(pm))
  pm = RunPCA(pm, features = VariableFeatures(object = pm), verbose = FALSE, weight.by.var = TRUE)
  
  pm <- FindNeighbors(object = pm, dims = 1:20)
  pm <- FindClusters(object = pm)
  
  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids 
  n.neighbors = 30;
  min.dist = 0.1
  pm <- RunUMAP(object = pm, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs), 
                spread = 1, n.neighbors = n.neighbors,
                min.dist = min.dist, verbose = TRUE)
  
  DimPlot(pm, group.by = 'annoted.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) + 
    NoLegend()
  
  
  pm <- RunVelocity(object = pm, spliced = 'spliced', unspliced = 'unspliced',
                    deltaT = 1, kCells = 25, fit.quantile = 0.02)
  
  Idents(pm) = pm$annoted.ids
  ident.colors <- (scales::hue_pal())(n = length(x = levels(x = pm)))
  names(x = ident.colors) <- levels(x = pm)
  cell.colors <- ident.colors[Idents(object = pm)]
  names(x = cell.colors) <- colnames(x = pm)
  
  show.velocity.on.embedding.cor(emb = Embeddings(object = pm, reduction = "umap"), vel = Tool(object = pm, slot = "RunVelocity"), 
                                 n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
                                 cex = 0.8, arrow.scale = 3, show.grid.flow = TRUE, min.grid.cell.mass = 0.5, grid.n = 40, arrow.lwd = 1, 
                                 do.par = FALSE, cell.border.alpha = 0.1)
  
  First.Test.Velocity.py.output = FALSE
  # original code found https://htmlpreview.github.io/?https://github.com/satijalab/seurat.wrappers/blob/master/docs/velocity.html
  if(First.Test.Velocity.py.output){
    # download an example for test
    #curl::curl_download(url = 'http://pklab.med.harvard.edu/velocyto/mouseBM/SCG71.loom', destfile= '~/Downloads/SCG71.loom')
    #ldat <- ReadVelocity(file = "~/Downloads/SCG71.loom")
    ldat <- ReadVelocity(file = "../data/velocyto_all_merged.loom")
    
    bm <- as.Seurat(x = ldat)
    bm <- SCTransform(object = bm, assay = "spliced")
    
    bm <- RunPCA(object = bm, verbose = FALSE)
    bm <- FindNeighbors(object = bm, dims = 1:20)
    bm <- FindClusters(object = bm)
    bm <- RunUMAP(object = bm, dims = 1:20)
    DimPlot(bm, reduction = "umap")
    
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
    pm <- RunVelocity(object = pm, deltaT = 1, kCells = 25, fit.quantile = 0.02)
    
    ident.colors <- (scales::hue_pal())(n = length(x = levels(x = pm)))
    names(x = ident.colors) <- levels(x = pm)
    cell.colors <- ident.colors[Idents(object = pm)]
    names(x = cell.colors) <- colnames(x = pm)
    show.velocity.on.embedding.cor(emb = Embeddings(object = pm, reduction = "umap"), vel = Tool(object = pm, 
                                                                                                 slot = "RunVelocity"), 
                                   n = 200, scale = "sqrt", cell.colors = ac(x = cell.colors, alpha = 0.5), 
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
    
  }
  
}

