##########################################################################
##########################################################################
# Project: scRNA-seq data analysis in MS lineage in c. elegans
# Script purpose: dissect MS terminal cells and manual annotation
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Oct  2 08:25:20 2020
##########################################################################
##########################################################################
verify.MSxpppp.MSxpppa.terminals = function(sub.obj)
{
  nfeatures = 3000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)
  
  nb.pcs = 30 # nb of pcs depends on the considered clusters or ids 
  n.neighbors = 50;
  min.dist = 0.01; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs), 
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()
  
  p2 = DimPlot(sub.obj, group.by = 'pred.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    ggtitle('pred.ids.seurat')
  p1 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()
  p1 + p2
  
  idntt = 'MSxpappp'
  cells_to_show <- list(colnames(sub.obj)[which(sub.obj$manual.annot.ids == idntt)])
  p1 = DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
               cells.highlight = cells_to_show, cols.highlight = 'blue', sizes.highlight = 1,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  cells_to_show <-  list(colnames(sub.obj)[which(sub.obj$predicted.ids.seurat.keep == idntt)])
  p2 = DimPlot(sub.obj, group.by = "pred.ids", reduction = 'umap', label = TRUE, repel = TRUE,
               cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1.5,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  p1 + p2
  
  p3 = FeaturePlot(sub.obj, reduction = 'umap', features = c('irx-1', 'sul-2'))
  
  (p1 + p2)/p3
  
  sub.obj$manual.annot.ids[which(sub.obj$pred.ids == 'MSxpppa')] = 'MSxpppa'
  
}

subset.terminal.mother.cells(sub.obj)
{
  ids.mothers = c('MSxpppa', 'MSxpppp', 'MSxppppx', 'MSxpppax', 'MSxpapp', 'MSxpapa', 'MSxppap')
  cells.mothers = unique(colnames(sub.obj)[!is.na(match(sub.obj$pred.ids, ids.mothers))| 
                                             !is.na(match(sub.obj$manual.annot.ids, ids.mothers))])
  ssub = subset(sub.obj, cells = cells.mothers)
  
  nfeatures = 3000;
  ssub <- FindVariableFeatures(ssub, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(ssub)), '\n')
  ssub = ScaleData(ssub, features = rownames(ssub))
  ssub <- RunPCA(object = ssub, features = VariableFeatures(ssub), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(ssub, ndims = 50)
  
  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids 
  n.neighbors = 20;
  min.dist = 0.01; spread = 1
  ssub <- RunUMAP(object = ssub, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs), 
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist)
  
  DimPlot(ssub, group.by = 'pred.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()
  
  DimPlot(ssub, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
}




