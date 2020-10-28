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


find.MSxpapap.in.iteration.32 = function()
{
  ids.sels = c('likely_MSxpapap', 'mixture_MSxppapp_MSxpappp', 'likely_MSxpaaap', 
               'MSxpapa', 'MSxppap', 'MSxpapp', 'mixture_MSxpppa_MSxpppp_MSxpapa_MSxpapp_MSxppap', 'likely_MSxppppx')
  
  ids.left = setdiff(ids.current, ids.sels)
  nchar(ids.left)
  
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])
  sub.obj = subset(seurat.obj, cells = cells.sels)
  
  xx = table(sub.obj$seurat_clusters)
  xx[which(xx > 0)]
  
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  sub.obj$pred.ids = sub.obj$predicted.ids.seurat.keep
  
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()
  
  barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)
  
  # sub.obj$manual.annot.ids = sub.obj$predicted.ids.seurat.keep
  ##########################################
  # check potential ids for selected clusters
  ##########################################
  # #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
  # threshold = 0.7
  # predicted.ids = sub.obj$scmap.pred.id.500
  # #predicted.ids[which(sub.obj$scmap.corr.500 < threshold)] = 'unassigned'
  # 
  # if(Refine.annotated.ids){
  #   counts = table(predicted.ids, sub.obj$manual.annot.ids)
  #   counts.seurat = table(as.character(sub.obj$seurat.pred.id), sub.obj$manual.annot.ids)
  # }else{
  #   counts = table(predicted.ids, as.character(sub.obj$seurat_clusters))
  #   counts.seurat = table(as.character(sub.obj$seurat.pred.id), as.character(sub.obj$seurat_clusters))
  # }
  counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters_split))
  barplot(counts, main="cluster compositions by scmap ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))
  # 
  # barplot(counts.seurat, main="cluster compositions by seurat ",
  #         xlab=NULL, col=c(1:nrow(counts)), las = 2,
  #         legend = rownames(counts))
  
  #counts[, match(c('31', '28', '52'), colnames(counts))]
  #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  Explore.umap.parameters.for.BWMcells = FALSE
  if(Explore.umap.parameters.for.BWMcells){
    
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    require(tictoc)
    tic()
    test.umap.params.for.BWM.cells(sub.obj, pdfname = 'BWM_middleCells_umap_param_terminalCell_dissect_mixtures_with_BWM_1_2.pdf',
                                   group.by = 'predicted.ids.seurat.keep', with_legend = TRUE,
                                   nfeatures.sampling = c(1000, 3000, 5000, 8000), nb.pcs.sampling = c(5, 10, 30, 50),
                                   n.neighbors.sampling = c(5, 10, 30, 50), 
                                   min.dist.sampling = c(0.01, 0.1)
    )
    toc()
    
  }
  
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
  
  p1 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat.keep', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
  p2 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()
  p1 + p2
  
  idntt = 'likely_MSxpapap'
  cells_to_show <- list(colnames(sub.obj)[which(sub.obj$manual.annot.ids == idntt)])
  p1 = DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
               cells.highlight = cells_to_show, cols.highlight = 'blue', sizes.highlight = 1,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  idntt = 'MSxpapap'
  cells_to_show <-  list(colnames(sub.obj)[which(sub.obj$predicted.ids.seurat.keep == idntt)])
  p2 = DimPlot(sub.obj, group.by = "pred.ids", reduction = 'umap', label = TRUE, repel = TRUE,
               cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1.5,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  p1 + p2
  
  FeaturePlot(sub.obj, reduction = 'umap', features = c('hnd-1', 'Y37E3.30', 'stn-2'))
  sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpapap')] = 'MSxpapap'
  
}



