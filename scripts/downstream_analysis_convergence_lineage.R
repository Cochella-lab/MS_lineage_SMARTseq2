##########################################################################
##########################################################################
# Project: MS lineage 
# Script purpose: analyze convergence lineages and identify potential regulators
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Dec 16 11:55:32 2020
##########################################################################
##########################################################################
library(pheatmap)
library(RColorBrewer)
library(grid)
library(Seurat)
library(scater)
library(SingleCellExperiment)
library(scran)

aggregate.cells.across.ids = function(seurat.obj)
{
  ids.current = names(table(seurat.obj$manual.annot.ids))
  cat(length(ids.current), ' annotated ids in total \n')
  
  ids.to.drop = ids.current[grep('doublets|like', ids.current)]
  if(length(ids.to.drop) >0) {
    ids.sel = setdiff(ids.current, ids.to.drop)
    cat('ids to drop : \n')
    print(ids.to.drop)
  }else{
    cat('no ids dropped \n')
    ids.sel = ids.current
  }
  
  cat(length(ids.sel), ' annotated ids were selected to aggregate \n')
  
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sel))])
  
  cat(length(cells.sels), ' cells selected to annotate \n')
  
  sub.obj = subset(seurat.obj, cells = cells.sels)
  
  # convert seurat object to SingleCellExperiment 
  sce = as.SingleCellExperiment(sub.obj)
  
  summed <- aggregateAcrossCells(sce, 
                                 id=colData(sce)[, c("manual.annot.ids")])
  summed
  
  # Creating up a DGEList object for use in edgeR:
  library(edgeR)
  y <- DGEList(counts(summed), samples=colData(summed))
  #y
  
  # normalize pseudo-bulk 
  y <- calcNormFactors(y)
  y$samples
  
  par(mfrow=c(2,2))
  for (i in seq_len(ncol(y))) {
    plotMD(y, column=i)
  }
  
  #plotMDS(cpm(y, log=TRUE), 
  #        col=ifelse(y$samples$tomato, "red", "blue"))
  
  # output is DGEList object from edgeR 
  return(y) 
  
}

compare.convergence.lineages.with.others = function(y, seurat.obj,  method = c('euclidean', 'correlation', 'jsd'))
{
  library("pheatmap")
  library("RColorBrewer")
  #library(philentropy)
  
  tfs = readxl::read_xlsx('../data/motifs_tfs/Table-S2-wTF-3.0-Fuxman-Bass-Mol-Sys-Biol-2016.xlsx', sheet = 1)
  
  ids.convergence = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
  ids.bwm = c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
  ids.phrx = c('MSxapa', 'MSxapap', 'MSxapapp', "MSpaaappp/MSxapappa")
  
  ##########################################
  # step 1: dimension reduction to visualize the convergence lineage, one bwm and one pharynx lineages    
  ##########################################
  ids.sel = unique(c('MSx', ids.bwm, ids.convergence, ids.phrx))
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sel))])
  sub.obj = subset(seurat.obj, cells = cells.sels)
  
  DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  USE.UMAP = FALSE
  if(USE.UMAP){
    nfeatures = 2000;
    sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
    sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
    sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
    ElbowPlot(sub.obj, ndims = 50)
    
    nb.pcs = 10 # nb of pcs depends on the considered clusters or ids 
    n.neighbors = 10;
    min.dist = 0.3; spread = 1
    sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs), 
                       spread = spread, n.neighbors = n.neighbors,
                       min.dist = min.dist, verbose = TRUE)
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) + 
      NoLegend()
  }
  
  USE.MDS = FALSE
  if(USE.MDS){
    library(pheatmap)
    library(RColorBrewer)
    library(grid)
    library(Seurat)
    library(scater)
    library(SingleCellExperiment)
    library(scran)
    sce = as.SingleCellExperiment(sub.obj)
    
    dec <- modelGeneVar(sce)
    plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
    curve(metadata(dec)$trend(x), col="blue", add=TRUE)
    
    top.hvgs <- getTopHVGs(dec, n=2000)
    
    sce <- runPCA(sce, subset_row=top.hvgs)
    reducedDimNames(sce)
    
    
    sce = runDiffusionMap(sce, ncomponents = 4, n_pcs = 50, k = 100)
    plotDiffusionMap(sce, ncomponents = c(1, 3), colour_by = "manual.annot.ids", text_by="manual.annot.ids")
    plotDiffusionMap(sce, ncomponents = c(1, 2), colour_by = "manual.annot.ids", text_by="manual.annot.ids")
    
    #install.packages("plot3D")
    library("plot3D")
    x = reducedDim(sce, 'DiffusionMap')[, 1]
    y = reducedDim(sce, 'DiffusionMap')[, 2]
    z = reducedDim(sce, 'DiffusionMap')[, 3]
    scatter3D(x, y, z, pch = 18,  theta = 20, phi = 20,
              main = "Iris data", xlab = "DC1",
              ylab ="DC2", zlab = "DC3")
    
    
    sce <- runMDS(sce, ncomponents=2, dimred = 'DiffusionMap', n_dimred = 3, scale_features = FALSE)
    plotMDS(sce, ncomponents = 2, colour_by = "manual.annot.ids", text_by="manual.annot.ids")
    
    #library(destiny)
    #ll.pca = reducedDim(sce, 'PCA')[, c(1:50)]
    #dm <- DiffusionMap(ll.pca, sigma = 'local', n_eigs = 3, k = 30, distance = 'euclidean', n_pcs = NA)
    #plot(dm)
    
  }
  
  ##########################################
  # step 2: comparer the distance or correlation between the convergence lineage, BWM and pharynx lienages
  ##########################################
  # here the input is the DGEList object from edgeR
  cpm = edgeR::cpm(y, log = TRUE, prior.count = 1)
  
  method = 'jsd'
  #colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  if(method == 'euclidean'){ sampleDists <- dist(t(cpm), method = 'euclidean'); xlim = c(0, 500)}
  
  if(method == 'correlation'){ sampleDists = cor(cpm, method = 'pearson'); xlim = c(0, 1)}
  
  if(method == 'jsd') {
    library(cummeRbund)
    cpm2 = log2(2^cpm + 1)
    probs <- cummeRbund::makeprobs(cpm2)
    sampleDists <-cummeRbund::JSdist(probs)
    #sampleDists = JSD(t(as.matrix(cpm)), unit = 'log2')
  }
  
  sampleDistMatrix <- as.matrix(sampleDists)
  
  
  pdfname = paste0(resDir, "/convergence_lineage_compared_with_one.BWM.lineage_one.Pharynx.lineage_", method, ".pdf")
  pdf(pdfname, width=12, height = 8)
  par(cex =0.7, mar = c(5,8,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  cc = matrix(NA, ncol = length(ids.convergence), nrow = length(c(ids.bwm, ids.phrx)))
  colnames(cc) = ids.convergence
  rownames(cc) = c(ids.bwm, ids.phrx)
  
  #ii1 = match(ids.convergence[-c(1:2)], colnames(sampleDistMatrix))
  #ii2 = match(ids.phrx, rownames(sampleDistMatrix))
  cc.vec = c()
  cc.vec.names = c()
  for(n in 1:ncol(cc))
  {
    # n = 1
    ii1 = match(colnames(cc)[n], colnames(sampleDistMatrix))
    ids2compare = c(ids.bwm[which(nchar(ids.bwm) <= (nchar(colnames(cc)[n]) + 1))], 
                    ids.phrx[which(nchar(ids.phrx) == nchar(colnames(cc)[n]))])
    ii2 = match(ids2compare, rownames(sampleDistMatrix))
    
    cc[match(ids2compare, rownames(cc)), n] = sampleDistMatrix[ii2, ii1]
    cc.vec = c(cc.vec, sampleDistMatrix[ii2, ii1])
    cc.dist = sampleDistMatrix[ii2, ii1]
    cc.vec.names = c(cc.vec.names, paste0(colnames(cc)[n], '_vs_', ids2compare))
    
    barplot(cc.dist, beside = TRUE, col=c(1:length(cc.dist)), main = colnames(cc)[n], 
            #legend.text = rownames(cc.vec), args.legend = c(x = 'topleft', bty = 'n'), 
            names.arg = names(cc.dist), las = 2, horiz = TRUE,
            xlim = xlim)
    
  }
  
  names(cc.vec) = cc.vec.names
  
  dev.off()
  
  ids.mothers = c('MSx', 'MSxa', 'MSxap')
   
  pdfname = paste0(resDir, "/convergence_lineage_compared_mothers_sisters_", method, ".pdf")
  pdf(pdfname, width=10, height = 6)
  par(cex =0.7, mar = c(4,10,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(n in 1:length(ids.mothers))
  {
    # n = 1
    ids2compare = c(ids.mothers[n], paste0(ids.mothers[n], c('p', 'a')))
    #ii1 = match(list1.ids, colnames(sampleDistMatrix))
    #ii2 = match(list2.ids, rownames(sampleDistMatrix))
    compares = sampleDistMatrix[match(ids2compare, rownames(sampleDistMatrix)), 
                                match(ids2compare, colnames(sampleDistMatrix))]
    compares[lower.tri(compares, diag = TRUE)] <- NA
    cc.names = expand.grid(colnames(compares), '_vs_', rownames(compares))
    cc.names = paste0(cc.names[, 3], cc.names[,2], cc.names[,1]) 
    cc = as.vector(compares)
    names(cc) = cc.names
    cc = cc[!is.na(cc)]
    
    barplot(cc, beside = TRUE, col = c(1:length(cc)), horiz = TRUE,
            names.arg = names(cc),
           las = 2, xlim = xlim)
    
  }
  
  dev.off()
  
}

find.regulators.for.convergence.lineage = function(y)
{
  library("pheatmap")
  library("RColorBrewer")
    
  # here the input is the DGEList object from edgeR
  cpm = edgeR::cpm(y, log = TRUE, prior.count = 1)
  
  ids.convergence = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
  ids.bwm = c('MSx', 'MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppx', 'MSxppppp')
  ids.phrx = c('MSxapa', 'MSxapap', 'MSxapapp', "MSpaaappp/MSxapappa")
  ids.phrx2 = c('MSxaa', 'MSxaaa', 'MSpaaap', 'MSpaaapp', 'MSpaaappp/MSxapappa')
  ids.bwm2 = c('MSxpa', 'MSxpaa', 'MSxpaaa', 'MSxpaaap')
  
  ids.sel = c(ids.bwm, ids.convergence, ids.phrx, ids.phrx2, ids.bwm2)
  #indexs = c(c(1:7), c(2:7), c(3:6))
  #cols = c('black', rep('darkblue', 6), rep('red', 6), rep('orange', 4))
  
  tfs = readxl::read_xlsx('../data/motifs_tfs/Table-S2-wTF-3.0-Fuxman-Bass-Mol-Sys-Biol-2016.xlsx', sheet = 1)
  gene2plot = unique(c('pha-4', 'hnd-1', 'hlh-1', 'unc-120', 'nhr-67', tfs$`Public name`))
  gene2plot = gene2plot[!is.na(match(gene2plot, rownames(cpm)))]
  
  cpmxx = cpm[match(gene2plot, rownames(cpm)), match(ids.sel, colnames(cpm))]
  kk = apply(as.matrix(cpmxx), 1, function(x) !all(as.numeric(x)<2))
  cpmxx = cpmxx[kk, ]
  
  pdfname = paste0(resDir, "/convergence_lineage_regulators_profiles_v3.pdf")
  pdf(pdfname, width=15, height = 10)
  par(cex =1.0, mar = c(4,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(n in 1:nrow(cpmxx))
  {
    # n = 1
    if(!all(is.na(cpmxx[n]))){
      cat(n, ' -- ', rownames(cpmxx)[n], '\n')
      ylims = range(cpmxx[n, ]) + c(-1, 1)
      plot(c(1:7), cpmxx[n, c(1:7)], ylim = ylims, col = 'darkblue', type='b', main = rownames(cpmxx)[n], pch = 15, lwd = 2.0,
           ylab = 'log2(cpm)', xlab = 'cell ids')
      #points(c(2:6), cpmxx[n, c(2, 23:26)], col = 'darkblue', type = 'b', pch =17, lty = 2, lwd = 2.0)
      points(c(1:7), cpmxx[n, c(1, 8:13)], col = 'red', type = 'b', lwd = 4.0, pch = 16)
      points(c(3:7), cpmxx[n, c(9, 14:17)], col = 'black', type = 'b', pch =17)
      points(c(2:7), cpmxx[n, c(8, 18:22)], col = 'black', type = 'b', pch =17, lty = 2)
      
      
      abline(v = c(4, 5), lwd =2.0, lty = 'dashed', col = 'gray')
      legend('topright', c('bwm', 'convergence', 'pharynx'),lwd = c(1, 2, 1), pch = c(15, 16, 17),
             col = c('darkblue', 'red', 'black'),  adj = c(0, 0.6), bty = 'n')
      # add text
      text(c(1:7), cpmxx[n, c(1:7)], ids.bwm, pos = 3, cex = 0.8)
      text(c(2:7), cpmxx[n, c(8:13)], ids.convergence, pos = 3, cex = 0.8)
      text(c(4:7), cpmxx[n, c(14:17)], ids.phrx, pos = 3, cex = 0.7)
      text(c(3:7), cpmxx[n, c(18:22)], ids.phrx2, pos = 3, cex = 0.5)
      #text(c(3:6), cpmxx[n, c(23:26)], ids.bwm2, pos = 3, cex = 0.5)
    }
   
  }
  
  dev.off()
  
  ##########################################
  # supporting argument for asymmetric cell division 
  ##########################################
  ms.early = c('MSx', 
               'MSxp', 'MSxa', 
               'MSxpp', 'MSxpa', 'MSxap', 'MSxaa', 
               'MSxppp', 'MSxppa', 'MSxpap', 'MSxpaa', 'MSxapp', 'MSxapa', 'MSxaap', 'MSxaaa')
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ms.early))])
  sub.obj = subset(seurat.obj, cells = cells.sels)
  
  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'manual.annot.ids')
  
  #FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  VlnPlot(sub.obj, features = c("par-1", 'par-3', 'par-2', 'par-5',  'par-6', 'pkc-3', 'num-1',
                                "mex-5", 'mex-6', 'ref-1', 'ref-2', 'glp-1' #'rnt-1', 'bro-1'
                                ), ncol = 3,
          group.by = 'manual.annot.ids') +
    ggsave(paste0(resDir, '/genes_in_asymmetric.cell.division_Notch.signaling.pdf'),  width = 18, height = 16)
  
  VlnPlot(sub.obj, features = c('hnd-1', 'pha-4' #"par-1", 'par-3', 'par-2', 'par-5',  'par-6', 'pkc-3', 'num-1',
                                #"mex-5", 'mex-6', 'ref-1', 'ref-2', 'glp-1' #'rnt-1', 'bro-1'
  ), ncol = 2,  group.by = 'manual.annot.ids')
  
  FeaturePlot(seurat.obj, reduction = 'umap', features = c('ref-1', 'ref-2', 'glp-1', 'pop-1', 'hnd-1', 'pha-4', 'sys-1', 'mom-4',
                                                           'mom-5'))
  
  FeaturePlot(seurat.obj, reduction = 'umap', features = c('hnd-1', 'pha-4', #'pal-1',
                                                           'sys-1', 'bar-1', 'wrm-1', 'hmp-2', 'pop-1', 'unc-37', 
                                                           'cwn-2', 'cwn-1', 'mom-2',
                                                           'mom-5', 'cam-1')) +
  ggsave(paste0(resDir, '/genes_in_asymmetric.cell.division_Wnt.signaling_cwn-2.pdf'),  width = 20, height = 16)
  
}




########################################################
########################################################
# Section : characterize the BWM terminal cells
# 
########################################################
########################################################
characterize.bwm.terminal.cells = function(seurat.obj)
{
  terminals = c('MSxppppp', 'MSxppppa', 'MSxpppap', 'MSxpppaa', 
                'MSxppapp', 
                'MSxpappp', 'MSxpappa', 
                'MSxpapap', 
                'MSxpaaap', 
                #'MSxapppp', 'MSxapppa',
                'MSxappppx', 'MSxapppax'
  )
  
  setdiff(terminals, ids.current)
  
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, terminals))])
  sub.obj = subset(seurat.obj, cells = cells.sels)
  
  Idents(sub.obj) = sub.obj$manual.annot.ids
  
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('egl-1'), label = FALSE)
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  
  markers <- FindAllMarkers(sub.obj, only.pos = FALSE, min.pct = 0.25, test.use = 'MAST')
  
  saveRDS(markers, file = paste0(resDir, '/bwm_terminal.cells_markers_FindAllMarkers.rds'))
  
  markers = readRDS(file = paste0(resDir, '/bwm_terminal.cells_markers_FindAllMarkers.rds'))
  
  qval.cutoff = 0.01;
  logfc.cutoff = 0.5
  ntop = 20
  
  markers = markers[which(abs(markers$avg_logFC) > logfc.cutoff & markers$p_val_adj < qval.cutoff), ]
  
  topGenes <- markers %>% group_by(cluster) %>% top_n(n = ntop, wt = avg_logFC)
  gene.sels = markers$gene[which(markers$p_val_adj < qval.cutoff)]
  gene.sels = topGenes$gene
  p1 = DoHeatmap(sub.obj, features = gene.sels) + NoLegend()
  
  plot(p1) + ggsave(paste0(resDir, '/heatmap_BWMterminalCells_avg_logFC.', logfc.cutoff, 
                           ' _qval.', qval.cutoff, '_top', ntop, '.pdf'),
                    width = 10, height = 16)
  
  #genesToshow = c('egl-1', rownames(topgenes)[1:30])
  #DoHeatmap(sub.obj, features = genesToshow) +
  
  write.csv(markers, file = paste0(resDir, 'BWM_terminalCell_markerGenes_avg_logFC.', logfc.cutoff, 
                                   ' _qval.', qval.cutoff, '.csv'), row.names = FALSE)
  
  
  
}

########################################################
########################################################
# Section : compare programmed cell death lineage
# 
########################################################
########################################################
heatmap.for.cell.death.lineage.MSxaap = function(seurat.obj)
{
  
  ms.sels = c('MSxaap', 'MSxaapa', 'MSaaapp', 'MSxaaa', 'MSxapa')
  
  #ms.sels = c('MSxaap', 'MSxaaa', 'MSxapa')
  setdiff(ms.sels, ids.current)
  
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ms.sels))])
  sub.obj = subset(seurat.obj, cells = cells.sels)
  
  Idents(sub.obj) = sub.obj$manual.annot.ids
  
  FeaturePlot(sub.obj, reduction = 'umap', features = c('egl-1'), label = FALSE)
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  
  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'manual.annot.ids') + NoLegend() + 
    ggsave(paste0(resDir, '/cellDeath.lineage_MSxaap_size.pdf'),  width = 12, height = 8)
  
  VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
          group.by = 'manual.annot.ids') + NoLegend() + 
    ggsave(paste0(resDir, '/cellDeath.lineage_MSxaap_estimatedTiming.pdf'),  width = 10, height = 8)
  
  ##########################################
  # 
  ##########################################
  test.uamp.params = FALSE
  if(test.umap.params){
    nfeatures = 1000;
    sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
    #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
    sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
    sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
    ElbowPlot(sub.obj, ndims = 50)
    
    nb.pcs = 10 # nb of pcs depends on the considered clusters or ids 
    n.neighbors = 30;
    min.dist = 0.01; spread = 1
    sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs), 
                       spread = spread, n.neighbors = n.neighbors,
                       min.dist = min.dist, verbose = TRUE)
    
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) + 
      NoLegend()
  }
  
  cluster1.markers <- FindMarkers(sub.obj, ident.1 = 'MSxaap', min.pct = 0.25, logfc.threshold = 0.5, only.pos = TRUE)
  head(cluster1.markers, n = 5)
  topgenes = cluster1.markers[which(cluster1.markers$p_val<0.001), ]
  topgenes <- topgenes[order(-topgenes$avg_logFC), ]
  
  markers <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  
  DoHeatmap(sub.obj, features = c('egl-1', top10$gene)) + NoLegend()
  
  genesToshow = c('egl-1', rownames(topgenes)[1:30])
  
  DoHeatmap(sub.obj, features = genesToshow) +
    ggsave(paste0(resDir, '/heatmap_MSxaap_and_daughters_cellDeath.lineage.pdf'),  width = 14, height = 8)
  
   
}





