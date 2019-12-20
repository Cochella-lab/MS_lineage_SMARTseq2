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
# Section : test trojectory methods
# The methods to test are following:
# 1) Monocle2
# 2) Slingshot
# 3) TSCAN
# 4) Celltree
# 5) PAGA
# 6) ELPiGrapha (possible)
# destiny is installed but not to test in the first time
# imputation method MAGIC is also installed
########################################################
########################################################
Split.data.per.timeWindow = function(ms)
{
  ## for the sake of memory empty some slots in seurat object
  ms@tools$RunFastMNN = list() 
  ms@reductions$umap_mnn = list()
  ms@commands = list()
  ms@tools$RunFastMNN = list() # reduce a little bit the size of seurat object
  #knn = data.frame(ms@graphs$RNA_nn)
  #snn = data.frame(ms@graphs$RNA_snn)
  
  ##########################################
  # look into the HVGs 
  ##########################################
  source.my.script('timingEst_functions.R')
  dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
  load(file = paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval_plus_timepoints.Rdata"))
  #head(timers)
  timerGenes.pval=0.001; timerGenes.ac=0.5;
  sels.timerGenes = which(timers$ac.max > timerGenes.ac & timers$pval.box < timerGenes.pval)
  timers = timers[sels.timerGenes, -c(1:4)]
  hvgs = VariableFeatures(ms)
  
  length(intersect(hvgs, rownames(timers)))/length(hvgs)
  
  ##########################################
  #  # fix time windows
  ##########################################
  timingEst = as.numeric(as.character(ms$timingEst))
  kk = which(timingEst <= 440)
  ms = subset(ms, cells = colnames(ms)[kk])
  
  #kk1 = which(timingEst <= 150)
  kk1 = which(timingEst <= 220)
  kk2 = which(timingEst > 220 & timingEst <= 280)
  kk3 = which(timingEst >280 & timingEst <=340)
  kk4 = which(timingEst > 340 & timingEst<=440)
  cat(length(kk1), length(kk2), length(kk3), length(kk4), '\n')
  
  ms$timingEst.group[kk1] = 1
  ms$timingEst.group[kk2] = 2
  ms$timingEst.group[kk3] = 3
  ms$timingEst.group[kk4] = 4
  #ms$timingEst.group[kk5] = 5
  
  timeWindowPlot = FALSE
  if(timeWindowPlot){
    p1 = DimPlot(ms, reduction = 'umap', cols.highlight = "red", cells.highlight = as.list(kk1)) + NoLegend()
    p2 = DimPlot(ms, reduction = 'umap', cols.highlight = "red", cells.highlight = as.list(kk2)) + NoLegend()
    p3 = DimPlot(ms, reduction = 'umap', cols.highlight = "red", cells.highlight = as.list(kk3)) + NoLegend()
    p4 = DimPlot(ms, reduction = 'umap', cols.highlight = "red", cells.highlight = as.list(kk4)) + NoLegend()
    #p5 = DimPlot(ms, reduction = 'umap', cols.highlight = "red", cells.highlight = as.list(kk5)) + NoLegend()
    
    CombinePlots(plots = list(p1, p2, p3, p4), ncol = 2)
  }
  
  HVG.using.scran = FALSE
  if(HVG.using.scran){
    library(scater)
    library(Seurat)
    library(scran)
    
    sce = as.SingleCellExperiment(ms)
    #rm(ms) # to save the memory
    #save(sce, file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_timingEst_Normed_bc_timeWidow_SCE.Rdata'))
    dec.sce = modelGeneVar(sce)
    length(getTopHVGs(dec.sce, var.threshold = 0))
    for(n in 1:4)
    {
      #n = 1
      jj = which(sce$timingEst.group == n)
      xx = sce[jj, ]
      
      dec.xx <- modelGeneVar(xx)
      dec.xx[order(dec.xx$bio, decreasing=TRUE), ]
      
      hvg.xx <- getTopHVGs(dec.xx, var.threshold = 0)
      cat('HVGs --', length(hvg.xx), '\n')
      
      write.table(match(hvg.xx, rownames(sce)), file = paste0(resDir, '/alex_graph/hvg_index_', n, '.txt'), row.names = FALSE, col.names = FALSE)
      write.table(logcounts(xx), file = paste0(resDir, '/alex_graph/X_norm_', n, '.txt'), sep='\t', row.names = FALSE, col.names = FALSE)
    }
    
    
  }else{
    
    saveTableForStitch = FALSE
    ms = FindVariableFeatures(ms, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    
    HVGs = cbind(VariableFeatures(ms), rep(0, length(VariableFeatures(ms))))
    nfeatures = c(500, 1000, 1500, 2000)
    
    for(n in 1:4){
      jj = which(ms$timingEst.group == n)
      ms1 = subset(ms, cells = colnames(ms)[jj])
      ms1 <- FindVariableFeatures(ms1, selection.method = "vst", nfeatures = nfeatures[n], verbose = FALSE)
      hvgs = VariableFeatures(ms1)
      
      HVGs = rbind(HVGs, cbind(hvgs, rep(n, length(hvgs))))
      #ms1 = ScaleData(ms1, features = rownames(ms1))
      #ms1 <- RunPCA(object = ms1, features = VariableFeatures(ms1), verbose = FALSE)
      #nb.pcs = 20; n.neighbors = 20; min.dist = 0.3;
      #ms1 <- RunUMAP(object = ms1, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
      #p1 = DimPlot(ms1, reduction = "umap", group.by = 'timingEst')
      #p0 = DimPlot(ms, reduction = 'umap', group.by = 'timingEst')
      #plot_grid(p0, p1, ncol = 2)
      hvg.index = match(VariableFeatures(ms1), rownames(ms1))
      cat(n, '-- ', length(hvg.index), 'hvg -- ', 
          length(intersect(hvgs, rownames(timers)))/length(hvgs),  ' timer genes --',
          length(setdiff(hvgs, VariableFeatures(ms))), 'specific for this timeWindow\n')
      X_norm = as.data.frame(ms1@assays$RNA@data)
      if(saveTableForStitch){
        write.table(hvg.index, file = paste0(resDir, '/alex_graph/hvg_index_', n, '.txt'), row.names = FALSE, col.names = FALSE)
        write.table(X_norm, file = paste0(resDir, '/alex_graph/X_norm_', n, '.txt'), sep='\t', row.names = FALSE, col.names = FALSE)
      }
    }
    colnames(HVGs) = c('hvg', 'timewindow')
    HVGs = data.frame(HVGs, stringsAsFactors = FALSE)
    
    genes2use = unique(HVGs$hvg)
    #genes2use = intersect(unique(HVGs$hvg), rownames(timers))
    #genes2use = setdiff(unique(HVGs$hvg), rownames(timers)) 
    cat('nb of genes to use --', length(genes2use), '\n')
    
    msx <- RunPCA(object = ms, features = genes2use, verbose = FALSE, npcs = 100, weight.by.var = FALSE)
    ElbowPlot(msx, ndims = 50)
    
    nb.pcs = 50; n.neighbors = 20; min.dist = 0.3; spread = 1;
    msx <- RunUMAP(object = msx, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist, spread = spread, metric = 'cosine')
    DimPlot(msx, reduction = "umap", group.by = 'timingEst') + ggtitle('HVGs test')
    
  }
  
}

########################################################
########################################################
# Section : test function (not used now)
# 
########################################################
########################################################

clean.NA.for.metadata = function(md)
{
  md = md[, c(1:6, 70:74, 121:130)]
  for(n in 1:ncol(md))
  {
    if(typeof(md[,n]) == 'character'){
      if(length(which(is.na(md[,n])==TRUE) >0 )) print(colnames(md)[n])
    } 
  }
}


Test.different.trajectory.methods = function()
{
  Test.Monocole = FALSE
  if(Test.Monocole){
    library(monocle)
    library(dplyr)
    library(reshape2)
    
    load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata')) 
    load(file = paste0("../results/R6875_scRNA_v1_20181106/Rdata/R6875_scRNA_v1_QCed_cells_genes_filtered_normalized_SCE.Rdata"))
    # importCDS(sce, import_all = FALSE)
    HSMM <- newCellDataSet(as.matrix(counts(sce)),
                           expressionFamily=VGAM::negbinomial.size())
    
    HSMM <- estimateSizeFactors(HSMM)
    HSMM <- estimateDispersions(HSMM)
    
    HSMM <- detectGenes(HSMM, min_expr = 0.1)
    print(head(fData(HSMM)))
    expressed_genes <- row.names(subset(fData(HSMM),
                                        num_cells_expressed >= 10))
    
    L <- log(exprs(HSMM[expressed_genes,]))
    melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))
    
    qplot(value, geom = "density", data = melted_dens_df) +
      stat_function(fun = dnorm, size = 0.5, color = 'red') +
      xlab("Standardized log(FPKM)") +
      ylab("Density")
    
    disp_table <- dispersionTable(HSMM)
    unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
    HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
    plot_ordering_genes(HSMM)
    
    plot_pc_variance_explained(HSMM, return_all = F, max_components = 40, norm_method = 'log') # norm_method='log'
    
    HSMM <- reduceDimension(HSMM, max_components = 2, norm_method = c("log"), 
                            reduction_method = 'ICA', verbose = T)
    
    HSMM <- clusterCells(HSMM, method = "DDRTree")
    plot_cell_clusters(HSMM, 1, 2, color = "Cluster")
    
    
    HSMM_myo = HSMM
    diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,], 
                                          fullModelFormulaStr = "~ Cluster", 
                                          reducedModelFormulaStr = "~1" 
    )
    ordering_genes <- row.names (diff_test_res[order(diff_test_res$pval), ])[1:200]
    
    HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
    plot_ordering_genes(HSMM_myo)
    
    HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2,
                                method = 'DDRTree')
    
    HSMM_myo <- orderCells(HSMM_myo)
    plot_cell_trajectory(HSMM_myo, color_by = "Cluster")
    
  }
  
  Test.cellTree = FALSE
  if(Test.cellTree){
    library(cellTree)
    library(scran)
    library(scater)
    #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata')) 
    load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
    
    lda.results = compute.lda(as.matrix(logcounts(sce.HVG.Brenneck)), k.topics=3:8, method="maptpx", log.scale = FALSE)
    # lda.results = compute.lda(as.matrix(logcounts(sce)), k.topics=6, method="Gibbs", log.scale = FALSE)
    
    mst.tree = compute.backbone.tree(HSMM_lda_model, only.mst=TRUE, rooting.method = "longest.path", width.scale.factor = 2.0)
    # plot the tree (showing topic distribution for each cell):
    mst.tree.with.layout = ct.plot.topics(mst.tree)
    
  }
  
  Test.Slingshot = FALSE
  if(Test.Slingshot){
    #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))
    #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))
    load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
    
    sim = sce.HVG.Brenneck
    
    library(slingshot, quietly = FALSE)
    
    assay(sim, "norm") <- exp(logcounts(sim))
    
    FQnorm <- function(counts){
      rk <- apply(counts,2,rank,ties.method='min')
      counts.sort <- apply(counts,2,sort)
      refdist <- apply(counts.sort,1,median)
      norm <- apply(rk,2,function(r){ refdist[r] })
      rownames(norm) <- rownames(counts)
      return(norm)
    }
    
    #assays(sim)$norm <- FQnorm(assays(sim)$counts)
    
    pca <- prcomp(t((assays(sim)$logcounts)), scale. = FALSE, center = TRUE)
    #xx = t((assays(sce)$logcounts))
    #vars = apply(xx, 2, sd)
    #ntop = 1000; 
    #xx = xx[, order(vars, decreasing = TRUE)]
    #xx = xx[, c(1:ntop)]
    #pca = prcomp(xx, scale. = TRUE, center = TRUE)
    
    rd1 <- pca$x[,1:2]
    plot(rd1, col = 'red', pch=16)
    
    
    library(destiny, quietly = TRUE)
    dm <- DiffusionMap(t(log(assays(sim)$norm)))
    rd2 <- cbind(DC1 = dm$DC1, DC2 = dm$DC2)
    
    plot(rd2, col = topo.colors(100), pch=16, asp = 1)
    
    reducedDims(sim) <- SimpleList(PCA = rd1, DiffMap = rd2)
    
    library(mclust, quietly = TRUE)
    
    cl1 <- Mclust(rd2)$classification
    colData(sim)$GMM <- cl1
    
    library(RColorBrewer)
    plot(rd2, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
    
    
    sce <- slingshot(sim, clusterLabels = 'GMM', reducedDim = 'DiffMap')
    
    summary(sce$slingPseudotime_1)
    
    colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
    plot(reducedDims(sce)$DiffMap, col = colors[cut(sce$slingPseudotime_1,breaks=100)], pch=16, asp = 1)
    lines(SlingshotDataSet(sce), lwd=2)
    
  }
  
  Test.ELPiGraph = FALSE
  if(Test.ELPiGraph){
    load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
    
    library(ElPiGraph.R)
    TreeEPG <- computeElasticPrincipalTree(X = t(logcounts(sce.HVG.Brenneck)), NumNodes = 7, NumEdges=6, Lambda = .03, Mu = .01)
    
  }
}

########################################################
########################################################
# Section : test velocyto.R with velocity.py output being input
# sctransform HGV from seurat and functio all from seurat too
########################################################
########################################################
test.velocity.in.R = function()
{
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
  
}

