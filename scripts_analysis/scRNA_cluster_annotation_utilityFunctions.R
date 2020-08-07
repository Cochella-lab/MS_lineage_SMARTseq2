##########################################################################
##########################################################################
# Project: single cell RNA-seq analysis
# Script purpose: clustering and annotation utility functions
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Aug  7 12:46:05 2020
##########################################################################
##########################################################################
########################################################
########################################################
# Section :
# test code and not used anymore 
########################################################
########################################################
test.clustering.method.scran = function(sce){
  if(Test.scran.clustering){
    ##########################################
    # test grapha method from scran
    ##########################################
    snn.gr <- buildSNNGraph(sce, use.dimred="MNN")
    clusters <- igraph::cluster_walktrap(snn.gr)
    table(clusters$membership, sce$Batch)
    
    sce$cluster <- factor(clusters$membership)
    plotTSNE(sce, colour_by="cluster") + ggtitle("scran -- graph based clustering")
    plotUMAP(sce, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
      fontsize + ggtitle("scran -- graph based clustering")
    
    ##########################################
    # test hierachy clustering from scran
    ##########################################
    pcs <- reducedDim(sce, "MNN")
    my.dist <- dist(pcs)
    my.tree <- hclust(my.dist, method="ward.D2")
    
    #hist(my.tree)
    library(dynamicTreeCut)
    my.clusters <- unname(cutreeDynamic(my.tree, cutHeight = 3, 
                                        distM=as.matrix(my.dist), 
                                        minClusterSize=5, verbose=0))
    
    table(my.clusters, sce$Batch)
    
    sce$cluster <- factor(my.clusters)
    
    #plotTSNE(sce, colour_by="cluster", size_by = "total_features_by_counts") + fontsize + ggtitle("scran -- hcluster")
    plotUMAP(sce, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
      fontsize + ggtitle("scran -- hcluster")
    
    #plotDiffusionMap(sce, colour_by="cluster", size_by = "total_features_by_counts") + fontsize
    library(cluster)
    clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
    sil <- silhouette(my.clusters, dist = my.dist)
    sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
    sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
    plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
         border=sil.cols, col=sil.cols, do.col.sort=FALSE) 
    
  }
}

##########################################
# test seurat for clustering 
# original codes in https://satijalab.org/seurat/pbmc3k_tutorial.html
##########################################
Test.Seurat.workflow = FALSE
if(Test.Seurat.workflow){
  library(Seurat)
  library(dplyr)
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))
  
  #pbmc = as.seurat(from = sce)
  #pbmc = Convert(from = sce, to = 'seurat', raw.data.slot = "counts",data.slot = "logcounts")
  pbmc <- CreateSeuratObject(raw.data = counts(sce), min.cells = 0, min.genes = 200, 
                             project = "seurat_test")
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
  
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.5, y.cutoff = 0.5)
  
  length(x = pbmc@var.genes)
  
  pbmc <- ScaleData(object = pbmc)
  
  pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
                 genes.print = 5)
  
  PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
  VizPCA(object = pbmc, pcs.use = 1:2)
  
  PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
  
  pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
  
  PCHeatmap(object = pbmc, pc.use = 1:4, cells.use = 80, do.balanced = TRUE, label.columns = FALSE)
  
  pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
  JackStrawPlot(object = pbmc, PCs = 1:12)
  
  PCElbowPlot(object = pbmc)
  
  pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                       resolution = 2, print.output = 0, save.SNN = TRUE)
  
  PrintFindClustersParams(object = pbmc)
  #PrintCalcParams(pbmc, calculation = "FindClusters")
  
  pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE, perplexity=10, eta=2000)
  TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 3.0)
  
  ## test phate plot
  pbmc_phate <- RunPHATE(object = pbmc, gamma=1, npca = 20, k=5, seed.use = 10, plot.optimal.t = TRUE, t=14)
  # Plot results
  DimPlot(object = pbmc_phate, reduction.use = 'phate', pt.size = 2.0)
  #DimPlot(object = pbmc_phate, reduction.use = 'pca')
}

##################################################
# test Hamberg's single-cell RNA seq analysis, clustering part
# original code https://hemberg-lab.github.io/scRNA.seq.course/biological-analysis.html
##################################################
### test SC3
TEST.Hamberg.workflow.clustering = FALSE
if(TEST.Hamberg.workflow.clustering)
{
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata')) 
  library(SC3)
  #library()
  expr_matrix =  exp(logcounts(sce))
  Brennecke_HVG <- BrenneckeGetVariableGenes(
    expr_mat = expr_matrix,
    spikes = NA,
    fdr = 0.2,
    minBiolDisp = 0.2
  )
  #HVG_genes <- Brennecke_HVG
  
  ## another method to identify by Kolodziejczyk AA, Kim JK, Tsang JCH et al. (2015)
  assay(sce, "normcounts") <- exp(logcounts(sce))
  means <- rowMeans(normcounts(sce))
  cv2 <- apply(normcounts(sce), 1, var)/means^2
  dm.stat <- DM(means, cv2)
  #head(dm.stat)
  
  DM_HVG = names(dm.stat)[which(dm.stat>0.3)]
  #DM_HVG = DM_HVG[which()]
  #dev.off()
  
  sce.HVG.Brenneck = sce[rownames(sce)%in%Brennecke_HVG, ]
  sce.HVG.DM = sce[rownames(sce)%in%DM_HVG, ]
  
  save(sce, sce.HVG.Brenneck, sce.HVG.DM, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  #HVG_genes <- Brennecke_HVG$Gene
  #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  
  sce.sels = sce.HVG.Brenneck
  #sce.sels = sce.HVG.DM
  sce = sce.sels
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  # remove features with duplicated names
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  
  plotPCA(
    sce,
    run_args = list(exprs_values = "logcounts"),
    colour_by = "total_counts",
    size_by = "total_features_by_counts"
  )
  
  ### test SC3 clustering method
  sce = sc3_estimate_k(sce)
  metadata(sce)$sc3$k_estimation
  
  sce <- sc3(sce.sels, gene_filter = FALSE, ks = 2:30, biology = TRUE, n_cores = 6)
  #rowData(sce)$feature_symbol
  #sce <- sc3(sce, ks = 2, gene_filter = TRUE, biology = TRUE)
  #deng <- plotTSNE(deng, rand_seed = 1, return_SCE = TRUE)
  
  col_data <- colData(sce)
  head(col_data[ , grep("sc3_", colnames(col_data))])
  plotPCA(
    sce, 
    colour_by = "sc3_10_clusters", 
    size_by = "sc3_10_log2_outlier_score"
  )
  
  plotPCA(
    sce, 
    colour_by = "sc3_3_clusters",
    shape_by = "celltypes",
    size_by = "total_features_by_counts"
  )
  
  row_data <- rowData(sce)
  head(row_data[ , grep("sc3_", colnames(row_data))])
  
  plotFeatureData(
    sce, 
    aes(
      x = sc3_3_markers_clusts, 
      y = sc3_3_markers_auroc, 
      colour = sc3_3_markers_padj
    )
  )
  
  set.seed(1)
  plotTSNE(sce, colour_by="sc3_6_clusters", size_by = "total_features_by_counts",
           run_args = list(perplexity=20)) 
  
  #sc3_plot_consensus(sce, k = 3)
  sc3_plot_consensus(
    sce, k = 3, 
    show_pdata = c(
      "celltypes", 
      "sc3_3_clusters", 
      "log10_total_features_by_counts",
      "log10_total_counts",
      
      "sc3_3_log2_outlier_score"
    )
  )
  
  sc3_plot_cluster_stability(sce, k = 6)
  
  sc3_plot_de_genes(sce, k = 3, p.val = 0.2)
  
  sc3_plot_markers(sce, k = 3)
  
}

##########################################
# # configuration of umap
# https://github.com/cran/umap/blob/master/R/umap.R
##########################################
config.umap = function(n_neighbors=15,
                       n_components=2,
                       metric="euclidean",
                       n_epochs=200,
                       input="data",
                       init="spectral",
                       min_dist=0.5,
                       set_op_mix_ratio=1,
                       local_connectivity=1,
                       bandwidth=1.0,
                       alpha=1,
                       gamma=1.0,
                       negative_sample_rate=5,
                       a=NA,
                       b=NA,
                       spread=1,
                       random_state=NA,
                       transform_state=NA,
                       knn_repeats=1,
                       verbose=TRUE,
                       umap_learn_args = NA)
{
  umap.defaults = list(
    n_neighbors = n_neighbors,
    n_components = n_components,
    metric = metric,
    n_epochs = n_epochs,
    input = input,
    init = init,
    min_dist = min_dist,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    bandwidth = bandwidth,
    alpha = alpha,
    gamma = gamma,
    negative_sample_rate = negative_sample_rate,
    a = a,
    b = b,
    spread = spread,
    random_state = random_state,
    transform_state = transform_state,
    knn_repeats = knn_repeats,
    verbose = verbose,
    umap_learn_args = umap_learn_args
  )
  
  class(umap.defaults) = "umap.config"
  
  return(umap.defaults)
  
}

