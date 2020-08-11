##########################################################################
##########################################################################
# Project: Aleks' scRNA-seq MS project
# Script purpose: functions for clustering  
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Nov 22 15:38:57 2019
##########################################################################
##########################################################################

########################################################
########################################################
# Section : Aleks' code for clustering and cell-identities assignment using correlation 
# 
########################################################
########################################################
cell.identities.assignment.correlation.based.aleks = function()
{
  dir.processed.data = "/Users/jiwang/workspace/imp/scRNAseq_MS_lineage_dev/results_aleks/results/all_batches_202005/Rdata/"
  load(paste0(dir.processed.data, "all_batches_QCed_cells_genes_filtered_timingEst_Normed_bc_Seurat.Rdata"))
  
  library(dplyr)
  #library(dplyr)
  #setwd("/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/results_aleks/results/all_batches_202005/Rdata/")
  #load("./all_batches_QCed_cells_genes_filtered_timingEst_Normed_bc_Seurat.Rdata")
  
  ms <- FindNeighbors(object = ms, reduction = "mnn", k.param = 20, dims = 1:20)
  ms <- FindClusters(ms, resolution = 12, algorithm = 3)
  
  DimPlot(object = ms, cells = colnames(ms), group.by = 'ident', label = TRUE, pt.size = 4)
  all_ms.markers <- FindAllMarkers(ms, min.pct = 0.25, logfc.threshold = 0.25)
  all_ms.top10 <- all_ms.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  
  lin_sc_expr_190602 <- readRDS("data/lin_sc_expr_190602.rds")
  MS_names <- lin_sc_expr_190602[,grepl(pattern = "MS.", x = colnames(lin_sc_expr_190602))]
  #lin_sc_expr_190602[,MS_names]
  #load(file = paste0(dir.processed.data, 'ms_correlation_idents.rds'))
  
  combined_MS_names <- data.frame(matrix(NA, dim(MS_names)[1], 1))
  used_cells <- c()
  
  for (cell in colnames(MS_names)){
    if (cell %in% used_cells){next}
    indx <- (MS_names[1:100,] == MS_names[1:100,cell])[1,]
    indx_name <- colnames(MS_names)[indx]
    tmp_df <- data.frame(MS_names[,indx_name[1]])
    colnames(tmp_df) <- paste(indx_name, collapse = "/")
    combined_MS_names <- cbind(combined_MS_names,tmp_df)
    used_cells <- c(used_cells, indx_name)
  }
  
  JM_data <- combined_MS_names[,-1]
  my_MS <- as.matrix(ms[["SCT"]]@data)
  
  indexes <- match(all_ms.top10$gene, rownames(my_MS))
  my_MS <- my_MS[indexes,]
  indexes <- match(rownames(my_MS), rownames(JM_data))
  mathed_gene_names <- rownames(JM_data)[indexes[!is.na(indexes)]]
  JM_data <- JM_data[mathed_gene_names,]
  my_MS <- my_MS[mathed_gene_names,]
  
  cor_vec <- c()
  
  combined_df <- data.frame(row.names = 1:(dim(JM_data)[2]))
  combined_df_corv <- data.frame(row.names = 1:(dim(JM_data)[2]))
  
  library(progress)
  pb <- progress_bar$new(
    format = " progress [:bar] :percent eta: :eta",
    total = length(colnames(my_MS)))
  
  for (a in colnames(my_MS)){
    for (i in 1:dim(JM_data)[2]){
      cor_vec[i] <- cor(my_MS[,a], JM_data[,i])
    }
    combined_df[,a] <- colnames(JM_data)[order(cor_vec, decreasing = T)]
    combined_df_corv[,a] <- cor_vec[order(cor_vec, decreasing = T)]
    pb$tick()
  }
  
  
  View(combined_df_corv[1:10,1:40])
  cell.names <- colnames(combined_df)
  
  Idents(ms, cells = cell.names) <- combined_df[1,]
  
  pb <- progress_bar$new(
    format = " progress [:bar] :percent eta: :eta",
    total = length(levels(Idents(ms))))
  
  pdfname = paste0(paste0(resDir, "/idents_correlation_aleks.pdf"))
  pdf(pdfname, width=16, height = 16)
  
  for (idntt in levels(Idents(ms))){
    
    #pdfname = paste0(paste0("~/Documents/plots/", which(levels(Idents(test_ms)) == idntt), ".pdf"))
    #pdf(pdfname, width=8, height = 8)
    
    cells_to_show <- list(c(WhichCells(ms, idents = idntt)))
    names(cells_to_show) <- idntt
    print(DimPlot(ms, reduction = "umap", label = F, cells.highlight = cells_to_show))
    #dev.off()
    pb$tick()
    
  }
  
  dev.off()
  
  ms_correlation_idents <- ms
  
  saveRDS(ms_correlation_idents, file = paste0(RdataDir,  'ms_correlation_idents.rds'))
  saveRDS(combined_df, file = paste0(RdataDir,  'cell.states_assignment_to_reference.JM_correlation_idents.rds'))
  saveRDS(combined_df_corv, file = paste0(RdataDir,  'cell.states.correlation_assignment_to_reference.JM_correlation_idents.rds'))
  
}


########################################################
########################################################
# Section : test umap parameters for best visulization
# 
########################################################
########################################################
test.umap.params = function(seurat.obj, pdfname = paste0(resDir, '/umap_params_test.pdf'), 
                            nb.pcs.sampling = seq(20, 50, by = 10),
                            n.neighbors.sampling = seq(10, 50, by = 10),
                            min.dist.sampling = c(0.01, 0.05, seq(0.1, 0.5, by = 0.1))
)
{
  #pdfname = paste0(resDir, "/cluster_annotations/test_umap_params_for_seurat.obj.pdf")
  pdf(pdfname, width=12, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #nb.pcs = ncol(seurat.obj[['pca']]);
  #n.neighbors = 50; min.dist = 0.01;
  
  for(nb.pcs in unique(c(nb.pcs.sampling, ncol(seurat.obj[['pca']]))))
  {
    if(nb.pcs <= ncol(seurat.obj[['pca']])){
      for(n.neighbors in n.neighbors.sampling)
      {
        for(min.dist in min.dist.sampling)
        {
          cat('nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, ', min.dist - ', min.dist, '\n')
          seurat.obj <- RunUMAP(object = seurat.obj, reduction = 'pca', 
                                dims = 1:nb.pcs, 
                                n.neighbors = n.neighbors, min.dist = min.dist)
          p1 =  DimPlot(seurat.obj, reduction = 'umap', label = TRUE, pt.size = 1, label.size = 6) + 
            NoLegend() + 
            ggtitle(paste0('nb.pcs - ', nb.pcs, '; n.neighbors - ', n.neighbors, ', min.dist - ', min.dist))
          plot(p1)
        }
      }
    }
  }
  
  dev.off()
  
}

########################################################
########################################################
# Section : label transferring from Murray dataset
# 1) seurat, scmap, svm, RF
# 2) utility functions
########################################################
########################################################
##########################################
# import Murray scRNA data and select the cell identities of interest (MS and some early stage cells)
##########################################
process.import.Murray.scRNA = function()
{
  library(VisCello.celegans)
  eset = readRDS(file = paste0('data/Parker_et_al_dataSet_afterFiltering_89701cell.rds'))
  pmeda = data.frame(pData(eset))
  
  ## select the cells for MS lineages
  kk = grep('^MS', pmeda$lineage)
  kk1 = which(pmeda$lineage == '28_cell_or_earlier'| pmeda$lineage == 'ABaxx'| pmeda$lineage == 'Cx'|
                pmeda$lineage == 'Dx'|pmeda$lineage == 'Dxa'|pmeda$lineage == 'Exx')
  
  kk = unique(c(kk, kk1))
  cat('nb of cell in reference -- ', length(kk), '\n')
  cat('nb of cell states in reference -- ', length(unique(pmeda$lineage[kk])), '\n')
  
  ee = CreateSeuratObject(counts = eset@assayData$exprs[,kk], assay = 'RNA', meta.data = pmeda[kk, ])
  ee@assays$RNA@data = eset@assayData$norm_exprs[,kk]
  
  return(ee)
  
}

reference.based.cluster.annotation = function(seurat.obj, redefine.clusters = TRUE,
                                              predict.unassignedCells = FALSE, threshold.svm = 0.5, threshold.rf = 0.5)
{
  # seurat.obj = ms; redefine.clusters = TRUE; predict.unassignedCells = FALSE;
  
  ##########################################
  # step 0) refine the clustering from seurat
  ##########################################
  if(redefine.clusters){
    nfeatures = 3000
    seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = nfeatures)
    
    seurat.obj = ScaleData(seurat.obj, features = rownames(seurat.obj))
    
    seurat.obj <- RunPCA(object = seurat.obj, features = VariableFeatures(seurat.obj), verbose = FALSE, weight.by.var = TRUE)
    ElbowPlot(seurat.obj, ndims = 50)
    
    seurat.obj <- FindNeighbors(object = seurat.obj, reduction = "pca", k.param = 10, dims = 1:20)
    
    #library(leiden)
    #cluster_leiden <- leiden(seurat.obj@graphs$RNA_snn)
    #library(reticulate)
    #use_python("/usr/local/bin/python", required = TRUE)
    #system('pip install leidenalg igraph')
    #system('python --version')
    #system('which python')
    
    seurat.obj <- FindClusters(seurat.obj, resolution = 3, algorithm = 3)
    
    cat(length(unique(seurat.obj$seurat_clusters)), 'clusters found \n')
    
    #nb.pcs = 30; n.neighbors = 40; min.dist = 0.25;
    #seurat.obj <- RunUMAP(object = seurat.obj, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
    
    DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
            na.value = "gray") + 
      ggtitle(paste0("Seurat_clustering_SLM_resolution_3_3000variableFeatures_20pca_k10")) +
      scale_colour_hue(drop = FALSE) + 
      NoLegend()
    
  }
  
  ##########################################
  # step 1): project aleks cells to the reference using scmap and seurat
  # transfer labels with stringent threshold
  ##########################################
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  
  ## import and process Murray data
  ee = process.import.Murray.scRNA()
  
  ## tranfer Murray labels with scmap
  seurat.obj = scmap.transfer.labels.from.Murray.scRNA(seurat.obj, ee)
  
  ## transfer Murray labels with seurat
  seurat.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA(seurat.obj, ee)
  
  rdsfile.saved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat.rds')
  saveRDS(seurat.obj, file = rdsfile.saved)
  
  ##########################################
  # step 2): predict unassigned cells using assgined cells with rf and svm
  ##########################################
  if(predict.unassignedCells){
    seurat.obj = prediction.unassinged.cells.rf.svm(seurat.obj) 
  }
  
  return(seurat.obj)
  
}


########################################################
########################################################
# Section : here manual annotate BWM lineages using various information:
# 1) clusters (splitting and merging if necessay)
# 2) predicted labels from seurat and scmap 
# 3) cell size info and estimated timing
# 4) cluster connection by PAGA or VarID 
# 5) RNA velocity (not sure ...)
########################################################
########################################################
annotate.clusters.using.predicted.id = function(seurat.obj, redefine.clusters = FALSE)
{
  ## current clusters were define using 3000 variable genes and resolution = 12
  ## those clusters are relatively small owe to the high value of resolution
  
  ## compare scmap and seurat reference-based annotation
  ## decide which methods and how many features to use or how to integrate different methods
  seurat.obj = readRDS(seurat.obj, file = rdsfile.saved)
  
  compare.reference.based.annotation.scmap.seurat(seurat.obj)
  
  ##########################################
  # step 3):  refine clusters (split or merge)
  ##########################################
  seurat.obj = annotate.clusters.using.predicted.id(seurat.obj)
  
  ##########################################
  # step 4): focus short list of cell identities and manual annotate with other information
  ##########################################
  
  
  # seurat.obj = readRDS(file = paste0(RdataDir, 'processed_5.4k.cells_scran.normalized_scmapProjection.rds'))
  
  seurat.obj$predicted.id.scmap[which(is.na(seurat.obj$predicted.id.scmap))] = 'unassigned'
  
  seurat.obj$seurat_clusters = seurat.obj$RNA_snn_res.5
  seurat.obj$cluster.idents = seurat.obj$seurat_clusters
  cluster.ids = unique(seurat.obj$seurat_clusters)
  cluster.ids = cluster.ids[order(cluster.ids)]
  
  pdfname = paste0(resDir, "/annotate_clusters_using_predicted.ids_split_merge_test1.pdf")
  pdf(pdfname, width=16, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(n in 1:10)
  #for(n in 1:length(cluster.ids))
  {
    # n = 6
    cat('cluster ', as.character(cluster.ids[n]), '\n')
    
    sels = which(seurat.obj$seurat_clusters == cluster.ids[n])
    
    p1 = DimPlot(seurat.obj, 
            cells.highlight = colnames(seurat.obj)[sels],
            group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 4,
            sizes.highlight = 2,
            na.value = "gray") + 
      ggtitle(paste0("cluster ", cluster.ids[n])) +
      NoLegend()
    
    p2 = DimPlot(seurat.obj, 
            cells = colnames(seurat.obj)[sels],
            group.by = "predicted.id.scmap", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3, label.size = 4,
            na.value = "gray") + 
      ggtitle(paste0("cluster ", cluster.ids[n])) 
    
    plot(CombinePlots(list(p1, p2), ncol = 2))
    
    pred.ids = seurat.obj$predicted.id.scmap[sels]
    pred.ids[which(is.na(pred.ids))] = 'unassigned'
    
    par(mfrow = c(1, 1))
    par(mar=c(5,8,4,2)) # increase y-axis margin.
    barplot(table(pred.ids)/length(pred.ids), las = 2, horiz = TRUE, xlim =c(0, 1), main = paste0("cluster ", cluster.ids[n]))
    
    par(mfrow = c(1, 2))
    hist(seurat.obj$BSC_log2[sels], breaks = 20)
    hist(seurat.obj$FSC_log2[sels], breaks = 20)
    
    ## refine clusters by splitting or outlier detection (k mean cluster and RaceID)
    clustering.splitting.kmean.outlier.detection(seurat.obj, sels, redefine.gene.to.use = TRUE, nfeatures = 200)
    cat('---------\n')
    
  }
  
  dev.off()
  
}

clustering.splitting.kmean.outlier.detection = function(seurat.obj, sels, redefine.gene.to.use = FALSE, nfeatures = 3000,
                                                        cluster.split.method = 'dbscan')
{
  # n = 6;  sels = which(seurat.obj$seurat_clusters == cluster.ids[n])
  subobj = subset(seurat.obj, cells = sels)
  
  p1 = DimPlot(subobj,
          group.by = "predicted.id.scmap", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3, label.size = 4,
          na.value = "gray") +
    ggtitle(paste0("cluster ", cluster.ids[n]))
  
  # plot(p1)
  
  if(redefine.gene.to.use){
    subobj <- FindVariableFeatures(subobj, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    # subobj = ScaleData(subobj, features = rownames(subobj))
    # subobj <- RunPCA(object = subobj, features = VariableFeatures(subobj), verbose = FALSE)
    # ElbowPlot(subobj, ndims = 50)
    # subobj <- FindNeighbors(object = subobj, reduction = "pca", k.param = 20, dims = 1:20)
    # subobj <- FindClusters(subobj, resolution = 1, algorithm = 3)
    ss = subset(subobj, features = VariableFeatures(subobj))
    sc = SCseq(ss@assays$RNA@counts)
    
  }else{
    sc = SCseq(subobj@assays$RNA@counts)
  }
  
  if(cluster.split.method == 'kmean')
  {
    library(RaceID) # refer to the vignett https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html
    library(Matrix)
    
    # test example from viggnette
    #sc <- SCseq(intestinalData)
    sc <- filterdata(sc, mintotal=2000, minexpr = 10, minnumber = 2)
    fdata <- getfdata(sc)
    #cat('nb of genes after filtering :', fdata@Dim[1], '\n')
    #cat('nb of cells after filtering :', fdata@Dim[2], '\n')
    
    #sc <- compdist(sc,metric="pearson", FSelect = FALSE)
    cat('use pca to calculate Pearson correlation and then distance and then k-mean\n')
    subobj.pca = subobj@reductions$pca@cell.embeddings[, c(1:30)]
    mat.dist = 1- cor(t(subobj.pca))
    sc@distances = mat.dist
    sc <- clustexp(sc, FUNcluster = 'kmedoids', verbose = FALSE)
    
    par(mfrow = c(1, 1))
    plotsaturation(sc,disp=FALSE)
    plotsaturation(sc,disp=TRUE)
    #plotjaccard(sc)
    
    nb.subcluster = sc@cluster$clb$nc
    cat('optimal subclusters found :', nb.subcluster, '\n')
    
    sc <- clustexp(sc,cln=nb.subcluster, sat=FALSE, verbose = FALSE)
    
    #sc <- findoutliers(sc, outminc = 10, outlg = 10)
    #plotbackground(sc)
    
    #plotoutlierprobs(sc)
    #clustheatmap(sc)
    
    #sc <- compumap(sc)
    #sc <- compfr(sc,knn=10)
    #plotmap(sc,um=TRUE)
    #plotmap(sc,fr=TRUE)
    
    subobj$seurat_clusters_split = sc@cluster$kpart
    
  }
  
  if(cluster.split.method == 'dbscan'){
    library('dbscan')
    library('lsa')
    #library("fpc")
    
    #sc <- compdist(sc,metric="pearson", FSelect = FALSE)
    cat('use pca to calculate Pearson correlation and then distance and then dbscan clustering \n')
    subobj.pca = subobj@reductions$pca@cell.embeddings[, c(1:20)]
    #mat.dist = dist(1- cor(t(subobj.pca)))
    mat.dist = dist(lsa::cosine(t(subobj.pca)))
    #mat.dist = dist(subobj.pca, method = 'euclidean')
    
    dbscan::kNNdistplot(mat.dist, k = 4)
    abline(h = 1, lty = 2)
    
    res.db <- dbscan(mat.dist, eps = 0.5, minPts = 4)
    cat('nb of clusters found by dbscan -- ', length(unique(res.db$cluster)), '\n')
    
    subobj$seurat_clusters_split = res.db$cluster
    
  }
  
  p2 = DimPlot(subobj,
          group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3, label.size = 4,
          na.value = "gray") +
    ggtitle(paste0("cluster ", cluster.ids[n]))
  
  plot(CombinePlots(list(p1, p2)))
  
  par(mfrow=c(1, 1))
  counts <- table(subobj$predicted.id.scmap, subobj$seurat_clusters_split)
  barplot(counts, main="composition of subclusters ",
          xlab="subcluster index", col=c(1:nrow(counts)),
          legend = rownames(counts))
  
  return(subobj)
  
}


########################################################
# Section : Clustering section by integrating various informations: 
# gene expression, fac info, estimated timing and velocity 
########################################################
########################################################
cluster.by.integration.gene.expression.fac.timingEst.velocity = function()
{
  library(scater)
  library(scran)
  library(scRNA.seq.funcs)
  library(matrixStats)
  #library(M3Drop)
  library(RColorBrewer)
  library(SingleCellExperiment)
  
  load(file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_bcMNN.Rdata')) 
  
  Seurat.clustering = TRUE
  Test.scran.clustering = FALSE
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
  
  ## double check the umap and tsnet visulization and 
  set.seed(100)
  sce <- runUMAP(sce, use_dimred="MNN", n_dimred = 15)
  p1 = plotUMAP(sce, ncomponents = 2, colour_by="batches", size_by = "FSC_log2", point_size= 0.01) + ggtitle("Corrected-umap") 
  plot(p1)
  
  sce = runTSNE(sce, use_dimred="MNN", perplexity = 20, n_dimred = 15)
  p2 = plotTSNE(sce, ncomponents = 2, colour_by="batches", size_by = "FSC_log2", point_size= 0.01) + ggtitle("Corrected-tSNE") 
  plot(p2)
  
  p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  multiplot(p1, p2, cols = 2)
  
  
  p1 = plotTSNE(sce, colour_by="pha-4", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  p2 = plotTSNE(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  multiplot(p1, p2, cols = 2)
  
  ##########################################
  # (test clustering method) but Seurat method is used at the end 
  # https://master.bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/work-1-reads.html/
  # DE analysis is tested with scran but will be considered to be replaced by MAST or Seurat
  # DE analysis (or marker gene finding) following the cluster analysis
  ##########################################
  if(Seurat.clustering)
  { 
    require(Seurat)
    library(cowplot)
    
    check.individualExample.geneMarker = FALSE
    ntops = 3 # nb of top gene markers
    resolutions = c(0.4, 0.6, 0.8, seq(1.0, 4.0, by = 0.5))
    
    for(rr in resolutions){
      
      rr = 1.2
      cat("--- resolution is :", rr, "---\n")
      
      pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_batchCorrected_clustering.Seurat_geneMarkers.scran_resolution_", 
                       rr, ".pdf")
      pdf(pdfname, width=22, height = 18)
      par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
      
      ##########################################
      # test graph-based Louvain algorithm 
      ##########################################
      pbmc = as.Seurat(sce)
      #pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
      #pbmc <- ScaleData(pbmc, features = rownames(pbmc))
      # pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
      pbmc = FindNeighbors(object = pbmc, reduction = "MNN", k.param = 20, dims = 1:20)
      pbmc = FindClusters(pbmc, resolution = rr, algorithm = 3)
      
      # note that you can set `label = TRUE` or use the LabelClusters function to help label
      # individual clusters
      #pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10, reduction = "MNN", reduction.key = "umap", n.neighbors = 15, repulsion.strength = 2)
      #DimPlot(pbmc, reduction = "umap")
      #FeaturePlot(pbmc, features = c("pha-4", "hnd-1"), reduction = "umap")
      # FeaturePlot(pbmc, features = c("pha-4", "hnd-1"), reduction = "UMAP")
      sce$cluster_seurat <- factor(pbmc@active.ident)
      sce$cluster <- factor(pbmc@active.ident)
      
      plotUMAP(sce, colour_by="cluster", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      
      plotTSNE(sce, colour_by="cluster", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      
      p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      multiplot(p1, p2, cols = 2)
      
      p1 = plotTSNE(sce, colour_by="pha-4", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      p2 = plotTSNE(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      multiplot(p1, p2, cols = 2)
      
      my.clusters = as.numeric(as.character(sce$cluster_seurat))
      cat(table(my.clusters), "\n")
      
      ## run the find markers and then collect markers for each clusters (scran) 
      # https://bioconductor.org/packages/3.10/workflows/vignettes/simpleSingleCell/inst/doc/de.html#blocking-on-uninteresting-factors-of-variation
      #design <- model.matrix( ~ sce$batches)
      #design <- design[,-1,drop=FALSE]
      markers <- findMarkers(sce, my.clusters, block=sce$batches, direction="up")
      top.markers = c()
      
      for(n in unique(my.clusters)){
        #n = 0
        marker.set <- markers[[as.character(n)]]
        #marker.set <- markers[["1"]]
        #head(marker.set, 5)
        top.markers <- c(top.markers, rownames(marker.set)[marker.set$Top <= ntops])  
      }
      
      top.markers = unique(top.markers)
      cat(length(top.markers), " gene markers found by scran \n")
      
      plotHeatmap(sce, features=top.markers,
                  columns=order(sce$cluster_seurat), 
                  colour_columns_by=c("cluster",  "FSC_log2", "batches"),
                  cluster_cols=FALSE, show_colnames = FALSE,
                  center=TRUE, symmetric=TRUE, zlim=c(-5, 5))
      
      if(check.individualExample.geneMarker){
        for(n in 1:length(top.markers)) {
          xx = plotTSNE(sce, colour_by = top.markers[n]) 
          plot(xx)
        }
      }
      
      dev.off()
      
    }
    
  }
  
  ##########################################
  # here select subset of whole dataset to redo clustering
  # and marker gene finding
  # this is the beginning of iterative clustering, 
  # the gene marker finding will be a criterion to stop iteration
  ##########################################
  Select.Early.timePoints = FALSE
  if(Select.Early.timePoints){
    
    pbmc = as.Seurat(sce)
    #Seurat::DimPlot(pbmc, dims = c(1, 2), reduction = "MNN")
    pbmc = FindNeighbors(object = pbmc, reduction = "MNN", k.param = 10, dims = 1:10)
    pbmc = FindClusters(pbmc, resolution = 2, algorithm = 3)
    sce$cluster_seurat <- factor(pbmc@active.ident)
    sce$cluster <- factor(pbmc@active.ident)
    
    xx = table(sce$cluster,sce$Batch)
    #colnames(xx) = sce$Batch
    cluster4early = rownames(xx)[which(xx[, 1]>=5|xx[,2]>=5)]
    
    mm = match(sce$cluster, factor(cluster4early))
    
    sels = which(!is.na(mm))
    
    cat("estimated cell in early stage -- ", length(cluster4early), 
        "clusters with", length(sels),  "cells\n")
    
    sce.sel = sce[, sels ]
    # here rerun the clustering for clusters or stages
    pcs <- reducedDim(sce.sel, "MNN")
    my.dist <- dist(pcs[, c(1:20)])
    my.tree <- hclust(my.dist, method="ward.D2")
    
    #hist(my.tree)
    library(dynamicTreeCut)
    my.clusters <- unname(cutreeDynamic(my.tree,
                                        distM=as.matrix(my.dist), 
                                        minClusterSize=5, verbose=0))
    
    table(my.clusters, sce.sel$Batch)
    
    sce.sel$cluster <- factor(my.clusters)
    
    set.seed(100)
    sce.sel <- runTSNE(sce.sel,  use_dimred = "MNN", perplexity = 20, n_dimred = 20)
    
    set.seed(100)
    sce.sel = runUMAP(sce.sel, use_dimred="MNN", perplexity = 20, n_dimred = 20)
    
    plotTSNE(sce.sel, colour_by="cluster", size_by = "total_features_by_counts") + 
      fontsize + ggtitle("scran -- hcluster")
    plotUMAP(sce.sel, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
      fontsize + ggtitle("scran -- hcluster")
    
    design <- model.matrix( ~ sce.sel$Batch)
    design <- design[,-1,drop=FALSE]
    
    markers <- findMarkers(sce.sel, sce.sel$cluster, design=design, direction = 'any')
    
    ntops = 5;
    top.markers = c()
    
    for(n in unique(my.clusters)){
      #n = 0
      marker.set <- markers[[as.character(n)]]
      #marker.set <- markers[["1"]]
      #head(marker.set, 5)
      top.markers <- c(top.markers, rownames(marker.set)[marker.set$Top <= ntops])  
    }
    top.markers = unique(top.markers)
    
    pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_batchCorrected_clustering_markerGenes_earlyTimepoint_tSNEexample.pdf")
    pdf(pdfname, width=12, height = 10)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    plotTSNE(sce.sel, colour_by="cluster", size_by = "total_features_by_counts") + 
      fontsize + ggtitle("scran -- hcluster (tSNE)")
    plotUMAP(sce.sel, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
      fontsize + ggtitle("scran -- hcluster (UMAP)")
    
    plotHeatmap(sce.sel, features=top.markers,
                columns=order(sce.sel$cluster), 
                colour_columns_by=c("cluster"),
                cluster_cols=FALSE, show_colnames = FALSE,
                center=TRUE, symmetric=TRUE, zlim = c(-5, 5))
    
    for(n in 1:length(top.markers)) {
      xx = plotTSNE(sce.sel, colour_by = top.markers[n]) 
      plot(xx)
    }
    
    dev.off()
    
  }
}




