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
  
  # transfer Tintori et al. labels 
  seurat.obj = scmap.transfer.labels.from.Tintor.scRNA(seurat.obj)
  
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
# 3.5) marker genes
# 4) cluster connection by PAGA or VarID 
# 5) RNA velocity (not sure ...)
########################################################
########################################################
##########################################
# here compare scmap and seurat, two reference-based cluster annotation
# the clusters were also given here in seurat.obj$seurat_clusters
##########################################
overview.and.compare.predicted.labels = function(seurat.obj)
{
  # seurat.obj = ms
  library("pheatmap")
  library("RColorBrewer")
  library(grid)
  
  # inital clusters
  p0 = DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, 
               label.size = 6,
               na.value = "gray") + 
    ggtitle(paste0("clusters ")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  plot(p0)
  
  ##########################################
  # check the predicted labels and cluster-label mapping without filtering
  ##########################################
  # here use seurat prediction as example
  seurat.obj$predicted.ids = seurat.obj$scmap.pred.id.500
  seurat.obj$predicted.scores = seurat.obj$scmap.corr.500
  threshold = 0.7
  
  seurat.obj$predicted.ids.filtered = seurat.obj$predicted.ids
  seurat.obj$predicted.ids.filtered[which(seurat.obj$predicted.ids.filtered == 'unassigned')] = NA
  seurat.obj$predicted.ids.filtered[which(seurat.obj$predicted.scores < threshold)] = NA
  
  
  pdfname = paste0(resDir, "/Overview_predictedLabels_seuratClusters_mapping_scmap.pdf")
  pdf(pdfname, width=16, height = 12)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  p1 = DimPlot(seurat.obj, group.by = "predicted.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, 
               label.size = 4,
               na.value = "gray") + 
    ggtitle(paste0("projection into Murray data with scmap (nfeature = 500)")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  plot(p1)
  
  par(mfrow = c(1, 1))
  hist(seurat.obj$predicted.scores, breaks = 100)
  abline(v= threshold, col = 'red')
  
  p2 = DimPlot(seurat.obj, group.by = "predicted.ids.filtered", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, 
               label.size = 4,
               na.value = "gray") + 
    ggtitle(paste0("filtered by threshold")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  plot(p2)
  
  #seurat.obj$predicted.ids[which(seurat.obj$predicted.scores < threshold)] = 'unassigned'
  # compositions of predicted labels for all clusters 
  counts <- table(seurat.obj$predicted.ids, seurat.obj$seurat_clusters)
  sels.bwm = c(match(c('MSx', 'MSxa', 'MSxap', 
                       'MSapaap', 'MSapaapp', 'MSappaaa', 
                       'MSpappa', 'MSpappax', 'MSppaap', 'MSppaapp', 'MSpppaaa'), rownames(counts)), 
               grep('MSxapp|MSxp', rownames(counts)))
  counts = counts[sels.bwm, ]
  counts = counts[!is.na(rownames(counts)), ]
  counts = counts[, apply(counts, 2, sum) >0]
  
  counts.norm = counts
  for(n in 1:nrow(counts)) counts.norm[n,] = counts.norm[n, ] /sum(counts.norm[n, ])
  cols = c(colorRampPalette((brewer.pal(n = 7, name="Reds")))(10))
  
  pheatmap(counts.norm, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, breaks = seq(0, 1, by = 0.1),
           cluster_cols=FALSE, main = paste0("cluster -- predicted labels mapping"), na_col = "white",
           color = cols, 
           #annotation_col = my_sample_col,
           #gaps_row = c(1:nrow(map)-1),
           fontsize_col = 10,
           height = 8,
           width = 30
  )
  
  # compositions of predicted labels for all clusters 
  #legend.text = rownames(counts)
  #cat('here \n')
  counts.norm = counts
  for(n in 1:ncol(counts)) {
    ss = sum(counts.norm[,n])
    if(ss>0) counts.norm[,n] = counts.norm[ ,n] /ss
  }
  
  ## heatmap of summarizing the mapping between cluster index and predicted cell identities
  #cols = c(colorRampPalette((brewer.pal(n = 7, name="Reds")))(10))
  pheatmap(counts.norm, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, breaks = seq(0, 1, by = 0.1),
           cluster_cols=FALSE, main = paste0("resolution -- 0.8"), na_col = "white",
           color = cols, 
           #annotation_col = my_sample_col,
           #gaps_row = c(1:nrow(map)-1),
           fontsize_col = 10,
           height = 8,
           width = 30
  )
  
  
  dev.off()
  
}

##########################################
# here is the main function for BWM cluster manual annotation
##########################################
manual.annotation.for.BWM.clusters = function(seurat.obj = ms, ids = c('MSx'))
{
  library(ggplot2)
  library(patchwork)
  library("pheatmap")
  library("RColorBrewer")
  library(grid)
  library(RaceID) # refer to the vignett https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html
  library(Matrix)
  library(lsa)
  
  ee = process.import.Murray.scRNA()
  murray.ids = unique(ee$lineage)
  bwms = unique(c('MSx', 'MSxa', 'MSxap', 
           'MSapaap', 'MSapaapp', 'MSappaaa', 
           'MSpappa', 'MSpappax', 'MSppaap', 'MSppaapp', 'MSpppaaa',
           murray.ids[grep('MSxapp|MSxp', murray.ids)]))
  
  markers = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet = 4, startRow = 8, colNames = TRUE)
  markers = markers[!is.na(match(markers$Lineage, bwms)), ]
  #write.csv(markers, file = paste0(tabDir, 'JM_marker_genes_BWM.csv'))
  
  dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
  load(file =paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval.Rdata"))
  timerGenes.pval=0.00001; timerGenes.ac=0.7
  sels.timerGenes = which(timers$ac.max > timerGenes.ac & timers$pval.box < timerGenes.pval)
  timers = rownames(timers[sels.timerGenes, -c(1:4)])
  
  seurat.obj = readRDS(file = paste0(RdataDir, 
                                     'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat.rds'))
  
  seurat.obj$manual.annot.ids = NA
  
  ##########################################
  # split the chosen clusters and manual annotate them
  ##########################################
  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxapp_4.pdf")
  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #cluster.sels = colnames(counts)
  #cluster.sels = c('29', '32',  '40', '42', '6', '14', '44', '4')
  #cluster.sels = c('29', '32', '35', '40', '42')
  cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4',  '22', '24') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4', '22', '42', '14', '6', '24', '3')
  
  sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  
  DimPlot(sub.obj, reduction = 'umap', label = TRUE)
  
  FeaturePlot(sub.obj, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67', 'pat-4'))
  
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
  ##########################################
  # check potential ids for selected clusters 
  ##########################################
  threshold = 0.7
  predicted.ids = sub.obj$scmap.pred.id.500
  selected.clusters = as.character(sub.obj$seurat_clusters)
  counts = table(predicted.ids, selected.clusters)
  
  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )
  
  predicted.ids[which(sub.obj$scmap.corr.500 <0.7)] = 'unassigned'
  counts = table(predicted.ids, selected.clusters)
  
  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )
  
  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
  
  #length(intersect(VariableFeatures(sub.obj), timers))
  #VariableFeatures(sub.obj) = setdiff(VariableFeatures(sub.obj), timers)
  cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 50)
  
  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids 
  n.neighbors = 20; min.dist = 0.2;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, n.neighbors = n.neighbors, 
                min.dist = min.dist)
  p0 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1
  
  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.6)
  
  p1  = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2)
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, 
               label.size = 5,
               na.value = "gray", combine = TRUE)
  p1 + p2
  
  manual.discovery.new.features = FALSE
  if(manual.discovery.new.features){
    Idents(sub.obj) = sub.obj$seurat_clusters_split
    markers <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
    
  }
  
  features.sels = c( 'nhr-67', 'pat-9', 'ham-1', 'tbx-8', 'cwn-2', 'unc-120', 'hnd-1', 'tab-1', 
                     'lin-39', 'rpm-1', 'hlh-1',  'Y48C3A.18', 'C05D9.4')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  p3 = VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split',
               idents = c('5', '4', '0', '1', '6'))
  
  p2 + p3
  
  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3, 
          group.by = 'seurat_clusters_split')
  
  counts = table(predicted.ids, as.character(sub.obj$seurat_clusters_split))
  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )
  
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  dev.off()
  
  ##########################################
  # update of manually annotated ids
  ##########################################
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '5')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxapp'
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '4')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSpappa'
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxappp'
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxappp/Mspappa.daugther'
  
    
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)
  
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'manual.annot.ids')
  
  #VlnPlot(seurat.obj, features = c('hnd-1', 'pha-4', 'ceh-76', 'fbxb-70'), group.by = 'manual.annot.ids')
  
  
  saveRDS(seurat.obj, 
          file = paste0(RdataDir, 
                        'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_2.rds'))
  
 
  
}

