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
source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

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
  library(dplyr)
  library(openxlsx)
  
  ee = process.import.Murray.scRNA()
  murray.ids = unique(ee$lineage)
  bwms = unique(c('MSx', 'MSxa', 'MSxap', 
           'MSapaap', 'MSapaapp', 'MSappaaa', 
           'MSpappa', 'MSpappax', 'MSppaap', 'MSppaapp', 'MSpppaaa',
           murray.ids[grep('MSxapp|MSxp', murray.ids)]))
  
  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
  #write.csv(markers, file = paste0(tabDir, 'JM_marker_genes_BWM.csv'))
  
  # dataDir.Hashimsholy = '../data/Hashimsholy_et_al'
  # load(file =paste0(dataDir.Hashimsholy, "/timer_genes_with_ac_pval.Rdata"))
  # timerGenes.pval=0.00001; timerGenes.ac=0.7
  # sels.timerGenes = which(timers$ac.max > timerGenes.ac & timers$pval.box < timerGenes.pval)
  # timers = rownames(timers[sels.timerGenes, -c(1:4)])
  
  # seurat.obj = readRDS(file = paste0(RdataDir, 
  #                                    'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat.rds'))
  # 
  # seurat.obj$manual.annot.ids = NA
  
  ########################################################
  ########################################################
  # Section : iteration 16
  # collect cells from BWM terminal cells () and mother cells () and also grand-mothers 
  # those cells are the BWM terminal and middle time points, the majority of BWM cells
  ## !!! HERE we got one of best projection from seurat for the terminal cells and transitions between middle and terminal cells 
  ## and keep the annotation as the manual.annot.ids
  # This is probably due to the 8000 variable genes used
  # the middle cells need to be refined.
  ########################################################
  ########################################################
  nb.iteration = 16
  Refine.annotated.ids = FALSE;
  
  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_', 
                    nb.iteration -1, '.rds')
  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
  RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_', 
                      nb.iteration, '.rds')
  
  seurat.obj = readRDS(file = RDSsaved)
  seurat.obj$predicted.ids.scmap = seurat.obj$scmap.pred.id.500
  seurat.obj$predicted.ids.seurat = seurat.obj$seurat.pred.id
  
  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}
  
  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  ##########################################
  # select subset of cells to annotate
  ##########################################
  cluster.index = '52'
  table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
  table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
  
  table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]
  
  xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSapaapp')])
  xx[which(xx > 0)]
  
  # select BWM terminal cells
  ##########################################
  # cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition 
  #                  '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
  #                  '25', # possible transition clusters
  #                  '24', # also transitions cells and many of them are not annotated
  #                  '2', '3', '16', '30' # all middle time points
  # )
  # cluster.sels = c('4',  '22', '24',  # middle time point for convergence branche
  #                  '3', '5', '16', '30', # middle time points for MSxp
  #                  '23', '43', '17', '41', '45', '34', '50', # well-separated clusters with clear mapped ids
  #                  '44', '28', '52', '31', '53',  # # well-separated clusters without clear mapped ids
  #                  '25', '36', '8', '39', '2', '19', '27', # BWM-1 cluster_25 transition to it and possibly cluster_28, 51 transition to BWM_2 
  #                  '24', '13', '1', '11', '33', '48', '18', '46', '15', '26' # BWM-2
  #                  )
  # 
  # cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition 
  #                  '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
  #                  '44', '31', '52', '28', '50', # we
  #                  '25', # possible transition clusters
  #                  '24' # many cells are not annotated
  #                 )
 
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))] 
  #                                   & is.na(seurat.obj$manual.annot.ids))]
  #cluster.sels = c('25', '31', '28', '52', '36', '8', '39', '2', '19', '27')
  #cluster.sels = c('25', '36', '8', '39', '2', '19', '27', '13', '1', '11', '33', '48', '18', '46', '15', '26')
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))] 
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))] 
  
  seurat.obj$BWM.cells[seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31'] = NA
  # select BWM terminal and middle time points cells
  ##########################################
  cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition 
                   '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                   '25', # possible transition clusters
                   '24', # also transitions cells and many of them are not annotated
                   #'44', '31', '52', '28', '50', # cluster '44', '31', '52', '28', '50' were not included here
                   '3', '5', '16', '30', '22', '4' # all middle time points
  )
  
  cells.sels = unique(colnames(seurat.obj)[seurat.obj$BWM.cells == 'BWM' &
    (!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
                                             seurat.obj$manual.annot.ids == 'MSxpppp'|
                                             seurat.obj$manual.annot.ids == 'MSxpppa'|
                                             seurat.obj$manual.annot.ids == 'MSxppap'|
                                             seurat.obj$manual.annot.ids == 'MSxppaa'|
                                             seurat.obj$manual.annot.ids == 'MSxpapp'|
                                             seurat.obj$manual.annot.ids == 'MSxpapa'|
                                             seurat.obj$manual.annot.ids == 'MSxpaap'|
                                             seurat.obj$manual.annot.ids == 'MSxpaaa'|
                                             
                                             seurat.obj$manual.annot.ids == 'MSxappp' |
                                             seurat.obj$manual.annot.ids == 'MSpappa' |
                                             
                                             seurat.obj$manual.annot.ids == 'MSxppp'|
                                             seurat.obj$manual.annot.ids == 'MSxppa'|
                                             seurat.obj$manual.annot.ids == 'MSxpap'|
                                             seurat.obj$manual.annot.ids == 'MSxpaa'|
                                             seurat.obj$manual.annot.ids == 'MSxapp'
                                             
                                             # seurat.obj$manual.annot.ids == 'MSxpp'|
                                             # seurat.obj$manual.annot.ids == 'MSxpa'|
                                             # seurat.obj$manual.annot.ids == 'MSxap'|
                                             # 
                                             # seurat.obj$manual.annot.ids == 'MSxp'|
                                             # seurat.obj$manual.annot.ids == 'MSxa'|
                                             # seurat.obj$manual.annot.ids == 'MSx'
                                           )])
  
  # seurat.obj$BWM.cells[!is.na(match(colnames(seurat.obj), cells.sels))] = 'BWM'
  # seurat.bwm = subset(seurat.obj, cells = cells.sels)
  # seurat.nonbwm = subset(seurat.obj, cells = setdiff(colnames(seurat.obj), cells.sels))
  # 
  # xx = table(seurat.nonbwm$seurat_clusters[which(seurat.nonbwm$manual.annot.ids == 'unknown.MSxpppaa')])
  # xx[which(xx > 0)]
  # sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[which(seurat.obj$BWM.cells == 'BWM')])
  
  sub.obj = subset(seurat.obj, cells = cells.sels)
  
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)
  barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)
  
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500', label = FALSE)
  #p2 = DimPlot(sub.obj, reduction = 'umap', group.by = 'seurat.pred.id', label = FALSE)  
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67', 'pat-4'))
  
  ##########################################
  # check potential ids for selected clusters
  ##########################################
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
  threshold = 0.7
  predicted.ids = sub.obj$scmap.pred.id.500
  #predicted.ids[which(sub.obj$scmap.corr.500 < threshold)] = 'unassigned'
  
  if(Refine.annotated.ids){
    counts = table(predicted.ids, sub.obj$manual.annot.ids)
    counts.seurat = table(as.character(sub.obj$seurat.pred.id), sub.obj$manual.annot.ids)
  }else{
    counts = table(predicted.ids, as.character(sub.obj$seurat_clusters))
    counts.seurat = table(as.character(sub.obj$seurat.pred.id), as.character(sub.obj$seurat_clusters))
  }
  barplot(counts, main="cluster compositions by scmap ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))
  
  barplot(counts.seurat, main="cluster compositions by seurat ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))
  
  #counts[, match(c('31', '28', '52'), colnames(counts))]
  #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  Explore.umap.parameters.for.BWMcells = FALSE
  if(Explore.umap.parameters.for.BWMcells){
    require(tictoc)
    tic()
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    test.umap.params.for.BWM.cells(sub.obj, pdfname = 'BWM_terminal_middle.timepoints_2.pdf',
                                   nfeatures.sampling = c(8000), nb.pcs.sampling = c(30, 50), n.neighbors.sampling = c(10, 30, 50), 
                                   min.dist.sampling = c(0.01, 0.05, 0.1, 0.2)
                                   )
    toc()
    
  }
  
  nfeatures = 8000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 50)
  
  nb.pcs = 30 # nb of pcs depends on the considered clusters or ids 
  n.neighbors = 50;
  min.dist = 0.01; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                     spread = spread, n.neighbors = n.neighbors, 
                     min.dist = min.dist)
  
  #DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
  
  ##########################################
  # rerun the seurat for label transferring
  ##########################################
  RErun.seurat.transferring.labels = TRUE
  if(RErun.seurat.transferring.labels){
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    terminals = c('MSxppppx', 'MSxpppax', 'MSxppppp','MSxppppa', 'MSxpppaa', 'MSxpppap',
                  'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 
                  'MSxpaaap', 'MSxapppp', 'MSxapppa', 
                  'MSxappppx', 'MSxapppax', 'MSpappax', 
                  'MSxppp', 'MSxppa', 'MSxppap', 'MSxpapp', 'MSxpapa', 'MSxpaaa', 'MSxappp', 'MSpappa',
                  'MSxpp', 'MSxpa', 'MSxap')
    
    #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
    sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30, 
                                                                              k.anchor = 5, k.filter = 200,
                                                                              terminals = terminals)
    sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
    sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA
    
  }
  
  p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
  p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) + NoLegend()
  #p2 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
  #  NoLegend()
  p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    NoLegend()
  
  p2 + p3
  
  DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
  
  p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1
  
  FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'pha-4', 'nhr-67'))  
  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:nb.pcs, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.2)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)
          
  p0 = DimPlot(sub.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)
  p1  = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) + 
    NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, 
               label.size = 6,
               na.value = "gray", combine = TRUE)
  p3 =   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
  (p0 + p3) / (p1 + p2)
  
  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3, 
          group.by = 'seurat_clusters_split')
  
  (p1 + p2) / p3
  
  manual.discovery.new.features = TRUE
  if(manual.discovery.new.features){
    Idents(sub.obj) = sub.obj$seurat_clusters_split
    markers <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
    
  }
  
  dev.off()
  
  #################################################################################################################################
  # update of manually annotated ids using marker genes and potential mapped labels from scmap or seurat
  # c('MSxppppx', 'MSxpppax', 'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 'MSxapppp', 'MSxapppa', 
  # 'MSxappppx', 'MSxapppax', 'MSpappax') 
  ##################################################################################################################################
  RErun.seurat.transferring.labels = FALSE
  if(RErun.seurat.transferring.labels){
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    terminals = c('MSxppppx', 'MSxpppax', 'MSxppppp','MSxppppa', 'MSxpppaa', 'MSxpppap',
                  'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 
                  'MSxpaaap', 'MSxapppp', 'MSxapppa', 
                  'MSxappppx', 'MSxapppax', 'MSpappax', 
                  'MSxppp', 'MSxppa', 'MSxppap', 'MSxpapp', 'MSxpapa', 'MSxpaaa', 'MSxappp', 'MSpappa')
    
    #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
    sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 5000,  npcs = 50, terminals = terminals)
    
  }
  
  sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
  sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
  sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
  sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA
  
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  counts.seurat = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
  counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
  
  #sub.obj$manual.annot.ids = sub.obj$predicted.ids
  p1 = DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
  p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()
  p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()
  p4 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
  (p4 + p3)/(p1 + p2)
  
  
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
  idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))
  
  idents.sel = c('23', '3', '13', '7', '2')
  
  ## chcek the reference-mapped ids for the rest of clusters
  #counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
  #counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
  counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
  counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
  counts.seurat.filter.sel = counts.seurat.filter[, !is.na(match(colnames(counts.seurat.filter), idents.sel))]
  counts.seurat.filter.sel = counts.seurat.filter.sel[apply(as.matrix(counts.seurat.filter.sel), 1, sum)>0, ]
  counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
  counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]
  
  #features.sels = c('unc-120', 'hnd-1', 'hlh-1', 'abts-1', 'ref-2', 'tbx-7', 'unc-39', 'cup-4', 'ins-2', 'F40H3.3', 'hot-1')
  features.sels = c('unc-120', 'ham-1', 'tab-1', 'pat-1')
  
  VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)
  
  # check info in JM data for specific lineage
  ids.sel = c('MSxppppx')
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  
  top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
  
  # to find new marker genes
  top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  top.markers[top.markers$cluster == '13',]
  
  FeaturePlot(sub.obj, reduction = 'umap', features = c('hnd-1', 'tbx-7', 'abts-1', 'Y66D12A.13', 'F41D9.2', 'unc-120', 'zig-8')) # MSxpaaa
  FeaturePlot(sub.obj, reduction = 'umap', features = c('tbx-8', 'F40H3.3', 'ins-2', 'Y42H9B.3', 'unc-120', 'ZK180.5')) # MSxpapa
  FeaturePlot(sub.obj, reduction = 'umap', features = c('zig-6', 'frpr-8','C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-4', 'maph-1.2', 'irx-1')) # MSxpapa
  
  FeaturePlot(sub.obj, reduction = 'umap', features = c('tnt-3', 'gana-1','F37H8.5','kvs-5'))
  
  
  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  #save first the seurat projection result, one of best projection prediciotn
  seurat.obj$predicted.ids.seurat.keep = NA
  mm = match(colnames(sub.obj), colnames(seurat.obj))
  seurat.obj$predicted.ids.seurat.keep[mm] = sub.obj$predicted.ids
  seurat.obj$manual.annot.ids.1 = seurat.obj$manual.annot.ids
  seurat.obj$manual.annot.ids = seurat.obj$predicted.ids.seurat.keep
  
  
  # cluster.assingment = list(c('23', 'MSxappp')
  #                           
  #                             
  # )
  # 
  # for(n in 1:length(cluster.assingment)){
  #   cluster.index = cluster.assingment[[n]][1]; 
  #   id2assign =  cluster.assingment[[n]][2]; 
  #   cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')    
  #   cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
  #   sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
  #   seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
  # }
  # 
  # #manual.assign.cluster.with.annotation(cluster.index = '17', id2assign = 'MSxppapp', sub.obj = sub.obj, seurat.obj = seurat.obj)
  # DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
  #         pt.size = 2)
  # 
  # 
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
     scale_colour_hue(drop = FALSE)

  #FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))  
  
  saveRDS(seurat.obj, file = RDS2save)
  
  
}

