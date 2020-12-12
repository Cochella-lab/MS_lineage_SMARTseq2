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
  
  ########################################################
  ########################################################
  # Section : iteration 36 (terminal cells) 
  ## we will always check the mother, current generation and daughter cells together
  ## After using the seurat prediction for middle time points, we focus on  
  ## all terminal cells and assocaited mother cells except convergence branch 
  ########################################################
  ########################################################
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
  
  # ee = process.import.Murray.scRNA()
  # murray.ids = unique(ee$lineage)
  # markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
  #markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
  load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
  
  ##########################################
  # manual.annot.ids.1 -- inital manual annotation for early and middle time points before using terminal cells annotation from seruat
  # manual.annot.ids.2 -- saved the manual annotation for non-BWM cells 
  # manual.annot.ids.3 -- saved the manual annotation before using seurat prediction for cells in middle time points
  # manual.annot.ids.6 -- saved manual annotation after iteration 34
  # predicted.ids.seurat.keep -- saved the seurat prediction for middle and terminal cells
  # pred.ids.seurat.keep.bwm.all -- saved the seurat prediction for all bwm cells (early, middle and temrinal cells)
  ##########################################
  nb.iteration = 37
  Refine.annotated.ids = TRUE;
  
  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_', 
                    nb.iteration -1, '.rds')
  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
  RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_', 
                      nb.iteration, '.rds')
  seurat.obj = readRDS(file = RDSsaved)
  
  #seurat.obj$predicted.ids.scmap = seurat.obj$scmap.pred.id.500
  #seurat.obj$predicted.ids.seurat = seurat.obj$seurat.pred.id
  #saveRDS(seurat.obj, file = RDS2save)
  
  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}
  
  ## 38 ids including MSx and redundant MSxppppx and MSxpppax; 
  # MSx missing and not necessary; 'MSapaap', 'MSppaap' for 'MSxpaap'
  bwms.all = c('MSxppppx', 'MSxpppax', # redundant
               'MSxppppp', 'MSxppppa', 'MSxpppaa', 'MSxpppap',
                'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 
                'MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax', 'MSpappax', # all terminal cells
                'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxppaa', 'MSxpapp', 'MSxpapa',  
                'MSapaap', 'MSppaap', #'MSxpaap', 
               'MSxpaaa', 'MSxappp', 'MSpappa',  
                'MSxppp', 'MSxppa', 'MSxpap', 'MSxpaa', 'MSxapp', 
                'MSxpp', 'MSxpa', 'MSxap',
                'MSxp', 'MSxa',
               'MSx')
  
  ## 16 terminal cells also including redudance
  terminals = c('MSxppppx', 'MSxpppax', # redundant
                'MSxppppp','MSxppppa', 'MSxpppaa', 'MSxpppap',
                'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 
                'MSxpaaap', 'MSxapppp', 'MSxapppa', 
                'MSxappppx', 'MSxapppax', 'MSpappax'
  )
  
  
  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  ##########################################
  # select subset of cells to annotate
  ##########################################
  # cluster.index = '25'
  # table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
  # 
  # table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
  # 
  # table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]
  # 
  # xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
  # xx[which(xx > 0)]
  # 
  # ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
  # xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
  # table(seurat.obj$manual.annot.ids[ii1])
  # #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
  
  # select BWM terminal cells
  ##########################################
  #cluster.sels = c('25', '31', '28', '52', '36', '8', '39', '2', '19', '27')
  #cluster.sels = c('25', '36', '8', '39', '2', '19', '27', '13', '1', '11', '33', '48', '18', '46', '15', '26')
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))] 
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))] 
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))] 
  #                                   & is.na(seurat.obj$manual.annot.ids))]
  
  # select BWM terminal and middle time points cells
  #' ##########################################
  #' cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition 
  #'                  '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
  #'                  '25', # possible transition clusters
  #'                  '24', # also transitions cells and many of them are not annotated
  #'                  #'44', '31', '52', '28', '50', # cluster '44', '31', '52', '28', '50' were not included here
  #'                  '3', '5', '16', '30', '22', '4' # all middle time points
  #' )
  #' 
  #cluster.sels = c('29', '32', '35', '40', '42')
  #sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  # cells.sels = unique(colnames(seurat.obj)[seurat.obj$BWM.cells == 'BWM' &
  #   (!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
  #    !is.na(match(seurat.obj$manual.annot.ids, ids.sels))                                       
  #     )])
  # cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
  #                                               !is.na(match(seurat.obj$manual.annot.ids, ids.sels))                                       
  #                                            ])
  #seurat.obj$manual.annot.ids.6 = seurat.obj$manual.annot.ids
  #seurat.obj$BWM.cells[which(seurat.obj$manual.annot.ids == 'unknown_MSxpppaa_MSxppppa_later')] = NA
  kk = which(seurat.obj$pred.ids.terminals.mothers.seurat == 'MSpappax')
  
  ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
  #ids.sels = ids.current[which(nchar(ids.current)>6)]
  #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')
  # ids.sels = c('MSxapp', 'MSxappp', 'MSpappa', 
  #              'MSxppp', 'MSxppa', 'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxppaa', 
  #              'MSxpap', 'MSxpaa', 'MSxpapp', 'MSxpapa', 'MSxpaap', 'MSxpaaa',
  #              "MSxpaaap/MSxppapp/MSxpapap/MSxppap", 
  #              "MSxpapp/MSxpppp",
  #              "MSxpapp/MSxpapa",
  #              "MSxppap/MSxpaaa/MSxpaaap",
  #              "MSxppap/MSxpaaa/MSxpapa/MSxpaaap",
  #              "MSxppppa/MSxppppp/MSxpppaa/MSxpappa/MSxappp",
  #              "MSxppppp/MSxappp" 
  #              #"MSxapppp/MSxapppa",
  #              #"MSxapppa", "MSxapppp" 
  #              )
  #ids.sels = c('MSxapp', 'MSxappp', 'MSpappa')
  #ids.sels = c('MSxpp', 'MSxppa', 'MSxppp', 'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxppaa')
  #ids.sels = c('MSxpa', 'MSxpaa', 'MSxpap', 'MSxpaaa', 'MSxpaap', 'MSxpapa', 'MSxpapp', 
  #             'MSxpp', 'MSxppa', 'MSxppp', 'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxppaa')
  #ids.sels = c('MSxpa', 'MSxpaa', 'MSxpaaa', 'MSxpap')
  #ids.sels = c('MSxpp', 'MSxppp', 'MSxppa', 'MSxppaa')
  #ids.sels = setdiff(ids.current, c('MSxa', 'MSxp', 'MSxap', 'MSxpp.early', 'MSxpa', 'MSxpa.early', 
  #                                  'MSxppp', 'MSxppa', 'MSxpap', 'MSxpaa', 'MSxapp', "MSxppaa"))
  # ids.sels = c(#'MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax', 
  #               'MSxppppx', 'MSxpppax', 'MSxppppa', 'MSxppppp', 'MSxpppaa', 'MSxpppap', 'MSxpappa', 
  #               'MSxppapp', 'MSxpapap', 'MSxpappp', 'MSxpaaap')
  #ids.sels = c('MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax')
  #ids.sels = c('MSxppppx', 'MSxpppax', 'MSxppppa', 'MSxppppp', 'MSxpppaa', 'MSxpppap', 'MSxpappa')
  #ids.sels = c('MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 'terminal.from.clusters.19.27', 
  #             'MSxppap',  "mixture_MSxpppp_MSxpppa_daughters")
  # ids.sels = c('MSxppppx', 'MSxpppax',
  #              'MSxppppp', 'MSxppppa', 'MSxpppaa', 'MSxpppap',  
  #              'MSxpppa', 'MSxpppp')
  # ids.sels = c('mixture_MSxppapp_MSxpappp', 'mixture_terminals', 'MSxpappa')
  # ids.sels = c('likely_MSxpapap', 'mixture_MSxppapp_MSxpappp', 'likely_MSxpaaap', 'MSxpaaap', 'MSxppapp', 
  #              'MSxpapa', 'MSxppap', 'MSxpapp', 'mixture_MSxpppa_MSxpppp_MSxpapa_MSxpapp_MSxppap', 'likely_MSxpappa', 'likely_MSxppppx',
  #              'MSxppppp', 'MSxpppaa', 'MSxppppa', 'MSxpppap', 
  #              'MSxpppp', 'MSxpppa')
  #' ids.sels = c('likely_MSxpappp', 'likely_MSxppapp',
  #'             #'MSxppap', 'MSxpapp', 
  #'              'MSxppapp'
  #'              )
  # ids.sels = c('MSxppppp', 'MSxpppaa', 'MSxppppa', 'MSxpppap', 'MSxpppp', 'MSxpppa',
  #              'likely_MSxpappa', 'likely_MSxppppx', 'mixture_MSxpppa_MSxpppp_MSxpapa_MSxpapp_MSxppap')
  #' ids.sels = c('likely_MSxppppx', 'mixture_MSxpppa_MSxpppp_MSxpapa_MSxpapp_MSxppap', 'mixture_terminal_mothers',
  #'            'MSxpappa') 
  #'            #''MSxpppa') 
  #'            #'MSxpppa', 'MSxpapa', 'MSxpapp', 'MSxppap')
  # ids.sels = c('mixture_terminal_MSxpppp_MSxpppa', 'MSxppppp', 'MSxppppa', 'MSxpppaa', 'MSxpppap', 'MSxpappa')
  # ids.sels = c('MSxppapp/MSxpappp', 'MSxppapp')
  # ids.sels = c('MSxpa.early', 'MSxpa')
  #' ids.sels = c('mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap', 
  #'              'MSxppapp/MSxpappp', 'MSxpapap', 'MSxpaaap', 'MSxppapp')
  #'              #'MSxpapa')
  ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)], 
                     c('MSxppaa'))
  
  ids.left = setdiff(ids.current, ids.sels)
  print(ids.left)
  nchar(ids.left)
  
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])
  sub.obj = subset(seurat.obj, cells = cells.sels)
  
  #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  sub.obj$pred.ids = sub.obj$predicted.ids.seurat.keep
  xx = table(sub.obj$predicted.ids.seurat.keep)
  xx[xx>10]
  sub.obj$pred.ids.filtered = sub.obj$pred.ids
  sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA
  
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
    test.umap.params.for.BWM.cells(sub.obj, 
                                   pdfname = 'UMAP_param_TEST_BWM_searching_for_MSx.pdf',
                                   group.by = 'manual.annot.ids', with_legend = TRUE,
                                   nfeatures.sampling = c(500, 1000), nb.pcs.sampling = c(10, 20, 30),
                                   n.neighbors.sampling = c(5, 10, 30, 50), 
                                   min.dist.sampling = c(0.01, 0.1)
                                   )
    toc()
    
  }
  
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
  
  #VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'manual.annot.ids') + NoLegend()
  #xx = table(sub.obj$predicted.ids.seurat.keep)
  #xx[xx>10]
  #sub.obj$pred.ids.filtered = sub.obj$pred.ids
  #sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA
  #jj2 = which(!is.na(match(sub.obj$predicted.ids.seurat.keep, c('MSxppppx', 'MSxpppax'))) == TRUE)
  p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + 
    NoLegend()
  p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
  #p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat.keep', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) 
  #p2 = DimPlot(sub.obj, group.by = 'pred.ids.terminals.mothers.svm', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) 
  p0 + p1
  #p0 + p2
  
  idntt = c('MSxpapap')
  cells_to_show <- list(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, idntt))])
  p1 = DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
               cells.highlight = cells_to_show, cols.highlight = 'blue', sizes.highlight = 1,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  #idntt = 'MSxpappp'
  cells_to_show <-  list(colnames(sub.obj)[which(sub.obj$predicted.ids.seurat.keep == idntt)])
  p2 = DimPlot(sub.obj, group.by = "pred.ids", reduction = 'umap', label = TRUE, repel = TRUE,
               cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1.5,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  p1 + p2
  
  idntt = c('likely_MSxpappp')
  cells_to_show <- list(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, idntt))])
  p3 = DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
               cells.highlight = cells_to_show, cols.highlight = 'blue', sizes.highlight = 1,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  idntt = 'MSxpappp'
  cells_to_show <-  list(colnames(sub.obj)[which(sub.obj$predicted.ids.seurat.keep == idntt)])
  p4 = DimPlot(sub.obj, group.by = "pred.ids", reduction = 'umap', label = TRUE, repel = TRUE,
               cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1.5,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  p3 + p4
  
  
  # MSx
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSx', group_2 = 'MSxp', ntop = 10)
  features.sels = c( 'sdz-1', 'hnd-1', 'pha-4', 'fbxb-70'
                    # 'clec-91', 'clec-87',  'clec-87','arrd-2','clec-88', 'otub-3'
  )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # MSxpa
  FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-39', 'R11A5.3', 'egl-43'))
  FeaturePlot(sub.obj, reduction = 'umap', features = c('zig-8','ham-1', 'F19C6.4',  'shc-2', 'F07C6.4', 
                                                         'unc-120', 'abts-1' # MSxpaaa
                                                        ))
  # MSxpappa
  #source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpappp', group_2 = 'MSxppapp', ntop = 20)
  features.sels = c('stn-2', 'Y9C2UA.1', 'hsp-12.1', 'unc-62', 
                    'ceh-34', 'zig-6',
                    'fbxb-88', 'fbxc-24', 'E01G4.5',  'igcm-4' 
                    )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpappp', group_2 = 'MSxppapp', ntop = 20)
  features.sels = c('maph-1.3', 'shc-2', 'twk-31', 'mec-2',   # MSxppapp
                    'D1086.12', 'tbx-2', 'K09G1.1', 'zig-6', 'stn-2', 'T25G12.11', 'unc-5') # MSxpappp   
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # MSxpapap (almost sure)
  features.sels = c('maph-1.2', 'maph-1.3', 'K09G1.1', 'T25G12.11', 'stn-2', 'hnd-1', 'Y37E3.30', 'zig-6', 'Y9C2UA.1')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  
  DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()
  
  p1 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
               pt.size = 2) + NoLegend() +
    ggtitle('manual ids')
  
  p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat.keep', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,  
               pt.size = 2) + 
    NoLegend() + ggtitle('seurat.pred.ids.keep') 
  
  p3 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,  
               pt.size = 2) + 
    NoLegend() + ggtitle('seurat.pred.ids.bwm.all')
  
  p1 /(p2 + p3)
  
  FeaturePlot(sub.obj, reduction = 'umap', features = 'maph-1.2')
 
  ##########################################
  # rerun the seurat for label transferring
  ##########################################
  RErun.seurat.transferring.labels = FALSE
  if(RErun.seurat.transferring.labels){
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    
    ids.refs = c(terminals, 
                 'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxpapp', 'MSxpapa', 'MSxpaaa', 'MSxappp', 'MSpappa')
    #ids.refs = terminals
    #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
    sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30, 
                                                                              reduction = 'cca',
                                                                              k.anchor = 5, k.filter = 200, max.features = 200,
                                                                              terminals = ids.refs)
    sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
    
    hist(sub.obj$predicted.ids.prob)
    abline(v = 0.5, col = 'red', lwd=2.0)
    
    seurat.obj$pred.ids.terminals.mothers.seurat = NA
    seurat.obj$pred.ids.terminals.mothers.seurat.prob = NA
    mm = match(colnames(sub.obj), colnames(seurat.obj))
    seurat.obj$pred.ids.terminals.mothers.seurat[mm] = sub.obj$predicted.ids.seurat.terminal
    seurat.obj$pred.ids.terminals.mothers.seurat.prob[mm] = sub.obj$predicted.ids.seurat.terminal.prob
    
    
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    sub.obj = reference.based.cell.projection.rf.svm(sub.obj, nfeatures = 3000, cost = 0.5, ntree = 50, scale.cell = TRUE,
                                                                              terminals = ids.refs)
    
    sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA
    sub.obj$pred.ids.svm.filter = sub.obj$pred.ids.svm
    sub.obj$pred.ids.svm.filter[sub.obj$pred.ids.svm.prob <0.5] = NA
    
    p1 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      ggtitle('manual ids')
    #p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
    p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) + 
     ggtitle('predicted ids. seurat')
    p3 = DimPlot(sub.obj, group.by = 'pred.ids.svm', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) + 
      ggtitle('predicted ids svm')
    p4 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
     ggtitle('seurat.pred.ids seurat. filtered')
    p5 = DimPlot(sub.obj, group.by = 'pred.ids.svm.filter', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) + 
      ggtitle('predicted ids svm filtered')
    
    (p2 + p3) / (p4 + p5)
    
    p2 + p3
    p1 + p2
    p1 + p3
    
    DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
    #DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE)
    
    # save the seurat prediction for all BWM cells
    seurat.obj$pred.ids.terminals.mothers.seurat = NA
    seurat.obj$pred.ids.terminals.mothers.svm = NA
    #sub.obj$predicted.ids.seurat.keep.bwm.all = NA
    #sub.obj$pred.ids.terminals.mothers.seurat = sub.obj$predicted.ids.seurat.terminal
    #seurat.obj$pred.ids.seurat.terminals.mothers = NA
    mm = match(colnames(sub.obj), colnames(seurat.obj))
    seurat.obj$pred.ids.terminals.mothers.seurat[mm] = sub.obj$predicted.ids.seurat.terminal
    seurat.obj$pred.ids.terminals.mothers.svm[mm] = sub.obj$pred.ids.svm
    
  }
  
  p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1
  
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'pha-4', 'hnd-1'))
  #FeaturePlot(seurat.obj, reduction = 'umap', features = c('unc-120', 'pha-4', 'hnd-1'))  
  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.5)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)
  
  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
     ggtitle('manual.annot.ids') + NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 6, na.value = "gray", combine = TRUE)
  p1 + p2
  
  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'seurat_clusters_split')
  
  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
               group.by = 'seurat_clusters_split') + NoLegend()
  
  VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
          group.by = 'manual.annot.ids') + NoLegend()
  
  p1 + p4
  p2 + p4
  plot(p3)
  
  for(idntt in unique(sub.obj$manual.annot.ids))
  {
    cells_to_show <- list(c(WhichCells(sub.obj, idents = idntt)))
    p1  = DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
                  cells.highlight = cells_to_show, cols.highlight = 'blue', sizes.highlight = 1,
                  pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
    plot(p1)
  }
  
  idntt = c('MSxpapap'); cells_to_show <- list(c(WhichCells(sub.obj, idents = idntt)))
  DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
                cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1,
                pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  
  # check info in JM data for specific lineage
  #' source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  #' extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppppp', group_2 = 'MSxpppaa', ntop = 10)
  #' features.sels = c(#'hnd-1',  'pha-4',
  #'   #'fbxb-70',
  #'   #'ceh-37', 'C45G7.4', # MSxap
  #'   #'tbx-8', 'hlh-1', 'hnd-1', 'unc-39' # MSxpap
  #'   #'ham-1',  'nhr-67', # MSxapp
  #'   #'unc-39', 'irx-1', 'egl-43'
  #'   #'unc-120', 'tab-1', 'rpm-1', 'F55C5.10', 'hil-3' # MSpappa
  #'   #'tbx-8', 'asic-2', 'skpo-1', 'pxn-2', 'D1086.12','tab-1'  # MSxappp
  #'   #'mnr-1', 'strl-1', 'zig-10', 'gana-1', 'fbxb-88', # MSxapppa
  #'   #'lin-39', 'F54D5.5'  # MSxapppp
  #'   #'gsnl-1', 'abts-1', 'stn-2', # MSxapppax
  #'   #'lin-39',  'tbx-2', 'B0379.1'  #MSxappppx
  #'   #'egl-43', 'R11A5.3' #MSxpa
  #'   #'zip-7', 'hlh-1', 'col-111' # MSxpap
  #'   #'tbx-7', 'unc-120', 'ref-2', 'hlh-16', 'unc-39', 'ttr-50', 'clec-266' # MSxpaa
  #'   #'F41D9.2', 'hnd-1', 'unc-120', 'tbx-7', 'abts-1', 'Y66D12A.13' # MSxpaaa
  #'   #'let-381', 'ins-2', 'F40H3.3', 'ZK183.5', # MSxpapa
  #'   #'unc-120', 'ceh-51', 'lag-2', 'fbxb-22', 'T02G6.11', 'C06A8.3', 'T05D4.2',   # MSxpp
  #'   #'zip-7', 'hnd-1', 'tbx-8', 'fkh-2', 'tbx-11', # MSxppp
  #'   #'hlh-16', 'rgs-7', 'tbx-8', 'unc-120', 'T24C2.2', 'T24C2.3'  #MSxppa
  #'   #'col-118', 'let-381', 'C03B1.1',  'T11B7.2' #MSxppaa
  #'   #'F40H3.3', 'Y42H9B.3', 'unc-120', 'ZK180.5', 'sfrp-1', 'ceh-34', 'irx-1', 'Y116A8C.3', 
  #'   #'let-381', 'ccg-1' # MSxpapa 
  #'   #'zig-8', 'pat-9', 'F19C6.4', 'tbx-8', 'her-1' # MSxpapp
  #'   #'sul-2', 'camt-1', 'irx-1', 'ctg-1',  # MSxpppa 
  #'   #'tbx-2', 'D1086.12' #MSxpppp
  #'   #'ceh-34', 'ham-1', 'F19C6.4', 'eya-1', 'shc-2', 'F07C6.4', 'ten-1', 'hil-7' # MSxppap
  #'   #'maph-1.3', 'fkh-2', 'ccg-1', 'zig-6',
  #'   #'ham-1', 'shc-2', 'lam-3', 'T25G12.11','F57F4.4', 'cnt-2', 'Y17G7B.23', 'Y37E3.30'  # MSxpaaap
  #' )
  #' 
  # MSxppppp
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppppp', group_2 = 'MSxpppaa', ntop = 10)
  features.sels = c('clec-264', 'ceh-13', 'gana-1','D1086.12', 'tbx-2', 'F37H8.5', 'B0379.1' 
  )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # MSxpppap
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpppap', group_2 = 'MSxppppp', ntop = 10)
  features.sels = c('ceh-13', 'fbxc-24', 'hot-1', 'pde-6',  'hmg-1.1' )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # MSxpppaa
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpppaa', group_2 = 'MSxppppa', ntop = 10)
  features.sels = c('bgal-1', 'sul-1', 'zig-7')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # MSxppppa
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppppa', group_2 = 'MSxpppaa', ntop = 10)
  features.sels = c('lntl-1', 'D1086.12', 'tbx-2', 'ZK512.1')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # MSxpppp and MSxpppa
  features.sels = c( 'sul-2', 'camt-1', 'irx-1', 'ctg-1',  # MSxpppa 
                     'tbx-2', 'D1086.12' #MSxpppp
  )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('zig-6', 'eya-1'))
  
  # MSxpaaap
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpaaap', group_2 = 'MSxpapap', ntop = 20, test.use = 'roc')
  
  features.sels = c('maph-1.3', 'maph-1.2', 'ZC449.5', 'shc-2',
                    'T25G12.11', 'zig-8', 'ceh-34', 'ham-1', 'lam-3', 'twk-31', 'stn-2')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # MSxppapp
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppapp', group_2 = 'MSxpappa', ntop = 10, test.use = 'roc')
  
  features.sels = c('K09G1.1', 'tbx-2', 'D1086.12', 'maph-1.3', 'shc-2', 'twk-31', 
                    'zig-6', 'mec-2')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # MSxpappp
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpappp', group_2 = 'MSxpapap', ntop = 10)
  
  features.sels = c('zig-8', 'D1086.12', 'K09G1.1', 'tbx-2', 'zig-6', 'stn-2',
                    'lipl-7')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # MSxpappa
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpappa', group_2 = 'MSxpppaa', ntop = 10)
  
  features.sels = c('Y9C2UA.1', 'hsp-12.1', 'fbxc-24',
                    'unc-62', 'cpn-3', 'lipl-7', 'fbxb-88', 'tnt-3', 'lbp-1', 'F37H8.5', 'E01G4.5',
                    'ceh-34', 'zig-6'
                    
  )
  
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # MSxpapap
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpapap', group_2 = 'MSxpappp', ntop = 10)
  
  features.sels = c('maph-1.2', 'maph-1.3', 'K09G1.1', 'T25G12.11', 'stn-2', 'hnd-1', 'Y37E3.30', 'zig-6', 'Y9C2UA.1')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  
  # MSxpapa and MSxpapp
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpapa', group_2 = 'MSxpapp', ntop = 10)
  
  features.sels = c('zig-8', 'tbx-8', 'hlh-1', 'pat-9','tag-196', 'fkh-2', 'spp-15', 'cnt-2',
                    'her-1','ceh-34', 'sfrp-1','let-381', 'irx-1', 'ins-2'
  )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  
  manual.discovery.new.features = FALSE
  if(manual.discovery.new.features){
    
  
  }
  
  dev.off()
   
  ##########################################
  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  ##########################################
  # save current manual.annot.ids, early and middle time pooints from myself and terminal cells annotated by seurat
  # seurat.obj$manual.annot.ids.3 = seurat.obj$manual.annot.ids 
  # jj1 = which(!is.na(sub.obj$predicted.ids.seurat.keep))
  # sub.obj$manual.annot.ids[jj1] = sub.obj$predicted.ids.seurat.keep[jj1] # use the seurat prediction for cells in the middle time points
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpaaap')] = 'MSxpaaap'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpapap')] = 'MSxpapap'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
  
  sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
  sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
  sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'
  
  sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
  seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
  #mm = match(colnames(sub.obj), colnames(seurat.obj))
  #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
  
  # sub.obj$ids.backup = sub.obj$manual.annot.ids
  # 
  cluster.assingment = list(#c('0', 'MSxp'),
                            c('5', 'MSx')
                            #c('1', 'mixture_BWM_terminal_2'),
                            #c('2', 'mixture_BWM_terminal_2'),
                            #c('3', 'mixture_BWM_terminal_2'),
                            #c('4', 'mixture_BWM_terminal_2')
  )
  
  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    if(is.na(id2assign)) {
      seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
    }else{
      seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
    }
  }
  
  #jj1 = which(sub.obj$ids.backup == 'MSxppap')
  #sub.obj$manual.annot.ids[jj1] = 'MSxppap'
  #sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxpppa_MSxpppp_MSxpapa_MSxpapp_MSxppap')] = 'mixture_BWM_terminal_2'
  #sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_terminal_mothers')] = 'mixture_BWM_terminal_2'
  seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
  #mm = match(colnames(sub.obj), colnames(seurat.obj))
  
  #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
           pt.size = 2) + NoLegend()
  
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
     scale_colour_hue(drop = FALSE) + NoLegend()
  
  # jj2 = which(seurat.obj$manual.annot.ids == 'mixture_BWM_terminal_2')
  # seurat.obj$manual.annot.ids[jj2] = 'mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap'
  # 
  # jj3 = which(seurat.obj$manual.annot.ids == 'mixture_terminal_MSxpppp_MSxpppa')
  # seurat.obj$manual.annot.ids[jj3] = 'mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa'
  # 
  # jj4 = which(seurat.obj$manual.annot.ids == 'MSxpa.early')
  # seurat.obj$manual.annot.ids[jj4] = 'MSxpa'
  # 
  # seurat.obj$manual.annot.ids[which(seurat.obj$manual.annot.ids == 'MSxpp.early')] = 'MSxpp'
  # 
  saveRDS(seurat.obj, file = RDS2save)
  
}

########################################################
########################################################
# Section : save and clean bwm annotation from most updated bwm annotation 
# prepare the seurat.obj for the pharynx annotation
########################################################
########################################################
Save.and.Clean.annotation.from.bwm = function(seurat.obj)
{
  # import the most updated annotation of bwm that is 37th iteration 
  nb.iteration = 37
  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_', 
                    nb.iteration, '.rds')
  seurat.obj = readRDS(file = RDSsaved)
  
  
  # 'manual.annot.ids.2' have the annotated non-BWM cell ids saved  
  DimPlot(seurat.obj, group.by = "manual.annot.ids.2", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_annotedIDs")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  seurat.obj$manual.annot.ids.backupBWM = seurat.obj$manual.annot.ids
  
  
  # put back non-BWM cluster annotations; but one can select bwm with metatdata 'BWM.cells'  
  kk = which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$manual.annot.ids.2)) 
  
  cat(length(kk), ' cells annotated as non-BWM cell ids; \nhowever due to early iteration some of them have BWM annotation  \n')
  
  ids.1 = names(table(seurat.obj$manual.annot.ids)) # lastest BWM annotations
  ids.2 = names(table(seurat.obj$manual.annot.ids.2[kk])) # cells with non-BWM annotations
  print(ids.1);
  print(ids.2)
  
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  
  ##########################################
  # here is what we wanna do:
  # if ids.2 have BWM annotations, it means that the annotations were outdated and put back to NA
  # if ids.2 are not annotated as BWM, they were put back to clean annotations
  ##########################################
  ids2clean = intersect(ids.2, ids.1)
  ids2add = setdiff(ids.2, ids.1)
  
  jj = which(!is.na(match(seurat.obj$manual.annot.ids.2, ids2add)))
  cells2add = unique(colnames(seurat.obj)[jj])
  cat(length(cells2add), 'cells will be added non-BWM annotatation\n')
  
  seurat.obj$manual.annot.ids[jj] = seurat.obj$manual.annot.ids.2[jj]
  
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  
  idntt = c("MSapaapp")
  cells_to_show <- list(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, idntt))])
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
          cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1,
          pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  
  idntt = c("MSppaapp")
  cells_to_show <- list(colnames(seurat.obj)[!is.na(match(seurat.obj$predicted.ids.scmap, idntt))])
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
          cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1,
          pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  
  saveRDS(seurat.obj, file = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot', 
                                 '_cleanedBWM_and_Pharynx_iteration_0.rds'))
  
  
}

########################################################
########################################################
# Section : main function for pharynx annotation
# here we will start the seurat.obj with iteration 37 from the BWM annotation
# the iteration will be added
# the newly annotated ids were added also in manual.annot.ids while manual.annot.ids.bwm.backup was to save the results from bwm annotaiton
########################################################
########################################################
manual.annotation.for.pharynx.clusters = function(seurat.obj = seurat.obj)
{
  library(ggplot2)
  library(patchwork)
  library("pheatmap")
  library("RColorBrewer")
  library(grid)
  #library(RaceID) # refer to the vignett https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html
  library(Matrix)
  #library(lsa)
  library(dplyr)
  library(openxlsx)
  
  ##########################################
  # Main aim:
  # here we solved the mixed terminal cells ids.sels = c("MSxppapp/MSxpappp") 
  # 
  # 
  # Notes:    
  # all mixture annotation were solved for the first time and requires comfirmation
  # 
  ##########################################
  GR.iteration = 7 # RG (revison global)
  Refine.annotated.ids = TRUE
  
  resDir = paste0("results/", version.analysis, '/annoted_pharynx')
  if(!dir.exists(resDir)){dir.create(resDir)}
  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}
  
  
  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot', 
                                 '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
  RDS2save =  paste0(RdataDir,
                     'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot', 
                     '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')
  
  seurat.obj = readRDS(file = RDSsaved)
  
  #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), 
  #    ' pharynx cells left to annotate \n')
  
  cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
  cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')
  
  #pdf(pdfname, width=18, height = 10)
  #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  ##########################################
  # select subset of cells 
  ##########################################
  
  # select cells with cluster index
  ##########################################
  jj = which(is.na(seurat.obj$manual.annot.ids))
  print(table(seurat.obj$seurat_clusters[jj]))
  
  #cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('6', '24', '20', '7')
  #cluster.sels = c('13', '17', '18', '11', '25', '10')
  # cluster.sels = c('1', '9', '12', '19')
  # #ids.sel = c('MSxaaa')
  # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
  #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
  cluster.sels = c('28', '52', '31')
  cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                   # & is.na(seurat.obj$manual.annot.ids)
                                    ]
  table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
          cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
          pt.size = 2, label.size = 5) + NoLegend() 
  
  #                                           
  # if(length(ids.sel)>0){
  #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
  #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
  #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
  # }else{
  #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  # }
  
  # select cells with ids
  ##########################################
  #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
  #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')
  
  clean.id.names = FALSE # clean id names
  if(clean.id.names){
    # jj = which(seurat.obj$manual.annot.ids == "MSpaaaap.sure")
    # seurat.obj$manual.annot.ids[jj] = "MSpaaaap"
    # jj = which(seurat.obj$manual.annot.ids == "MSaaaappp.MSxapaapp")
    # seurat.obj$manual.annot.ids[jj] = "MSaaaappp/MSxapaapp"
    # 
    # jj = which(seurat.obj$manual.annot.ids == "MSpaaappp.MSxapappa")
    # seurat.obj$manual.annot.ids[jj] = "MSpaaappp/MSxapappa"
    # 
    # jj = which(seurat.obj$manual.annot.ids == 'unknown_MSxpppaa_MSxppppa_later')
    # seurat.obj$manual.annot.ids[jj] = "MSxpppaa_MSxppppa_later_unknown"
    # 
    # jj = which(seurat.obj$manual.annot.ids == 'unknown.MSxpppaa')
    # seurat.obj$manual.annot.ids[jj] = NA
    # 
    # jj = which(seurat.obj$manual.annot.ids == 'MSxapapp.pharynx')
    # seurat.obj$manual.annot.ids[jj] = NA
    # 
    # jj = which(seurat.obj$manual.annot.ids == 'likely.nonBWM_origCluster_17')
    # seurat.obj$manual.annot.ids[jj] = NA
    # 
    # jj = which(seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31')
    # seurat.obj$manual.annot.ids[jj] = NA
    
    # jj = which(seurat.obj$manual.annot.ids == 'MSxpppaa_MSxppppa_later_unknown')
    # seurat.obj$manual.annot.ids[jj] = 'MSxpppaa_MSxppppa_later_like'
    # 
    # jj = which(seurat.obj$manual.annot.ids == 'MSxpaap.MSppaapp.likely')
    # seurat.obj$manual.annot.ids[jj] = 'MSxpaap.MSppaapp.like'
    # 
    # jj = which(seurat.obj$manual.annot.ids == 'MSxapappp.likely')
    # seurat.obj$manual.annot.ids[jj] = 'MSxapappp.like'
    
    #jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")
    #seurat.obj$manual.annot.ids[jj] = 'mixture_transitionToTerminal'
    
  }
  
  
  ids.current = names(table(seurat.obj$manual.annot.ids))
  #ids.sels = ids.current
  ids.sels = setdiff(ids.current, 
                     c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp', 
                       'MSxaaa', 'MSpaaapp', 'MSpaaapa', 'MSaappa', 'MSpaaaap.sure', 
                        "MSaaaaa",  "MSaaaap", "MSpaaaa", "MSpaaap"))
   
  #ids.sels = c("mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap", "mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa",
  #             "mixture_transitionToTerminal")
  ids.sels = c("MSxppapp/MSxpappp")
  
  ids.left = setdiff(ids.current, ids.sels)
  # print(ids.left)
  # nchar(ids.left)
  
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])
  
  ##########################################
  # subset seurat object with selected cells
  ##########################################
  cat(length(cells.sels), ' cells selected to annotate \n')
 
  sub.obj = subset(seurat.obj, cells = cells.sels)
  
  #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
    
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()
  
  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  Explore.umap.parameters.for.BWMcells = FALSE
  if(Explore.umap.parameters.for.BWMcells){
    
    source.my.script('scRNA_functions.R')
    require(tictoc)
    tic()
    test.umap.params.for.BWM.cells(sub.obj, 
                                   pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                   group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                   nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                   n.neighbors.sampling = c(10, 30, 50), 
                                   min.dist.sampling = c(0.01, 0.1)
    )
    toc()
    
  }
  
  nfeatures = 5000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)
  
  nb.pcs = 5 # nb of pcs depends on the considered clusters or ids 
  n.neighbors = 10;
  min.dist = 0.01; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs), 
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist, verbose = TRUE)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) + 
    NoLegend()
  
  
  if(Refine.annotated.ids){
    
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) + 
      NoLegend()
    
  }else{
    DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) + 
      NoLegend()
    
    p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + 
      NoLegend()
    
    #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
    p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + 
      NoLegend() + ggtitle('scamp prediction')
    
    p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
      NoLegend() + ggtitle('seurat prediction')
    
    p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + 
      NoLegend()
    
    p0 + p2
    p0 + p1
    p1 + p2
    p2 + p3
    
  }
  
  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 1.0)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:5, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)
  
  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
    ggtitle('manual.ids') + NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
               label.size = 8, na.value = "gray", combine = TRUE)
  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()
  
  (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration, '.pdf'), 
                           width = 18, height = 16)
  
  #jj.missing = which(is.na(sub.obj$manual.annot.ids)) 
  #table(sub.obj$seurat_clusters_split[jj.missing])
  
  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'seurat_clusters_split')
  #dev.off()
  ##########################################
  # check the counts of predicted ids for newly split clusters
  ##########################################
  idents.sel = c('0', '1', '2', '3', '4')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
  counts.seurat = table(sub.obj$predicted.ids.seurat.keep, sub.obj$seurat_clusters_split)
  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
  if(exists('idents.sel')){
    if(length(setdiff(colnames(counts), idents.sel)) > 0){
      counts = counts[, !is.na(match(colnames(counts), idents.sel))]
      counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
      counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
      counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
      counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
      counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
    }
  }
  
  # features of early stage
  features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                    'ceh-36', # MSxaap
                    'tbx-2', 'ceh-27') # MSxaaa
  # features of lineage MSxapa  
  features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                    'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                    'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                    'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                    'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                    'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                    'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                    'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa
                    
  )
  # features of MSxaap lineage
  features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                    'sptf-1', 'ceh-22',  #MSaaapaa
                    'lin-12', 'pax-1', # MSxaapap
                    'ngn-1', 'ces-1', # MSpaapaa
                    'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                    'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                    'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                    'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
  )
  # features of MSxaaa
  features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                    'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                    #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                    'ceh-13', # MSaaaap not MSpaaap
                    'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                    'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                    'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                    'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                    'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                    'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                    'igcm-4', 'F54E2.2', 'K04G2.12', 
                    'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                    'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                    'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                    )
  
  features.sels = c('unc-62','gana-1','clec-264', 'maph-1.3', 'zig-6', 'maph-1.2', 
                    'mec-2', 'twk-31', 'stn-2', 'ZC449.5',  'shc-2', 'ham-1', 'ceh-34',
                    'K09G1.1', 'T25G12.11', 'shc-2',
                    'ceh-13',  'tbx-2', 'D1086.12', 'sul-1'
                    
                    #'gsnl-2', 'abts-1'
                    )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  # update the annotation with marker genes
  cluster.assingment = list(
    c('3', 'MSxpappp_new2'),
    c('1', 'MSxpappp_new2'),
    c('0', 'MSxppapp_new2'), 
    c('2', 'MSxppapp_new2'), 
    c('4', 'MSxppapp_new2')
  )
  
  cat(length(cluster.assingment), 'clusters assigned \n')
  cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
  nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
  if(nb.unassigned >= 1){
    cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
  } 
  
  # check info in JM data for specific lineage
  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
  #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  
  ids.sel = c('MSxpaaap'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)
  
  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  # add.missing.annotation.for.pharynx.left.over = FALSE
  # if(add.missing.annotation.for.pharynx.left.over){
  #   jj = which(is.na(sub.obj$manual.annot.ids))
  #   
  #   for(j in jj){
  #     #j = 1;
  #     cat('cell : ',  j, '\n')
  #     index.cluster = sub.obj$seurat_clusters_split[j]
  #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
  #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
  #   }
  #   
  #   mm = match(colnames(sub.obj), colnames(seurat.obj))
  #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
  #   
  # }
  
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'
  
  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    # if(is.na(id2assign)) {
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
    # }else{
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
    # }
  }
  
  #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
  
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
          pt.size = 2) + NoLegend()
  
  
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE) + NoLegend()
  
  
  saveRDS(seurat.obj, file = RDS2save)
  
  
}




