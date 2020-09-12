manual.annotation.for.BWM.clusters = function(seurat.obj = ms, ids = c('MSx'))
{
  ########################################################
  ########################################################
  # Section : initial test for manual cluster annotation of BDW 
  # 
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
  
  
  pdfname = paste0(resDir, "/Overview_predictedLabels_seuratClusters_mapping_seurat.pdf")
  pdf(pdfname, width=16, height = 12)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #ids = c('MSx')
  ids = c('MSx', 'MSxa', 'MSxap', 'MSxapp','MSxappa', 'MSxappp', 'MSxp', 'MSxpa', 'MSxpp')
  ids = c('MSxa', 'MSxap', 'MSxp', 'MSxpa', 'MSxpp')
  
  # given the fact that scmap seems to be working better than seurat 
  # for the manual annotation, the scmap predicted labels will be mainly considered, 
  # sstart with scmap projection with nb.features = 500 and threshold = 0.7
  # seurat.obj$predicted.ids = seurat.obj$scmap.tintori.id
  # seurat.obj$predicted.scores = seurat.obj$scmap.tintori.cor
  # threshold = 0.7
  seurat.obj$predicted.ids = seurat.obj$scmap.pred.id.500
  seurat.obj$predicted.scores = seurat.obj$scmap.corr.500
  threshold = 0.7
  
  seurat.obj$predicted.ids[which(seurat.obj$predicted.ids == 'unassigned')] = NA
  seurat.obj$predicted.ids.fitered = seurat.obj$predicted.ids
  seurat.obj$predicted.ids.fitered[which(seurat.obj$predicted.scores < threshold)] = NA;
  
  ##########################################
  # check the predicted labels and cluster-label mapping without filtering
  ##########################################
  cells.sels = which(!is.na(match(seurat.obj$predicted.ids, ids)))
  cells.fitered = cells.sels[which(seurat.obj$predicted.scores[cells.sels] > threshold)]
  #cat(length(cells.sels), ' cells with predicted labels \n')
  cat(length(cells.fitered), ' cells after filtering with predicted labels \n')
  #cat('clusters involved : '); print(as.character(unique(seurat.obj$seurat_clusters[cells.sels])));
  cluster.sels =  unique(seurat.obj$seurat_clusters[cells.fitered])
  cat('clusters involved after fitering : '); print(as.character(cluster.sels));
  
  p0 = DimPlot(seurat.obj,
               cells.highlight = colnames(seurat.obj)[cells.fitered],
               group.by = "predicted.ids.fitered", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 4,
               sizes.highlight = 1,
               na.value = "gray") + 
    ggtitle(paste0("cells predicted by scmap with 500 features and threshold 0.7 ")) +
    NoLegend()
  
  p00 = DimPlot(seurat.obj,
                cells.highlight = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))],
                group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 4,
                sizes.highlight = 1,
                na.value = "gray") +
    #ggtitle(paste0("cells predicted by scmap with 500 features and threshold 0.7 ")) +
    NoLegend()
  
  p0 + p00
  
  ##########################################
  # found involved clusters to consider in the following 
  ##########################################
  # make the count table for ids considered and involved cluster index
  predicted.ids = seurat.obj$predicted.ids.fitered
  predicted.ids[is.na(predicted.ids)] = 'unassigned'
  counts <- table(predicted.ids, seurat.obj$seurat_clusters)
  counts = counts[grep('MS|unassigned', rownames(counts)), ]
  
  nb.cells.clusters = apply(counts[rownames(counts) != 'unassigned', ], 2, sum)
  
  counts = counts[!is.na(match(rownames(counts), c(bwms, 'unassigned'))), ]
  nb.cells.bwm.clusters = apply(counts[rownames(counts) != 'unassigned',], 2, sum)
  
  counts = counts[!is.na(match(rownames(counts), c(ids, 'unassigned'))), ]
  nb.cells.ids.clusters = apply(counts[rownames(counts) != 'unassigned',], 2, sum)
  
  cols = c(colorRampPalette((brewer.pal(n = 7, name="Reds")))(max(counts)))
  pheatmap(counts, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, breaks = seq(0, max(counts), by = 1),
           cluster_cols=FALSE, main = paste0("cluster -- predicted labels mapping"), na_col = "white",
           color = cols, 
           #display_numbers = TRUE,
           #number_format = "%.0f",
           #annotation_col = my_sample_col,
           #gaps_row = c(1:nrow(map)-1),
           fontsize_col = 10,
           height = 8,
           width = 30
  )
  
  ## filter clusters that are not likely for the considered ids due to noise
  kk = which(nb.cells.ids.clusters >0 ) # > 0 cell in the cluster
  cat(length(kk), 'clusters for considered ids \n')
  print(colnames(counts)[kk])
  
  rr.bwm = (nb.cells.bwm.clusters/nb.cells.clusters)[kk]
  barplot(rr.bwm, main = '% of predited cells found in BWM');abline(h=0.5, col ='red')
  rr.ids = (nb.cells.ids.clusters/nb.cells.clusters)[kk]
  
  #kk = 
  kk = which(nb.cells.ids.clusters > 0 
             & nb.cells.ids.clusters/nb.cells.bwm.clusters > 0.5 # >50% of bwm cells have the ids considered 
             & nb.cells.bwm.clusters/nb.cells.clusters > 0.5 # > 50% of cell predicted to be in bwm
  ) 
  
  counts = counts[, kk]
  
  par(mfrow = c(1, 2))
  barplot(t(counts), main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = colnames(counts)
  )
  
  barplot((counts), main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))
  
  
  #counts.norm1 = counts
  #for(n in 1:nrow(counts)) counts.norm[n,] = counts.norm[n, ] /sum(counts.norm[n, ])
  
  #dev.off()
  
  #cluster.sels = colnames(counts)
  cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('29', '32', '35')
  
  sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
  
  #length(intersect(VariableFeatures(sub.obj), timers))
  #VariableFeatures(sub.obj) = setdiff(VariableFeatures(sub.obj), timers)
  cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 50)
  
  
  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids 
  n.neighbors = 10; min.dist = 0.2;
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
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.4)
  
  #DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5)
  
  p1  = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2)
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, 
               label.size = 5,
               na.value = "gray", combine = TRUE)
  p1 + p2
  
  VlnPlot(sub.obj, features = c('hnd-1', 'pha-4', 'sdz-1', 'sdz-31', 'ceh-51', 'ceh-36', 'unc-39', 'unc-120'),
          group.by = 'seurat_clusters_split')
  
  VlnPlot(sub.obj, features = c('timingEst', "FSC_log2", "BSC_log2"), ncol = 3, 
          group.by = 'seurat_clusters_split')
  
  # sc = SCseq(sub.obj@assays$RNA@counts)
  # sc <- filterdata(sc, mintotal=2000, minexpr = 10, minnumber = 2)
  # #sc <- compdist(sc,metric="pearson", FSelect = FALSE)
  # cat('use pca to calculate Pearson correlation and then distance and then k-mean\n')
  # sub.obj.pca = sub.obj@reductions$pca@cell.embeddings[, c(1:nb.pcs)]
  # #mat.dist = 1- cor(t(sub.obj.pca))
  # mat.dist = 1 - lsa::cosine(t(sub.obj.pca))
  # sc@distances = mat.dist
  # sc <- clustexp(sc, FUNcluster = 'kmedoids', verbose = FALSE)
  # 
  # par(mfrow = c(1, 2))
  # plotsaturation(sc,disp=FALSE)
  # plotsaturation(sc,disp=TRUE)
  # #plotjaccard(sc)
  # 
  # nb.subcluster = sc@cluster$clb$nc
  # cat('optimal subclusters found :', nb.subcluster, '\n')
  # 
  # sc <- clustexp(sc,cln=nb.subcluster, sat=FALSE, verbose = FALSE)
  # sub.obj$seurat_clusters_split = sc@cluster$kpart
  # 
  p1  = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 3)
  
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3, 
               label.size = 6,
               na.value = "gray", combine = TRUE)
  
  p1 + p2
  
  VlnPlot(sub.obj, features = c('hnd-1', 'pha-4', 'sdz-1', 'sdz-31', 'ceh-51', 'ceh-36', 'unc-39', 'unc-120'),
          group.by = 'seurat_clusters_split')
  
  #sub.obj$timingEst = as.numeric(sub.obj$timingEst)
  VlnPlot(sub.obj, features = c('timingEst', "FSC_log2", "BSC_log2"), ncol = 3, 
          group.by = 'seurat_clusters_split')
  
  # p3 = DimPlot(sub.obj, group.by = "timingEst", reduction = 'umap', label = FALSE, repel = FALSE, pt.size = 3, 
  #         label.size = 6,
  #         na.value = "gray") 
  # 
  # par(mfrow=c(1, 1))
  # counts <- table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
  # counts = counts[!is.na(match(rownames(counts), c(ids, 'unassigned'))), ]
  # 
  # barplot(counts, main="composition of subclusters ",
  #         xlab="subcluster index", col=c(1:nrow(counts)),
  #         legend = rownames(counts))
  # 
  # marker genes
  Idents(sub.obj) = sub.obj$seurat_clusters_split 
  markers <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
  
  top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  
  #Idents(seurat.cistopic) = $lineage
  DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend() 
  
  
  ##########################################
  # check the update of manually annotated ids
  ##########################################
  seurat.obj$manual.annot.ids[1] = 'unknown'
  
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  
  
  ########################################################
  ########################################################
  # Section : manual annotation for c('29', '32', '35', '40', '42')
  # 
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
  
  #ids = c('MSx')
  ids = c('MSx', 'MSxa', 'MSxap', 'MSxapp','MSxappa', 'MSxappp', 'MSxp', 'MSxpa', 'MSxpp')
  ids = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSpappa', 'MSxp', 'MSxpa', 'MSxpp')
  
  # given the fact that scmap seems to be working better than seurat 
  # for the manual annotation, the scmap predicted labels will be mainly considered, 
  # sstart with scmap projection with nb.features = 500 and threshold = 0.7
  # seurat.obj$predicted.ids = seurat.obj$scmap.tintori.id
  # seurat.obj$predicted.scores = seurat.obj$scmap.tintori.cor
  # threshold = 0.7
  seurat.obj$predicted.ids = seurat.obj$scmap.pred.id.500
  seurat.obj$predicted.scores = seurat.obj$scmap.corr.500
  threshold = 0.7
  
  seurat.obj$predicted.ids[which(seurat.obj$predicted.ids == 'unassigned')] = NA
  seurat.obj$predicted.ids.fitered = seurat.obj$predicted.ids
  seurat.obj$predicted.ids.fitered[which(seurat.obj$predicted.scores < threshold)] = NA;
  
  ##########################################
  # check the predicted labels and cluster-label mapping without filtering
  ##########################################
  cells.sels = which(!is.na(match(seurat.obj$predicted.ids, ids)))
  cells.fitered = cells.sels[which(seurat.obj$predicted.scores[cells.sels] > threshold)]
  #cat(length(cells.sels), ' cells with predicted labels \n')
  cat(length(cells.fitered), ' cells after filtering with predicted labels \n')
  #cat('clusters involved : '); print(as.character(unique(seurat.obj$seurat_clusters[cells.sels])));
  cluster.sels =  unique(seurat.obj$seurat_clusters[cells.fitered])
  cat('clusters involved after fitering : '); print(as.character(cluster.sels));
  
  p0 = DimPlot(seurat.obj,
               cells.highlight = colnames(seurat.obj)[cells.fitered],
               group.by = "predicted.ids.fitered", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 4,
               sizes.highlight = 1,
               na.value = "gray") + 
    ggtitle(paste0("cells predicted by scmap with 500 features and threshold 0.7 ")) +
    NoLegend()
  
  p00 = DimPlot(seurat.obj,
                cells.highlight = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))],
                group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 4,
                sizes.highlight = 1,
                na.value = "gray") +
    #ggtitle(paste0("cells predicted by scmap with 500 features and threshold 0.7 ")) +
    NoLegend()
  
  p0 + p00
  
  ##########################################
  # found involved clusters to consider in the following 
  ##########################################
  # make the count table for ids considered and involved cluster index
  predicted.ids = seurat.obj$predicted.ids.fitered
  predicted.ids[is.na(predicted.ids)] = 'unassigned'
  counts <- table(predicted.ids, seurat.obj$seurat_clusters)
  counts = counts[grep('MS|unassigned', rownames(counts)), ]
  
  nb.cells.clusters = apply(counts[rownames(counts) != 'unassigned', ], 2, sum)
  
  counts = counts[!is.na(match(rownames(counts), c(bwms, 'unassigned'))), ]
  nb.cells.bwm.clusters = apply(counts[rownames(counts) != 'unassigned',], 2, sum)
  
  counts = counts[!is.na(match(rownames(counts), c(ids, 'unassigned'))), ]
  nb.cells.ids.clusters = apply(counts[rownames(counts) != 'unassigned',], 2, sum)
  
  kk = which(nb.cells.ids.clusters >0 ) # > 0 cell in the cluster
  cat(length(kk), 'clusters for considered ids \n')
  print(colnames(counts)[kk])
  
  cols = c(colorRampPalette((brewer.pal(n = 7, name="Reds")))(max(counts)))
  pheatmap(counts, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, breaks = seq(0, max(counts), by = 1),
           cluster_cols=FALSE, main = paste0("cluster -- predicted labels mapping"), na_col = "white",
           color = cols, 
           #display_numbers = TRUE,
           #number_format = "%.0f",
           #annotation_col = my_sample_col,
           #gaps_row = c(1:nrow(map)-1),
           fontsize_col = 10,
           height = 8,
           width = 30
  )
  
  ## filter clusters that are not likely for the considered ids due to noise
  rr.bwm = (nb.cells.bwm.clusters/nb.cells.clusters)[kk]
  barplot(rr.bwm, main = '% of predited cells found in BWM');abline(h=0.5, col ='red')
  rr.ids = (nb.cells.ids.clusters/nb.cells.clusters)[kk]
  
  kk = which(nb.cells.ids.clusters > 0 
             & nb.cells.ids.clusters/nb.cells.bwm.clusters > 0.5 # >50% of bwm cells have the ids considered 
             & nb.cells.bwm.clusters/nb.cells.clusters > 0.5 # > 50% of cell predicted to be in bwm
  ) 
  
  counts = counts[, kk]
  
  par(mfrow = c(1, 2))
  barplot(t(counts), main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = colnames(counts)
  )
  
  barplot((counts), main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))
  
  ##########################################
  # split the chosen clusters and manual annotate them
  ##########################################
  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_3.pdf")
  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #cluster.sels = colnames(counts)
  #cluster.sels = c('29', '32',  '40', '42', '6', '14', '44', '4')
  cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('29')
  features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67',
                    'unc-120', 'unc-39', 'irx-1', 'egl-43')
  
  sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  
  DimPlot(sub.obj, reduction = 'umap', label = TRUE)
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = 'manual.annot.ids')
  
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
  
  #length(intersect(VariableFeatures(sub.obj), timers))
  #VariableFeatures(sub.obj) = setdiff(VariableFeatures(sub.obj), timers)
  cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 50)
  
  
  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids 
  n.neighbors = 10; min.dist = 0.1;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, n.neighbors = n.neighbors, 
                     min.dist = min.dist)
  p0 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1
  
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
  
  p1  = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2)
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, 
               label.size = 5,
               na.value = "gray", combine = TRUE)
  p1 + p2
  
  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3, 
          group.by = 'seurat_clusters_split')
  
  VlnPlot(sub.obj, features = features.sels,
          group.by = 'seurat_clusters_split')
  
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  dev.off()
  
  ##########################################
  # update of manually annotated ids
  ##########################################
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '7')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSx'
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxa/p'
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '4')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxap'
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '5')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxap'
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '6')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpp'
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '2')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpa'
  
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  seurat.obj = split.cluster.with.specific.gene.exrepssion(seurat.obj)
  
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'manual.annot.ids')
  
  #VlnPlot(seurat.obj, features = c('hnd-1', 'pha-4', 'ceh-76', 'fbxb-70'), group.by = 'manual.annot.ids')
  
  
  saveRDS(seurat.obj, 
          file = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_1.rds'))
  
  
  ########################################################
  ########################################################
  # Section : annoate parts of cluster 4 and 22 (convergence branches)
  # 
  ########################################################
  ########################################################
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
  
  
  ########################################################
  ########################################################
  # Section : have more MSapp cells from cluster 3, 24, 6, 14
  # 
  ########################################################
  ########################################################
  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxapp_4.pdf")
  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  seurat.obj = readRDS(
    file = paste0(RdataDir, 
                  'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_2.rds'))
  
  #cluster.sels = colnames(counts)
  #cluster.sels = c('29', '32',  '40', '42', '6', '14', '44', '4')
  #cluster.sels = c('29', '32', '35', '40', '42')
  cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4',  '22', '24') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('14', '6', '3')
  cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) | 
                                      seurat.obj$manual.annot.ids == 'MSxapp'| seurat.obj$manual.annot.ids == 'MSxap']
  
  sub.obj = subset(seurat.obj, cells = cells.sels)
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  
  DimPlot(sub.obj, reduction = 'umap', label = TRUE)
  
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67', 'pat-4'))
  
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
  n.neighbors = 10; min.dist = 0.2;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, n.neighbors = n.neighbors, 
                     min.dist = min.dist)
  p0 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1
  
  DimPlot(sub.obj, reduction = 'umap', group.by = 'manual.annot.ids')
  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.4)
  
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
  
  features.sels = c( 'nhr-67', 'pat-9', 'ceh-37', 'C45G7.4', 'irx-1', 'egl-43')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  p3 = VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = c('1', '6', '3')
  )
  
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
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '3')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxap'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxapp'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '6')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxapp'
  
  # cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '4')]
  # seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSpappa'
  # 
  # cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  # seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxappp'
  # 
  # cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  # seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxappp/Mspappa.daugther'
  
  
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)
  
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'manual.annot.ids')
  
  #VlnPlot(seurat.obj, features = c('hnd-1', 'pha-4', 'ceh-76', 'fbxb-70'), group.by = 'manual.annot.ids')
  
  
  saveRDS(seurat.obj, 
          file = paste0(RdataDir, 
                        'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_3.rds'))
  
  
  ########################################################
  ########################################################
  # Section : confirm manually anotated MSxa and MSxap form good clusters in our data 
  # 
  ########################################################
  ########################################################
  seurat.obj = readRDS(
    file = paste0(RdataDir, 
                  'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_3.rds'))
  
  
  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxapp_5.pdf")
  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #cluster.sels = colnames(counts)
  #cluster.sels = c('29', '32',  '40', '42', '6', '14', '44', '4')
  #cluster.sels = c('29', '32', '35', '40', '42')
  cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4',  '22', '24') # cluster 4 and 22 have some cells in Msxp lineage
  #cluster.sels = c('14', '6', '3')
  cells.sels = colnames(seurat.obj)[seurat.obj$manual.annot.ids == 'MSxapp' |
                                      seurat.obj$manual.annot.ids == 'MSxappp'|
                                      seurat.obj$manual.annot.ids == 'MSpappa'|
                                      seurat.obj$manual.annot.ids == 'MSxappp/Mspappa.daugther']
  
  sub.obj = subset(seurat.obj, cells = cells.sels)
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = 'manual.annot.ids')
  
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67', 'pat-4'))
  
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
  ##########################################
  # check potential ids for selected clusters 
  ##########################################
  threshold = 0.7
  predicted.ids = sub.obj$scmap.pred.id.500
  selected.clusters = as.character(sub.obj$manual.annot.ids)
  counts = table(predicted.ids, selected.clusters)
  
  predicted.ids[which(sub.obj$scmap.corr.500 <0.7)] = 'unassigned'
  counts = table(predicted.ids, selected.clusters)
  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )
  
  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 2000)
  
  #length(intersect(VariableFeatures(sub.obj), timers))
  #VariableFeatures(sub.obj) = setdiff(VariableFeatures(sub.obj), timers)
  cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 50)
  
  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids 
  n.neighbors = 20; min.dist = 0.1;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, n.neighbors = n.neighbors, 
                     min.dist = min.dist)
  p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5)
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
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
  
  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
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
  
  features.sels = c( 'nhr-67', 'pat-9', 'ceh-37', 'hlh-1', 'tbx-8', 'lin-39')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  p3 = VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split'
  )
  
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
  
  
  
  ########################################################
  ########################################################
  # Section : clean some MSxapp cells 
  # 
  ########################################################
  ########################################################
  seurat.obj = readRDS(
    file = paste0(RdataDir, 
                  'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_3.rds'))
  
  
  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxapp_5.pdf")
  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #cluster.sels = colnames(counts)
  #cluster.sels = c('29', '32',  '40', '42', '6', '14', '44', '4')
  #cluster.sels = c('29', '32', '35', '40', '42')
  cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4',  '22', '24') # cluster 4 and 22 have some cells in Msxp lineage
  #cluster.sels = c('14', '6', '3')
  cells.sels = colnames(seurat.obj)[seurat.obj$manual.annot.ids == 'MSxapp' |
                                      seurat.obj$manual.annot.ids == 'MSxappp'|
                                      seurat.obj$manual.annot.ids == 'MSpappa'|
                                      seurat.obj$manual.annot.ids == 'MSxappp/Mspappa.daugther']
  
  sub.obj = subset(seurat.obj, cells = cells.sels)
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = 'manual.annot.ids')
  
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67', 'pat-4'))
  
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
  ##########################################
  # check potential ids for selected clusters 
  ##########################################
  threshold = 0.7
  predicted.ids = sub.obj$scmap.pred.id.500
  selected.clusters = as.character(sub.obj$manual.annot.ids)
  counts = table(predicted.ids, selected.clusters)
  
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
  n.neighbors = 10; min.dist = 0.1;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, n.neighbors = n.neighbors, 
                     min.dist = min.dist)
  p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1
  
  FeaturePlot(sub.obj, reduction = 'umap', features = c('nhr-67', 'hnd-1', 'pat-9'))
  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
  
  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
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
  
  features.sels = c( 'nhr-67', 'pat-9', 'ceh-37', 'hlh-1', 'tbx-8', 'lin-39')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  p3 = VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split'
  )
  
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
  ##########################################
  # update of manually annotated ids
  ##########################################
  ##########################################
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '3')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = NA
  #cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  #seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxapp'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '6')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = NA
  
  # cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '4')]
  # seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSpappa'
  # 
  # cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  # seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxappp'
  # 
  # cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  # seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxappp/Mspappa.daugther'
  
  
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)
  
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'manual.annot.ids')
  
  #VlnPlot(seurat.obj, features = c('hnd-1', 'pha-4', 'ceh-76', 'fbxb-70'), group.by = 'manual.annot.ids')
  
  
  saveRDS(seurat.obj, 
          file = paste0(RdataDir, 
                        'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_4.rds'))
  
  
  
  ########################################################
  ########################################################
  # Section : refine MSxapp convergence branche
  # 
  ########################################################
  ########################################################
  
  seurat.obj = readRDS(
    file = paste0(RdataDir, 
                  'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_4.rds'))
  
  
  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxapp_6.pdf")
  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #cluster.sels = colnames(counts)
  #cluster.sels = c('29', '32',  '40', '42', '6', '14', '44', '4')
  #cluster.sels = c('29', '32', '35', '40', '42')
  cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4',  '22', '24') # cluster 4 and 22 have some cells in Msxp lineage
  #cluster.sels = c('14', '6', '3')
  cells.sels = colnames(seurat.obj)[seurat.obj$manual.annot.ids == 'MSxapp' |
                                      seurat.obj$manual.annot.ids == 'MSxappp'|
                                      seurat.obj$manual.annot.ids == 'MSpappa'|
                                      seurat.obj$manual.annot.ids == 'MSxappp/Mspappa.daugther']
  
  sub.obj = subset(seurat.obj, cells = cells.sels)
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = 'manual.annot.ids')
  
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67', 'pat-4'))
  
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
  ##########################################
  # check potential ids for selected clusters 
  ##########################################
  threshold = 0.7
  predicted.ids = sub.obj$scmap.pred.id.500
  selected.clusters = as.character(sub.obj$manual.annot.ids)
  counts = table(predicted.ids, selected.clusters)
  
  predicted.ids[which(sub.obj$scmap.corr.500 <0.7)] = 'unassigned'
  counts = table(predicted.ids, selected.clusters)
  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )
  
  counts.seurat = table(sub.obj$seurat.pred.id, sub.obj$manual.annot.ids)
  
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
  n.neighbors = 10; min.dist = 0.05;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, n.neighbors = n.neighbors, 
                     min.dist = min.dist)
  p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1
  
  FeaturePlot(sub.obj, reduction = 'umap', features = c('nhr-67', 'hnd-1', 'pat-9'))
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
  
  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2)
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, 
               label.size = 5,
               na.value = "gray", combine = TRUE)
  p1 + p2
  
  manual.discovery.new.features = TRUE
  if(manual.discovery.new.features){
    Idents(sub.obj) = sub.obj$seurat_clusters_split
    markers <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
    
  }
  
  features.sels = c( 'nhr-67', 'pat-9', 'hlh-1', 'tbx-8', 'rpm-1', 'tab-1', 'gana-1', 'spp-15', 'F37H8.5', 'tnt-3',
                     'asic-2', 'skpo-1')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  p3 = VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split'
  )
  
  p2 + p3
  
  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3, 
          group.by = 'seurat_clusters_split')
  
  counts = table(predicted.ids, as.character(sub.obj$seurat_clusters_split))
  
  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )
  
  counts = table(sub.obj$seurat.pred.id, sub.obj$manual.annot.ids)
  
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  
  dev.off()
  
  
  ##########################################
  ##########################################
  # update of manually annotated ids
  ##########################################
  ##########################################
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '2')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxapppa'
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '3')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxappp'
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '4')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxapp'
  
  
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSpappa'
  
  ## cluster 5 and 1 are likely to be other cell ids, but unclear so keep them as MSxapp
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '5')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxapp'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxapp'
  
  # cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '4')]
  # seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSpappa'
  # 
  # cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  # seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxappp'
  # 
  # cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  # seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxappp/Mspappa.daugther'
  
  
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)
  
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'manual.annot.ids')
  
  #VlnPlot(seurat.obj, features = c('hnd-1', 'pha-4', 'ceh-76', 'fbxb-70'), group.by = 'manual.annot.ids')
  
  
  saveRDS(seurat.obj, 
          file = paste0(RdataDir, 
                        'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_5.rds'))
  
  
  
}





