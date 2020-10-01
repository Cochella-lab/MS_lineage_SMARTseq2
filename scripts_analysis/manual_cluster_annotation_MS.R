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


  ########################################################
  ########################################################
  # Section : iteration 6
  # first annotation of MSxp for middle time points
  ########################################################
  ########################################################
  seurat.obj = readRDS(
    file = paste0(RdataDir,
                  'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_5.rds'))


  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_7.pdf")
  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

  #cluster.sels = colnames(counts)
  #cluster.sels = c('29', '32',  '40', '42', '6', '14', '44', '4')
  #cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4', '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage

  cells.sels = colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                     & is.na(seurat.obj$manual.annot.ids))]

  sub.obj = subset(seurat.obj, cells = cells.sels)
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = 'seurat_clusters')

  #FeaturePlot(sub.obj, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67', 'pat-4'))

  #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
  ##########################################
  # check potential ids for selected clusters
  ##########################################
  threshold = 0.7
  predicted.ids = sub.obj$scmap.pred.id.500
  selected.clusters = as.character(sub.obj$seurat_clusters)
  counts = table(predicted.ids, selected.clusters)

  predicted.ids[which(sub.obj$scmap.corr.500 <0.7)] = 'unassigned'
  counts = table(predicted.ids, selected.clusters)
  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )

  counts.seurat = table(sub.obj$seurat.pred.id, sub.obj$seurat_clusters)

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
  p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1

  FeaturePlot(sub.obj, reduction = 'umap', features = c('hnd-1'))
  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }

  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)

  p1  = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
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

  counts = table(sub.obj$seurat.pred.id, sub.obj$seurat_clusters_split)

  FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'hnd-1', 'hlh-1', 'abts-1', 'ref-2', 'tbx-7',
                                                        'unc-39', 'cup-4', 'ins-2', 'F40H3.3', 'hot-1'))

  dev.off()

  ##########################################
  ##########################################
  # update of manually annotated ids
  ##########################################
  ##########################################

  ##########################################
  # repeat marker genes to annote obtained clusters
  ##########################################
  #features.sels = c('unc-120', 'hnd-1', 'hlh-1', 'abts-1', 'ref-2', 'tbx-7', 'unc-39', 'cup-4', 'ins-2', 'F40H3.3', 'hot-1')
  features.sels = c('hnd-1', 'tbx-7', 'unc-120', 'abts-1', 'Y66D12A.13', 'F41D9.2') # markers for MSxpaaa
  features.sels = c('unc-39', 'ref-2', 'tbx-7', 'unc-120') # markers for MSxpaa
  features.sels = c('hlh-1', 'tbx-8', 'hnd-1', 'unc-39', 'ref-2', 'unc-120')
  features.sels = c('tbx-8', 'hlh-1', 'ZK180.5', 'F19C6.4', 'zig-8', 'hnd-1', 'unc-39', 'unc-120', 'F40H3.3', 'ins-2', 'Y42H9B.3')

  features.sels = c('col-118', 'unc-120', 'C03B1.1', 'F53C3.7')

  features.sels = c('hlh-1', 'tbx-8', 'unc-120', 'F19C6.4')
  features.sels = c('hlh-1', 'ins-2', 'sul-2', 'camt-1', 'tbx-8')

  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(c(0:12))
  idents.sel = setdiff(idents.sel, c('1', '10', '4', '12', '3', '6', '11', '0', '2', '5', '9'))
  VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '12')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxapp'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '4')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpaaa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpaa/MSxpaap'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '10')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpaa/MSxpaap'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '3')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpap'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '6')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpapp'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '11')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpapa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '9')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '5')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppp'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '2')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpppa/MSxpppp'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '7')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpppa/MSxpppp daughters'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '8')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpppa/MSxpppp daughters'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppaa/MSxppap'



  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)

  #DimPlot(sub.obj, reduction = 'umap', group.by = 'manual.annot.ids')

  #VlnPlot(seurat.obj, features = c('hnd-1', 'pha-4', 'ceh-76', 'fbxb-70'), group.by = 'manual.annot.ids')

  saveRDS(seurat.obj,
          file = paste0(RdataDir,
                        'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_6.rds'))



  ########################################################
  ########################################################
  # Section : iteration 7
  #  refine MSxpaa, MSxpaaa and MSxpaap and MSxpap and its daughters
  ########################################################
  ########################################################
  seurat.obj = readRDS(
    file = paste0(RdataDir,
                  'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_6.rds'))


  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_7.pdf")
  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

  #cluster.sels = colnames(counts)
  #cluster.sels = c('29', '32',  '40', '42', '6', '14', '44', '4')
  #cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  #cluster.sels = c('4', '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  #cells.sels = colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
  #                                   & is.na(seurat.obj$manual.annot.ids))]

  cells.sels = colnames(seurat.obj)[seurat.obj$manual.annot.ids == 'MSxpaaa' |
                                      seurat.obj$manual.annot.ids == 'MSxpaa/MSxpaap' |
                                      seurat.obj$manual.annot.ids == 'MSxpap'|
                                      seurat.obj$manual.annot.ids == 'MSxpapa'|
                                      seurat.obj$manual.annot.ids == 'MSxpapp']
  Refine.annotated.ids = TRUE;
  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}

  sub.obj = subset(seurat.obj, cells = cells.sels)
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

  #FeaturePlot(sub.obj, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67', 'pat-4'))

  ##########################################
  # check potential ids for selected clusters
  ##########################################
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
  threshold = 0.7
  predicted.ids = sub.obj$scmap.pred.id.500
  predicted.ids[which(sub.obj$scmap.corr.500 <0.7)] = 'unassigned'
  if(Refine.annotated.ids){
    counts = table(predicted.ids, sub.obj$manual.annot.ids)
    counts.seurat = table(sub.obj$seurat.pred.id, sub.obj$manual.annot.ids)
  }else{
    counts = table(predicted.ids, as.character(sub.obj$seurat_clusters))
    counts.seurat = table(sub.obj$seurat.pred.id, sub.obj$seurat_clusters)
  }
  barplot(counts, main="cluster compositions by scmap ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))

  barplot(counts.seurat, main="cluster compositions by seurat ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))

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
  p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1

  FeaturePlot(sub.obj, reduction = 'umap', features = c('hnd-1'))
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

  p1  = DimPlot(sub.obj, group.by = by.group, reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
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

  #p3 = VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split')
  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3,
               group.by = 'seurat_clusters_split')

  (p1 + p2) / p3

  counts = table(predicted.ids, as.character(sub.obj$seurat_clusters_split))
  counts.seurat = table(sub.obj$seurat.pred.id, sub.obj$seurat_clusters_split)

  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )

  dev.off()

  ##########################################
  ##########################################
  # update of manually annotated ids
  ##########################################
  ##########################################

  ##########################################
  # repeat marker genes to annote obtained clusters
  ##########################################
  #features.sels = c('unc-120', 'hnd-1', 'hlh-1', 'abts-1', 'ref-2', 'tbx-7', 'unc-39', 'cup-4', 'ins-2', 'F40H3.3', 'hot-1')
  features.sels = c('hnd-1', 'tbx-7', 'unc-120', 'abts-1', 'Y66D12A.13', 'F41D9.2') # markers for MSxpaaa
  features.sels = c('unc-39', 'ref-2', 'tbx-7', 'unc-120') # markers for MSxpaa
  features.sels = c('hlh-1', 'tbx-8', 'hnd-1', 'unc-39', 'ref-2', 'unc-120')
  features.sels = c('tbx-8', 'hlh-1', 'ZK180.5', 'F19C6.4', 'zig-8', 'hnd-1', 'unc-39', 'unc-120', 'F40H3.3', 'ins-2', 'Y42H9B.3')

  features.sels = c('col-118', 'unc-120', 'C03B1.1', 'F53C3.7')

  features.sels = c('hlh-1', 'tbx-8', 'unc-120', 'F19C6.4')
  features.sels = c('hlh-1', 'ins-2', 'sul-2', 'camt-1', 'tbx-8')

  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(c(0:12))
  idents.sel = setdiff(idents.sel, c('1', '10', '4', '12', '3', '6', '11', '0', '2', '5', '9'))
  VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)


  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '4')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpaaa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '5')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpapa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpapp'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpaap'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '2')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpaa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '3')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpap'

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)

  #DimPlot(sub.obj, reduction = 'umap', group.by = 'manual.annot.ids')

  #VlnPlot(seurat.obj, features = c('hnd-1', 'pha-4', 'ceh-76', 'fbxb-70'), group.by = 'manual.annot.ids')

  saveRDS(seurat.obj,
          file = paste0(RdataDir,
                        'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_7.rds'))


  ########################################################
  ########################################################
  # Section : iteration 9
  #  try to find MSppaap and MSapaap instead of symetric MSxpaap
  ########################################################
  ########################################################
  nb.iteration = 9
  Refine.annotated.ids = TRUE;

  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_',
                    nb.iteration -1, '.rds')
  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
  RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_',
                     nb.iteration, '.rds')

  seurat.obj = readRDS(file = RDSsaved)
  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}

  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  #cluster.sels = colnames(counts)
  #cluster.sels = c('29', '32',  '40', '42', '6', '14', '44', '4')
  #cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  #cluster.sels = c('4', '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  #cells.sels = colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
  #                                   & is.na(seurat.obj$manual.annot.ids))]

  cells.sels = colnames(seurat.obj)[seurat.obj$manual.annot.ids == 'MSxppp' |
                                      seurat.obj$manual.annot.ids == 'MSxppa' |
                                      seurat.obj$manual.annot.ids == 'MSxpppa/MSxpppp'|
                                      seurat.obj$manual.annot.ids == 'MSxppaa'|
                                      seurat.obj$manual.annot.ids == 'MSxppaa/MSxppap' |
                                      seurat.obj$manual.annot.ids == 'MSxpppa/MSxpppp daughters']

  sub.obj = subset(seurat.obj, cells = cells.sels)
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

  #FeaturePlot(sub.obj, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67', 'pat-4'))

  ##########################################
  # check potential ids for selected clusters
  ##########################################
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
  threshold = 0.7
  predicted.ids = sub.obj$scmap.pred.id.500
  predicted.ids[which(sub.obj$scmap.corr.500 <0.7)] = 'unassigned'
  if(Refine.annotated.ids){
    counts = table(predicted.ids, sub.obj$manual.annot.ids)
    counts.seurat = table(sub.obj$seurat.pred.id, sub.obj$manual.annot.ids)
  }else{
    counts = table(predicted.ids, as.character(sub.obj$seurat_clusters))
    counts.seurat = table(sub.obj$seurat.pred.id, sub.obj$seurat_clusters)
  }
  barplot(counts, main="cluster compositions by scmap ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))

  barplot(counts.seurat, main="cluster compositions by seurat ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))

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
  p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1

  #FeaturePlot(sub.obj, reduction = 'umap', features = c('hnd-1'))
  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.8)

  p1  = DimPlot(sub.obj, group.by = by.group, reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
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

  #p3 = VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split')
  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3,
               group.by = 'seurat_clusters_split')

  (p1 + p2) / p3

  counts = table(predicted.ids, as.character(sub.obj$seurat_clusters_split))
  counts.seurat = table(sub.obj$seurat.pred.id, sub.obj$seurat_clusters_split)

  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )

  dev.off()

  ##########################################
  ##########################################
  # update of manually annotated ids
  ##########################################
  ##########################################

  ##########################################
  # repeat marker genes to annote obtained clusters
  ##########################################
  #features.sels = c('unc-120', 'hnd-1', 'hlh-1', 'abts-1', 'ref-2', 'tbx-7', 'unc-39', 'cup-4', 'ins-2', 'F40H3.3', 'hot-1')
  features.sels = c('hnd-1', 'tbx-7', 'unc-120', 'abts-1', 'Y66D12A.13', 'F41D9.2') # markers for MSxpaaa
  features.sels = c('unc-39', 'ref-2', 'tbx-7', 'unc-120') # markers for MSxpaa
  features.sels = c('hlh-1', 'tbx-8', 'hnd-1', 'unc-39', 'ref-2', 'unc-120')
  features.sels = c('tbx-8', 'hlh-1', 'ZK180.5', 'F19C6.4', 'zig-8', 'hnd-1', 'unc-39', 'unc-120', 'F40H3.3', 'ins-2', 'Y42H9B.3')

  features.sels = c('col-118', 'unc-120', 'C03B1.1', 'F53C3.7')

  features.sels = c('hlh-1', 'ins-2', 'sul-2', 'camt-1', 'tbx-8')

  features.sels = c('hlh-1', 'unc-120', 'F19C6.4', 'col-118', 'F53C3.7', 'C03B1.1')


  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(c(0:6))
  idents.sel = setdiff(idents.sel, c('0', '4', '1', '6'))
  VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)


  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppp'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '4')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpppa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '6')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpppp'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '5')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppaa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '2')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppap'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '2')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppap'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '3')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppppx/MSxpppax'

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)


  saveRDS(seurat.obj, file = RDS2save)


  ########################################################
  ########################################################
  # Section : iteration 9
  # Continue with terminal cells that are in well separated clusters
  ########################################################
  ########################################################
  nb.iteration = 9
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

  #cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4',  '22', '24',  # middle time point for convergence branche
                   '3', '5', '16', '30', # middle time points for MSxp
                   '23', '43', '17', '41', '45', '31', '34', '50', '53', # well-separated clusters with clear mapped ids
                   #'23', '43', '17', '41', '45', '34', '31', '34', '50', '53', # well-separated clusters without clear mapped ids
                   '44', '28', '52', # unknow yet
                   '25', '36', '8', '39', '2', '19', '27', # BWM-1 cluster_25 transition to it and possibly cluster_28, 51 transition to BWM_2
                   '24', '13', '1', '11', '33', '48', '18', '46', '15', '26' # BWM-2
  )

  cluster.sels = c('23', '43', '17', '41', '45', '34', '31', '50', '53' # well-separated clusters with clear mapped ids
  )
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]

  cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
  #                                   & is.na(seurat.obj$manual.annot.ids))]
  # cells.sels = colnames(seurat.obj)[seurat.obj$manual.annot.ids == 'MSxpaaa' |
  #                                     seurat.obj$manual.annot.ids == 'MSxpaa' |
  #                                     seurat.obj$manual.annot.ids == 'MSxpaap' |
  #                                     seurat.obj$manual.annot.ids == 'MSxpap'|
  #                                     seurat.obj$manual.annot.ids == 'MSxpapa'|
  #                                     seurat.obj$manual.annot.ids == 'MSxpapp']
  #
  sub.obj = subset(seurat.obj, cells = cells.sels)
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

  DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500', label = FALSE)
  #p2 = DimPlot(sub.obj, reduction = 'umap', group.by = 'seurat.pred.id', label = FALSE)
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67', 'pat-4'))

  ##########################################
  # check potential ids for selected clusters
  ##########################################
  #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
  threshold = 0.7
  predicted.ids = sub.obj$scmap.pred.id.500
  predicted.ids[which(sub.obj$scmap.corr.500 <0.7)] = 'unassigned'
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
  DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)

  DimPlot(sub.obj, group.by = 'scmap.pred.id.500', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) +
    NoLegend()

  p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1

  #FeaturePlot(sub.obj, reduction = 'umap', features = c('hnd-1'))
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

  p1  = DimPlot(sub.obj, group.by = by.group, reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2)
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 5,
               na.value = "gray", combine = TRUE)
  p1 + p2


  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3,
               group.by = 'seurat_clusters_split')

  (p1 + p2) / p3

  counts = table(sub.obj$scmap.pred.id.500, as.character(sub.obj$seurat_clusters_split))
  counts.seurat = table(sub.obj$seurat.pred.id, sub.obj$seurat_clusters_split)

  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )

  manual.discovery.new.features = TRUE
  if(manual.discovery.new.features){
    Idents(sub.obj) = sub.obj$seurat_clusters_split
    markers <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

  }


  dev.off()

  ##########################################
  ##########################################
  # update of manually annotated ids
  ##########################################
  ##########################################

  ##########################################
  # repeat marker genes to annote obtained clusters
  ##########################################
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(c(0:9))
  idents.sel = setdiff(idents.sel, c('5', '0', '1', '4', '3', '2', '9', '8'))

  ## chcek the reference-mapped ids for the rest of clusters
  counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
  counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
  counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
  counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]

  ids.sel = c('MSxapapp', 'MSxppaaa', 'MSxapppax')
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

  #features.sels = c('unc-120', 'hnd-1', 'hlh-1', 'abts-1', 'ref-2', 'tbx-7', 'unc-39', 'cup-4', 'ins-2', 'F40H3.3', 'hot-1')
  features.sels = c('pha-4', 'unc-62', 'C45G7.4', 'asp-4', 'cup-4', 'lgc-23', 'hlh-1', 'lin-39', 'tnt-3', 'gana-1', 'spp-15', 'F37H8.5')
  VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)

  # to find new marker genes
  DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '5')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppaap'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppaaa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpapaa'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '4')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpapaa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '3')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSapaapp'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '2')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpaaaa/MSxpaaaax'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '9')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'weirdos.possible.doublets'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '8')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxapapp.pharynx'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '6')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'unknown.MSxapppax.MSxpppaa.MSxpppaa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '7')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'unknown.MSxpppaa'


  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)

  FeaturePlot(seurat.obj, reduction = 'umap', features = 'pha-4')

  saveRDS(seurat.obj, file = RDS2save)

  ########################################################
  ########################################################
  # Section : iteration 10.1
  # manual annotate terminal cells BWM_1 without transition
  # but fail to find common marker genes between JM data and our data, possibly due to the sparsity of terminal cells
  ########################################################
  ########################################################
  nb.iteration = 10
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

  #cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4',  '22', '24',  # middle time point for convergence branche
                   '3', '5', '16', '30', # middle time points for MSxp
                   '23', '43', '17', '41', '45', '34', '50', # well-separated clusters with clear mapped ids
                   '44', '28', '52', '31', '53',  # # well-separated clusters without clear mapped ids
                   '25', '36', '8', '39', '2', '19', '27', # BWM-1 cluster_25 transition to it and possibly cluster_28, 51 transition to BWM_2
                   '24', '13', '1', '11', '33', '48', '18', '46', '15', '26' # BWM-2
  )

  cluster.sels = c('22','30', '25', '36', '8', '39', '2', '19', '27', # BWM-1 without transition cluster_25
                   '24', '13', '1', '11', '33', '48', '18', '46', '15', '26' # BWM-2
  )

  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]

  cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
  #                                   & is.na(seurat.obj$manual.annot.ids))]
  # cells.sels = colnames(seurat.obj)[seurat.obj$manual.annot.ids == 'MSxpaaa' |
  #                                     seurat.obj$manual.annot.ids == 'MSxpaa' |
  #                                     seurat.obj$manual.annot.ids == 'MSxpaap' |
  #                                     seurat.obj$manual.annot.ids == 'MSxpap'|
  #                                     seurat.obj$manual.annot.ids == 'MSxpapa'|
  #                                     seurat.obj$manual.annot.ids == 'MSxpapp']
  #
  sub.obj = subset(seurat.obj, cells = cells.sels)
  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

  DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500', label = FALSE)
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
  n.neighbors = 30; min.dist = 0.1;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, n.neighbors = n.neighbors,
                     min.dist = min.dist)
  DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)

  DimPlot(sub.obj, group.by = 'scmap.pred.id.500', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) +
    NoLegend()

  p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1

  #FeaturePlot(sub.obj, reduction = 'umap', features = c('hnd-1'))
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

  p1  = DimPlot(sub.obj, group.by = by.group, reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2)
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 5,
               na.value = "gray", combine = TRUE)
  p1 + p2


  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3,
               group.by = 'seurat_clusters_split')

  (p1 + p2) / p3

  counts = table(sub.obj$scmap.pred.id.500, as.character(sub.obj$seurat_clusters_split))
  counts.seurat = table(sub.obj$seurat.pred.id, sub.obj$seurat_clusters_split)

  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts)
  )

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
  ##################################################################################################################################
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
  #idents.sel = setdiff(idents.sel, c('5', '0', '1', '4', '3', '2', '9', '8'))

  ## chcek the reference-mapped ids for the rest of clusters
  counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
  counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
  counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
  counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]

  ids.sel = c('MSxpppax')
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)

  #features.sels = c('unc-120', 'hnd-1', 'hlh-1', 'abts-1', 'ref-2', 'tbx-7', 'unc-39', 'cup-4', 'ins-2', 'F40H3.3', 'hot-1')
  features.sels = c('unc-120', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'kvs-5')

  VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)

  # to find new marker genes
  DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)


  ########################################################
  ########################################################
  # Section : iteration 10.5
  # first try of manual annotate terminal cells by considering all BDW terminal and their mothers
  # we learnt that to select specific lineages in reference and to use variable gene for that time points help a lot
  # for reference-based annotation, in particular seurat works well
  # at the end of the this section, we quickly annotated potential BWM temrinal cells
  ########################################################
  ########################################################
  nb.iteration = 11
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

  #cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4',  '22', '24',  # middle time point for convergence branche
                   '3', '5', '16', '30', # middle time points for MSxp
                   '23', '43', '17', '41', '45', '34', '50', # well-separated clusters with clear mapped ids
                   '44', '28', '52', '31', '53',  # # well-separated clusters without clear mapped ids
                   '25', '36', '8', '39', '2', '19', '27', # BWM-1 cluster_25 transition to it and possibly cluster_28, 51 transition to BWM_2
                   '24', '13', '1', '11', '33', '48', '18', '46', '15', '26' # BWM-2
  )

  cluster.sels = c('36', '8', '39', '2', '19', '27', # terminal BWM-1 without transition
                   '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM-2
                   '44', '31', '52', '28',
                   '25', # possible transition clusters
                   '24' # many cells are not annotated
  )

  table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == '50'], useNA = 'ifany')

  table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

  xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSapaapp')])
  xx[which(xx > 0)]

  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
  #                                   & is.na(seurat.obj$manual.annot.ids))]
  cluster.sels = c('25', '31', '28', '52', '36', '8', '39', '2', '19', '27')

  cluster.sels = c('25', '36', '8', '39', '2', '19', '27', '13', '1', '11', '33', '48', '18', '46', '15', '26')
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
  cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
                                      seurat.obj$manual.annot.ids == 'MSpappa' |
                                      seurat.obj$manual.annot.ids == 'MSxappp' |
                                      seurat.obj$manual.annot.ids == 'MSxapppa' |
                                      seurat.obj$manual.annot.ids == 'MSxppppx/MSxpppax'|

                                      seurat.obj$manual.annot.ids == 'MSxpppp'|
                                      seurat.obj$manual.annot.ids == 'MSxpppa'|
                                      seurat.obj$manual.annot.ids == 'MSxppap'|
                                      seurat.obj$manual.annot.ids == 'MSxppaa'|
                                      seurat.obj$manual.annot.ids == 'MSxpapp'|
                                      seurat.obj$manual.annot.ids == 'MSxpapa'|
                                      seurat.obj$manual.annot.ids == 'MSxpaap'|
                                      #seurat.obj$manual.annot.ids == 'MSapaapp'|
                                      seurat.obj$manual.annot.ids == 'MSxpaaa'
                                    ]

  sub.obj = subset(seurat.obj, cells = cells.sels)

  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

  DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500', label = FALSE)
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
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 1000)
  #length(intersect(VariableFeatures(sub.obj), timers))
  #VariableFeatures(sub.obj) = setdiff(VariableFeatures(sub.obj), timers)
  cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')

  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20; min.dist = 0.1; spread = 1;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist)

  #DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label = TRUE, label.size = 6)
  DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE)

  p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1

  FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'pha-4'))

  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }

  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1)

  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2)
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 5,
               na.value = "gray", combine = TRUE)
  p1 + p2

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
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  sub.obj = find.reference.mapped.ids.variableGenes.scamp(sub.obj, nfeatures = 2000, filter.ref.MS.terminal = TRUE)
  sub.obj$predicted.ids = sub.obj$predicted.ids.scmap.newfeatures
  #sub.obj$predicted.ids[sub.obj$predicted.ids.scmap.newfeatures.cor<0.5] = NA
  counts = table(sub.obj$predicted.ids.scmap.newfeatures, as.character(sub.obj$seurat_clusters_split))
  counts.seurat = table(sub.obj$predicted.ids.seurat.terminal, as.character(sub.obj$seurat_clusters_split))
  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)

  p1 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()
  p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat.terminal', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()

  p1 + p2

  legend.txt = NULL
  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = legend.txt
  )

  barplot(counts.seurat, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = legend.txt
  )

  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
  idents.sel = setdiff(idents.sel, c('0', '1', '4', '5'))
  idents.sel = c('18', '14', '5', '17', '16', '4')
  idents.sel = c('1', '10', '3', '8', '18', '5', '16')
  idents.sel = c('0', '1', '2', '3', '5', '6', '7', '8', '0', '14', '15', '18' )
  ## chcek the reference-mapped ids for the rest of clusters
  counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
  counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
  counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
  counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
  counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
  counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

  #features.sels = c('unc-120', 'hnd-1', 'hlh-1', 'abts-1', 'ref-2', 'tbx-7', 'unc-39', 'cup-4', 'ins-2', 'F40H3.3', 'hot-1')
  features.sels = c('unc-120', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'kvs-5', 'hlh-1', 'tnt-3', 'spp-15', 'gana-1', 'lin-39', 'rpm-1')
  #features.sels = c('unc-120', 'clec-264', 'T08H10.1', 'kvs-5', 'pha-4')
  features.terminal = c('gana-1', 'spp-15', 'tnt-3', 'F37H8.5', 'lin-39', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'far-5',
                        'clec-264', 'T08H10.1', 'kvs-5', 'B0379.1', 'zig-6', 'frpr-8', 'C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-1')

  features.sels = c('unc-120', 'pha-4',
                    'hot-1', 'wago-1', 'pde-6','clec-264', 'T08H10.1',  'B0379.1',
                    'tnt-3', 'gana-1','F37H8.5','kvs-5', 'spp-15'
  )

  VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)

  # check info in JM data for specific lineage
  ids.sel = c('MSxppppx')
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)

  # to find new marker genes
  top.markers[top.markers$cluster == '1',]

  DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

  FeaturePlot(sub.obj, reduction = 'umap', features = c('hnd-1', 'tbx-7', 'abts-1', 'Y66D12A.13', 'F41D9.2', 'unc-120', 'zig-8')) # MSxpaaa
  FeaturePlot(sub.obj, reduction = 'umap', features = c('tbx-8', 'F40H3.3', 'ins-2', 'Y42H9B.3', 'unc-120', 'ZK180.5')) # MSxpapa
  FeaturePlot(sub.obj, reduction = 'umap', features = c('zig-6', 'frpr-8','C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-4', 'maph-1.2', 'irx-1')) # MSxpapa


  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '13')]
  sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == '13')] = 'MSxpaaa'
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpaaa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '17')]
  sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == '17')] = 'MSxpapa'
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpapa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '17')]
  sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == '17')] = 'MSxpapa'
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxpapa'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '0')]
  sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == '0')] = 'MSxppapp'
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppapp'

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '6')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '7')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '14')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '2')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '15')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '5')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '18')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '8')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '3')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '10')]
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'BWD_terminal'

  #manual.assign.cluster.with.annotation(cluster.index = '17', id2assign = 'MSxppapp', sub.obj = sub.obj, seurat.obj = seurat.obj)

  #cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == '1')]
  #seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = 'MSxppapp'

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)

  FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

  saveRDS(seurat.obj, file = RDS2save)


  ########################################################
  ########################################################
  # Section : iteration 11
  # first try of manual annotate terminal cells by considering all BDW terminal and their mothers
  # we learnt that:
  # 1) select specific lineages in reference and to use variable gene for that time points help a lot
  # for reference-based annotation, in particular seurat works well for the terminal cells;
  # using the subset reference and relevant variable genes, seurat works quite well to transfer labels for terminal cells
  # 2) one characteristic of terminal cells are more likely smooth trajectory, which poses probolem for the clustering
  # 3) long list of newly identified cluster-specific genes to annotation terminal cells works a little better than a few existing marker genes
  # but less well than seurat. The reason is the trajectory shape of cell ids.
  #
  ########################################################
  ########################################################
  nb.iteration = 11
  Refine.annotated.ids = FALSE;

  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_',
                    nb.iteration -1, '.rds')
  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, "_redo.pdf")
  RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_',
                     nb.iteration, '.rds')

  seurat.obj = readRDS(file = RDSsaved)
  seurat.obj$predicted.ids.scmap = seurat.obj$scmap.pred.id.500
  seurat.obj$predicted.ids.seurat = seurat.obj$seurat.pred.id
  seurat.obj$previous.iteration.clusters = NA
  seurat.obj$BWM.cells = NA

  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}


  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

  #cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4',  '22', '24',  # middle time point for convergence branche
                   '3', '5', '16', '30', # middle time points for MSxp
                   '23', '43', '17', '41', '45', '34', '50', # well-separated clusters with clear mapped ids
                   '44', '28', '52', '31', '53',  # # well-separated clusters without clear mapped ids
                   '25', '36', '8', '39', '2', '19', '27', # BWM-1 cluster_25 transition to it and possibly cluster_28, 51 transition to BWM_2
                   '24', '13', '1', '11', '33', '48', '18', '46', '15', '26' # BWM-2
  )

  cluster.sels = c('36', '8', '39', '2', '19', '27', # terminal BWM-1 without transition
                   '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM-2
                   '44', '31', '52', '28',
                   '25', # possible transition clusters
                   '24' # many cells are not annotated
  )

  table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == '50'], useNA = 'ifany')

  table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

  xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxppapp')])
  xx[which(xx > 0)]

  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
  #                                   & is.na(seurat.obj$manual.annot.ids))]
  cluster.sels = c('25', '31', '28', '52', '36', '8', '39', '2', '19', '27')
  cluster.sels = c('25', '36', '8', '39', '2', '19', '27', '13', '1', '11', '33', '48', '18', '46', '15', '26')
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
  # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
  # seurat.obj$manual.annot.ids == 'MSpappa' |
  # seurat.obj$manual.annot.ids == 'MSxappp' |
  # seurat.obj$manual.annot.ids == 'MSxapppa' |
  # seurat.obj$manual.annot.ids == 'MSxppppx/MSxpppax'|
  #
  # seurat.obj$manual.annot.ids == 'MSxpppp'|
  # seurat.obj$manual.annot.ids == 'MSxpppa'|
  # seurat.obj$manual.annot.ids == 'MSxppap'|
  # seurat.obj$manual.annot.ids == 'MSxppaa'|
  # seurat.obj$manual.annot.ids == 'MSxpapp'|
  # seurat.obj$manual.annot.ids == 'MSxpapa'|
  # seurat.obj$manual.annot.ids == 'MSxpaap'|
  # #seurat.obj$manual.annot.ids == 'MSapaapp'|
  # seurat.obj$manual.annot.ids == 'MSxpaaa'
  #]
  ids.annotated.sel = c('BWD_terminal',
                    'MSxppapp', # from cluster 2 19 and 39
                    'MSxapppa', # this annotation mainly from cluster 24 and 22
                    'MSxppppx/MSxpppax', # mainly from cluster 24
                    'unknown.MSxapppax.MSxpppaa.MSxpppaa' # mainly from cluster 31 and cluster 17
                    #'unknown.MSxpppaa' # mainly from cluster 50
                    )
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.annotated.sel))])

  sub.obj = subset(seurat.obj, cells = cells.sels)

  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

  DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500', label = FALSE)
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
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 5000)
  #length(intersect(VariableFeatures(sub.obj), timers))
  #VariableFeatures(sub.obj) = setdiff(VariableFeatures(sub.obj), timers)
  cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')

  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, npcs = 50)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 50 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20; min.dist = 0.1; spread = 1;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist, verbose = FALSE)

  #DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label = TRUE, label.size = 6)
  DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE)

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE)

  p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1

  FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'pha-4'))

  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }

  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = n.neighbors, dims = 1:nb.pcs, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 2.5)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  ##########################################
  # redo the seurat reference-projection
  ##########################################
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  terminals = c('MSxppppx', 'MSxpppax', 'MSxppppp','MSxppppa', 'MSxpppaa', 'MSxpppap',
                'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap',
                'MSxpaaap', 'MSxapppp', 'MSxapppa',
                'MSxappppx', 'MSxapppax', 'MSpappax',
                'MSxppp', 'MSxppa', 'MSxppap', 'MSxpapp', 'MSxpapa', 'MSxpaaa', 'MSxappp', 'MSpappa')

  #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
  #sub.obj$predicted.ids[sub.obj$predicted.ids.scmap.newfeatures.cor<0.5] = NA
  #counts = table(sub.obj$predicted.ids.scmap.newfeatures, as.character(sub.obj$seurat_clusters_split))
  sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 5000, terminals = terminals)

  sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
  sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
  sub.obj$predicted.ids.filtered = sub.obj$predicted.ids.seurat.terminal
  sub.obj$predicted.ids.filtered[sub.obj$predicted.ids.prob < 0.5] = NA


  p1  = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE, pt.size = 2)
  p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE) + NoLegend()
  p22 = DimPlot(sub.obj, group.by = 'predicted.ids.filtered', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE) + NoLegend()

  p1 + p2

  p1 + p22


  p3 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 6)
  p2 + p3

  p4 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3,
               group.by = 'seurat_clusters_split')

  (p1 + p3) / p4


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
  counts.seurat = table(sub.obj$predicted.ids, as.character(sub.obj$seurat_clusters_split))
  counts.seurat.filter = table(sub.obj$predicted.ids.fitered, as.character(sub.obj$seurat_clusters_split))
  #counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)


  p1 = DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
  p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()
  # DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()

  p1 + p2


  legend.txt = NULL
  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = legend.txt
  )

  barplot(counts.seurat, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = legend.txt
  )


  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
  idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

  idents.sel = c('4', '10', '8', '6', '5', '15', '9')
  idents.sel = c('1', '12', '2', '17', '18', '14', '16')
  idents.sel = c('13', '3', '0', '11', '7')
  ## chcek the reference-mapped ids for the rest of clusters
  counts = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))];
  counts.filtered = counts.seurat.filter[, !is.na(match(colnames(counts.seurat), idents.sel))]
  counts = counts[apply(as.matrix(counts), 1, sum)>0, ];
  counts.filtered = counts.filtered[apply(as.matrix(counts.filtered), 1, sum)>0, ]

  #counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
  #counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
  #counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
  #counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

  features.sels = c('unc-120', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'kvs-5', 'hlh-1', 'tnt-3', 'spp-15', 'gana-1', 'lin-39', 'rpm-1')
  features.terminal = c('gana-1', 'spp-15', 'tnt-3', 'F37H8.5', 'lin-39', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'far-5',
                        'clec-264', 'T08H10.1', 'kvs-5', 'B0379.1', 'zig-6', 'frpr-8', 'C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-1')

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

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  cluster.assingment = list(c('4', 'MSxappppx'), # orignal_clsuter 33
                            c('5', 'MSxpppaa'),
                            c('6', 'MSxpppaa'),
                            c('8', 'MSxppppp'),
                            c('10', 'MSxappppx/MSxapppp/MSxapppax/MSxapppa'),
                            c('15', 'MSxppppx'),
                            c('9', 'MSxppppp/MSxappp'),
                            c('0', 'MSxpaaap/MSxppapp'),
                            c('11', 'MSxppapp'),
                            c('13', 'MSxpapap'),
                            c('3', 'MSxpaaap/MSxpapap'),
                            c('7', 'MSxpaaap/MSxppapp/MSxpapap/MSxppap'),
                            c('1', 'MSxppppa/MSxppppp/MSxpppaa/MSxpappa'),
                            c('2', 'MSxppap/MSxpaaa/MSxpaaap'),
                            c('12', 'MSxpapp'),
                            c('14', 'MSxappp'),
                            c('16', 'likely_nonBWM_origCluster_31'),
                            c('17', 'BWM_terminal_origCluster_48'),
                            c('18', 'likely.nonBWM_origCluster_17')
  )
  length(cluster.assingment)

  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];

    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    seurat.obj$previous.iteration.clusters[match(cells, colnames(seurat.obj))] = cluster.index
  }

  #manual.assign.cluster.with.annotation(cluster.index = '17', id2assign = 'MSxppapp', sub.obj = sub.obj, seurat.obj = seurat.obj)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
          pt.size = 2)


  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + NoLegend()
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)

  FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

  saveRDS(seurat.obj, file = RDS2save)



  ########################################################
  ########################################################
  # Section : iteration 12
  # second try of manually annotate BWM-terminal cells by
  # selecting cells from original clusters c('26', '15', '18', '1', '11', '46','48', '25',  '33', '13', '24', '22') # MSxappp lineage
  # mainly for one branche of BWM terminal cells mainly from MSxppppx (MSxppppa/p), MSxpppax, MSxapppp(x) and MSxapppa(x)
  ########################################################
  ########################################################
  nb.iteration = 12
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

  #cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
  cluster.sels = c('4',  '22', '24',  # middle time point for convergence branche
                   '3', '5', '16', '30', # middle time points for MSxp
                   '23', '43', '17', '41', '45', '34', '50', # well-separated clusters with clear mapped ids
                   '44', '28', '52', '31', '53',  # # well-separated clusters without clear mapped ids
                   '25', '36', '8', '39', '2', '19', '27', # BWM-1 cluster_25 transition to it and possibly cluster_28, 51 transition to BWM_2
                   '24', '13', '1', '11', '33', '48', '18', '46', '15', '26' # BWM-2
  )

  cluster.sels = c('36', '8', '39', '2', '19', '27', # terminal BWM-1 without transition
                   '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM-2
                   '44', '31', '52', '28',
                   '25', # possible transition clusters
                   '24' # many cells are not annotated
  )

  table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == '24'], useNA = 'ifany')

  table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

  xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSapaapp')])
  xx[which(xx > 0)]

  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
  #                                   & is.na(seurat.obj$manual.annot.ids))]
  #cluster.sels = c('25', '31', '28', '52', '36', '8', '39', '2', '19', '27')
  #cluster.sels = c('25', '36', '8', '39', '2', '19', '27', '13', '1', '11', '33', '48', '18', '46', '15', '26')
  # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
  # seurat.obj$manual.annot.ids == 'MSpappa' |
  # seurat.obj$manual.annot.ids == 'MSxappp' |
  # seurat.obj$manual.annot.ids == 'MSxapppa' |
  # seurat.obj$manual.annot.ids == 'MSxppppx/MSxpppax'|
  #
  # seurat.obj$manual.annot.ids == 'MSxpppp'|
  # seurat.obj$manual.annot.ids == 'MSxpppa'|
  # seurat.obj$manual.annot.ids == 'MSxppap'|
  # seurat.obj$manual.annot.ids == 'MSxppaa'|
  # seurat.obj$manual.annot.ids == 'MSxpapp'|
  # seurat.obj$manual.annot.ids == 'MSxpapa'|
  # seurat.obj$manual.annot.ids == 'MSxpaap'|
  # #seurat.obj$manual.annot.ids == 'MSapaapp'|
  # seurat.obj$manual.annot.ids == 'MSxpaaa'
  #]
  # cells.sels = colnames(seurat.obj)[seurat.obj$manual.annot.ids == 'MSxpappa/MSxapppax/MSxppppx' | # cluster 0 and cluster 5
  #                                     seurat.obj$manual.annot.ids == 'MSxppppx/MSxappppx/MSxapppa/MSxapppp' | # cluster 3 'MSxppppx/MSxappppx/MSxapppa/MSxapppp'
  #                                     seurat.obj$manual.annot.ids == 'MSxpappa/MSxapppax' | # cluster 1
  #                                     seurat.obj$manual.annot.ids == 'MSxappppx/MSxapppax' ] # cluster 6
  #

  cluster.sels = c('26', '15', '18', '1',
                  '11', '46','48', '25',
                   '33', '13', '24', '22') # MSxappp lineage
  cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
  sub.obj = subset(seurat.obj, cells = cells.sels)

  barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 0.7, las=2)

  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

  DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500', label = FALSE)
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
  # rerun the seurat for label transferring
  ##########################################
  RErun.seurat.transferring.labels = TRUE
  if(RErun.seurat.transferring.labels){
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    terminals = c('MSxppppx', 'MSxpppax', 'MSxppppp','MSxppppa', 'MSxpppaa', 'MSxpppap',
                  'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap',
                  'MSxpaaap', 'MSxapppp', 'MSxapppa',
                  'MSxappppx', 'MSxapppax', 'MSpappax',
                  'MSxppp', 'MSxppa', 'MSxppap', 'MSxpapp', 'MSxpapa', 'MSxpaaa', 'MSxappp', 'MSpappa')

    #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
    sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 5000,  npcs = 50, terminals = terminals)

    sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
    sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA

  }

  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = 3000)
  #length(intersect(VariableFeatures(sub.obj), timers))
  #VariableFeatures(sub.obj) = setdiff(VariableFeatures(sub.obj), timers)
  cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')

  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 30 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 10;
  min.dist = 0.01; spread = 2;
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist)

  #DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label = TRUE, label.size = 6)
  #DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
  DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

  p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
  p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
  p0 + p1

  FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'pha-4'))

  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }

  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:npcs, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)

  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p0 = DimPlot(sub.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)
  p1  = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE,  pt.size = 2) +
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

  p1 = DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
  p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()
  # DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()

  p1 + p2


  legend.txt = NULL
  barplot(counts, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = legend.txt
  )

  barplot(counts.seurat, main="cluster compositions for predicted labels ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = legend.txt
  )


  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
  idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

  idents.sel = c('4', '10', '8', '6', '5', '15', '9')
  idents.sel = c('1', '12', '2', '17', '18', '14', '16')
  idents.sel = c('13', '3', '0', '11', '7')
  ## chcek the reference-mapped ids for the rest of clusters
  counts = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))];
  counts.filtered = counts.seurat.filter[, !is.na(match(colnames(counts.seurat), idents.sel))]
  counts = counts[apply(as.matrix(counts), 1, sum)>0, ];
  counts.filtered = counts.filtered[apply(as.matrix(counts.filtered), 1, sum)>0, ]

  #counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
  #counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
  #counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
  #counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

  features.sels = c('unc-120', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'kvs-5', 'hlh-1', 'tnt-3', 'spp-15', 'gana-1', 'lin-39', 'rpm-1')
  features.terminal = c('gana-1', 'spp-15', 'tnt-3', 'F37H8.5', 'lin-39', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'far-5',
                        'clec-264', 'T08H10.1', 'kvs-5', 'B0379.1', 'zig-6', 'frpr-8', 'C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-1')

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

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  #n = 5; cbind(counts.annot[, n], counts.seurat.filter[,n])
  cluster.assingment = list(c('0', 'MSxpppaa'), # orignal_clsuter 33
                            c('1', 'MSxppppa'),
                            c('2', 'MSxppap/MSxpaaa/MSxpapa/MSxpaaap'),
                            c('3', 'MSxappppx'),
                            c('4', 'MSxppppp/MSxpppap'),
                            c('5', 'MSxappp'),
                            c('6', 'MSxapppp/MSxapppa'),
                            c('7', 'MSxpapp/MSxpppp'),
                            c('8', 'MSxppppx/MSxpppaa/MSxpappa'),
                            c('9', 'MSxppppa/MSxppppp/MSxpppaa/MSxpappa/MSxappp')
                            #c('10', 'MSxppppp/MSxpppaa')

  )
  length(cluster.assingment)

  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];

    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    seurat.obj$previous.iteration.clusters[match(cells, colnames(seurat.obj))] = cluster.index
  }

  #manual.assign.cluster.with.annotation(cluster.index = '17', id2assign = 'MSxppapp', sub.obj = sub.obj, seurat.obj = seurat.obj)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
          pt.size = 2)


  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") + NoLegend()
  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE)

  FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

  saveRDS(seurat.obj, file = RDS2save)


}


########################################################
########################################################
# Section : iteration 13
# further annotate the BWM terminal cells
# selecting cells from original clusters 27, 19, 2, 39, 30, 16, 36, 8, 25 48, 46
# mainly for one branche of BWM terminal cells mainly from MSxppapp, MSxpapa, MSxpappp, MSxpaaap
# However, it turns out that this does not work, because seurat yields confusing predicted annotation, MSxapppx MSxppppx lineages
# were found here
# at the end, we were not able to improve too much the annotation for this Branch of BWM-terminal cells
########################################################
########################################################
nb.iteration = 13
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

#cluster.sels = c('4',  '22', '24', '3', '5', '16', '30') # cluster 4 and 22 have some cells in Msxp lineage
cluster.sels = c('4',  '22', '24',  # middle time point for convergence branche
                 '3', '5', '16', '30', # middle time points for MSxp
                 '23', '43', '17', '41', '45', '34', '50', # well-separated clusters with clear mapped ids
                 '44', '28', '52', '31', '53',  # # well-separated clusters without clear mapped ids
                 '25', '36', '8', '39', '2', '19', '27', # BWM-1 cluster_25 transition to it and possibly cluster_28, 51 transition to BWM_2
                 '24', '13', '1', '11', '33', '48', '18', '46', '15', '26' # BWM-2
)

cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                 '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                 '44', '31', '52', '28', '50', # we
                 '25', # possible transition clusters
                 '24' # many cells are not annotated
)

cluster.index = '52'
table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSapaapp')])
xx[which(xx > 0)]

#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
#                                   & is.na(seurat.obj$manual.annot.ids))]
cluster.sels = c('25', '31', '28', '52', '36', '8', '39', '2', '19', '27')
cluster.sels = c('25', '36', '8', '39', '2', '19', '27', '13', '1', '11', '33', '48', '18', '46', '15', '26')
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]

# BWM_terminal_2
# cells.sels = colnames(seurat.obj)[seurat.obj$manual.annot.ids == 'MSxpappa/MSxapppax/MSxppppx' | # cluster 0 and cluster 5 in iteration 11
#                                     seurat.obj$manual.annot.ids == 'MSxppppx/MSxappppx/MSxapppa/MSxapppp' | # cluster 3
#                                     seurat.obj$manual.annot.ids == 'MSxpappa/MSxapppax' | # cluster 1
#                                     seurat.obj$manual.annot.ids == 'MSxappppx/MSxapppax' ] # cluster 6
# BWM_terminal_1
# cells.sels = colnames(seurat.obj)[seurat.obj$manual.annot.ids == 'MSxpaaap/MSxppapp/MSxpappp'| # cluster 2 in iteration 11
#                                     seurat.obj$manual.annot.ids == 'MSxpaaap/MSxppapp/MSxpappp' | # cluster 9
#                                     seurat.obj$manual.annot.ids == 'MSxpaaap'| # cluster 4
#                                     seurat.obj$manual.annot.ids == 'MSxppapp' | # cluster 14
#                                     seurat.obj$manual.annot.ids == 'MSxpaaap/MSxppapp/MSxpapap'] # cluster 10 and 11

##########################################
# try to select all BWM_cells
##########################################
cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                 '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                 '25', # possible transition clusters
                 '24', # also transitions cells and many of them are not annotated
                 '2', '3', '16', '30' # all middle time points
)

cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
                                           seurat.obj$manual.annot.ids == 'MSxppppx/MSxpppax'|
                                           seurat.obj$manual.annot.ids == 'MSxapppa' |
                                           seurat.obj$manual.annot.ids == 'unknown.MSxapppax.MSxpppaa.MSxpppaa' |

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
                                           seurat.obj$manual.annot.ids == 'MSxapp'|

                                           seurat.obj$manual.annot.ids == 'MSxpp'|
                                           seurat.obj$manual.annot.ids == 'MSxpa'|
                                           seurat.obj$manual.annot.ids == 'MSxap'|

                                           seurat.obj$manual.annot.ids == 'MSxp'|
                                           seurat.obj$manual.annot.ids == 'MSxa'|
                                           seurat.obj$manual.annot.ids == 'MSx'
                                         ])


cluster.sels = c('27', '19', '2', '39', '30', '36', '8', '25')
# cluster.sels = c('27', '19', '2', '39', '36', '8')
cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
sub.obj = subset(seurat.obj, cells = cells.sels)
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
nfeatures = 3000;
nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
n.neighbors = 5;
min.dist = 0.05; spread = 1;

sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
#cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 50)
sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                   spread = spread, n.neighbors = n.neighbors,
                   min.dist = min.dist)

DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)

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
                'MSxppp', 'MSxppa', 'MSxppap', 'MSxpapp', 'MSxpapa', 'MSxpaaa', 'MSxappp', 'MSpappa')

  #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
  sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000,  npcs = 30,
                                                                            k.anchor = 10, k.filter = 50,
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

p1 + p2 + p3

DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
p0 + p1

FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'pha-4'))

##########################################
# redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
##########################################
FindClusters_subclusters = function(sub.obj, resolution = 0.4)
{
  sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
  return(sub.obj$seurat_clusters)
}

sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:nb.pcs, compute.SNN = TRUE)
sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.8)
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

p1 = DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()
p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()

(p1 + p2)/p3

legend.txt = NULL
barplot(counts, main="cluster compositions for predicted labels ",
        xlab=NULL, col=c(1:nrow(counts)), las = 2,
        legend = legend.txt
)

barplot(counts.seurat, main="cluster compositions for predicted labels ",
        xlab=NULL, col=c(1:nrow(counts)), las = 2,
        legend = legend.txt
)


Idents(sub.obj) = sub.obj$seurat_clusters_split
idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

#idents.sel = c('0', '1', '2', '3', '5', '6', '7', '8', '0', '14', '15', '18' )

## chcek the reference-mapped ids for the rest of clusters
counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

#features.sels = c('unc-120', 'hnd-1', 'hlh-1', 'abts-1', 'ref-2', 'tbx-7', 'unc-39', 'cup-4', 'ins-2', 'F40H3.3', 'hot-1')
features.sels = c('unc-120', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'kvs-5', 'hlh-1', 'tnt-3', 'spp-15', 'gana-1', 'lin-39', 'rpm-1')
#features.sels = c('unc-120', 'clec-264', 'T08H10.1', 'kvs-5', 'pha-4')
features.terminal = c('gana-1', 'spp-15', 'tnt-3', 'F37H8.5', 'lin-39', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'far-5',
                      'clec-264', 'T08H10.1', 'kvs-5', 'B0379.1', 'zig-6', 'frpr-8', 'C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-1')

features.sels = c('unc-120', 'pha-4',
                  'hot-1', 'wago-1', 'pde-6','clec-264', 'T08H10.1',
                  'gana-1','kvs-5', 'far-5','F37H8.5'
)

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
cluster.assingment = list(c('0', 'MSxpaaap/MSxppapp.300min'),
                          c('1', 'MSxpaaap.300min'),
                          #c('2', 'MSxppap/MSxpaaa/MSxpaaap, MSxppap/MSxpaaa/MSxpapa/MSxpaaap'),
                          #c('3', ''),
                          c('4', 'MSxpapp/MSxpapa'),
                          #c('5', ''),
                          #c('6', ''),
                          #c('7', ''),
                          #c('8', ''),
                          #c('9', ''),
                          #c('10', ''),
                          c('11', 'MSxpapa')
                          #c('12', '')
)

for(n in 1:length(cluster.assingment)){
  cluster.index = cluster.assingment[[n]][1];
  id2assign =  cluster.assingment[[n]][2];

  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
  sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
}

#manual.assign.cluster.with.annotation(cluster.index = '17', id2assign = 'MSxppapp', sub.obj = sub.obj, seurat.obj = seurat.obj)
DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
        pt.size = 2)


DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
        na.value = "gray") +
  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
  scale_colour_hue(drop = FALSE)

#FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

saveRDS(seurat.obj, file = RDS2save)


########################################################
########################################################
# Section : iteration 14
# start to integrate BWM terminal and mother cells
# cells were selected from original cluster 36, 39, 2, 19, 27, and also cluster 24, 13, 11, 1, 46, 18, 33, 15, 26, 48
#
########################################################
########################################################
nb.iteration = 14
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
# select cells of interest
##########################################
cluster.index = '24'
table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxppppx/MSxpppax')])
xx[which(xx > 0)]

cluster.sels = c('4',  '22', '24',  # middle time point for convergence branche
                 '3', '5', '16', '30', # middle time points for MSxp
                 '23', '43', '17', '41', '45', '34', '50', # well-separated clusters with clear mapped ids
                 '44', '28', '52', '31', '53',  # # well-separated clusters without clear mapped ids
                 '25', '36', '8', '39', '2', '19', '27', # BWM-1 cluster_25 transition to it and possibly cluster_28, 51 transition to BWM_2
                 '24', '13', '1', '11', '33', '48', '18', '46', '15', '26' # BWM-2
)

cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                 '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                 '44', '31', '52', '28', '50', # we
                 '25', # possible transition clusters
                 '24' # many cells are not annotated
)

### BWM terminal cells from cluster 36, 39, 2, 19, 27, 8 and also cluster 24, 13, 11, 1, 46, 18, 33, 15, 26, 48
cluster.sels = c('27', '19', '2', '39', '36', '8', '24', '13', '11', '1', '46', '18', '33', '15', '26', '48')
cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
#                                   & is.na(seurat.obj$manual.annot.ids))]
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]

### try to select all BWM_cells;
### cluster '44', '31', '52', '28', '50' were not included here
cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                 '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                 '25', # possible transition clusters
                 '24', # also transitions cells and many of them are not annotated
                 '3', '5', '16', '30', '22', '4' # all middle time points
)
cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |

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
                                           seurat.obj$manual.annot.ids == 'MSxapp'|

                                           seurat.obj$manual.annot.ids == 'MSxpp'|
                                           seurat.obj$manual.annot.ids == 'MSxpa'|
                                           seurat.obj$manual.annot.ids == 'MSxap'|

                                           seurat.obj$manual.annot.ids == 'MSxp'|
                                           seurat.obj$manual.annot.ids == 'MSxa'|
                                           seurat.obj$manual.annot.ids == 'MSx'
                                         ])

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
nfeatures = 5000;
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
#cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE)
ElbowPlot(sub.obj, ndims = 50)
nb.pcs = 50 # nb of pcs depends on the considered clusters or ids
n.neighbors = 50;
min.dist = 0.1; spread = 1;
sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                   spread = spread, n.neighbors = n.neighbors,
                   min.dist = min.dist)

DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)

DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4) + NoLegend()
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
                'MSxppp', 'MSxppa', 'MSxppap', 'MSxpapp', 'MSxpapa', 'MSxpaaa', 'MSxappp', 'MSpappa')

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

p1 + p2 + p3

DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
p0 + p1

FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'pha-4'))
##########################################
# redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
##########################################
FindClusters_subclusters = function(sub.obj, resolution = 0.4)
{
  sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
  return(sub.obj$seurat_clusters)
}

sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:nb.pcs, compute.SNN = TRUE)
sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.5)
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

p1 = DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()
p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()

(p1 + p2)/p3

legend.txt = NULL
barplot(counts, main="cluster compositions for predicted labels ",
        xlab=NULL, col=c(1:nrow(counts)), las = 2,
        legend = legend.txt
)

barplot(counts.seurat, main="cluster compositions for predicted labels ",
        xlab=NULL, col=c(1:nrow(counts)), las = 2,
        legend = legend.txt
)


Idents(sub.obj) = sub.obj$seurat_clusters_split
idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

idents.sel = c('0', '1', '2', '3', '5', '6', '7', '8' )

## chcek the reference-mapped ids for the rest of clusters
counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

#features.sels = c('unc-120', 'hnd-1', 'hlh-1', 'abts-1', 'ref-2', 'tbx-7', 'unc-39', 'cup-4', 'ins-2', 'F40H3.3', 'hot-1')
features.sels = c('unc-120', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'kvs-5', 'hlh-1', 'tnt-3', 'spp-15', 'gana-1', 'lin-39', 'rpm-1')
#features.sels = c('unc-120', 'clec-264', 'T08H10.1', 'kvs-5', 'pha-4')
features.terminal = c('gana-1', 'spp-15', 'tnt-3', 'F37H8.5', 'lin-39', 'hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2', 'far-5',
                      'clec-264', 'T08H10.1', 'kvs-5', 'B0379.1', 'zig-6', 'frpr-8', 'C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-1')

features.sels = c('unc-120', 'pha-4',
                  'hot-1', 'wago-1', 'pde-6','clec-264', 'T08H10.1',
                  'gana-1','kvs-5', 'far-5','F37H8.5'
)

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
# cluster.assingment = list(c('0', 'MSxppppa, MSxppppp/MSxpppap, MSxppppx/MSxpppaa/MSxpappa, MSxppppa/MSxppppp/MSxpppaa/MSxpappa/MSxappp '),
#                           #c('1', ''),
#                           #c('2', 'MSxppap/MSxpaaa/MSxpaaap, MSxppap/MSxpaaa/MSxpapa/MSxpaaap'),
#                           #c('3', ''),
#                           c('4', 'MSxpapp/MSxpapa'),
#                           #c('5', ''),
#                           #c('6', ''),
#                           #c('7', ''),
#                           #c('8', ''),
#                           #c('9', ''),
#                           #c('10', ''),
#                           c('11', 'MSxpapa'),
#                           #c('12', ''),
#                           c('13', ''),
#                           c('14', ''),
#                           c('15', ''),
#                           c('16', ''),
#                           c('17', ''),
#                           c('18', '')
#
# )
#
# for(n in 1:length(cluster.assingment)){
#   cluster.index = cluster.assingment[[n]][1];
#   id2assign =  cluster.assingment[[n]][2];
#
#   cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
#   sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
#   seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
# }

#manual.assign.cluster.with.annotation(cluster.index = '17', id2assign = 'MSxppapp', sub.obj = sub.obj, seurat.obj = seurat.obj)
DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
        pt.size = 2)


DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
        na.value = "gray") +
  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
  scale_colour_hue(drop = FALSE)

#FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

saveRDS(seurat.obj, file = RDS2save)


########################################################
########################################################
# Section : iteration 15
# start to integrate BWM terminal and mother cells
# cells were selected from original cluster 36, 39, 2, 19, 27, and also cluster 24, 13, 11, 1, 46, 18, 33, 15, 26, 48
#
########################################################
########################################################
nb.iteration = 15
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
# select cells of interest
##########################################
cluster.index = '24'
table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxppppx/MSxpppax')])
xx[which(xx > 0)]

cluster.sels = c('4',  '22', '24',  # middle time point for convergence branche
                 '3', '5', '16', '30', # middle time points for MSxp
                 '23', '43', '17', '41', '45', '34', '50', # well-separated clusters with clear mapped ids
                 '44', '28', '52', '31', '53',  # # well-separated clusters without clear mapped ids
                 '25', '36', '8', '39', '2', '19', '27', # BWM-1 cluster_25 transition to it and possibly cluster_28, 51 transition to BWM_2
                 '24', '13', '1', '11', '33', '48', '18', '46', '15', '26' # BWM-2
)

cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                 '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                 '44', '31', '52', '28', '50', # we
                 '25', # possible transition clusters
                 '24' # many cells are not annotated
)

### BWM terminal cells from cluster 36, 39, 2, 19, 27, 8 and also cluster 24, 13, 11, 1, 46, 18, 33, 15, 26, 48
cluster.sels = c('27', '19', '2', '39', '36', '8', '24', '13', '11', '1', '46', '18', '33', '15', '26', '48')
cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
#                                   & is.na(seurat.obj$manual.annot.ids))]
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]

### try to select all BWM_cells;
### cluster '44', '31', '52', '28', '50' were not included here
cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                 '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                 '25', # possible transition clusters
                 '24', # also transitions cells and many of them are not annotated
                 '3', '5', '16', '30', '22', '4' # all middle time points
)
cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |

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
                                           seurat.obj$manual.annot.ids == 'MSxapp'|

                                           seurat.obj$manual.annot.ids == 'MSxpp'|
                                           seurat.obj$manual.annot.ids == 'MSxpa'|
                                           seurat.obj$manual.annot.ids == 'MSxap'|

                                           seurat.obj$manual.annot.ids == 'MSxp'|
                                           seurat.obj$manual.annot.ids == 'MSxa'|
                                           seurat.obj$manual.annot.ids == 'MSx'
                                         ])

seurat.obj$BWM.cells[!is.na(match(colnames(seurat.obj), cells.sels))] = 'BWM'
seurat.bwm = subset(seurat.obj, cells = cells.sels)
seurat.nonbwm = subset(seurat.obj, cells = setdiff(colnames(seurat.obj), cells.sels))

xx = table(seurat.nonbwm$seurat_clusters[which(seurat.nonbwm$manual.annot.ids == 'unknown.MSxpppaa')])
xx[which(xx > 0)]
sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[which(seurat.obj$BWM.cells == 'BWM')])

sub.obj = subset(seurat.obj, cells = cells.sels)

sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)
barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)


saveRDS(seurat.obj, file = RDS2save)



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
#jj = which(seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31')
#seurat.obj$BWM.cells[] = NA
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
DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
        pt.size = 2) + NoLegend()

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

sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:nb.pcs, compute.SNN = TRUE)
sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
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

sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.keep
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


########################################################
########################################################
# Section : iteration 17
# the second round start
## first of all, we need to clean up the seurat.obj
########################################################
########################################################
nb.iteration = 17
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

## 37 ids including MSx and redundant MSxppppx and MSxpppax
bwms.all = c('MSxppppx', 'MSxpppax', # redundant
             'MSxppppp', 'MSxppppa', 'MSxpppaa', 'MSxpppap',
             'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap',
             'MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax', 'MSpappax', # all terminal cells
             'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxppaa', 'MSxpapp', 'MSxpapa', 'MSxpaap', 'MSxpaaa', 'MSxappp', 'MSpappa',
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
cluster.index = '30'
table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids.1 == 'MSxa')])
xx[which(xx > 0)]

#FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
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

#seurat.obj$BWM.cells[seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31'] = NA

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

nfeatures = 5000;
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
#cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(sub.obj, ndims = 50)

nb.pcs = 30 # nb of pcs depends on the considered clusters or ids
n.neighbors = 30;
min.dist = 0.1; spread = 1
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

  #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
  sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30,
                                                                            k.anchor = 5, k.filter = 200,
                                                                            terminals = bwms.all)
  sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
  sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
  sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
  sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA


}

p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
  NoLegend() + ggtitle('predicted ids')
#p2 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
#  NoLegend()
p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
  NoLegend() + ggtitle('manual ids')

p2 + p3

DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
p0 + p1

FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'pha-4', 'hnd-1'))
FeaturePlot(seurat.obj, reduction = 'umap', features = c('unc-120', 'pha-4', 'hnd-1'))

##########################################
# update the manual annotation if good marker genes or mapped ids were found
##########################################
#save first the seurat projection result, one of best projection prediciotn
seurat.obj$predicted.ids.seurat.keep = NA
mm = match(colnames(sub.obj), colnames(seurat.obj))
seurat.obj$predicted.ids.seurat.keep[mm] = sub.obj$predicted.ids
#seurat.obj$manual.annot.ids.1 = seurat.obj$manual.annot.ids
seurat.obj$manual.annot.ids = seurat.obj$manual.annot.ids.1

# integrate the good prediction from seurat
jj = which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$predicted.ids.seurat.keep))
seurat.obj$manual.annot.ids[jj] = seurat.obj$predicted.ids.seurat.keep[jj]
kk = which(!is.na(match(seurat.obj$predicted.ids.seurat.keep, terminals)))
seurat.obj$manual.annot.ids[kk] = seurat.obj$predicted.ids.seurat.keep[kk]

seurat.obj$manual.annot.ids.2 = seurat.obj$manual.annot.ids # save the current annotation
seurat.obj$manual.annot.ids[is.na(seurat.obj$BWM.cells)] = NA # keep only BWM cell annotation

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
  scale_colour_hue(drop = FALSE) + NoLegend()

#FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

saveRDS(seurat.obj, file = RDS2save)



########################################################
########################################################
# Section : iteration 18
# the second round
## we will always check the mother, current generation and daughter cells together
## after cleaning up and we focus on : MSx, MSxa, MSxp, MSxpp, MSxpa, MSxap
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

ee = process.import.Murray.scRNA()

murray.ids = unique(ee$lineage)
markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
# bwms = unique(c('MSx', 'MSxa', 'MSxap',
#          'MSapaap', 'MSapaapp', 'MSappaaa',
#          'MSpappa', 'MSpappax', 'MSppaap', 'MSppaapp', 'MSpppaaa',
#          murray.ids[grep('MSxapp|MSxp', murray.ids)]))

#markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
#markers = markers[!is.na(match(markers$Lineage, bwms)), ]
# seurat.obj = readRDS(file = paste0(RdataDir,
#                                    'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat.rds'))
#
# seurat.obj$manual.annot.ids = NA


nb.iteration = 18
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

## 37 ids including MSx and redundant MSxppppx and MSxpppax
bwms.all = c('MSxppppx', 'MSxpppax', # redundant
             'MSxppppp', 'MSxppppa', 'MSxpppaa', 'MSxpppap',
             'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap',
             'MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax', 'MSpappax', # all terminal cells
             'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxppaa', 'MSxpapp', 'MSxpapa', 'MSxpaap', 'MSxpaaa', 'MSxappp', 'MSpappa',
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
cluster.index = '14'
table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids.1 == 'MSxpp')])
xx[which(xx > 0)]

#FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
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
##########################################
cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                 '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                 '25', # possible transition clusters
                 '24', # also transitions cells and many of them are not annotated
                 #'44', '31', '52', '28', '50', # cluster '44', '31', '52', '28', '50' were not included here
                 '3', '5', '16', '30', '22', '4' # all middle time points
)

cluster.sels = c('29', '32', '35', '40', '42')
ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')

cells.sels = unique(colnames(seurat.obj)[seurat.obj$BWM.cells == 'BWM' &
                                           (!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
                                              !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
                                           )])

#cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
#                                           !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
#                                         ])

sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])

#seurat.obj$BWM.cells[seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31'] = NA
#sub.obj = subset(seurat.obj, cells = cells.sels)

table(sub.obj$seurat_clusters)[table(sub.obj$seurat_clusters)>0]

sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)

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
  test.umap.params.for.BWM.cells(sub.obj, pdfname = 'BWM_early_cells_umap_param_test_2.pdf',
                                 nfeatures.sampling = c(500, 1000, 3000), nb.pcs.sampling = c(5, 10, 20),
                                 n.neighbors.sampling = c(5, 10, 30),
                                 min.dist.sampling = c(0.01, 0.05, 0.1)
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
min.dist = 0.05; spread = 1
sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                   spread = spread, n.neighbors = n.neighbors,
                   min.dist = min.dist)
DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5)

features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)


p1 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend() +
  ggtitle('manual ids')
p2 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = , repel = TRUE,  pt.size = 2) +
  NoLegend() + ggtitle('predicted ids scmap')
p3 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,  pt.size = 2) +
  NoLegend() + ggtitle('predicted ids scmap')

p1 + p2 + p3
##########################################
# rerun the seurat for label transferring
##########################################
RErun.seurat.transferring.labels = FALSE
if(RErun.seurat.transferring.labels){
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

  #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
  sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30,
                                                                            k.anchor = 5, k.filter = 200,
                                                                            terminals = bwms.all)
  sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
  sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
  sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
  sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA

  p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
  p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    NoLegend() + ggtitle('predicted ids')
  #p2 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
  #  NoLegend()
  p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    NoLegend() + ggtitle('manual ids')

  p2 + p3

  DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()


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

sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:10, compute.SNN = TRUE)
sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

p0 = DimPlot(sub.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)
p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
  NoLegend()
p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
             label.size = 6,
             na.value = "gray", combine = TRUE)
p3 =   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
(p0 + p3) / (p1 + p2)

p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3,
             group.by = 'seurat_clusters_split')

(p1 + p2) / p3

features.sels = c('ceh-76', # MSx
                  'hnd-1',  # MSxp
                  'pha-4', 'fbxb-70', # MSxa
                  'unc-120', 'fbxa-81', # MSxpp
                  'unc-39', 'irx-1', 'egl-43', # MSxpa
                  'ceh-37', 'C45G7.4', 'nhr-67', # MSxap
                  'cft-1'
)
FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

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

# check info in JM data for specific lineage
ids.sel = c('MSxppppx')
source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
#find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
extrack.markers.from.JM(markers = markers.JM, id = 'MSx', ntop = 5)

top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

# to find new marker genes
top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
top.markers[top.markers$cluster == '13',]


features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)


VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split')

##########################################
# update the manual annotation if good marker genes or mapped ids were found
##########################################
#save first the seurat projection result, one of best projection prediciotn
seurat.obj$predicted.ids.seurat.keep = NA
mm = match(colnames(sub.obj), colnames(seurat.obj))
seurat.obj$predicted.ids.seurat.keep[mm] = sub.obj$predicted.ids
#seurat.obj$manual.annot.ids.1 = seurat.obj$manual.annot.ids
seurat.obj$manual.annot.ids = seurat.obj$manual.annot.ids.1

# integrate the good prediction from seurat
jj = which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$predicted.ids.seurat.keep))
seurat.obj$manual.annot.ids[jj] = seurat.obj$predicted.ids.seurat.keep[jj]
kk = which(!is.na(match(seurat.obj$predicted.ids.seurat.keep, terminals)))
seurat.obj$manual.annot.ids[kk] = seurat.obj$predicted.ids.seurat.keep[kk]

seurat.obj$manual.annot.ids.2 = seurat.obj$manual.annot.ids # save the current annotation
seurat.obj$manual.annot.ids[is.na(seurat.obj$BWM.cells)] = NA # keep only BWM cell annotation

cluster.assingment = list(c('0', 'MSxp'),
                          #c('1', ''), # don't change cluster 1
                          c('2', 'MSxpa'),
                          c('3', 'MSxa'),
                          c('4', 'MSxa'),
                          c('5', 'MSxap'),
                          c('6', 'MSxap'),
                          c('7', NA),
                          c('8', 'MSxpp'),
                          c('9', NA),
                          c('10', NA),
                          c('11', NA),
                          c('12', 'MSxa'),
                          c('13', 'MSxa')

)

for(n in 1:length(cluster.assingment)){

  cluster.index = cluster.assingment[[n]][1];
  id2assign =  cluster.assingment[[n]][2];
  cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
  sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
  if(is.na(id2assign)) seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA

}


# #manual.assign.cluster.with.annotation(cluster.index = '17', id2assign = 'MSxppapp', sub.obj = sub.obj, seurat.obj = seurat.obj)
DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
        pt.size = 2)

#
DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
        na.value = "gray") +
  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
  scale_colour_hue(drop = FALSE) + NoLegend()

#FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

saveRDS(seurat.obj, file = RDS2save)


########################################################
  ########################################################
  # Section : iteration 19 (Second round)
  ## we will always check the mother, current generation and daughter cells together
  ## after cleaning up and we focus on :
  # MSxa, MSxp, MSxpp, MSxpa, MSxap, MSxapp, MSxappp, MSpappa
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
  markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
  # bwms = unique(c('MSx', 'MSxa', 'MSxap',
  #          'MSapaap', 'MSapaapp', 'MSappaaa',
  #          'MSpappa', 'MSpappax', 'MSppaap', 'MSppaapp', 'MSpppaaa',
  #          murray.ids[grep('MSxapp|MSxp', murray.ids)]))

  #markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
  # seurat.obj = readRDS(file = paste0(RdataDir,
  #                                    'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat.rds'))
  #
  # seurat.obj$manual.annot.ids = NA


  nb.iteration = 19
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

  ## 36 ids including MSx and redundant MSxppppx and MSxpppax; MSx missing and not necessary
  bwms.all = c('MSxppppx', 'MSxpppax', # redundant
               'MSxppppp', 'MSxppppa', 'MSxpppaa', 'MSxpppap',
                'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap',
                'MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax', 'MSpappax', # all terminal cells
                'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxppaa', 'MSxpapp', 'MSxpapa', 'MSxpaap', 'MSxpaaa', 'MSxappp', 'MSpappa',
                'MSxppp', 'MSxppa', 'MSxpap', 'MSxpaa', 'MSxapp',
                'MSxpp', 'MSxpa', 'MSxap',
                'MSxp', 'MSxa')

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
  cluster.index = '22'
  table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

  table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

  table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

  xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxapp')])
  xx[which(xx > 0)]

  #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
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
  ##########################################
  cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                   '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                   '25', # possible transition clusters
                   '24', # also transitions cells and many of them are not annotated
                   #'44', '31', '52', '28', '50', # cluster '44', '31', '52', '28', '50' were not included here
                   '3', '5', '16', '30', '22', '4' # all middle time points
  )

  #cluster.sels = c('29', '32', '35', '40', '42')
  #sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])

  #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')
  ids.sels = c('MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap', 'MSxapp', 'MSxappp', 'MSpappa')

  # cells.sels = unique(colnames(seurat.obj)[seurat.obj$BWM.cells == 'BWM' &
  #   (!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
  #    !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
  #     )])
  # cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
  #                                               !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
  #                                            ])

  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

  #seurat.obj$BWM.cells[seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31'] = NA

  sub.obj = subset(seurat.obj, cells = cells.sels)
  xx = table(sub.obj$seurat_clusters)
  xx[which(xx > 0)]

  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

  barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)

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
    test.umap.params.for.BWM.cells(sub.obj, pdfname = 'BWM_early_cells_umap_param_test_2.pdf',
                                   nfeatures.sampling = c(500, 1000, 3000), nb.pcs.sampling = c(5, 10, 20),
                                   n.neighbors.sampling = c(5, 10, 30),
                                   min.dist.sampling = c(0.01, 0.05, 0.1)
                                   )
    toc()

  }

  nfeatures = 500;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 10;
  min.dist = 0.01; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5)

  #features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
  features.sels = c('fbxb-70', 'ceh-37', 'C45G7.4', 'nhr-67')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)


  p1 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend() +
    ggtitle('manual ids')
  p2 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = , repel = TRUE,  pt.size = 2) +
    NoLegend() + ggtitle('predicted ids scmap')
  p3 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,  pt.size = 2) +
    NoLegend() + ggtitle('predicted ids scmap')

  p1 + p2 + p3
  ##########################################
  # rerun the seurat for label transferring
  ##########################################
  RErun.seurat.transferring.labels = FALSE
  if(RErun.seurat.transferring.labels){
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

    #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
    sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30,
                                                                              k.anchor = 5, k.filter = 200,
                                                                              terminals = bwms.all)
    sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
    sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA

    p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
    p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      NoLegend() + ggtitle('predicted ids')
    #p2 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    #  NoLegend()
    p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      NoLegend() + ggtitle('manual ids')

    p2 + p3

    DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()


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
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p0 = DimPlot(sub.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)
  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 6,
               na.value = "gray", combine = TRUE)
  p3 =   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
  (p0 + p3) / (p1 + p2)

  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3,
          group.by = 'seurat_clusters_split')

  (p1 + p2) / p3

  features.sels = c('hnd-1',  # MSxp
                    'pha-4', 'fbxb-70', # MSxa
                    'ceh-37', 'C45G7.4', 'nhr-67', # MSxap
                    'cft-1', # MSxaa
                    'pat-9', 'ham-1'
                    )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  manual.discovery.new.features = TRUE
  if(manual.discovery.new.features){
    Idents(sub.obj) = sub.obj$seurat_clusters_split
    markers <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

  }

  dev.off()

  FeaturePlot(sub.obj, reduction = 'umap', features = c('cft-1', 'fbxb-70', 'ceh-37'))

  #################################################################################################################################
  # update of manually annotated ids using marker genes and potential mapped labels from scmap or seurat
  # c('MSxppppx', 'MSxpppax', 'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 'MSxapppp', 'MSxapppa',
  # 'MSxappppx', 'MSxapppax', 'MSpappax')
  ##################################################################################################################################
  RErun.seurat.transferring.labels = FALSE
  if(RErun.seurat.transferring.labels){
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
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

  # check info in JM data for specific lineage
  ids.sel = c('MSxppppx')
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  extrack.markers.from.JM(markers = markers.JM, id = 'MSx', ntop = 5)

  top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

  # to find new marker genes
  top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  top.markers[top.markers$cluster == '13',]


  features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)


  VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split')

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  #save first the seurat projection result, one of best projection prediciotn
  seurat.obj$predicted.ids.seurat.keep = NA
  mm = match(colnames(sub.obj), colnames(seurat.obj))
  seurat.obj$predicted.ids.seurat.keep[mm] = sub.obj$predicted.ids
  #seurat.obj$manual.annot.ids.1 = seurat.obj$manual.annot.ids
  seurat.obj$manual.annot.ids = seurat.obj$manual.annot.ids.1

  # integrate the good prediction from seurat
  jj = which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$predicted.ids.seurat.keep))
  seurat.obj$manual.annot.ids[jj] = seurat.obj$predicted.ids.seurat.keep[jj]
  kk = which(!is.na(match(seurat.obj$predicted.ids.seurat.keep, terminals)))
  seurat.obj$manual.annot.ids[kk] = seurat.obj$predicted.ids.seurat.keep[kk]

  seurat.obj$manual.annot.ids.2 = seurat.obj$manual.annot.ids # save the current annotation
  seurat.obj$manual.annot.ids[is.na(seurat.obj$BWM.cells)] = NA # keep only BWM cell annotation

  cluster.assingment = list(#c('0', 'MSxp'),
                            #c('1', ''), # don't change cluster 1
                            c('8', NA),
                            c('9', NA),
                            c('6', 'MSxapp'),
                            c('3', 'MSxap')
  )

  for(n in 1:length(cluster.assingment)){

    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    if(is.na(id2assign)) seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA

  }

  # #manual.assign.cluster.with.annotation(cluster.index = '17', id2assign = 'MSxppapp', sub.obj = sub.obj, seurat.obj = seurat.obj)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
           pt.size = 2)

  #
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
     scale_colour_hue(drop = FALSE) + NoLegend()

  #FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

  saveRDS(seurat.obj, file = RDS2save)



  ########################################################
    ########################################################
    # Section : iteration 20 (Second round)
    ## we will always check the mother, current generation and daughter cells together
    ## after cleaning up and we focus on :
    # MSxa, MSxp, MSxpp, MSxpa, MSxap, MSxapp, MSxappp, MSpappa
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
    markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
    # bwms = unique(c('MSx', 'MSxa', 'MSxap',
    #          'MSapaap', 'MSapaapp', 'MSappaaa',
    #          'MSpappa', 'MSpappax', 'MSppaap', 'MSppaapp', 'MSpppaaa',
    #          murray.ids[grep('MSxapp|MSxp', murray.ids)]))

    #markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
    #markers = markers[!is.na(match(markers$Lineage, bwms)), ]

    ##########################################
    # manual.annot.ids.1 -- inital manual annotation for early and middle time points before using terminal cells annotation from seruat
    # manual.annot.ids.2 -- save the manual annotation for non-BWM cells
    # manual.annot.ids.3 -- save the manual annotation before using seurat prediction for cells in middle time points
    ##########################################
    nb.iteration = 20
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

    ## 36 ids including MSx and redundant MSxppppx and MSxpppax; MSx missing and not necessary
    bwms.all = c('MSxppppx', 'MSxpppax', # redundant
                 'MSxppppp', 'MSxppppa', 'MSxpppaa', 'MSxpppap',
                  'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap',
                  'MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax', 'MSpappax', # all terminal cells
                  'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxppaa', 'MSxpapp', 'MSxpapa', 'MSxpaap', 'MSxpaaa', 'MSxappp', 'MSpappa',
                  'MSxppp', 'MSxppa', 'MSxpap', 'MSxpaa', 'MSxapp',
                  'MSxpp', 'MSxpa', 'MSxap',
                  'MSxp', 'MSxa')

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
    cluster.index = '25'
    table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

    table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

    table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

    xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxapp')])
    xx[which(xx > 0)]

    #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
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
    ##########################################
    cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                     '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                     '25', # possible transition clusters
                     '24', # also transitions cells and many of them are not annotated
                     #'44', '31', '52', '28', '50', # cluster '44', '31', '52', '28', '50' were not included here
                     '3', '5', '16', '30', '22', '4' # all middle time points
    )

    #cluster.sels = c('29', '32', '35', '40', '42')
    #sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
    ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))

    #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')
    ids.sels = c('MSxapp', 'MSxappp', 'MSpappa',
                 'MSxppp', 'MSxppa', 'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxppaa',
                 'MSxpap', 'MSxpaa', 'MSxpapp', 'MSxpapa', 'MSxpaap', 'MSxpaaa',
                 "MSxpaaap/MSxppapp/MSxpapap/MSxppap",
                 "MSxpapp/MSxpppp",
                 "MSxpapp/MSxpapa",
                 "MSxppap/MSxpaaa/MSxpaaap",
                 "MSxppap/MSxpaaa/MSxpapa/MSxpaaap",
                 "MSxppppa/MSxppppp/MSxpppaa/MSxpappa/MSxappp",
                 "MSxppppp/MSxappp"
                 #"MSxapppp/MSxapppa",
                 #"MSxapppa", "MSxapppp"
                 )
    ids.left = setdiff(ids.current, ids.sels)
    nchar(ids.left)
    # cells.sels = unique(colnames(seurat.obj)[seurat.obj$BWM.cells == 'BWM' &
    #   (!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
    #    !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
    #     )])
    # cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
    #                                               !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
    #                                            ])

    cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

    #seurat.obj$BWM.cells[seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31'] = NA

    sub.obj = subset(seurat.obj, cells = cells.sels)
    xx = table(sub.obj$seurat_clusters)
    xx[which(xx > 0)]

    sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
    sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
    DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

    barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)

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

      sub.obj$manual.annot.ids = sub.obj$predicted.ids.seurat.keep
      test.umap.params.for.BWM.cells(sub.obj, pdfname = 'BWM_middleCells_umap_param_test_3.pdf',
                                     nfeatures.sampling = c(500, 1000), nb.pcs.sampling = c(5, 10, 30),
                                     n.neighbors.sampling = c(10, 30),
                                     min.dist.sampling = c(0.01, 0.05, 0.1)
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
    min.dist = 0.1; spread = 1
    sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                       spread = spread, n.neighbors = n.neighbors,
                       min.dist = min.dist)
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()

    #features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
    features.sels = c(#'hnd-1',  'pha-4',
                      #'ceh-37', 'C45G7.4', # MSxap
                      #'ham-1',  'nhr-67', # MSxapp
                      'unc-120',
                      'unc-39', 'irx-1', 'egl-43'
                      #'unc-120', 'tab-1', 'rpm-1', # MSpappa
                      #'tbx-8', 'asic-2', 'skpo-1', # MSxappp
                      #'tbx-8', 'hlh-1', 'hnd-1', 'unc-39' # MSxpap
    )

    FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)


    p1 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                 pt.size = 2) + NoLegend() +
      ggtitle('manual ids')
    #p2 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = , repel = TRUE,  pt.size = 2) +
    #  NoLegend() + ggtitle('predicted ids scmap')
    p3 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat.keep', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                 pt.size = 2) +
      NoLegend() + ggtitle('predicted ids seurat')

    p1 + p3


    ##########################################
    # rerun the seurat for label transferring
    ##########################################
    RErun.seurat.transferring.labels = FALSE
    if(RErun.seurat.transferring.labels){
      source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

      #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
      sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30,
                                                                                k.anchor = 5, k.filter = 200,
                                                                                terminals = bwms.all)
      sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
      sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
      sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
      sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA

      p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
      p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
        NoLegend() + ggtitle('predicted ids')
      #p2 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      #  NoLegend()
      p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
        NoLegend() + ggtitle('manual ids')

      p2 + p3

      DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
      DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()


    }

    p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
    p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
    p0 + p1


    #FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'pha-4', 'hnd-1'))
    #FeaturePlot(seurat.obj, reduction = 'umap', features = c('unc-120', 'pha-4', 'hnd-1'))
    ##########################################
    # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
    ##########################################
    # FindClusters_subclusters = function(sub.obj, resolution = 0.4)
    # {
    #   sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    #   return(sub.obj$seurat_clusters)
    # }
    #
    # sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
    # sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.7)
    # DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)
    #
    # p0 = DimPlot(sub.obj, group.by = "predicted.ids.seurat.keep", reduction = 'umap', label = TRUE, repel = TRUE,
    #              pt.size = 2, label.size = 5) + NoLegend() + ggtitle('seurat.pred.ids')
    # p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    #   NoLegend() + ggtitle('manual.annot.ids')
    # p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
    #              label.size = 6,
    #              na.value = "gray", combine = TRUE)
    #
    # p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3,
    #         group.by = 'seurat_clusters_split')
    #
    # p0 + p1
    # p2 / p3
    #
    # manual.discovery.new.features = TRUE
    # if(manual.discovery.new.features){
    #   Idents(sub.obj) = sub.obj$seurat_clusters_split
    #   markers <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    #   top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    #   DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
    #
    # }
    #
    # dev.off()
    #
    # FeaturePlot(sub.obj, reduction = 'umap', features = c('cft-1', 'fbxb-70', 'ceh-37'))
    #
    # #################################################################################################################################
    # # check marker genes and potential mapped labels from scmap or seurat
    # # c('MSxppppx', 'MSxpppax', 'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 'MSxapppp', 'MSxapppa',
    # # 'MSxappppx', 'MSxapppax', 'MSpappax')
    # ##################################################################################################################################
    #
    # # check info in JM data for specific lineage
    # ids.sel = c('MSxppppx')
    # source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    # #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
    # extrack.markers.from.JM(markers = markers.JM, id = 'MSx', ntop = 5)
    #
    # top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    # DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
    #
    # # to find new marker genes
    # top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
    # top.markers[top.markers$cluster == '13',]
    #
    #
    # #features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
    # #FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
    #
    # VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split')

    ##########################################
    ##########################################
    # update the manual annotation if good marker genes or mapped ids were found
    ##########################################
    ##########################################
    # save current manual.annot.ids, early and middle time pooints from myself and terminal cells annotated by seurat
    seurat.obj$manual.annot.ids.3 = seurat.obj$manual.annot.ids
    jj1 = which(!is.na(sub.obj$predicted.ids.seurat.keep))
    sub.obj$manual.annot.ids[jj1] = sub.obj$predicted.ids.seurat.keep[jj1] # use the seurat prediction for cells in the middle time points
    sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxp')] = 'MSxpa'
    sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxap')] = 'MSxpp'

    seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

    # cluster.assingment = list(#c('0', 'MSxp'),
    #                           #c('1', ''), # don't change cluster 1
    #                           c('8', 'MSxpap'),
    #                           c('9', NA),
    #                           c('6', 'MSxapp'),
    #                           c('3', 'MSxap')
    # )
    #
    # for(n in 1:length(cluster.assingment)){
    #
    #   cluster.index = cluster.assingment[[n]][1];
    #   id2assign =  cluster.assingment[[n]][2];
    #   cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    #   cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    #   sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    #   seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    #   if(is.na(id2assign)) seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
    #
    # }

    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
             pt.size = 2)

    DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
             na.value = "gray") +
       ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
       scale_colour_hue(drop = FALSE) + NoLegend()

    #FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

    saveRDS(seurat.obj, file = RDS2save)



    ########################################################
  ########################################################
  # Section : iteration 21 (Second round)
  ## we will always check the mother, current generation and daughter cells together
  ## After using the seurat prediction for middle time points, we focus on
  # MSxap, MSxapp, MSxappp, MSpappa
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
  # bwms = unique(c('MSx', 'MSxa', 'MSxap',
  #          'MSapaap', 'MSapaapp', 'MSappaaa',
  #          'MSpappa', 'MSpappax', 'MSppaap', 'MSppaapp', 'MSpppaaa',
  #          murray.ids[grep('MSxapp|MSxp', murray.ids)]))

  #markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
  load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
  ##########################################
  # manual.annot.ids.1 -- inital manual annotation for early and middle time points before using terminal cells annotation from seruat
  # manual.annot.ids.2 -- save the manual annotation for non-BWM cells
  # manual.annot.ids.3 -- save the manual annotation before using seurat prediction for cells in middle time points
  ##########################################
  nb.iteration = 21
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

  ## 38 ids including MSx and redundant MSxppppx and MSxpppax;
  # MSx missing and not necessary; 'MSapaap', 'MSppaap' for 'MSxpaap'
  bwms.all = c('MSxppppx', 'MSxpppax', # redundant
               'MSxppppp', 'MSxppppa', 'MSxpppaa', 'MSxpppap',
                'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap',
                'MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax', 'MSpappax', # all terminal cells
                'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxppaa', 'MSxpapp', 'MSxpapa',  'MSapaap', 'MSppaap', #'MSxpaap',
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
  cluster.index = '25'
  table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

  table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

  table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

  xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSpappax')])
  xx[which(xx > 0)]

  xx = table(seurat.obj$seurat_clusters[which(seurat.obj$pred.ids.seurat.keep.bwm.all == 'MSpappax')])
  xx[which(xx > 0)]
  #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
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
  ##########################################
  cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                   '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                   '25', # possible transition clusters
                   '24', # also transitions cells and many of them are not annotated
                   #'44', '31', '52', '28', '50', # cluster '44', '31', '52', '28', '50' were not included here
                   '3', '5', '16', '30', '22', '4' # all middle time points
  )

  #cluster.sels = c('29', '32', '35', '40', '42')
  #sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))

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

  ids.sels = c('MSxapp', 'MSxappp', 'MSpappa')
  ids.left = setdiff(ids.current, ids.sels)
  nchar(ids.left)
  # cells.sels = unique(colnames(seurat.obj)[seurat.obj$BWM.cells == 'BWM' &
  #   (!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
  #    !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
  #     )])
  # cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
  #                                               !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
  #                                            ])

  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

  #seurat.obj$BWM.cells[seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31'] = NA

  sub.obj = subset(seurat.obj, cells = cells.sels)
  xx = table(sub.obj$seurat_clusters)
  xx[which(xx > 0)]

  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

  barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)

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

    #sub.obj$manual.annot.ids = sub.obj$predicted.ids.seurat.keep
    test.umap.params.for.BWM.cells(sub.obj, pdfname = 'BWM_middleCells_umap_param_convergence_branche_test_2.pdf',
                                   nfeatures.sampling = c(500, 1000), nb.pcs.sampling = c(5, 10, 30),
                                   n.neighbors.sampling = c(5, 10,  30),
                                   min.dist.sampling = c(0.01, 0.05, 0.1)
                                   )
    toc()

  }


  nfeatures = 1000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 5 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 5;
  min.dist = 0.1; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()


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

  # check info in JM data for specific lineage
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSpappa', group_2 = 'MSxappp',  ntop = 5)

  #features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
  features.sels = c('hnd-1',  'pha-4',
                    #'fbxb-70',
                    #'ceh-37', 'C45G7.4', # MSxap
                    #'tbx-8', 'hlh-1', 'hnd-1', 'unc-39' # MSxpap
                    'ham-1',  'nhr-67', # MSxapp
                    #'unc-120',
                    #'unc-39', 'irx-1', 'egl-43'
                    'unc-120', 'tab-1', 'rpm-1', 'F55C5.10', 'hil-3', # MSpappa
                    'tbx-8', 'asic-2', 'skpo-1', 'pxn-2', 'ost-1', 'pck-2', 'fkb-4', 'D1086.12' # MSxappp

  )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  ##########################################
  # rerun the seurat for label transferring
  ##########################################
  RErun.seurat.transferring.labels = FALSE
  if(RErun.seurat.transferring.labels){
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

    #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
    sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30,
                                                                              k.anchor = 5, k.filter = 200,
                                                                              terminals = bwms.all)
    sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
    sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA

    p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
    p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      NoLegend() + ggtitle('predicted ids')
    #p2 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    #  NoLegend()
    p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      NoLegend() + ggtitle('manual ids')

    p3 + p2

    DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

    # save the seurat prediction for all BWM cells
    #sub.obj$predicted.ids.seurat.keep.bwm.all = NA
    sub.obj$pred.ids.seurat.keep.bwm.all = sub.obj$predicted.ids.seurat.terminal
    seurat.obj$pred.ids.seurat.keep.bwm.all = NA
    seurat.obj$pred.ids.seurat.keep.bwm.all[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$pred.ids.seurat.keep.bwm.all

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

  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:5, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.2)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p0 = DimPlot(sub.obj, group.by = "predicted.ids.seurat.keep", reduction = 'umap', label = TRUE, repel = TRUE,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle('seurat.pred.ids')
  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    NoLegend() + ggtitle('manual.annot.ids')
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 6,
               na.value = "gray", combine = TRUE)

  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2", 'timingEst'), ncol = 3,
          group.by = 'seurat_clusters_split')

  p0 + p1
  # p2 / p3
  #
  # manual.discovery.new.features = TRUE
  # if(manual.discovery.new.features){
  #   Idents(sub.obj) = sub.obj$seurat_clusters_split
  #   markers <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
  #   top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  #   DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
  #
  # }
  #
  # dev.off()
  #
  # FeaturePlot(sub.obj, reduction = 'umap', features = c('cft-1', 'fbxb-70', 'ceh-37'))
  #
  # #################################################################################################################################
  # # check marker genes and potential mapped labels from scmap or seurat
  # # c('MSxppppx', 'MSxpppax', 'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 'MSxapppp', 'MSxapppa',
  # # 'MSxappppx', 'MSxapppax', 'MSpappax')
  # ##################################################################################################################################
  #
  # # check info in JM data for specific lineage
  # ids.sel = c('MSxppppx')
  # source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  # #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  # extrack.markers.from.JM(markers = markers.JM, id = 'MSx', ntop = 5)
  #
  # top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
  #
  # # to find new marker genes
  # top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  # top.markers[top.markers$cluster == '13',]
  #
  #
  # #features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
  # #FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  #
  # VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split')

  ##########################################
  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  ##########################################
  # save current manual.annot.ids, early and middle time pooints from myself and terminal cells annotated by seurat
  # seurat.obj$manual.annot.ids.3 = seurat.obj$manual.annot.ids
  # jj1 = which(!is.na(sub.obj$predicted.ids.seurat.keep))
  # sub.obj$manual.annot.ids[jj1] = sub.obj$predicted.ids.seurat.keep[jj1] # use the seurat prediction for cells in the middle time points
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxp')] = 'MSxpa'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxap')] = 'MSxpp'
  #
  # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

  cluster.assingment = list(#c('0', 'MSxp'),
                            #c('1', ''), # don't change cluster 1
                            c('0', 'MSpappa'),
                            c('1', 'MSxappp'),
                            c('2', 'MSxappp'),
                            c('3', 'MSxapp'),
                            c('4', 'MSxappp'),
                            c('5', 'MSxappp'),
                            c('6', 'MSxapp')
  )

  for(n in 1:length(cluster.assingment)){

    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    if(is.na(id2assign)) seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA

  }

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
           pt.size = 2)

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
     scale_colour_hue(drop = FALSE) + NoLegend()

  #FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

  saveRDS(seurat.obj, file = RDS2save)

  ##########################################
  # manual.annot.ids.1 -- inital manual annotation for early and middle time points before using terminal cells annotation from seruat
  # manual.annot.ids.2 -- saved the manual annotation for non-BWM cells
  # manual.annot.ids.3 -- saved the manual annotation before using seurat prediction for cells in middle time points
  # predicted.ids.seurat.keep -- saved the seurat prediction for middle and terminal cells
  # pred.ids.seurat.keep.bwm.all -- saved the seurat prediction for all bwm cells (early, middle and temrinal cells)
  ##########################################
  nb.iteration = 24
  Refine.annotated.ids = TRUE;

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
  cluster.index = '25'
  table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

  table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

  table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

  xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
  xx[which(xx > 0)]

  ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
  xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
  table(seurat.obj$manual.annot.ids[ii1])
  #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
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
  ##########################################
  cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                   '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                   '25', # possible transition clusters
                   '24', # also transitions cells and many of them are not annotated
                   #'44', '31', '52', '28', '50', # cluster '44', '31', '52', '28', '50' were not included here
                   '3', '5', '16', '30', '22', '4' # all middle time points
  )

  #cluster.sels = c('29', '32', '35', '40', '42')
  #sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  # cells.sels = unique(colnames(seurat.obj)[seurat.obj$BWM.cells == 'BWM' &
  #   (!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
  #    !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
  #     )])
  # cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
  #                                               !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
  #                                            ])

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
  ids.sels = c('MSxpp', 'MSxppp', 'MSxppa', 'MSxppaa')

  ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
  ids.left = setdiff(ids.current, ids.sels)
  nchar(ids.left)


  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])
  #seurat.obj$BWM.cells[seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31'] = NA

  sub.obj = subset(seurat.obj, cells = cells.sels)
  xx = table(sub.obj$seurat_clusters)
  xx[which(xx > 0)]

  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

  barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)

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
  # barplot(counts, main="cluster compositions by scmap ",
  #         xlab=NULL, col=c(1:nrow(counts)), las = 2,
  #         legend = rownames(counts))
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
    require(tictoc)
    tic()
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

    #sub.obj$manual.annot.ids = sub.obj$predicted.ids.seurat.keep
    test.umap.params.for.BWM.cells(sub.obj, pdfname = 'BWM_middleCells_umap_param_MSxpa_branch_test_1.pdf',
                                   nfeatures.sampling = c(500, 1000, 3000), nb.pcs.sampling = c(10, 30),
                                   n.neighbors.sampling = c(5, 10,  30),
                                   min.dist.sampling = c(0.01, 0.05, 0.1)
                                   )
    toc()

  }


  nfeatures = 500;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 5 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 30;
  min.dist = 0.1; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()

  #DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()


  p1 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
               pt.size = 2) + NoLegend() +
    ggtitle('manual ids')
  p2 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
               pt.size = 2) +
    NoLegend() + ggtitle('seurat.pred.ids.keep')

  p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids.3', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
               pt.size = 2) +
    NoLegend() + ggtitle('manual.annot.middle')

  p1 /(p2 + p3)

  # check info in JM data for specific lineage
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppaa', group_2 = NULL,   ntop = 5)

  #features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
  features.sels = c(#'hnd-1',  'pha-4',
                    #'fbxb-70',
                    #'ceh-37', 'C45G7.4', # MSxap
                    #'tbx-8', 'hlh-1', 'hnd-1', 'unc-39' # MSxpap
                    #'ham-1',  'nhr-67', # MSxapp
                    #'unc-39', 'irx-1', 'egl-43'
                    #'unc-120', 'tab-1', 'rpm-1', 'F55C5.10', 'hil-3', # MSpappa
                    #'tbx-8', 'asic-2', 'skpo-1', 'pxn-2', 'ost-1', 'pck-2', 'fkb-4', 'D1086.12' # MSxappp
                    #'egl-43', 'R11A5.3' #MSxpa
                    #'zip-7', 'hlh-1', 'col-111' # MSxpap
                    #'tbx-7', 'unc-120', 'ref-2', 'hlh-16', 'unc-39', 'ttr-50', 'clec-266' # MSxpaa
                    #'F41D9.2', 'hnd-1', 'unc-120', 'tbx-7', 'abts-1', 'Y66D12A.13' # MSxpaaa
                    #'let-381', 'ins-2', 'F40H3.3', 'ZK183.5', # MSxpapa
                    #'unc-120', 'ceh-51', 'lag-2', 'fbxb-22', 'T02G6.11', 'C06A8.3', 'T05D4.2',   # MSxpp
                    #'zip-7', 'hnd-1', 'tbx-8', 'fkh-2', 'tbx-11', # MSxppp
                    #'hlh-16', 'rgs-7', 'tbx-8', 'unc-120', 'T24C2.2', 'T24C2.3'  #MSxppa
                    'col-118', 'let-381', 'C03B1.1',  'T11B7.2' #MSxppaa

  )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  ##########################################
  # rerun the seurat for label transferring
  ##########################################
  RErun.seurat.transferring.labels = FALSE
  if(RErun.seurat.transferring.labels){
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

    #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
    sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30,
                                                                              k.anchor = 5, k.filter = 200,
                                                                              terminals = bwms.all)
    sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
    sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA

    p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
    p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      NoLegend() + ggtitle('predicted ids')
    #p2 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    #  NoLegend()
    p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      NoLegend() + ggtitle('manual ids')

    p3 + p2

    DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

    # save the seurat prediction for all BWM cells
    #sub.obj$predicted.ids.seurat.keep.bwm.all = NA
    sub.obj$pred.ids.seurat.keep.bwm.all = sub.obj$predicted.ids.seurat.terminal
    seurat.obj$pred.ids.seurat.keep.bwm.all = NA
    seurat.obj$pred.ids.seurat.keep.bwm.all[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$pred.ids.seurat.keep.bwm.all

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

  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:5, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p0 = DimPlot(sub.obj, group.by = "predicted.ids.seurat.keep", reduction = 'umap', label = TRUE, repel = TRUE,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle('seurat.pred.ids')
  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    NoLegend() + ggtitle('manual.annot.ids')

  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 6,
               na.value = "gray", combine = TRUE)

  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'seurat_clusters_split')

  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
               group.by = 'seurat_clusters_split') + NoLegend()

  p1 + p2
  p1 + p4
  p2 + p4

  plot(p3)

  manual.discovery.new.features = TRUE
  if(manual.discovery.new.features){
    Idents(sub.obj) = sub.obj$seurat_clusters_split
    markers.new <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    top.markers <- markers.new %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

  }

  dev.off()

  # FeaturePlot(sub.obj, reduction = 'umap', features = c('cft-1', 'fbxb-70', 'ceh-37'))
  #
  # #################################################################################################################################
  # # check marker genes and potential mapped labels from scmap or seurat
  # # c('MSxppppx', 'MSxpppax', 'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 'MSxapppp', 'MSxapppa',
  # # 'MSxappppx', 'MSxapppax', 'MSpappax')
  # ##################################################################################################################################
  #
  # # check info in JM data for specific lineage
  # ids.sel = c('MSxppppx')
  # source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  # #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  # extrack.markers.from.JM(markers = markers.JM, id = 'MSx', ntop = 5)
  #
  # top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
  #
  # # to find new marker genes
  # top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  # top.markers[top.markers$cluster == '13',]
  #
  #
  # #features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
  # #FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  #
  # VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split')

  ##########################################
  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  ##########################################
  # save current manual.annot.ids, early and middle time pooints from myself and terminal cells annotated by seurat
  # seurat.obj$manual.annot.ids.3 = seurat.obj$manual.annot.ids
  # jj1 = which(!is.na(sub.obj$predicted.ids.seurat.keep))
  # sub.obj$manual.annot.ids[jj1] = sub.obj$predicted.ids.seurat.keep[jj1] # use the seurat prediction for cells in the middle time points
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxp')] = 'MSxpa'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxap')] = 'MSxpp'
  #
  # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

  cluster.assingment = list(#c('0', 'MSxp'),
                             c('0', 'MSxppa'),
                             c('1', 'MSxppaa'),
                             c('2', 'MSxppp'),
                             c('3', 'MSxppaa'),
                             c('4', 'MSxppaa'),
                             c('5', 'MSxppa')
  )

  for(n in 1:length(cluster.assingment)){

    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    if(is.na(id2assign)) seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA

  }

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
           pt.size = 2)

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
     scale_colour_hue(drop = FALSE) + NoLegend()

  #FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

  saveRDS(seurat.obj, file = RDS2save)

  ########################################################
    ########################################################
    # Section : iteration 25 (Second round)
    ## we will always check the mother, current generation and daughter cells together
    ## After using the seurat prediction for middle time points, we focus on
    # MSxap, MSxapp, MSxappp, MSpappa
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
    # predicted.ids.seurat.keep -- saved the seurat prediction for middle and terminal cells
    # pred.ids.seurat.keep.bwm.all -- saved the seurat prediction for all bwm cells (early, middle and temrinal cells)
    ##########################################
    nb.iteration = 25
    Refine.annotated.ids = TRUE;

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
    cluster.index = '25'
    table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

    table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

    table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

    xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
    xx[which(xx > 0)]

    ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
    xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
    table(seurat.obj$manual.annot.ids[ii1])
    #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
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
    ##########################################
    cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
                     '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
                     '25', # possible transition clusters
                     '24', # also transitions cells and many of them are not annotated
                     #'44', '31', '52', '28', '50', # cluster '44', '31', '52', '28', '50' were not included here
                     '3', '5', '16', '30', '22', '4' # all middle time points
    )

    #cluster.sels = c('29', '32', '35', '40', '42')
    #sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
    # cells.sels = unique(colnames(seurat.obj)[seurat.obj$BWM.cells == 'BWM' &
    #   (!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
    #    !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
    #     )])
    # cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
    #                                               !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
    #                                            ])

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
    ids.sels = c('MSxpppa', 'MSxpppp', 'MSxpapa', 'MSxpapp', 'MSxppap')

    ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
    ids.left = setdiff(ids.current, ids.sels)
    nchar(ids.left)


    cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])
    #seurat.obj$BWM.cells[seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31'] = NA

    sub.obj = subset(seurat.obj, cells = cells.sels)
    xx = table(sub.obj$seurat_clusters)
    xx[which(xx > 0)]

    sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
    sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
    DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group)

    barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)

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
    # barplot(counts, main="cluster compositions by scmap ",
    #         xlab=NULL, col=c(1:nrow(counts)), las = 2,
    #         legend = rownames(counts))
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
      require(tictoc)
      tic()
      source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

      #sub.obj$manual.annot.ids = sub.obj$predicted.ids.seurat.keep
      test.umap.params.for.BWM.cells(sub.obj, pdfname = 'BWM_middleCells_umap_param_MSxpppp_MSxpppa_MSxpapa_MSxpapp_test_2.pdf',
                                     nfeatures.sampling = c(500, 1000, 3000, 5000), nb.pcs.sampling = c(10, 20),
                                     n.neighbors.sampling = c(5, 10, 20, 30),
                                     min.dist.sampling = c(0.01, 0.05, 0.1)
                                     )
      toc()

    }


    nfeatures = 500;
    sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
    #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
    sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
    sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
    ElbowPlot(sub.obj, ndims = 50)

    nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
    n.neighbors = 10;
    min.dist = 0.01; spread = 1
    sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                       spread = spread, n.neighbors = n.neighbors,
                       min.dist = min.dist)
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()

    #DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()


    p1 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                 pt.size = 2) + NoLegend() +
      ggtitle('manual ids')
    p2 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                 pt.size = 2) +
      NoLegend() + ggtitle('seurat.pred.ids.keep')

    p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids.3', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                 pt.size = 2) +
      NoLegend() + ggtitle('manual.annot.middle')

    p1 /(p2 + p3)


    ##########################################
    # rerun the seurat for label transferring
    ##########################################
    RErun.seurat.transferring.labels = FALSE
    if(RErun.seurat.transferring.labels){
      source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

      #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
      sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30,
                                                                                k.anchor = 5, k.filter = 200,
                                                                                terminals = bwms.all)
      sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
      sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
      sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
      sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA

      p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
      p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
        NoLegend() + ggtitle('predicted ids')
      #p2 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      #  NoLegend()
      p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
        NoLegend() + ggtitle('manual ids')

      p3 + p2

      DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
      DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

      # save the seurat prediction for all BWM cells
      #sub.obj$predicted.ids.seurat.keep.bwm.all = NA
      sub.obj$pred.ids.seurat.keep.bwm.all = sub.obj$predicted.ids.seurat.terminal
      seurat.obj$pred.ids.seurat.keep.bwm.all = NA
      seurat.obj$pred.ids.seurat.keep.bwm.all[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$pred.ids.seurat.keep.bwm.all

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

    sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 20, dims = 1:10, compute.SNN = TRUE)
    sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
    DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

    p0 = DimPlot(sub.obj, group.by = "predicted.ids.seurat.keep", reduction = 'umap', label = TRUE, repel = TRUE,
                 pt.size = 2, label.size = 5) + NoLegend() + ggtitle('seurat.pred.ids')
    p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      NoLegend() + ggtitle('manual.annot.ids')

    p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
                 label.size = 6,
                 na.value = "gray", combine = TRUE)

    p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
            group.by = 'seurat_clusters_split')

    p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                 group.by = 'seurat_clusters_split') + NoLegend()

    p1 + p2

    p1 + p4
    p2 + p4
    plot(p3)

    # check info in JM data for specific lineage
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppap', group_2 = 'MSxpppp',   ntop = 10)

    features.sels = c(#'hnd-1',  'pha-4',
      #'fbxb-70',
      #'ceh-37', 'C45G7.4', # MSxap
      #'tbx-8', 'hlh-1', 'hnd-1', 'unc-39' # MSxpap
      #'ham-1',  'nhr-67', # MSxapp
      #'unc-39', 'irx-1', 'egl-43'
      #'unc-120', 'tab-1', 'rpm-1', 'F55C5.10', 'hil-3', # MSpappa
      #'tbx-8', 'asic-2', 'skpo-1', 'pxn-2', 'ost-1', 'pck-2', 'fkb-4', 'D1086.12' # MSxappp
      #'egl-43', 'R11A5.3' #MSxpa
      #'zip-7', 'hlh-1', 'col-111' # MSxpap
      #'tbx-7', 'unc-120', 'ref-2', 'hlh-16', 'unc-39', 'ttr-50', 'clec-266' # MSxpaa
      #'F41D9.2', 'hnd-1', 'unc-120', 'tbx-7', 'abts-1', 'Y66D12A.13' # MSxpaaa
      #'let-381', 'ins-2', 'F40H3.3', 'ZK183.5', # MSxpapa
      #'unc-120', 'ceh-51', 'lag-2', 'fbxb-22', 'T02G6.11', 'C06A8.3', 'T05D4.2',   # MSxpp
      #'zip-7', 'hnd-1', 'tbx-8', 'fkh-2', 'tbx-11', # MSxppp
      #'hlh-16', 'rgs-7', 'tbx-8', 'unc-120', 'T24C2.2', 'T24C2.3'  #MSxppa
      #'col-118', 'let-381', 'C03B1.1',  'T11B7.2' #MSxppaa
      #'F40H3.3', 'Y42H9B.3', 'unc-120', 'ZK180.5', 'sfrp-1', 'ceh-34', 'irx-1', 'Y116A8C.3',
      #'let-381', 'ccg-1' # MSxpapa
      #'zig-8', 'pat-9', 'F19C6.4', 'tbx-8', 'her-1' # MSxpapp
      #'sul-2', 'camt-1', 'irx-1', 'ctg-1',  # MSxpppa
      #'tbx-2', 'D1086.12' #MSxpppp
      'ceh-34', 'ham-1', 'F19C6.4', 'eya-1', 'shc-2', 'F07C6.4', 'ten-1', 'hil-7' # MSxppap
    )
    FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)



    manual.discovery.new.features = TRUE
    if(manual.discovery.new.features){
      Idents(sub.obj) = sub.obj$seurat_clusters_split
      markers.new <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
      top.markers <- markers.new %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
      DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

    }

    dev.off()

    # FeaturePlot(sub.obj, reduction = 'umap', features = c('cft-1', 'fbxb-70', 'ceh-37'))
    #
    # #################################################################################################################################
    # # check marker genes and potential mapped labels from scmap or seurat
    # # c('MSxppppx', 'MSxpppax', 'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 'MSxapppp', 'MSxapppa',
    # # 'MSxappppx', 'MSxapppax', 'MSpappax')
    # ##################################################################################################################################
    #
    # # check info in JM data for specific lineage
    # ids.sel = c('MSxppppx')
    # source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    # #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
    # extrack.markers.from.JM(markers = markers.JM, id = 'MSx', ntop = 5)
    #
    # top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    # DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
    #
    # # to find new marker genes
    # top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
    # top.markers[top.markers$cluster == '13',]
    #
    #
    # #features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
    # #FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
    #
    # VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split')

    ##########################################
    ##########################################
    # update the manual annotation if good marker genes or mapped ids were found
    ##########################################
    ##########################################
    # save current manual.annot.ids, early and middle time pooints from myself and terminal cells annotated by seurat
    # seurat.obj$manual.annot.ids.3 = seurat.obj$manual.annot.ids
    # jj1 = which(!is.na(sub.obj$predicted.ids.seurat.keep))
    # sub.obj$manual.annot.ids[jj1] = sub.obj$predicted.ids.seurat.keep[jj1] # use the seurat prediction for cells in the middle time points
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxp')] = 'MSxpa'
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxap')] = 'MSxpp'
    #
    # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

    cluster.assingment = list(#c('0', 'MSxp'),
                               #c('0', 'MSxpppa'),
                               c('1', 'MSxppap')
                               # c('2', 'MSxpppp'),
                               # c('3', 'MSxpapa'),
                               # c('4', 'MSxpppp'),
                               # c('5', 'MSxpppp'),
                               # c('6', 'MSxpapa')
    )

    for(n in 1:length(cluster.assingment)){

      cluster.index = cluster.assingment[[n]][1];
      id2assign =  cluster.assingment[[n]][2];
      cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
      cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
      sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
      seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
      if(is.na(id2assign)) seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA

    }

    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
             pt.size = 2)

    DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
             na.value = "gray") +
       ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
       scale_colour_hue(drop = FALSE) + NoLegend()

    #FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))

    saveRDS(seurat.obj, file = RDS2save)





########################################################
########################################################
# Section : iteration 26 (Second round)
## we will always check the mother, current generation and daughter cells together
## After using the seurat prediction for middle time points, we focus on
# MSxap, MSxapp, MSxappp, MSpappa
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
# predicted.ids.seurat.keep -- saved the seurat prediction for middle and terminal cells
# pred.ids.seurat.keep.bwm.all -- saved the seurat prediction for all bwm cells (early, middle and temrinal cells)
##########################################
nb.iteration = 26
Refine.annotated.ids = TRUE;

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
cluster.index = '25'
table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
xx[which(xx > 0)]

ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
table(seurat.obj$manual.annot.ids[ii1])
#FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
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
ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
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
# ids.sels = c('MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax',
#              'MSxppppx', 'MSxpppax', 'MSxppppa', 'MSxppppp', 'MSxpppaa', 'MSxpppap', 'MSxpappa',
#              'MSxppapp', 'MSxpapap', 'MSxpappp', 'MSxpaaap')

ids.sels = c('MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax')
ids.left = setdiff(ids.current, ids.sels)
nchar(ids.left)


cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])
#seurat.obj$BWM.cells[seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31'] = NA

sub.obj = subset(seurat.obj, cells = cells.sels)
xx = table(sub.obj$seurat_clusters)
xx[which(xx > 0)]

sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)

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

  # jj1 = grep('likely.nonBWM_origCluster_17|MSxapppp/MSxapppa|MSxpaaap/MSxpapap|MSxppppx/MSxpppaa/MSxpappa',
  #            sub.obj$manual.annot.ids)
  # sub.obj$manual.annot.ids[jj1] = sub.obj$predicted.ids.seurat.keep[jj1]
  #
  require(tictoc)
  tic()
  test.umap.params.for.BWM.cells(sub.obj, pdfname = 'BWM_middleCells_umap_param_terminalCell_and_mothers_MSxappp_MSxpppp_MSxpppa.pdf',
                                 nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(5, 10, 30),
                                 n.neighbors.sampling = c(5, 10, 30),
                                 min.dist.sampling = c(0.01, 0.05, 0.1)
                                 )
  toc()

}

nfeatures = 1000;
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
#cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(sub.obj, ndims = 50)

nb.pcs = 5 # nb of pcs depends on the considered clusters or ids
n.neighbors = 10;
min.dist = 0.01; spread = 1
sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs,
                   spread = spread, n.neighbors = n.neighbors,
                   min.dist = min.dist)
DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)

#DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()

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


##########################################
# rerun the seurat for label transferring
##########################################
RErun.seurat.transferring.labels = FALSE
if(RErun.seurat.transferring.labels){
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

  #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
  sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30,
                                                                            k.anchor = 5, k.filter = 100,
                                                                            terminals = terminals)
  sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
  sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob
  sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
  sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA

  p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
  p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
   ggtitle('predicted ids')
  #p2 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
  #  NoLegend()
  p3 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
   ggtitle('manual ids')

  p3 + p2

  DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()

  # save the seurat prediction for all BWM cells
  #sub.obj$predicted.ids.seurat.keep.bwm.all = NA
  sub.obj$pred.ids.seurat.keep.bwm.all = sub.obj$predicted.ids.seurat.terminal
  seurat.obj$pred.ids.seurat.keep.bwm.all = NA
  seurat.obj$pred.ids.seurat.keep.bwm.all[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$pred.ids.seurat.keep.bwm.all

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
sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

p0 = DimPlot(sub.obj, group.by = "predicted.ids.seurat.keep", reduction = 'umap', label = TRUE, repel = TRUE,
             pt.size = 2, label.size = 5) + NoLegend() + ggtitle('seurat.pred.ids')
p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
   ggtitle('manual.annot.ids')

p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
             label.size = 6,
             na.value = "gray", combine = TRUE)

p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
        group.by = 'seurat_clusters_split')

p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
             group.by = 'seurat_clusters_split') + NoLegend()

p1 + p2

p1 + p4
p2 + p4
plot(p3)

# check info in JM data for specific lineage
source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxapppp', group_2 = NULL, ntop = 10)

features.sels = c(#'hnd-1',  'pha-4',
  #'fbxb-70',
  #'ceh-37', 'C45G7.4', # MSxap
  #'tbx-8', 'hlh-1', 'hnd-1', 'unc-39' # MSxpap
  #'ham-1',  'nhr-67', # MSxapp
  #'unc-39', 'irx-1', 'egl-43'
  #'unc-120', 'tab-1', 'rpm-1', 'F55C5.10', 'hil-3' # MSpappa
  #'tbx-8', 'asic-2', 'skpo-1', 'pxn-2', 'D1086.12','tab-1'  # MSxappp
  'mnr-1', 'strl-1', 'zig-10', 'gana-1', 'fbxb-88', # MSxapppa
  'lin-39', 'F54D5.5'  # MSxapppp
  #'gsnl-1', 'abts-1', 'stn-2', # MSxapppax
  #'lin-39',  'tbx-2', 'B0379.1'  #MSxappppx
  #'egl-43', 'R11A5.3' #MSxpa
  #'zip-7', 'hlh-1', 'col-111' # MSxpap
  #'tbx-7', 'unc-120', 'ref-2', 'hlh-16', 'unc-39', 'ttr-50', 'clec-266' # MSxpaa
  #'F41D9.2', 'hnd-1', 'unc-120', 'tbx-7', 'abts-1', 'Y66D12A.13' # MSxpaaa
  #'let-381', 'ins-2', 'F40H3.3', 'ZK183.5', # MSxpapa
  #'unc-120', 'ceh-51', 'lag-2', 'fbxb-22', 'T02G6.11', 'C06A8.3', 'T05D4.2',   # MSxpp
  #'zip-7', 'hnd-1', 'tbx-8', 'fkh-2', 'tbx-11', # MSxppp
  #'hlh-16', 'rgs-7', 'tbx-8', 'unc-120', 'T24C2.2', 'T24C2.3'  #MSxppa
  #'col-118', 'let-381', 'C03B1.1',  'T11B7.2' #MSxppaa
  #'F40H3.3', 'Y42H9B.3', 'unc-120', 'ZK180.5', 'sfrp-1', 'ceh-34', 'irx-1', 'Y116A8C.3',
  #'let-381', 'ccg-1' # MSxpapa
  #'zig-8', 'pat-9', 'F19C6.4', 'tbx-8', 'her-1' # MSxpapp
  #'sul-2', 'camt-1', 'irx-1', 'ctg-1',  # MSxpppa
  #'tbx-2', 'D1086.12' #MSxpppp
  #'ceh-34', 'ham-1', 'F19C6.4', 'eya-1', 'shc-2', 'F07C6.4', 'ten-1', 'hil-7' # MSxppap
)

FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)


manual.discovery.new.features = TRUE
if(manual.discovery.new.features){
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  markers.new <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
  top.markers <- markers.new %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

}

dev.off()

# FeaturePlot(sub.obj, reduction = 'umap', features = c('cft-1', 'fbxb-70', 'ceh-37'))
#
# #################################################################################################################################
# # check marker genes and potential mapped labels from scmap or seurat
# # c('MSxppppx', 'MSxpppax', 'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 'MSxapppp', 'MSxapppa',
# # 'MSxappppx', 'MSxapppax', 'MSpappax')
# ##################################################################################################################################
#
# # check info in JM data for specific lineage
# ids.sel = c('MSxppppx')
# source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
# #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
# extrack.markers.from.JM(markers = markers.JM, id = 'MSx', ntop = 5)
#
# top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
# DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
#
# # to find new marker genes
# top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
# top.markers[top.markers$cluster == '13',]
#
#
# #features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
# #FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
#
# VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split')

##########################################
##########################################
# update the manual annotation if good marker genes or mapped ids were found
##########################################
##########################################
# save current manual.annot.ids, early and middle time pooints from myself and terminal cells annotated by seurat
# seurat.obj$manual.annot.ids.3 = seurat.obj$manual.annot.ids
# jj1 = which(!is.na(sub.obj$predicted.ids.seurat.keep))
# sub.obj$manual.annot.ids[jj1] = sub.obj$predicted.ids.seurat.keep[jj1] # use the seurat prediction for cells in the middle time points
# sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxp')] = 'MSxpa'
# sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxap')] = 'MSxpp'
#
# seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

#mm = match(colnames(sub.obj), colnames(seurat.obj))
#seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids


cluster.assingment = list(#c('0', 'MSxp'),
                           c('0', 'MSxappppx'),
                           c('1', 'MSxapppax'),
                           c('2', 'MSxapppp'),
                           c('3', 'MSxapppa'),
                           c('4', 'MSxapppp'),
                           c('5', 'MSxapppax')
                            #c('6', 'MSxpapa')
)

for(n in 1:length(cluster.assingment)){

  cluster.index = cluster.assingment[[n]][1];
  id2assign =  cluster.assingment[[n]][2];
  cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
  sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
  if(is.na(id2assign)) seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA

}

DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
         pt.size = 2)

DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
   scale_colour_hue(drop = FALSE) + NoLegend()

#FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))
saveRDS(seurat.obj, file = RDS2save)


########################################################
  ########################################################
  # Section : iteration 27 (Second round)
  ## we will always check the mother, current generation and daughter cells together
  ## After using the seurat prediction for middle time points, we focus on

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
  # predicted.ids.seurat.keep -- saved the seurat prediction for middle and terminal cells
  # pred.ids.seurat.keep.bwm.all -- saved the seurat prediction for all bwm cells (early, middle and temrinal cells)
  ##########################################
  nb.iteration = 27
  Refine.annotated.ids = TRUE;

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
  cluster.index = '25'
  table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

  table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')

  table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]

  xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
  xx[which(xx > 0)]

  ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
  xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
  table(seurat.obj$manual.annot.ids[ii1])
  #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
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
  ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
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
  ids.sels = c(#'MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax',
                'MSxppppx', 'MSxpppax', 'MSxppppa', 'MSxppppp', 'MSxpppaa', 'MSxpppap', 'MSxpappa',
                'MSxppapp', 'MSxpapap', 'MSxpappp', 'MSxpaaap')
  #ids.sels = c('MSxapppp', 'MSxapppa', 'MSxappppx', 'MSxapppax')
  #ids.sels = c('MSxppppx', 'MSxpppax', 'MSxppppa', 'MSxppppp', 'MSxpppaa', 'MSxpppap', 'MSxpappa')
  #ids.sels = c('MSxppapp', 'MSxpapap', 'MSxpappp', 'MSxpaaap', 'terminal.from.clusters.19.27',
  #             'MSxppap')

  ids.left = setdiff(ids.current, ids.sels)
  nchar(ids.left)

  #jj1 = which(!is.na(match(seurat.obj$manual.annot.ids, ids.sels)) & (seurat.obj$seurat_clusters == '19'|
  #                                                                    seurat.obj$seurat_clusters == '27'))
  #seurat.obj$manual.annot.ids[jj1] = 'terminal.from.clusters.19.27'
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])
  #seurat.obj$BWM.cells[seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31'] = NA

  sub.obj = subset(seurat.obj, cells = cells.sels)
  xx = table(sub.obj$seurat_clusters)
  xx[which(xx > 0)]

  sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

  barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)

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

    # jj1 = grep('likely.nonBWM_origCluster_17|MSxapppp/MSxapppa|MSxpaaap/MSxpapap|MSxppppx/MSxpppaa/MSxpappa',
    #            sub.obj$manual.annot.ids)
    # sub.obj$manual.annot.ids[jj1] = sub.obj$predicted.ids.seurat.keep[jj1]
    #
    require(tictoc)
    tic()
    test.umap.params.for.BWM.cells(sub.obj, pdfname = 'BWM_middleCells_umap_param_terminalCell_and_mothers_MSxpp_MSxpa_1.pdf',
                                   nfeatures.sampling = c(1000, 3000, 5000, 80000), nb.pcs.sampling = c(10, 20, 30, 50),
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

  nb.pcs = 30 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 10;
  min.dist = 0.01; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist)

  # idntt = c('MSxpappp'); cells_to_show <- list(c(WhichCells(sub.obj, idents = idntt)))
  # DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
  #         cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1,
  #         pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
  #
  # idntt = c('MSxpapap'); cells_to_show <- list(c(WhichCells(sub.obj, idents = idntt)))
  # DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
  #         cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1.5,
  #         pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)

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

    #ids.refs = c(terminals, 'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxpapp', 'MSxpapa', 'Msxpaaa')
    ids.refs = terminals
    #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
    sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30,
                                                                              reduction = 'pcaproject',
                                                                              k.anchor = 30, k.filter = 200, max.features = 200,
                                                                              terminals = ids.refs)
    sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
    sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob

    hist(sub.obj$predicted.ids.prob)
    abline(v = 0.5, col = 'red', lwd=2.0)

    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    sub.obj = reference.based.cell.projection.rf.svm(sub.obj, nfeatures = 3000, cost = 0.1, ntree = 50, scale.cell = TRUE,
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

    DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
    #DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE)

    # save the seurat prediction for all BWM cells
    #sub.obj$predicted.ids.seurat.keep.bwm.all = NA
    sub.obj$pred.ids.seurat.keep.bwm.all = sub.obj$predicted.ids.seurat.terminal
    seurat.obj$pred.ids.seurat.keep.bwm.all = NA
    seurat.obj$pred.ids.seurat.keep.bwm.all[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$pred.ids.seurat.keep.bwm.all

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
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:30, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 2.0)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p0 = DimPlot(sub.obj, group.by = "predicted.ids.seurat.keep", reduction = 'umap', label = TRUE, repel = TRUE,
               pt.size = 2, label.size = 5) + NoLegend() + ggtitle('seurat.pred.ids')
  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
     ggtitle('manual.annot.ids')

  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 6,
               na.value = "gray", combine = TRUE)

  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'seurat_clusters_split')

  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
               group.by = 'seurat_clusters_split') + NoLegend()

  p1 + p2

  p1 + p4
  p2 + p4
  plot(p3)

  Idents(sub.obj) = sub.obj$manual.annot.ids
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

  # MSxpaaap
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpaaap', group_2 = 'MSxppppp', ntop = 10)

  features.sels = c('maph-1.3', 'maph-1.2', 'ZC449.5', 'shc-2')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  # MSxppapp
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppapp', group_2 = 'MSxpaaap', ntop = 10)

  features.sels = c('K09G1.1', 'tbx-2', 'D1086.12')
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)



  manual.discovery.new.features = FALSE
  if(manual.discovery.new.features){
    Idents(sub.obj) = sub.obj$seurat_clusters_split
    markers.new <- FindAllMarkers(sub.obj, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
    top.markers <- markers.new %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

  }

  dev.off()

  # FeaturePlot(sub.obj, reduction = 'umap', features = c('cft-1', 'fbxb-70', 'ceh-37'))
  #
  # #################################################################################################################################
  # # check marker genes and potential mapped labels from scmap or seurat
  # # c('MSxppppx', 'MSxpppax', 'MSxppapp', 'MSxpappp', 'MSxpappa', 'MSxpapap', 'MSxpaaap', 'MSxapppp', 'MSxapppa',
  # # 'MSxappppx', 'MSxapppax', 'MSpappax')
  # ##################################################################################################################################
  #
  # # check info in JM data for specific lineage
  # ids.sel = c('MSxppppx')
  # source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  # #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  # extrack.markers.from.JM(markers = markers.JM, id = 'MSx', ntop = 5)
  #
  # top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  # DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
  #
  # # to find new marker genes
  # top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  # top.markers[top.markers$cluster == '13',]
  #
  #
  # #features.sels = c('hnd-1', 'pha-4', 'fbxb-70', 'ceh-37', 'C45G7.4', 'pat-9', 'nhr-67', 'unc-120', 'unc-39', 'irx-1', 'egl-43')
  # #FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
  #
  # VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split')

  ##########################################
  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  ##########################################
  # save current manual.annot.ids, early and middle time pooints from myself and terminal cells annotated by seurat
  # seurat.obj$manual.annot.ids.3 = seurat.obj$manual.annot.ids
  # jj1 = which(!is.na(sub.obj$predicted.ids.seurat.keep))
  # sub.obj$manual.annot.ids[jj1] = sub.obj$predicted.ids.seurat.keep[jj1] # use the seurat prediction for cells in the middle time points
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxp')] = 'MSxpa'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'MSxap')] = 'MSxpp'
  #
  # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

  #mm = match(colnames(sub.obj), colnames(seurat.obj))
  #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids


  cluster.assingment = list(#c('0', 'MSxp'),
                             c('0', 'MSxappppx'),
                             c('1', 'MSxapppax'),
                             c('2', 'MSxapppp'),
                             c('3', 'MSxapppa'),
                             c('4', 'MSxapppp'),
                             c('5', 'MSxapppax')
                              #c('6', 'MSxpapa')
  )

  for(n in 1:length(cluster.assingment)){

    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    if(is.na(id2assign)) seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA

  }

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
           pt.size = 2)

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
     scale_colour_hue(drop = FALSE) + NoLegend()

  #FeaturePlot(seurat.obj, reduction = 'umap', features = c('hot-1', 'wago-1', 'pde-6', 'rrc-1', 'maph-1.2'))
  saveRDS(seurat.obj, file = RDS2save)

  
