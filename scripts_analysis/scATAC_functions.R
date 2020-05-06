##########################################################################
##########################################################################
# Project:
# Script purpose: main functions for scATAC-seq data analysis
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Mar  6 10:24:10 2020
##########################################################################
##########################################################################

########################################################
########################################################
# Section :
# functions from Andrew Hill
# andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html
# Andrew's cistopic analysis was done with cisTopic 0.2 
########################################################
########################################################
options(stringsAsFactors = FALSE)
library(Matrix)
library(data.table)
library(tictoc)
library(DelayedArray)
library(Seurat)
library(ggplot2)
library(irlba)
library(patchwork)
library(plyr)
library(dplyr)
library(stringr)
#library(SnapATAC)
library(GenomicRanges)
library(cisTopic)

########################################################
########################################################
# Section : scATAC-seq peak annotation
# here using ChIPseeker
########################################################
########################################################
run_scATAC_peak_annotation = function(peak.file)
{
  ## loading packages
  library(ChIPseeker)
  library(TxDb.Celegans.UCSC.ce11.ensGene)
  #library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Celegans.UCSC.ce11.ensGene
  library(clusterProfiler)
  
  # peak.file = paste0(filtered_mtx_dir, '/peaks.bed')
  peak <- readPeakFile(peakfile = peak.file, as = 'GRanges')
  
  upstream = 4000; downstream = 4000
  promoter <- getPromoters(TxDb=txdb, upstream=upstream, downstream=downstream, by = 'gene')
  tagMatrix <- getTagMatrix(peak, windows=promoter)
  
  plotAvgProf(tagMatrix, xlim=c(-upstream, downstream),
              xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
  
  tagHeatmap(tagMatrix, xlim=c(-upstream, downstream), color="red")
  
  peakAnno <- annotatePeak(peak = peak, tssRegion=c(-5000, 5000), level = 'gene', TxDb=txdb)
  
  #par(mfrow=c(1,2))
  plotAnnoPie(peakAnno)
  #plotDistToTSS(peakAnno, title="Distribution of transcription factor-binding loci\nrelative to TSS")
  #p1 + p2
  ## to speed up the compilation of this vignettes, we use a precalculated tagMatrix
  #data("tagMatrixList")
  #tagMatrix <- tagMatrixList[[4]]
  #promoter <- getPromoters(TxDb=txdb, upstream=5000, downstream=5000) #ChIPseeker method
  #tagMatrixList <- lapply(allPeaks, getTagMatrix, windows=promoter)
  
  #tagMatrixList <- tagMatrixList[!sapply(lapply(tagMatrixList, nrow), is.null)]
  #tagMatrixList[sapply(tagMatrixList, nrow) == 0] <- NULL
  
  #promoter <- promoters(genes(txdb), upstream=5000, downstream=5000) #GRanges
  #tagMatrixList <- lapply(allPeaks, getTagMatrix, windows=promoter)
  #plotAvgProf(tagMatrixList, xlim=c(-5000, 5000))
  #tagHeatmap(tagMatrixList, xlim=c(-5000, 5000), color=NULL)
  
  #peakAnnoList <- lapply(allPeaks, annotatePeak, TxDb=txdb, tssRegion=c(-2000, 200), verbose=FALSE)
  
  # for (n in names(peakAnnoList))
  # {
  #   par(mfrow=c(1,1))
  #   vennpie(peakAnnoList[[n]])
  #   text(x=0, y=-1, n)
  #   par(mfrow=c(1,1))
  #   plotAnnoPie(peakAnnoList[[n]])
  #   #        par(mfrow=c(1,1))
  #   #        upsetplot(peakAnnoList[[n]]) #only in TB-3.2.1-dev #  vennpie=TRUE,
  # }
  # 
}

########################################
# Utility functions
########################################
# Loads 10x dataset into a sparse matrix object
# Args:
#   matrix_fn (string): name of mtx or mtx.gz file coontaining MM formatted matrix
#   peaks_fn (string): name of peaks BED file (matches rows of matrix)
#   barcodes_fn (string): name of file containing cell barcodes (matches columns of matrix)
# Returns:
#   sparse matrix: matrix represented by files above with row and column names set (chr_start_stop format used for rownames)
load_tenx_atac = function(matrix_fn, peaks_fn, barcodes_fn) {
  atac_matrix = readMM(matrix_fn)
  colnames(atac_matrix) = read.delim(barcodes_fn, header=FALSE)$V1
  peaks = read.delim(peaks_fn, header=FALSE)
  peaks = paste(peaks$V1, peaks$V2, peaks$V3, sep = '_')
  rownames(atac_matrix) = peaks
  return(atac_matrix)
}

# Allows filtering of sites measured as non-zero in less than a given number of cells
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   cells (int): filter sites if they have less than this number of cells above zero
# Returns:
#   sparse matrix: filtered sparse matrix
filter_features = function(bmat, cells=10) {
  bmat = bmat[Matrix::rowSums(bmat) >= cells, ]
  return(bmat)
}

# Allows filtering of cells with below a given number of non-zero features
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   features_above_zero (int): filter cells if they have less than this number of features above zero
# Returns:
#   sparse matrix: filtered sparse matrix
filter_cells = function(bmat, features_above_zero=100) {
  bmat = bmat[, Matrix::colSums(bmat > 0) >= features_above_zero]
  return(bmat)
}

# Takes sparse matrix object and downsamples to a given fraction of entries remaining.
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
#   fraction_remaining (float): float (0, 1) that indicates fraction of non-zero entries to retain
#   cells_per_site_min (int): min cells a site must be measured in to retain the site in matrix
#   sites_per_cell_min (int): min sites a cell must have non-zero entries in to retain the cell in matrix
# Returns:
#   sparse matrix: downsampled sparse matrix
downsample_matrix = function(bmat, fraction_remaining=0.5, cells_per_site_min=1, sites_per_cell_min=1) {
  set.seed(2019)
  non_zero_entries = which(bmat@x > 0)
  indices_to_zero = sample(non_zero_entries, size=ceiling(length(non_zero_entries) * (1 - fraction_remaining)))
  bmat@x[indices_to_zero] = 0
  
  # Make sure to get rid of stuff that has gone to ~ 0 after downsampling
  bmat = filter_features(bmat, cells=cells_per_site_min)
  bmat = filter_cells(bmat, features_above_zero=sites_per_cell_min)
  return(bmat)
}

########################################
# Functions for LSI
########################################
# Helper function to do fast version of row scaling of sparse TF matrix by IDF vector.
# Exploits the way that data is stored within sparseMatrix object. Seems to be much more memory efficient than tf * idf and faster than DelayedArray.
# Args:
#   tf (sparse matrix): term frequency matrix
#   idf (vector): IDF vector
# Returns:
#   sparse matrix: TF-IDF matrix
safe_tfidf_multiply = function(tf, idf) {
  tf = t(tf)
  tf@x <- tf@x * rep.int(idf, diff(tf@p))
  tf = t(tf)
  return(tf)
}

# Perform TF-IDF on binary matrix
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
#   frequencies (bool): divide bmat by colSums (if FALSE simply use bmat for TF matrix)
#   log_scale_tf (bool): log scale TF matrix if TRUE
#   scale_factor (float): multiply terms in TF matrix by scale_factor prior to log1p. Equivalent to adding small pseudocount but doesn't cast to dense matrix at any point.
# Returns:
#   sparse matrix: TF-IDF matrix
tfidf = function(bmat, frequencies=TRUE, log_scale_tf=TRUE, scale_factor=100000) {
  # Use either raw counts or divide by total counts in each cell
  if (frequencies) {
    # "term frequency" method
    tf = t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf = bmat
  }
  
  # Either TF method can optionally be log scaled
  if (log_scale_tf) {
    if (frequencies) {
      tf@x = log1p(tf@x * scale_factor)
    } else {
      tf@x = log1p(tf@x * 1)
    }
  }
  
  # IDF w/ "inverse document frequency smooth" method
  idf = log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  
  # TF-IDF
  tf_idf_counts = safe_tfidf_multiply(tf, idf)
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  return(tf_idf_counts)
}

# Perform current version of TF-IDF used by 10x on binary matrix
# Args:
#   bmat (sparse matrix): sparse matrix to downsample
# Returns:
#   sparse matrix: TF-IDF matrix
tenx_tfidf = function(bmat) {
  idf = log(ncol(bmat) + 1) - log(1 + Matrix::rowSums(bmat))
  tf_idf_counts = safe_tfidf_multiply(bmat, idf)
  
  rownames(tf_idf_counts) = rownames(bmat)
  colnames(tf_idf_counts) = colnames(bmat)
  tf_idf_counts = as(tf_idf_counts, "sparseMatrix")
  return(tf_idf_counts)
}

# Perform fast PCA (irlba) on matrix, retaining observation names
# Args:
#   mat (sparse matrix): matrix to use for PCA (no further scaling or centering done)
#   dims (int): number of PCs to calculate
# Returns:
#   sparse matrix: TF-IDF matrix
do_pca = function(mat, dims=50) {
  pca.results = irlba(t(mat), nv=dims)
  final_result = pca.results$u %*% diag(pca.results$d)
  rownames(final_result) = colnames(mat)
  colnames(final_result) = paste0('PC_', 1:dims)
  return(final_result)
}

########################################
# Helper functions for dim reduction
########################################
# Wrapper for performing further dim reduction (tSNE/UMAP) and clustering given PCA space via Seurat.
# Args:
#   atac_matrix (sparse matrix): matrix to store in Seurat object (not used in computations)
#   cell_embeddings (matrix): typically PCA coordinates of cells but could be any set of reduced dim coordinates
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   metadata (dataframe): dataframe of metadata (rowonames are cell names) to add to Seurat object
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
# Returns:
#   Seurat object: seurat object
run_dim_reduction = function(atac_matrix, cell_embeddings, dims, metadata=NULL, reduction='pca.l2') {
  if (is.null(metadata)) {
    seurat_obj = Seurat::CreateSeuratObject(atac_matrix, assay = 'peaks')
  } else {
    seurat_obj = Seurat::CreateSeuratObject(atac_matrix, meta.data = metadata, assay = 'peaks')
  }
  
  seurat_obj[['pca']] = Seurat::CreateDimReducObject(embeddings=cell_embeddings, key='PC_', assay='peaks')
  seurat_obj = seurat_obj %>%
    Seurat::L2Dim(reduction='pca') %>%
    Seurat::RunUMAP(reduction = reduction, dims = dims) %>%
    Seurat::RunTSNE(reduction = reduction, dims = dims) %>%
    Seurat::FindNeighbors(reduction=reduction, nn.eps=0.25, dims=dims)
  return(seurat_obj)
}

# Helper function for plotting Spearman correlations of a given metadata column with all dimensions in a reduced space.
# Args:
#   seurat_obj (seurat object): Seurat object to use
#   reduction (string): name of reduction to use
#   column (string): name of column in metadata to use
# Returns:
#   ggplot object: plot object
plot_pc_correlation = function(seurat_obj, reduction, column='nCount_RNA') {
  coords = Seurat::Embeddings(seurat_obj, reduction=reduction)
  column_value = seurat_obj@meta.data[, column]
  correlations = apply(coords, 2, function(x) {cor(x, column_value, method='spearman')})
  correlations_df = data.frame(correlation=correlations, PC=1:ncol(coords))
  
  plot_obj = ggplot(correlations_df, aes(PC, correlation)) +
    geom_point() +
    theme_classic() +
    geom_hline(yintercept = 0, linetype='dashed', color='red')
  
  return(plot_obj)
}

########################################
# Wrapper functions for workflows
########################################
# Wrapper for full LSI workflow (TF-IDF and PCA + clustering + further dim reduction)
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   metadata: dataframe of metadata (rowonames are cell names) to add to Seurat object
#   log_scale_tf (bool): log scale TF matrix if TRUE
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#   resolution (float): resolution parameter to Seurat Louvain clustering
# Returns:
#   Seurat object: Seurat object. clustering + tSNE + UMAP done on PCA results from TF-IDF matrix.
lsi_workflow = function(bmat, dims, metadata=NULL, log_scale_tf=TRUE, reduction='pca.l2', resolution=0.3) {
  tfidf_mat = tfidf(bmat, log_scale_tf=log_scale_tf)
  pca_mat = do_pca(tfidf_mat, dims=max(dims))
  
  seurat_obj = run_dim_reduction(bmat, pca_mat, dims, metadata) %>%
    Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)
  return(seurat_obj)
}

# Wrapper for 10x version of full LSI workflow (TF-IDF and PCA + clustering + further dim reduction). Only TF-IDF step is modified.
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   metadata: dataframe of metadata (rownames are cell names) to add to Seurat object
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#   resolution (float): resolution parameter to Seurat Louvain clustering
# Returns:
#   Seurat object: Seurat object. clustering + tSNE + UMAP done on PCA results from TF-IDF matrix.
tenx_lsi_workflow = function(bmat, dims, metadata=NULL, reduction='pca.l2', resolution=0.3) {
  tfidf_mat = tenx_tfidf(bmat)
  pca_mat = do_pca(tfidf_mat, dims=max(dims))
  
  seurat_obj = run_dim_reduction(bmat, pca_mat, dims, metadata) %>%
    Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)
  return(seurat_obj)
}

# Runs cisTopic on binary matrix.
# Args:
#   bmat (sparse matrix): sparse matrix (binarized)
#   topic (vector of int): topic numbers to generate models for
# Returns:
#   cisTopic object: cisTopic object with models generated
cistopic_workflow = function(bmat, topic=seq(20, 120, by=10)) {
  coords = str_split_fixed(rownames(bmat), '_', 3)
  new_coords = paste0(coords[, 1], ':', coords[, 2], '-', coords[, 3])
  rownames(bmat) = new_coords
  
  cisTopicObject = cisTopic::createcisTopicObject(bmat, project.name='mouse_atac')
  
  #cisTopicObject = cisTopic::runModels(cisTopicObject, topic, seed=2019, nCores=4, burnin = 250, iterations = 500)
  cisTopicObject = cisTopic::runWarpLDAModels(cisTopicObject, topic = topic, seed=2019, nCores=4, iterations = 500, addModels = FALSE)
  #cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=nb.topcis, seed=987, nCores=6, iterations = 200, addModels=FALSE)
  
  return(cisTopicObject)
}

# Wrapper for SnapATAC workflow up until PCA step.
# Args:
#   snap_file (string): path to snap file
#   promooter.df (dataframe): dataframe with promoter definitions as shown in SnapATAC tutorials (dataframe not file name; see examples in markdown)
#   blacklist.df (dataframe): dataframe with blacklist region definitions as shown in SnapATAC tutorials (dataframe not file name; see examples in markdown)
#   fragment_number_threshold (int): threshold for number of unique fragments per cell
#   promoter_ratio_range (c(float, float)): vector with lower and upper bound of acceptable fraction of reads in promoter regions for cell filtering as used in SnapATAC tutorial.
#   window_z_range (c(float, float)): vector with lower and upper bound of acceptable window z-scores for non-zero entries for site filtering as used in SnapATAC tutorial.
#   sample_name (string): sample_name provided to SnapATAC
#   pc.num (int): total PCs to compute
# Returns:
#   snap object: SnapATAC object
# Notes:
#   This uses single core because multithreaded implementation interferes with Knitr. In running this code, any do.par=FALSE and num.cores=1 could be changed as needed.
snapatac_workflow = function(snap_file, promoter.df=NULL, blacklist.df=NULL, fragment_number_threshold=500, promoter_ratio_range=c(0.2, 0.8), window_z_range=c(-1.5, 1.5), sample_name='default', pc.num=50) {
  x.sp = createSnap(
    file=snap_file,
    sample=sample_name,
    do.par=FALSE,
    num.cores=1
  )
  
  plotBarcode(
    obj=x.sp, 
    pdf.file.name=NULL, 
    pdf.width=7, 
    pdf.height=7, 
    col="grey",
    border="grey",
    breaks=50
  )
  
  x.sp = filterCells(
    obj=x.sp, 
    subset.names=c("fragment.num", "UMI"),
    low.thresholds=c(fragment_number_threshold, fragment_number_threshold),
    high.thresholds=c(Inf, Inf)
  )
  
  x.sp = addBmatToSnap(x.sp, bin.size=5000, num.cores=1)
  
  # Optionally filter cells based on ratio of reads in promoters
  if (!is.null(promoter.df)) {
    promoter.gr = GRanges(promoter.df[,1], IRanges(promoter.df[,2], promoter.df[,3]))
    ov = findOverlaps(x.sp@feature, promoter.gr)
    idy = queryHits(ov)
    promoter_ratio = SnapATAC::rowSums(x.sp[,idy, mat="bmat"], mat="bmat") / SnapATAC::rowSums(x.sp, mat="bmat")
    plot(log(SnapATAC::rowSums(x.sp, mat="bmat") + 1,10), promoter_ratio, cex=0.5, col="grey", xlab="log(count)", ylab="FIP Ratio", ylim=c(0,1 ))
    idx = which(promoter_ratio > promoter_ratio_range[1] & promoter_ratio < promoter_ratio_range[2])
    x.sp = x.sp[idx,]
  }
  
  x.sp = makeBinary(x.sp, mat="bmat");
  
  # Filter out non-standard contigs if present
  idy2 = grep("chrM|random", x.sp@feature)
  
  if (!is.null(blacklist.df)) {
    black_list.gr = GRanges(
      blacklist.df[,1], 
      IRanges(blacklist.df[,2], blacklist.df[,3])
    )
    idy1 = queryHits(findOverlaps(x.sp@feature, black_list.gr))
    
  } else {
    # No blacklist provided, so just ignore
    idy1 = c()
  }
  
  idy = unique(c(idy1, idy2))
  x.sp = x.sp[,-idy, mat="bmat"]
  
  # Filter based on frequency
  x.sp = filterBins(
    x.sp,
    low.threshold=window_z_range[1],
    high.threshold=window_z_range[2],
    mat="bmat"
  )
  
  plotBinCoverage(
    x.sp,
    pdf.file.name=NULL,
    col="grey",
    border="grey",
    breaks=10,
    xlim=c(-6,6)
  )
  
  x.sp = runJaccard(
    obj = x.sp,
    tmp.folder=tempdir(),
    mat = "bmat",
    max.var=2000,
    ncell.chunk=1000,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )
  
  x.sp = runNormJaccard(
    obj = x.sp,
    tmp.folder=tempdir(),
    ncell.chunk=1000,
    method="normOVE",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )
  
  x.sp = runDimReduct(
    x.sp,
    pc.num=pc.num,
    input.mat="jmat",
    method="svd",
    center=TRUE,
    scale=FALSE,
    seed.use=10
  )
  
  rownames(x.sp@bmat) = x.sp@barcode
  colnames(x.sp@bmat) = as.character(1:ncol(x.sp@bmat))
  
  return(x.sp)
}

# Reperform Jaccard + PCA on SnapATAC object (used for redoing after modifying matrix)
# Args:
#   snao_obj (snap object): snap object
#   pc.num (int): total PCs to compute
# Returns:
#   snap object: SnapATAC object
# Notes:
#   This uses single core because multithreaded implementation interferes with Knitr. In running this code, any do.par=FALSE and num.cores=1 could be changed as needed.
snapatac_rerun_jaccard = function(snap_obj, pc.num=50) {
  
  snap_obj = runJaccard(
    obj = snap_obj,
    tmp.folder=tempdir(),
    mat = "bmat",
    max.var=2000,
    ncell.chunk=1000,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )
  
  snap_obj = runNormJaccard(
    obj = snap_obj,
    tmp.folder=tempdir(),
    ncell.chunk=1000,
    method="normOVE",
    row.center=TRUE,
    row.scale=TRUE,
    low.threshold=-5,
    high.threshold=5,
    do.par=FALSE,
    num.cores=1,
    seed.use=10
  )
  
  snap_obj = runDimReduct(
    snap_obj,
    pc.num=pc.num,
    input.mat="jmat",
    method="svd",
    center=TRUE,
    scale=FALSE,
    seed.use=10
  )
}

# Wrapper for cisTopic workflow (choose from models that have been run + clustering + further dim reduction)
# Args:
#   cistopicObject (cistopic object): cisTopic object with runModels having already been run
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   method (string): method argument to modelMatSelection function in cisTopic
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#   resolution (float): resolution parameter to Seurat Louvain clustering
# Returns:
#   Seurat object: Seurat object. clustering + tSNE + UMAP done on topic matrix from cisTopic
seurat_workflow_on_cistopic = function(cistopicObject, method='Z-score', reduction='pca', resolution=0.3) {
  
  # cistopicObject = mouse.cistopic
  cistopicObject = cisTopic::selectModel(cistopicObject)
  
  cistopicObject.reduced_space = t(cisTopic::modelMatSelection(cistopicObject, target='cell', method=method))
  colnames(cistopicObject.reduced_space) = paste0('PC_', 1:ncol(cistopicObject.reduced_space))
  dimensions = ncol(cistopicObject.reduced_space)
  
  cistopicObject.seurat = run_dim_reduction(cistopicObject@binary.count.matrix, cistopicObject.reduced_space, dims=1:dimensions, reduction='pca')
  
  cistopicObject.seurat = cistopicObject.seurat %>% 
    Seurat::FindNeighbors(reduction=reduction, nn.eps=0.25, dims=1:dimensions) %>%
    Seurat::FindClusters(reduction=reduction, n.start=20, resolution=resolution)
  
  return(cistopicObject.seurat)
}

# Wrapper for running downstream Seurat workflow (clustering + further dim reduction) on PCA from Jaccard matrix generated by SnapATAC
# Args:
#   snao_obj (snap object): snap object with runDimReduct already run
#   dims (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   reduction (string): reduction to use for downstream steps. Can be 'pca' (cell_embeddings) or 'pca.l2' (L2 normalized cell_embeddings)
#   resolution (float): resolution parameter to Seurat Louvain clustering
# Returns:
#   Seurat object: Seurat object. clustering + tSNE + UMAP done on PCA matrix from SnapATAC (note PCA is weighted by variance explained)
seurat_workflow_on_jaccard_pca = function(snap_obj, dims, reduction='pca', resolution=0.3) {
  pca_results.jaccard = snap_obj@smat@dmat %*% diag(snap_obj@smat@sdev)
  colnames(pca_results.jaccard) = paste0('PC_', 1:ncol(pca_results.jaccard))
  rownames(pca_results.jaccard) = snap_obj@barcode
  
  seurat_obj.jaccard = run_dim_reduction(t(snap_obj@bmat), pca_results.jaccard, dims, reduction=reduction)
  seurat_obj.jaccard = seurat_obj.jaccard %>%
    Seurat::FindNeighbors(nn.eps=0.25, dims=dims, reduction=reduction) %>%
    Seurat::FindClusters(n.start=20, resolution=resolution, dims=dims, reduction=reduction)
}

# Helper function to plot embeddings corresponding to same set of cells with one another's cluster assignments for comparison.
# Args:
#   seurat_obj1 (snap object): snap object with runDimReduct already run
#   seurat_obj2 (vector of int): vector of dims to use from cell_embeddings in downstream analysis
#   reduction (string): reduction to use for plot (can be tsne or umap)
#   description1 (string): title for seurat_obj1 (used in plot)
#   description2 (string): title for seurat_obj1 (used in plot)
#   cluster_column1 (string): column from metadata of seurat_obj1 to use for coloring plot
#   cluster_column2 (string): column from metadata of seurat_obj2 to use for coloring plot
# Returns:
#   ggolot object: ggplot object where each embedding is plotted colored by its own clusters and then again with the opposite object's clusters assigned for comparison. Four total panels shown.
plot_clustering_comparison = function(seurat_obj1, seurat_obj2, reduction, description1='', description2 = '', cluster_column1='RNA_snn_res.0.3', cluster_column2='RNA_snn_res.0.3') {
  # Clusters as called on each dataset
  seurat_obj1 = SetIdent(seurat_obj1, value=cluster_column1)
  seurat_obj2 = SetIdent(seurat_obj2, value=cluster_column2)
  
  p1 = DimPlot(seurat_obj1, reduction = reduction, pt.size=0.25) +
    ggtitle(description1)
  
  p2 = DimPlot(seurat_obj2, reduction = reduction, pt.size=0.25) +
    ggtitle(description2)
  
  # Now swap the labels to verify they are finding the same groups
  seurat_obj1@meta.data$cluster_seurat_obj2 = seurat_obj2@meta.data[, cluster_column2]
  seurat_obj2@meta.data$cluster_seurat_obj1 = seurat_obj1@meta.data[, cluster_column1]
  
  seurat_obj1 = SetIdent(seurat_obj1, value='cluster_seurat_obj2')
  seurat_obj2 = SetIdent(seurat_obj2, value='cluster_seurat_obj1')
  
  p3 = DimPlot(seurat_obj1, reduction = reduction, pt.size=0.25) +
    ggtitle(paste0(description1, ' (', description2, ' clusters)'))
  
  p4 = DimPlot(seurat_obj2, reduction = reduction, pt.size=0.25) +
    ggtitle(paste0(description2, ' (', description1, ' clusters)'))
  
  (p1 + p3) / (p2 + p4)
}

########################################################
########################################################
# Section : my own functions for cluster annotations
# 1) compute gene activity scores
# 2) motif enrichment analysis
# 3) integrate scRNA-seq data 
########################################################
########################################################
compare.methods.for.gene.activity.matrix = function()
{
  seurat.cistopic = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object.rds'))
  
  DefaultAssay(seurat.cistopic) <- 'peaks'
  Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
  new.cluster.ids <- c("P0", "1.ABaxx", "2.MS", "3.AB", "4.ABa/p/MS/E", "5.E",
                       '6.E', '7.AB/EMS','8.ABp', '9.ABaxx', '10.MS', 
                       '11.ABp', '12.AB', '13.MS', '14.MS', '15.ABp',
                       '16.MS', '17.ABp', '18.AB','19.C/Ea', '20.ABa/p/MS/E', '21.P4')
  names(new.cluster.ids) <- levels(seurat.cistopic)
  seurat.cistopic <- RenameIdents(seurat.cistopic, new.cluster.ids)
  
  DimPlot(seurat.cistopic, reduction = "umap", label = TRUE, pt.size = 2,  label.size = 10, repel = TRUE) + NoLegend()
  
  ##########################################
  # computing gene body and promoter activities using seurat
  ##########################################
  source.my.script('scATAC_functions.R')
  fragment.file = '../output_cellranger.ce11_scATACpro/cellranger_atac_ce11/fragments.tsv.gz'
  
  for (promoter.size in c(2000))
  {
    for (regions in c('promoter'))
    {
      #promoter.size = 5000
      #regions = 'promoter'
      compute.gene.acitivity.scores(seurat.cistopic, 
                                    method = 'seurat',
                                    fragment.file = fragment.file,
                                    regions = regions,
                                    promoter.upstream = promoter.size,
                                    promoter.downstream = 500,
                                    saveDir = RdataDir)
    }
  }
  
  ##########################################
  # compare promoters of different sizes with peaks 
  ##########################################
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneBody.promoter.activityscores.rds')
  
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.2000bpDownstream.rds')
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_promoterOnly.activityscores.rds')
  
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.200bpDownstream.rds')
  
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.2000bpDownstream.rds')
  
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.500bpDownstream.rds')
  
  seurat.cistopic = readRDS(file =  gamat)
  #seurat.cistopic = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object_geneBody.promoter.activityscores.rds'))
  
  DefaultAssay(seurat.cistopic) <- 'RNA'
  
  hist(log2(as.numeric(GetAssayData(seurat.cistopic, slot = 'counts'))+1), breaks = 10000)
  
  pdfname = paste0(resDir, "/UMAP_geneActivityMatrix_promoter.genebody_test_nbHVGs.pdf")
  pdf(pdfname, width=16, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(nb.variableFeatures in c(seq(1000, 5000, by = 1000), 8000, 10000, 12000, 15000, nrow(seurat.cistopic)))
  {
    # nb.variableFeatures = 20000
    cat('nb of HVGs -- ', nb.variableFeatures, '\n')
    seurat.cistopic <- FindVariableFeatures(seurat.cistopic, nfeatures = nb.variableFeatures)
    seurat.cistopic <- NormalizeData(seurat.cistopic, 
                                     #normalization.method = 'CLR',
                                     scale.factor = median(Matrix::colSums(seurat.cistopic@assays$RNA@counts)))
    seurat.cistopic <- ScaleData(seurat.cistopic)
    
    labels = Idents(seurat.cistopic)
    Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
    
    seurat.cistopic <- RunPCA(seurat.cistopic, npcs = 50, verbose = FALSE, reduction.name = 'pca.ga')
    
    nb.pcs = 50; n.neighbors = 20; min.dist = 0.1;
    seurat.cistopic <- RunUMAP(seurat.cistopic, reduction = "pca.ga", dims = 1:nb.pcs,
                               n.neighbors = n.neighbors, min.dist = min.dist, reduction.name = "umap.ga")
    p1 = DimPlot(seurat.cistopic, reduction = "umap.ga", label = TRUE, pt.size = 2, label.size = 5, repel = FALSE) + 
      ggtitle(paste0('nb of HVGs - ', nb.variableFeatures))
    plot(p1)
  }
  
  dev.off()
  
  ##########################################
  # compute gene activity matrix using cicero
  ##########################################
  #source.my.script('scATAC_functions.R')
  #compute.cicero.gene.activity.scores(seuratObj, output_dir, tss_file, genome_size_file)
  
  
}

compute.gene.acitivity.scores = function(seurat.cistopic, 
                                         method = 'seurat',
                                         fragment.file = 'fragments.tsv.gz',
                                         regions = 'promoter.geneBody',
                                         promoter.upstream = 2000,
                                         promoter.downstream = 2000, 
                                         saveDir = "results/scATAC_earlyEmbryo_20200302/Rdata/")
{
  if(method == 'cicero'){
    seurat.cistopic  = compute.cicero.gene.activity.scores(seuratObj = seurat.cistopic)
  }
  
  if(method == 'seurat')
  {
    library(Signac)
    library(Seurat)
    library(GenomeInfoDb)
    library(TxDb.Celegans.UCSC.ce11.ensGene)
    library(TxDb.Celegans.UCSC.ce11.refGene)
    library(ggplot2)
    set.seed(1234)
    
    cat('method : ', method, '- regions: ', regions, '- promoter.upstream: ', promoter.upstream, '- promoter.downstream: ', 
        promoter.downstream, '\n')
    #fragment.file = '/Volumes/groups/cochella/jiwang/Projects/Aleks/R8898_scATAC/cellranger_atac_wbcel235/outs/fragments.tsv.gz'
    #fragment.file = '../output_cellranger.ce11_scATACpro/cellranger_atac_ce11/fragments.tsv.gz'
    
    annot = read.csv(file = "data/annotations/BioMart_WBcel235_noFilters.csv", header = TRUE)
    annot = annot[which((annot$Gene.type == 'miRNA' | annot$Gene.type == 'protein_coding') 
                        & annot$Chromosome.scaffold.name != 'MtDNA'), ]
    
    #extract gene coordinates from Ensembl, and ensure name formatting is consistent with  Seurat object 
    if(regions == 'promoter.geneBody'){
      gene.coords <- GenomicFeatures::genes(TxDb.Celegans.UCSC.ce11.ensGene)
      seqlevelsStyle(gene.coords) <- 'UCSC'
      
      mm = match(gene.coords$gene_id, annot$Gene.stable.ID)
      gene.coords = gene.coords[which(!is.na(mm))]
      
      genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
      genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = promoter.upstream, downstream = 0)
      #features.coords = genebodyandpromoter.coords
      
      seqends <- seqlengths(genebodyandpromoter.coords)[as.character(seqnames(genebodyandpromoter.coords))]
      kk = which(start(genebodyandpromoter.coords) > 1L & end(genebodyandpromoter.coords) <= seqends)
      genebodyandpromoter.coords = genebodyandpromoter.coords[kk, ]
      
    }else{
      promoter.ens = GenomicRanges::promoters(TxDb.Celegans.UCSC.ce11.ensGene, 
                                              upstream = promoter.upstream, downstream = promoter.downstream)
      
      # only select promoter overlapped with refseq annotation
      select.only.promoter.overlapped.by.Refseq = FALSE
      if(select.only.promoter.overlapped.by.Refseq){
        promoter.refseq = GenomicFeatures::promoters(TxDb.Celegans.UCSC.ce11.refGene, 
                                                     upstream = promoter.upstream, downstream = promoter.downstream)
        ii = findOverlaps(promoter.ens, promoter.refseq,
                          maxgap=-1L, minoverlap=0L,
                          type=c("equal"),
                          select=c("all"),
                          ignore.strand=FALSE)
        gene.coords = promoter.ens[unique(ii@from), ] 
      }else{
        gene.coords = promoter.ens
      }
      
      # select promoters for protein coding genes and miRNAs
      mm = match(gene.coords$tx_name, annot$Transcript.stable.ID)
      gene.coords = gene.coords[which(!is.na(mm))]
      
      # random choose one promoter for 
      #mm = match(rownames(gene.activities), annot$Transcript.stable.ID)
      #mapNames = data.frame(c(1:nrow(gene.activities)), mm, geneNames = annot$Gene.name[mm], stringsAsFactors = FALSE)
      #genes.names = unique(mapNames$geneNames[!is.na(mapNames$geneNames)])
      #kk = match(genes.names, mapNames$geneNames)
      
      seqends <- seqlengths(gene.coords)[as.character(seqnames(gene.coords))]
      gene.coords = gene.coords[which(start(gene.coords) > 1L & end(gene.coords) <= seqends), ]
      #which(end(gene.coords) > seqlengths(gene.coords)[as.character(seqnames(gene.coords))])
      
      genebody.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
      genebodyandpromoter.coords <- Extend(x = gene.coords, upstream = 0, downstream = 0)
      #features.coords = genebodyandpromoter.coords
    }
    
    # build a gene by cell matrix
    tic('counting fragments for gene body and promoter')
    gene.activities <- FeatureMatrix(
      fragments = fragment.file,
      features = genebodyandpromoter.coords,
      cells = colnames(seurat.cistopic),
      chunk = 20
    )
    toc()
    
    if(regions == 'promoter.geneBody'){
      # convert rownames from chromsomal coordinates into gene names
      gene.key <- genebodyandpromoter.coords$gene_id
      names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
      rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])
      gene.activities <- gene.activities[rownames(gene.activities)!="",]
      
      mm = match(rownames(gene.activities), annot$Gene.stable.ID)
      gene.activities = gene.activities[which(!is.na(mm)), ]
      rownames(gene.activities) = annot$Gene.name[match(rownames(gene.activities), annot$Gene.stable.ID) ]
      
      #Add the gene activity matrix to the Seurat object as a new assay, and normalize it
      seurat.cistopic[['RNA']] <- CreateAssayObject(counts = gene.activities)
      saveRDS(seurat.cistopic, file =  paste0(saveDir, 'atac_LDA_seurat_object_geneActivityMatrix_', method, '_',
                                              regions,  '_', promoter.upstream, 'bpUpstream.rds'))
      
    }else{
      #saveRDS(gene.activities, file =  paste0(RdataDir, 'gene_activities_matrix_Promoter.rds'))
      
      # convert rownames from chromsomal coordinates into gene names
      gene.key <- genebodyandpromoter.coords$tx_name
      names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
      rownames(gene.activities) <- make.unique(gene.key[rownames(gene.activities)])
      gene.activities <- gene.activities[rownames(gene.activities)!="",]
      
      mm = match(rownames(gene.activities), annot$Transcript.stable.ID)
      mapNames = data.frame(c(1:nrow(gene.activities)), mm, geneNames = annot$Gene.name[mm], stringsAsFactors = FALSE)
      genes.names = unique(mapNames$geneNames[!is.na(mapNames$geneNames)])
      kk = match(genes.names, mapNames$geneNames)
      
      gene.activities = gene.activities[mapNames[kk, 1], ]
      rownames(gene.activities) = genes.names
      
      #Add the gene activity matrix to the Seurat object as a new assay, and normalize it
      seurat.cistopic[['RNA']] <- CreateAssayObject(counts = gene.activities)
      saveRDS(seurat.cistopic, file =  paste0(saveDir, 'atac_LDA_seurat_object_geneActivityMatrix_', method, '_',
                                              regions,  '_', promoter.upstream, 'bpUpstream.', promoter.downstream,'bpDownstream.rds'))
    }
    
  }
  
}

##########################################
# do cicero given a Seurat object, output gene activity score
##########################################
compute.cicero.gene.activity.scores = function(seuratObj, output_dir, tss_file, genome_size_file)
{
  library(cicero)
  #args = commandArgs(T)
  #seuratObj_file = args[1]
  #output_dir = args[2]
  #tss_file = args[3]
  #genome_size_file = args[4]
  source.my.script('scATAC_functions.R')
  
  seurat.obj = seurat.cistopic
  output_dir = RdataDir
  reduction = 'umap'
  npc = 30
  
  tss_file = 'data/annotations/ce11_tss.bed'
  genome_size_file = 'data/annotations/ce11.chrom.sizes'
  chr_sizes = genome_size_file
    
  # process tss annotation
  tss_ann <- fread(tss_file, header = F)
  names(tss_ann)[c(1:4)] <- c('chr', 'start', 'end', 'gene_name')
  annot = read.csv(file = "data/annotations/BioMart_WBcel235_noFilters.csv", header = TRUE)
  annot = annot[which((annot$Gene.type == 'miRNA' | annot$Gene.type == 'protein_coding') 
                      & annot$Chromosome.scaffold.name != 'MtDNA'), ]
  mm = match(tss_ann$gene_name, annot$Gene.stable.ID)
  tss_ann$gene_name = annot$Gene.name[mm]
  tss_ann = tss_ann[which(!is.na(tss_ann$gene_name)), ]
  #tss_ann <- tss_ann[gene_type %in% c('miRNA', 'lincRNA', 'protein_coding'), ]
  
  seurat.obj$active_clusters = as.character(seurat.obj$seurat_clusters)
  
  # main function calling Cicero
  res = doCicero_gascore(seurat.obj, reduction = 'umap', genome_size_file, tss_ann, npc = 30)
  
  ga_score = log1p(res$ga_score * 10000)
  
  conns = res$conns
  
  saveRDS(ga_score, file = paste0(output_dir, '/cicero_gene_activity.rds'))
  
  #conns$Peak1 = assignGene2Peak_coords(conns$Peak1, tss_ann)
  #conns$Peak2 = assignGene2Peak_coords(conns$Peak2, tss_ann)
  write.table(conns, file = paste0(output_dir, '/cicero_interactions.txt'), row.names = F,
              sep = '\t', quote = F)
  
}

doCicero_gascore <- function(seurat.obj, reduction = 'tsne', chr_sizes,
                             gene_ann, npc = 30, coaccess_thr = 0.25){
  ## gene_ann: the first four columns: chr, start, end, gene name
  set.seed(2019)
  
  mtx = GetAssayData(seurat.obj, slot = 'counts')
  # change rownames using _ to delimited
  rnames = rownames(mtx)
  new.rnames = sapply(rnames, function(x) unlist(strsplit(x, ','))[1])
  new.rnames = sapply(new.rnames, function(x) gsub('-', '_', x))
  rownames(mtx) <- new.rnames
  
  dt = reshape2::melt(as.matrix(mtx), value.name = 'count')
  rm(mtx)
  dt = dt[dt$count > 0, ]
  names(dt) = c('Peak', 'Cell', 'Count')
  dt$Cell = as.character(dt$Cell)
  dt$Peak = as.character(dt$Peak)
  
  # The output of make_atac_cds is a valid CDS object ready to be input into downstream Cicero functions.
  input_cds <- make_atac_cds(dt, binarize = TRUE)
  rm(dt)
  
  input_cds <- detectGenes(input_cds)
  input_cds <- estimateSizeFactors(input_cds)
  
  # Next, we access the tsne or umap coordinates from the input CDS object where they are stored by Monocle and run make_cicero_cds:
  if(reduction == 'tsne') {
    if(is.null(seurat.obj@reductions$tsne))
      seurat.obj <- RunTSNE(seurat.obj, dims = 1:npc, check_duplicates = F)
    redu.coords = seurat.obj@reductions$tsne@cell.embeddings
  }
  if(reduction == 'umap') {
    if(is.null(seurat.obj@reductions$umap))
      seurat.obj <- RunUMAP(seurat.obj, dims = 1:npc)
    redu.coords = seurat.obj@reductions$umap@cell.embeddings
  }
  
  
  #make the cell id consistet
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = redu.coords, k = 100)
  
  ## get connections
  #conns <- run_cicero(cicero_cds, chr_sizes, sample_num = 2)
  distance_parameters = estimate_distance_parameter(cicero_cds, genomic_coords = chr_sizes, 
                                                    window = 100000, s = 0.8, distance_constraint = 50000, sample_num = 100)
  
  ## get cicero gene activity score
  names(gene_ann)[4] <- "gene"
  
  input_cds <- annotate_cds_by_site(input_cds, gene_ann)
  
  # generate unnormalized gene activity matrix
  unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
  
  # make a list of num_genes_expressed
  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))
  
  # normalize
  cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
  
  # if you had two datasets to normalize, you would pass both:
  # num_genes should then include all cells from both sets
  #unnorm_ga2 <- unnorm_ga
  #cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga, unnorm_ga2), num_genes)
  conns = data.table(conns)
  conns = conns[coaccess > coaccess_thr, ]
  res = list('conns' = conns, 'ga_score' = cicero_gene_activities)
  return(res)  
}

##########################################
# motif analysis 
##########################################
compute.motif.enrichment = function(seurat.cistopic)
{
  library(Signac)
  library(Seurat)
  library(JASPAR2018)
  library(TFBSTools)
  library(BSgenome.Celegans.UCSC.ce11)
  set.seed(1234)
  
  ##########################################
  # add motif class in assay 'peaks' for seurat object
  ##########################################
  Convert.MEME.to.JASPAR = FALSE
  if(Convert.MEME.to.JASPAR){
    library(universalmotif) 
    library(MotifDb)
    pwm.meme = read_meme(file = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015.meme')
    write_jaspar(pwm.meme, file = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015.pfm')
  }
  
  pfm = readJASPARMatrix(fn = 'data/motifs/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015.pfm', 
                         matrixClass = 'PFM')
  
  # Scan the DNA sequence of each peak for the presence of each motif
  tic('scan the motif for peak regions')
  motif.matrix <- CreateMotifMatrix(
    features = StringToGRanges(rownames(seurat.cistopic), sep = c(":", "-")),
    pwm = pfm,
    genome = BSgenome.Celegans.UCSC.ce11,
    sep = c(":", "-"),
    use.counts = TRUE
  )
  #saveRDS(motif.matrix, file =  paste0(RdataDir, 'motifOccurence_inPeaks.rds'))
  toc()
  
  # Create a new Mofif object to store the results
  motif <- CreateMotifObject(
    data = motif.matrix,
    pwm = pfm
  )
  
  # Add the Motif object to the assay
  seurat.cistopic[['peaks']] <- AddMotifObject(
    object = seurat.cistopic[['peaks']],
    motif.object = motif
  )
  
  seurat.cistopic <- RegionStats(
    object = seurat.cistopic,
    genome = BSgenome.Celegans.UCSC.ce11,
    sep = c(":", "-")
  )
  
  ##########################################
  # motif analysis 1) :Finding overrepresented motifs for differentially accessible peaks
  ##########################################
  da_peaks <- FindMarkers(
    object = seurat.cistopic,
    ident.1 = '10', 
    ident.2 = '0',
    min.pct = 0.1,
    test.use = 'LR',
    latent.vars = 'nFeature_peaks'
  )
  
  head(da_peaks)
  
  enriched.motifs <- FindMotifs(
    object = seurat.cistopic,
    #background = 2000,
    features = head(rownames(da_peaks), 2000)
  )
  
  head(enriched.motifs, 20)
  
  ##########################################
  # motif analysis 2): run ChromVAR
  ##########################################
  seurat.cistopic <- RunChromVAR(
    object = seurat.cistopic,
    genome = BSgenome.Celegans.UCSC.ce11
    
  )
  
  return(seurat.cistopic)
  
}

##########################################
# integrate scRNA-seq data from published dataset in https://github.com/qinzhu/VisCello.celegans
# Packer, J. S J. I. Murray (2019). A lineage-resolved molecular atlas of C. elegans embryogenesis at single-cell resolution. Science: eaax1971.
# here we are using Seurat label transferring method

##########################################
transfer.labels.from.scRNA.to.scATAC = function(seurat.cistopic, tintori, method = c('seurat', 'liger'))
{
  library("pheatmap")
  library("RColorBrewer")
  
  if(method == 'seurat'){
    
    # process.tintori.et.al()
    # processed raw data
    tintori = readRDS(file = paste0(RdataDir, 'Tintori.et.al_rawCounts_processed_sizefactorNormalization.rds')) 
    #xx <- readRDS(paste0(RdataDir, 'Tintori_et_al_highQualtiyCells.rds')) # processed data from cello 
    tintori = tintori[!is.na(match(rownames(tintori), rownames(seurat.cistopic)))] # select only protein-coding genes and miRNAs
    
    tintori <- FindVariableFeatures(
      object = tintori,
      nfeatures = 5000
    )
    
    tintori <- ScaleData(object = tintori)
    tintori <- RunPCA(object = tintori, npcs = 50, verbose = FALSE)
    
    Idents(tintori) = tintori$lineage
    tintori <- RunUMAP(object = tintori, dims = 1:30, n.neighbors = 10, min.dist = 0.2)
    DimPlot(tintori, reduction = "umap", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
    
    Merge.cell.types.tintori = FALSE
    if(Merge.cell.types.tintori){
      #Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
      new.cluster.ids <- c("P0/AB/P1", "P0/AB/P1", "P0/AB/P1","ABa/p/EMS", "ABa/p/EMS", "ABa/p/EMS",
                           "P2/P3/C","ABa/pr/l",  "ABa/pr/l",  "ABa/pr/l", "ABa/pr/l",
                           "P2/P3/C", "MS/E","P2/P3/C", "MS/E",
                           "ABa/pr/lx", "ABa/pr/lx", "ABa/pr/lx", "ABa/pr/lx",
                           "MSx1/2", "MSx1/2",  "Cx1/2",   
                           "Ea/p",    "Ea/p",    "Cx1/2",   "P4/D",    "P4/D" )   
      names(new.cluster.ids) <- levels(tintori)
      tintori <- RenameIdents(tintori, new.cluster.ids)
      
      DimPlot(tintori, reduction = "umap", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
    }
    
    gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.500bpDownstream.rds')
    #gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.2000bpDownstream.rds')
    seurat.cistopic = readRDS(file =  gamat)
    
    resolution = 0.8;
    
    pdfname = paste0(resDir, "/label_transferring_seurat_resolution_rawTintoriLineage_", resolution,  ".pdf")
    pdf(pdfname, width=12, height = 8)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    DefaultAssay(seurat.cistopic) <- 'peaks'
    seurat.cistopic = FindClusters(object = seurat.cistopic, 
                                   n.start=20, resolution=resolution,
                                   algorithm = 3)
    #seurat.cistopic <- RunPCA(seurat.cistopic, npcs = 50, verbose = FALSE, reduction.name = 'pca.ga',  features = tops$gene)
    
    #nb.pcs = 30; n.neighbors = 20; min.dist = 0.2;
    #seurat.cistopic <- RunUMAP(seurat.cistopic, reduction = "pca", dims = 1:nb.pcs,
    #                           n.neighbors = n.neighbors, min.dist = min.dist, reduction.name = "umap.ga")
    DimPlot(seurat.cistopic, reduction = "umap", label = TRUE, pt.size = 2,  label.size = 8, repel = FALSE) + NoLegend()
    
    # Here, we process the gene activity matrix 
    # in order to find anchors between cells in the scATAC-seq dataset 
    # and the scRNA-seq dataset.
    DefaultAssay(seurat.cistopic) <- 'RNA'
    nb.variableFeatures = 15000
    
    seurat.cistopic <- FindVariableFeatures(seurat.cistopic, nfeatures = nb.variableFeatures)
    seurat.cistopic <- NormalizeData(seurat.cistopic)
    seurat.cistopic <- ScaleData(seurat.cistopic)
    
    # features.to.use = unique(c('erm-1', # AB enriched
    #                     'cpg-2',
    #                     'med-2', # marker for EMS
    #                     'era-1', #'W02F12.3',
    #                     'tbx-35',
    #                     'end-3',
    #                     'mex-5',
    #                     'pigv-1', #'T09B4.1'
    #                     'ceh-51',
    #                     'end-1',
    #                     'pal-1',
    #                     'hlh-27',
    #                     'ref-1',
    #                     'tbx-38',
    #                     'sdz-38',
    #                     'lsy-6', # marker for ABaxx
    #                     'tbx-37','tbx-38', 'ceh-27',
    #                     'ceh-36', 'elt-1',
    #                     'dmd-4', 
    #                     'dpy-31', 'elt-7', 'end-1', 'nrh-79', 
    #                       "mab-5","med-1","med-2", 'pal-1', 'pha-4', 
    #                        'pes-1', 'nos-1', 'nos-2'
    # ))
    
    transfer.anchors <- FindTransferAnchors(
      reference = seurat.cistopic,
      query = tintori,
      features = unique(c(VariableFeatures(object = tintori), VariableFeatures(seurat.cistopic))),
      #features = features.to.use,
      reference.assay = 'RNA',
      query.assay = 'RNA',
      reduction = 'cca',
      k.anchor = 50, # k.anchor is neighborhood size for MNN big k.anchor, the bigger, the more anchors found
      k.filter = 150, # retain the anchor (cell from one dataset to annother) if within k.filter neighbors, the bigger, the more retained  
      max.features = 5000, # max nb of features used for anchor filtering
      k.score = 50, 
      npcs = 50, 
      dims = 1:50
    )
    
    cat('nb of cells in scATAC and in tintori as anchors : ', 
        length(unique(transfer.anchors@anchors[, 1])), '--',  length(unique(transfer.anchors@anchors[, 2])), '\n')
    
    predicted.labels <- TransferData(
      anchorset = transfer.anchors,
      #refdata = Idents(tintori),
      refdata = Idents(seurat.cistopic),
      #refdata = as.vector(Idents(seurat.cistopic)),
      #weight.reduction = seurat.cistopic[['pca']],
      weight.reduction = tintori[['pca']],
      dims = 1:50,
      k.weight = 50
    )
    
    tintori <- AddMetaData(object = tintori, metadata = predicted.labels)
    
    p1 = DimPlot(tintori, group.by = "predicted.id",reduction = 'umap', pt.size = 2,  label.size = 6,
                 label = TRUE, repel = FALSE) + ggtitle("scRNA-seq cells") + 
      scale_colour_hue(drop = FALSE)
    p2 = DimPlot(tintori,reduction = 'umap', label = TRUE, repel = FALSE, pt.size = 2,  label.size = 6) + 
      ggtitle("scRNA-seq cells") + 
      scale_colour_hue(drop = FALSE)
    
    p1 + p2
    
    res = data.frame(lineage = Idents(tintori), predicted.labels, stringsAsFactors = FALSE)
    
    map = matrix(0, nrow = length(levels(seurat.cistopic)),
                 ncol = length(levels(tintori)))
    rownames(map) = levels(seurat.cistopic)
    colnames(map) = levels(tintori)
    
    for(n in 1:ncol(map))
    {
      # n = 1
      kk = which(res$lineage == colnames(map)[n])
      #cat('nb of cells in the tintori et al.: ' ,length(kk), 'for lineage ', colnames(map)[n],  '\n')
      
      if(length(kk)>1){
        map[,n] = apply(res[kk, match(paste0('prediction.score.', rownames(map)), colnames(res))], 2, median)
        map[,n] = map[,n]/sum(map[,n])
      }else{
        cat('ERROR \n')
      }
    }
    
    cols = c('black',colorRampPalette((brewer.pal(n = 7, name="Reds")))(9))
    pheatmap(map, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, breaks = seq(0, 1, by = 0.1),
             cluster_cols=FALSE, main = paste0("resolution -- ",  resolution), na_col = "white",
             color = cols)
    
    map.fitered = map
    map.fitered[which(map.fitered<0.1)] = 0
    pheatmap(map.fitered, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
             cluster_cols=FALSE, main = paste0("fitlered < 0.1 with resolution -- ",  resolution), na_col = "white",
             color = cols)
    
    dev.off()
    
    
    Manual.modify.cluster.annotation = FALSE
    if(Manual.modify.cluster.annotation){
      DefaultAssay(seurat.cistopic) <- 'peaks'
      Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
      new.cluster.ids <- c("0.?", "1.ABaxx", "2.MS", "3.ABxx", "4.EMS", "5.D/P4",
                           '6.Ea/p', '7.P0/AB/P1/ABa/p/P2/P3','8.ABp', '9.ABaxx', '10.MS', 
                           '11.ABp', '12.AB', '13.MS', '14.MS', '15.ABp',
                           '16.MS', '17.ABp', '18.AB','19.C/Ea', '20.ABa/p/MS/E', '21.P4')
      names(new.cluster.ids) <- levels(seurat.cistopic)
      seurat.cistopic <- RenameIdents(seurat.cistopic, new.cluster.ids)
      
      DimPlot(seurat.cistopic, reduction = "umap", label = TRUE, pt.size = 2,  label.size = 5, repel = FALSE) + NoLegend()
      
      
    }
    
    # xx = as.matrix(map)
    # pheatmap(xx, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
    #          cluster_cols=FALSE, main = paste0("alpha = ",  " -- "), na_col = "white",
    #          color = cols)
    #table(seurat.cistopic$prediction.score.max > 0.5)
    
    # 
    # seurat.cistopic <- AddMetaData(object = seurat.cistopic, metadata = predicted.labels)
    # DimPlot(seurat.cistopic, group.by = "predicted.id",reduction = 'umap', label = TRUE, repel = FALSE) + ggtitle("scATAC-seq cells") + 
    #   scale_colour_hue(drop = FALSE)
    # 
    # table(seurat.cistopic$prediction.score.max > 0.5)
    
    hist(seurat.cistopic$prediction.score.max)
    abline(v = 0.5, col = "red")
    
    seurat.cistopic.filtered <- subset(seurat.cistopic, subset = prediction.score.max > 0.5)
    # to make the colors match
    seurat.cistopic.filtered$predicted.id <- factor(seurat.cistopic.filtered$predicted.id, levels = levels(tintori))  
    DimPlot(seurat.cistopic.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
      NoLegend() + scale_colour_hue(drop = FALSE)
    
      
  }
  
  if(method == 'liger'){
    transfer.labels.from.scRNA.to.scATAC.with.liger()
  }
    
}

process.tintori.et.al = function()
{
  # DefaultAssay(seurat.cistopic) <- 'peaks'
  # seurat.cistopic <- RunTFIDF(seurat.cistopic)
  # seurat.cistopic <- FindTopFeatures(seurat.cistopic, min.cutoff = 'q25')
  # seurat.cistopic <- RunSVD(
  #   object = seurat.cistopic,
  #   assay = 'peaks',
  #   reduction.key = 'LSI_',
  #   reduction.name = 'lsi'
  # )
  
  #VariableFeatures(seurat.cistopic) <- names(which(Matrix::rowSums(seurat.cistopic) > 100))
  #seurat.cistopic <- RunLSI(seurat.cistopic, n = 50, scale.max = NULL)
  #seurat.cistopic <- RunUMAP(seurat.cistopic, reduction = "lsi", dims = 1:50)
  
  # Here, we process the gene activity matrix 
  # in order to find anchors between cells in the scATAC-seq dataset 
  # and the scRNA-seq dataset.
  DefaultAssay(seurat.cistopic) <- 'RNA'
  nb.variableFeatures = nrow(seurat.cistopic)
  seurat.cistopic <- FindVariableFeatures(seurat.cistopic, nfeatures = nb.variableFeatures)
  seurat.cistopic <- NormalizeData(seurat.cistopic)
  seurat.cistopic <- ScaleData(seurat.cistopic)
  
  Refine.HGV.Gene.Activity = FALSE
  if(Refine.HGV.Gene.Activity){
    
    labels = Idents(seurat.cistopic)
    Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
    
    seurat.cistopic <- RunPCA(seurat.cistopic, npcs = 50, verbose = FALSE, reduction.name = 'pca.ga')
    
    nb.pcs = 30; n.neighbors = 20; min.dist = 0.3;
    seurat.cistopic <- RunUMAP(seurat.cistopic, reduction = "pca.ga", dims = 1:nb.pcs,
                               n.neighbors = n.neighbors, min.dist = min.dist, reduction.name = "umap.ga")
    DimPlot(seurat.cistopic, reduction = "umap.ga", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
    
    
    # find all markers of cluster 1
    markers.ga <- FindAllMarkers(seurat.cistopic, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    tops <- markers.ga %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
    #length(intersect(tops$gene, VariableFeatures(seurat.cistopic)))
    #VariableFeatures(seurat.cistopic) = tops$gene
    seurat.cistopic <- RunPCA(seurat.cistopic, npcs = 50, verbose = FALSE, reduction.name = 'pca.ga',  features = tops$gene)
    
    nb.pcs = 30; n.neighbors = 20; min.dist = 0.3;
    seurat.cistopic <- RunUMAP(seurat.cistopic, reduction = "pca.ga", dims = 1:nb.pcs,
                               n.neighbors = n.neighbors, min.dist = min.dist, reduction.name = "umap.ga")
    DimPlot(seurat.cistopic, reduction = "umap.ga", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
    
  }
  
  # Load the pre-processed scRNA-seq data 
  tintori <- readRDS(paste0(RdataDir, 'Tintori_et_al_highQualtiyCells.rds'))
  tintori = tintori[!is.na(match(rownames(tintori), rownames(seurat.cistopic)))] # select only protein-coding genes and miRNAs
  
  tintori <- FindVariableFeatures(
    object = tintori,
    nfeatures = 5000
  )
  
  DimPlot(tintori, reduction = "umap", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
  
  #Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
  new.cluster.ids <- c("P0/AB/P1", "P0/AB/P1", "P0/AB/P1","ABa/p/EMS", "ABa/p/EMS", "ABa/p/EMS",
                       "P2/P3/C","ABa/pr/l",  "ABa/pr/l",  "ABa/pr/l", "ABa/pr/l",
                       "P2/P3/C", "MS/E","P2/P3/C", "MS/E",
                       "ABa/pr/lx", "ABa/pr/lx", "ABa/pr/lx", "ABa/pr/lx",
                       "MSx1/2", "MSx1/2",  "Cx1/2",   
                       "Ea/p",    "Ea/p",    "Cx1/2",   "P4/D",    "P4/D" )   
  names(new.cluster.ids) <- levels(tintori)
  tintori <- RenameIdents(tintori, new.cluster.ids)
  
  DimPlot(tintori, reduction = "umap", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
  
}

##########################################
# test liger to transfer lables 
# see details in https://macoskolab.github.io/liger/walkthrough_rna_atac.html
##########################################
transfer.labels.from.scRNA.to.scATAC.with.liger = function(seurat.cistopic, tintori)
{
  
  # rna_clusts = readRDS("../liger-rna-atac-vignette/rna_cluster_assignments.RDS")
  # atac_clusts = readRDS("../liger-rna-atac-vignette/atac_cluster_assignments.RDS")
  # pbmc.atac <- readRDS('../liger-rna-atac-vignette/pbmc.atac.expression.mat.RDS')
  # pbmc.rna <- readRDS('../liger-rna-atac-vignette/pbmc.rna.expression.mat.RDS')
  
  # sels = sample(c(1:length(rna_clusts)), 200, replace = FALSE)
  # rna_clusts = rna_clusts[sels]
  #pbmc.rna = pbmc.rna[, sels]
  
  atac_clusts = seurat.cistopic$peaks_snn_res.0.8
  rna_clusts = Idents(tintori)
  pbmc.atac = GetAssayData(seurat.cistopic, slot = 'counts', assay = 'RNA')
  pbmc.rna = GetAssayData(tintori, slot = 'counts', assay = 'RNA')
  
  # data preprocessing
  library(liger)
  
  convertSeurt_toLiger = FALSE
  if(convertSeurt_toLiger){
    ligerex2 <- seuratToLiger(list(seurat.cistopic, tintori), assays.use = c('RNA', 'RNA'), num.hvg.info = 2500, renormalize = TRUE)
    
    cat('nb of HVGs : ', length(ligerex2@var.genes), '\n')
    int.pbmc = ligerex2
    
  }else{
    ggplot2::theme_set(theme_cowplot())
    #xx = pbmc.atac[,names(atac_clusts)]
    pbmc.data = list(atac=pbmc.atac[,names(atac_clusts)], rna=pbmc.rna[,names(rna_clusts)])
    int.pbmc <- createLiger(pbmc.data)
    
    int.pbmc <- liger::normalize(int.pbmc)
    
    par(mfrow=c(1,2))
    int.pbmc <- selectGenes(int.pbmc, datasets.use = c(1, 2), 
                            num.genes = c(20000, 5000), do.plot = TRUE)
    cat('nb of HVGs : ', length(unique(int.pbmc@var.genes)), '\n')
    
  }
  
  int.keep = int.pbmc
  
  pdfname = paste0(resDir, "/Test_liger_lambda.pdf")
  pdf(pdfname, width=16, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(lambda in c(1:4, seq(5, 100, by = 5)))
  {
    lambda = 50
    int.pbmc <- scaleNotCenter(int.keep)
    
    Suggest.k.lambda = FALSE
    if(Suggest.k.lambda){
      suggestK(int.pbmc) # plot entropy metric to find an elbow that can be used to select the number of factors
      suggestLambda(int.pbmc, k = 20, num.cores = 5) # plot alignment metric to find an elbow that can be used to select the value of lambda
    }
    
    # Factorization and Quantile Normalization
    int.pbmc <- optimizeALS(int.pbmc, k=30, lambda = lambda, nrep = 1)
    
    #int.pbmc <- liger::runUMAP(int.pbmc, use.raw = T)
    int.pbmc <- liger::quantile_norm(int.pbmc)
    
    # visualization
    int.pbmc <- liger::runUMAP(int.pbmc, use.raw = FALSE,  dims.use = 1:30, distance = 'euclidean', n_neighbors = 20, min_dist = 0.2)
    
    plots1 <- plotByDatasetAndCluster(int.pbmc, return.plots = T, clusters=rna_clusts) 
    #print(plots1[[1]])
    p1 =  print(plots1[[2]])
    plots2 <- plotByDatasetAndCluster(int.pbmc, return.plots = T, clusters=atac_clusts) 
    #print(plots2[[1]])
    p2 = print(plots2[[2]])
    
    plot(p1 + p2) + ggtitle(paste0('lambda -', lambda))
    
  }
  
  dev.off()
  
  calcAlignment(int.pbmc)
  
  xx = ligerToSeurat(int.pbmc, use.liger.genes = T)
  
  # NKG7 = plotGene(int.pbmc,gene="NKG7",return.plots=T)
  # MS4A1 = plotGene(int.pbmc,gene="MS4A1",return.plots=T)
  # plot_grid(NKG7[[2]],MS4A1[[2]],NKG7[[1]],MS4A1[[1]],ncol=2)
  
  #seurat_liger = ligerToSeurat(int.pbmc, use.liger.genes = T)
  
  
}


########################################################
########################################################
# Section : functions for scATAC-seq cluster annotation
# 
########################################################
########################################################
annotate.scATAC.clusters.with.marker.genes = function(seurat.cistopic)
{
  
}

manually.edit.table_marker.genes_tf.motif_peaks = function(seurat.cistopic)
{
  
  # processed from raw data
  tintori = readRDS(file = paste0(RdataDir, 'Tintori.et.al_rawCounts_processed_sizefactorNormalization.rds')) 
  tintori = tintori[!is.na(match(rownames(tintori), rownames(seurat.cistopic)))] # select only protein-coding genes and miRNAs
  
  tintori <- FindVariableFeatures(
    object = tintori,
    nfeatures = 5000
  )
  
  tintori <- ScaleData(object = tintori)
  tintori <- RunPCA(object = tintori, npcs = 50, verbose = FALSE)
  
  Idents(tintori) = tintori$lineage
  tintori <- RunUMAP(object = tintori, dims = 1:30, n.neighbors = 10, min.dist = 0.2)
  
  DimPlot(tintori, reduction = "umap", label = TRUE, pt.size = 2, label.size = 5, repel = FALSE)
  
  ## import peak, RNA (gene activity matrix) as Assays
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter.geneBody_2000bp.rds')
  #gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.500bpDownstream.rds')
  seurat.cistopic = readRDS(file = gamat)
  
  # add motif activities into seurat.cistopic 
  xx = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object_motifClassChromVar.rds'))
  
  seurat.cistopic[['chromvar']] = CreateAssayObject(data = xx@assays$chromvar@data)
  
  DefaultAssay(seurat.cistopic) <- 'peaks'
  
  FeaturePlot(
    object = seurat.cistopic,
    features = "chrV:10647260-10647485",
    pt.size = 0.1,
    #max.cutoff = ,
    #min.cutoff = 'q20',
    ncol = 1
  )
  
  
  #seurat.cistopic <- NormalizeData(seurat.cistopic, 
  #                                 #normalization.method = 'CLR'
  #                                 scale.factor = median(Matrix::colSums(seurat.cistopic@assays$RNA@counts))
  #                                 )
  #seurat.cistopic <- ScaleData(seurat.cistopic)
  
  features.to.plot = c('erm-1', # AB enriched
    'cpg-2',
    'med-2', # marker for EMS
    'era-1', #'W02F12.3',
    'tbx-35',
    'end-3',
    'mex-5',
    'pigv-1', #'T09B4.1'
    'ceh-51',
    'end-1',
    'pal-1',
    'hlh-27',
    'ref-1',
    'tbx-38',
    'sdz-38',
    'pha-4', 'hnd-1',
    'lsy-6' # marker for ABaxx
    #'tbx-37','tbx-38', 'ceh-27',
    #'ceh-36', 'elt-1'
    #'dmd-4', 
    #'dpy-31', 'elt-7', 'end-1', 'nrh-79', 
    #  "mab-5","med-1","med-2", 'pal-1', 'pha-4', 
    #   'pes-1', 'nos-1', 'nos-2'
  )
  
  FeaturePlot(
    object = tintori,
    features = features.to.plot,
    pt.size = 0.1,
    #max.cutoff = 2,
    #min.cutoff = 'q5',
    ncol = 1
  )
  
  features.to.plot = 'pal-1'
  
  DefaultAssay(seurat.cistopic) <- 'RNA'
  p1 = FeaturePlot(
    object = seurat.cistopic,
    features = features.to.plot,
    pt.size = 0.1,
    max.cutoff = 'q90',
    min.cutoff = 'q25',
    ncol = 1
  )
  
  DefaultAssay(seurat.cistopic) <- 'chromvar'
  motifs2check = rownames(seurat.cistopic)[grep('tbx-35|tbx-38|MA0542.1|MA0924.1 pal-1|MA0546.1|end-1|MA0547.1 skn-1|
                                                hnd|unc-130|mom-2|med-2|ref', 
                                                rownames(seurat.cistopic))]
  motifs2check = rownames(seurat.cistopic)[grep(features.to.plot, rownames(seurat.cistopic))]
  p2 = FeaturePlot(
    object = seurat.cistopic,
    features = motifs2check,
    pt.size = 0.1,
    #cols = c('red', 'blue'),
    #max.cutoff = 1,
    min.cutoff = 0.5,
    ncol = 1
  )
  
  p1 + p2
  
  
  ##########################################
  # test tf normalization and clustering 
  ##########################################
  # seurat.lsi = tenx.seurat.lsi
  # gene.activity = GetAssayData(object = seurat.cistopic, slot = "counts", assay = 'RNA')  
  # seurat.lsi[['RNA']] = CreateAssayObject(counts = gene.activity)
  # 
  # 
  # FeaturePlot(
  #   object = seurat.lsi,
  #   features = features.to.plot,
  #   pt.size = 0.1,
  #   max.cutoff = 1.,
  #   #min.cutoff = 'q5',
  #   ncol = 4
  # )
  
  
}

##########################################
# first test cisTopic and Seurat
##########################################
# library(data.table)
# library(Matrix)
# library(tictoc)
# DownSample.mtx = FALSE
# 
# #loadRDS(filter.out, file = paste0(filtered_mtx_dir, '/EmptyDrop_obj.rds'))
# raw_mtx_dir = '/Volumes/groups/cochella/jiwang/Projects/Aleks/R8898_scATAC/output/filtered_matrix'
# #metrics <- paste0(pathTo10X, 'atac_v1_pbmc_5k_singlecell.csv')
# 
# mat = readMM(paste0(raw_mtx_dir, "/matrix.mtx"))
# features = fread(paste0(raw_mtx_dir, '/features.txt'), header = F)
# barcodes = fread(paste0(raw_mtx_dir, '/barcodes.txt'), header = F)
# rownames(mat) = features$V1
# colnames(mat) = barcodes$V1
# 
# if(DownSample.mtx){
#   ## downsample
#   set.seed(1)
#   input = mat[, sample(c(1:ncol(mat)), 2000, replace = FALSE)]
# }else{
#   input = mat
# }
# 
# peaknames = rownames(input)
# peaknames = sapply(peaknames, function(x) {xx = unlist(strsplit(x, '-')); return(paste0(xx[1], ':', xx[2], "-", xx[3]))} )
# rownames(input) = peaknames
# 
# 
# library(cisTopic)
# cisTopicObject <- createcisTopicObject(count.matrix = input,  project.name='earlyEmbro')
# 
# tic('run model for cisTopic')
# #toc()
# nb.topcis = c(seq(10, 100, 10))
# cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=nb.topcis, seed=987, nCores=6, iterations = 200, addModels=FALSE)
# 
# save(cisTopicObject, file = paste0(RdataDir, 'cisTopicOject_runWarpLDAmodels_localrunning.Rdata'))
# 
# toc()
# #logLikelihoodByIter(cisTopicObject, select=nb.topcis)
# 
# ##########################################
# # model selection and interpretation
# ##########################################
# load(file = 'results/scATAC_earlyEmbryo_20200302/Rdata/cisTopicOject_runWarpLDAmodels.Rdata')
# 
# par(mfrow=c(3,3))
# cisTopicObject <- selectModel(cisTopicObject, type='maximum')
# cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
# cisTopicObject <- selectModel(cisTopicObject, type='derivative')
# 
# cisTopicObject <- runPCA(cisTopicObject, target='cell', seed=123, method='Probability')
# #cisTopicObject <- runtSNE(cisTopicObject, target='cell', seed=123, pca = TRUE, method='Probability')
# cisTopicObject <- runDM(cisTopicObject, target='cell', seed=123, pca=TRUE, method='Probability')
# 
# nb.pcs = 20; n.neighbors = 20; min.dist = 0.2;
# #ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
# cisTopicObject <- runUmap(cisTopicObject, target='cell', seed=123, pca = TRUE, method='Probability', n.neighbors = n.neighbors, min.dist = min.dist)
# 
# 
# cellassign <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
# dim(cellassign)
# cellassign[1:5,1:5]
# 
# #sum(colnames(cellassign) == rownames(metadata))
# ### make quick clustering using tSNE and densityClust
# set.seed(123)
# library(Rtsne)
# DR <- Rtsne(t(cellassign), pca=F)
# DRdist <- dist(DR$Y)
# library(densityClust)
# dclust <- densityClust(DRdist,gaussian=T)
# 
# cutoff.rho = 100
# cutoff.delta = 5
# dclust <- findClusters(dclust, rho = cutoff.rho, delta = cutoff.delta)
# par(mfrow=c(1,1))
# options(repr.plot.width=6, repr.plot.height=6)
# plot(dclust$rho,dclust$delta, pch=20,cex=0.6,xlab='rho', ylab='delta')
# points(dclust$rho[dclust$peaks],dclust$delta[dclust$peaks],col="red",pch=20,cex=0.8)
# text(dclust$rho[dclust$peaks]-0.2,dclust$delta[dclust$peaks]+0.2,labels=dclust$clusters[dclust$peaks])
# abline(v=cutoff.rho)
# abline(h=cutoff.delta)
# 
# densityClust <- dclust$clusters
# densityClust <- as.data.frame(densityClust)
# rownames(densityClust) <- cisTopicObject@cell.names
# colnames(densityClust) <- 'densityClust'
# densityClust[,1] <- as.factor(densityClust[,1])
# cisTopicObject <- addCellMetadata(cisTopicObject, densityClust)
# 
# 
# par(mfrow=c(1,1))
# plotFeatures(cisTopicObject, method='Umap', target='cell', topic_contr=NULL, 
#              colorBy=c('densityClust'), 
#              cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, 
#              col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
# 
# plotFeatures(cisTopicObject, method='tSNE', target='cell', topic_contr=NULL, 
#              colorBy=c('densityClust'), 
#              cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, 
#              col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)
# 
# plotFeatures(cisTopicObject, method='DM', target='cell', topic_contr=NULL, 
#              colorBy=c('densityClust'), 
#              cex.legend = 0.8, factor.max=.75, dim=2, legend=TRUE, 
#              col.low='darkgreen', col.mid='yellow', col.high='brown1', intervals=10)


# library(Signac)
# library(Seurat)
# library(ggplot2)
# set.seed(1234)
# 
# pbmc <- CreateSeuratObject(
#   counts = mat,
#   assay = 'peaks',
#   project = 'scATAC_earlyEmbro',
#   min.cells = 10000,
#   meta.data = NULL
# )
# 
# pbmc <- RunTFIDF(pbmc)
# pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q75')
# pbmc <- RunSVD(object = pbmc, n = 100,
#                assay = 'peaks',
#                reduction.key = 'LSI_',
#                reduction.name = 'lsi',
#                irlba.work=500
# )
# 
# pbmc = ScaleData(pbmc, features = VariableFeatures(pbmc))
# #ms <- RunPCA(object = ms, features = VariableFeatures(ms), verbose = FALSE)
# pbmc = RunPCA(pbmc, npcs = 100, assay = 'peaks', features = VariableFeatures(pbmc), verbose = FALSE)
# 
# pbmc <- FindNeighbors(object = pbmc, reduction = 'pca', dims = 2:nb.pcs)
# pbmc <- FindClusters(object = pbmc, verbose = FALSE, algorithm = 3, resolution = 1.0)
# 
# nb.pcs = 60; n.neighbors = 20; min.dist = 0.1;
# pbmc <- RunUMAP(object = pbmc, reduction = 'pca', dims = 2:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
# DimPlot(object = pbmc, label = TRUE, reduction = 'umap') + NoLegend()
# 
# #pbmc = RunTSNE(pbmc, reduction = 'lsi', dims = 2:nb.pcs)
# #DimPlot(object = pbmc, label = TRUE, reduction = 'tsne') + NoLegend()
# save(pbmc, file = paste0(RdataDir, 'pbmc_Seurat_TFIDtransform_pca_umap.Rdata'))
# 

