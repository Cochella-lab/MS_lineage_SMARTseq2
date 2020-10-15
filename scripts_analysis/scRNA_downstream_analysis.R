########################################################
########################################################
# Section : clustering and DE analysis or gene markers discovery
# feature selection and dimension reduction was done already in the batch correction part
# So we start with the PCs from fastMNN in scran, with which we define distance for clustering
# 1) different clustering methods will be tested  
# 2) special design matrix will be used for DE analysis 
# http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/
# doc/de.html#2_blocking_on_uninteresting_factors_of_variation
########################################################
########################################################
source.my.script <- function(name.of.function){
  tryCatch(path <- rstudioapi::getSourceEditorContext()$path,
           error = function(e){
             install.packages("rstudioapi")
             path <-  rstudioapi::getSourceEditorContext()$path})
  source.path <- sub(basename(path), "", path)
  source(paste0(source.path,name.of.function))
}

## set up the paths for the data and results
tryCatch(path <- rstudioapi::getSourceEditorContext()$path, 
         error = function(e){
           install.packages("rstudioapi")
           path <-  rstudioapi::getSourceEditorContext()$path})
source.path <- sub(basename(path), "", path)

user <- "results_jiwang/"
setwd(paste0("../", user))
#setwd(paste0("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/",user))

version.DATA = 'all_batches'
version.analysis =  paste0(version.DATA, '_202008')
dataDir = paste0("../data/gene_counts/")
resDir = paste0("results/", version.analysis)
tabDir = paste0("results/", version.analysis, "/tables/")

RdataDir = paste0("results/", version.analysis, "/Rdata/")
if(!dir.exists("results/")){dir.create("results/")}
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

library(Seurat)
library(ggplot2)
library(dplyr)
#source.my.script('scRNA_cluster_annotation_functions.R')

########################################################
########################################################
# Section : annotate scRNA clusters by mapping to reference
# 
########################################################
########################################################
# import processed data by Aleks
dir.processed.data = "/Users/jiwang/workspace/imp/scRNAseq_MS_lineage_dev/results_aleks/results/all_batches_202008_6.5k_cells/Rdata"
load(paste0(dir.processed.data, "/6.5k_cells_QCed_cells_genes_filtered_timingEst_Normed_bc_Seurat.Rdata"))

ms <- FindNeighbors(object = ms, reduction = "mnn", k.param = 20, dims = 1:20)
ms <- FindClusters(ms, resolution = 12, algorithm = 3)

DimPlot(object = ms, cells = colnames(ms), group.by = 'ident', label = TRUE, pt.size = 1)

#all_ms.markers <- FindAllMarkers(ms, min.pct = 0.25, logfc.threshold = 0.25)
#all_ms.top10 <- all_ms.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

#sce = ms
#rm(ms)
#DefaultAssay(ms) = 'RNA'
##########################################
# redo the normalization
##########################################
library(scater)
library(scran)
options(stringsAsFactors = FALSE)

sce = as.SingleCellExperiment(ms, assay = 'RNA')

plotColData(sce,
            x = "log10_total_counts",
            y = "log10_total_features_by_counts",
            #colour_by = "percent_mapped",
            colour_by = "request",
            size_by = "pct_counts_Mt"
) + scale_x_continuous(limits=c(4, 7)) +
  scale_y_continuous(limits = c(2.5, 4.1)) +
  geom_hline(yintercept=log10(c(500, 1000, 5000)) , linetype="dashed", color = "darkgray", size=0.5) +
  geom_vline(xintercept = c(4:6), linetype="dotted", color = "black", size=0.5)

##########################################
# here estimat the timing with timer genes
### first test 5 lineages from Hashimsholy et al. paper
##########################################
reEstimate.timing.using.timer.genes.using.cpmNorm = FALSE
if(reEstimate.timing.using.timer.genes.using.cpmNorm){
  # test the timingEst main function using the 5 lineages from Hashimshony et al. 
  #Test.timingEstimate()
  
  ## Here we are sampling a range of parameters and timing estimation were done with each of them
  ## Whereby we assess the sensibility of our timingEst
  ## this will take some time to finish
  source.my.script('timingEst_functions.R')
  sce = estimate.timing.and.variance.with.timer.genes(sce, lineageCorrs = NA)
  
  save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_timingEst_SCE.Rdata'))
  
}


par(mfrow = c(1, 3))
plot(sce$FSC_log2, as.numeric(as.character(sce$timingEst)), type='p', cex = 0.5)
plot(sce$BSC_log2, as.numeric(as.character(sce$timingEst)), type='p', cex = 0.5)
plot(sce$FSC_log2, sce$BSC_log2, type = 'p', cex = 0.5)

plotColData(sce,
            x = "FSC_log2",
            y = "BSC_log2",
            colour_by = "timingEst",
            point_size = 1
)

########################################################
########################################################
# Section : scRNA-seq data processing steps
# 1) normalization 
# 2) cell cycle correction (optional)
# 3) batch correction
########################################################
########################################################
#load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_timingEst_SCE.Rdata'))
library(scater)
library(SingleCellExperiment)
library(scran)
reducedDim(sce) <- NULL
#endog_genes <- !rowData(sce)$is_feature_control

Normalization.Testing = FALSE
source.my.script("normalization_HVGs_cellCycle_batchCorrection_functions.R")

if(Normalization.Testing){
  pdfname = paste0(resDir, "/scRNAseq_filtered_normalization_testing.pdf")
  pdf(pdfname, width=14, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  compare.scran.seurat.sctransform(sce, using.HVGs = TRUE)
  
  dev.off()
}

sce$library.size = apply(counts(sce), 2, sum)
qclust <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = qclust)
sce <- logNormCounts(sce, log = TRUE, pseudo_count = 1)

par(mfrow = c(1, 1))
plot(sce$library.size/1e6, sizeFactors(sce), log="xy", xlab="Library size (millions)", ylab="Size factor", cex = 0.5)


ms = as.Seurat(sce, counts = 'counts', data = 'logcounts', assay = "RNA") # scran normalized data were kept in Seurat

##########################################
# HVGs were performed 
##########################################
nfeatures = 3000
ms <- FindVariableFeatures(ms, selection.method = "vst", nfeatures = nfeatures)

top10 <- head(VariableFeatures(ms), 10) # Identify the 10 most highly variable genes

plot1 <- VariableFeaturePlot(ms)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2)) # plot variable features with and without labels

ms = ScaleData(ms, features = rownames(ms))

# new normalization from Seurat
# tried regress out the pct_counts_Mt but works less well
#ms <- SCTransform(object = ms, variable.features.n = nfeatures) 
ms <- RunPCA(object = ms, features = VariableFeatures(ms), verbose = FALSE)
ElbowPlot(ms, ndims = 50)

ms <- FindNeighbors(object = ms, reduction = "pca", k.param = 20, dims = 1:20)
ms <- FindClusters(ms, resolution = 12, algorithm = 3)

nb.pcs = 30; n.neighbors = 30; min.dist = 0.4;
ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)

DimPlot(ms, reduction = "umap", group.by = 'seurat_clusters') + ggtitle('scran normalization')

DimPlot(ms, reduction = "umap", group.by = 'timingEst') + ggtitle(paste0(nfeatures, ' HVGs'))

source.my.script("scRNA_cluster_annotation_functions.R")
test.umap.params(seurat.obj = ms)

# save(ms, file=paste0(RdataDir, version.DATA, '_QCleaned_sctransformNorm.Rdata'))
##########################################
# (Optional!!) correct the cell cycle confounder using Seurat
# !!! not used, because there is no clear cell cycle pattern when trying to correct the cell cycle
##########################################
if(correct.cellCycle){
  source.my.script("normalization_HVGs_cellCycle_batchCorrection_functions.R")
  # cellCycle.correction(sce, method = "seurat")
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2.Rdata'))
}

##########################################
# Batch correction using fastMNN from scran
# here we are calling fastMNN from Seurat 
##########################################
Correction.Batch.using.fastMNN = FALSE
if(Correction.Batch.using.fastMNN){
  library(Seurat)
  library(SeuratWrappers)
  library(cowplot)
  
  ## set parameters for fastMNN, most important the merging order
  bcs = unique(ms$request)
  ## the ordered used before
  # c('R7130', 'R8729', 'R8612', 'R8526', 'R7926', # 3 plates for each request
  # 'R6875','R8613','R8348') # 1 plate for each request
  for(n in 1:length(bcs)) 
    eval(parse(text= paste0('p', n, '=  DimPlot(ms, cols.highlight = "red", cells.highlight = as.list(which(ms$request =="', 
                            bcs[n], '"))) + NoLegend() + ggtitle("', bcs[n], '")')))
  
  CombinePlots(plots = list(p1, p2, p3, p4, p5, p6, p7, p8), ncol = 4)
  
  #order2correct = list(list(5, 8, 4), list(list(6, 7), list(2, 1, 3)))
  order2correct = list(2, 6, list(5, 8), 3, 4, 7, 1)
  # HVGs = find.HVGs(sce, Norm.Vars.per.batch = Norm.Vars.per.batch, method = "scran", ntop = 2000) # if use batch-specific HGVs
  # gene.chosen = match(HVGs, rownames(sce))
  # cat("nb of HGV : ", length(gene.chosen), "\n")
  
  ## note that after calling fastMNN, many features stored in Seurat object get lost
  msc <- RunFastMNN(object.list = SplitObject(ms, split.by = "request"), assay = "RNA", 
                    features = VariableFeatures(ms), reduction.name = 'mnn', 
                    cos.norm = TRUE, merge.order = order2correct, min.batch.skip = 0.6)
  metadata(msc@tools$RunFastMNN)$merge.info # check the mergint thresholds
  
  ms@reductions$mnn = Reductions(msc, slot = 'mnn')
  ms@tools$RunFastMNN = msc@tools$RunFastMNN
  
  nb.pcs = 30; n.neighbors = 30; min.dist = 0.3;
  ms <- RunUMAP(object = ms, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                n.neighbors = n.neighbors, min.dist = min.dist)
  ms <- RunUMAP(ms, reduction = "mnn", reduction.name = "umap_mnn", reduction.key = 'umap_mnn_',
                dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  
  p0 =DimPlot(ms, reduction = 'umap',  group.by = c("request"))
  p1 = DimPlot(ms, reduction = 'umap_mnn', group.by = c("request"))
  
  plot_grid(p0, p1)
  
  save(ms, file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_timingEst_Normed_bc_Seurat.Rdata'))
  
}

########################################################
########################################################
# Section : label transfer from John Murray's scRNA-seq data
# 
########################################################
########################################################
#load(file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_timingEst_Normed_bc_Seurat.Rdata'))

#ms <- FindNeighbors(object = ms, reduction = "pca", k.param = 20, dims = 1:20)
#ms <- FindClusters(ms, resolution = 12, algorithm = 3)
library(Seurat)
library(SeuratWrappers)
library(cowplot)
library(ggplot2)

#ElbowPlot(ms)
#ms <- FindNeighbors(object = ms, reduction = 'pca', dims = 1:20)
nb.pcs = 50; n.neighbors = 50; min.dist = 0.4;

ms <- RunUMAP(object = ms, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, n.neighbors = n.neighbors, 
              min.dist = min.dist)
DimPlot(ms, reduction = "umap", group.by = 'SCT_snn_res.12', label = TRUE) 

Idents(ms) = ms$SCT_snn_res.12

DimPlot(ms, reduction = "umap", label = TRUE) 
ms[['umap_sct']] = ms[['UMAP']]
ms[['UMAP']] = NULL

FeaturePlot(ms, reduction = 'umap', features = c('pha-4', 'hnd-1', 'nhr-67'))

saveRDS(ms, file = paste0(RdataDir, 'processed_6.5k.cells_scran.normalized.rds'))


##########################################
# mapping reference with scmap, seurat, svm, rf
# alternatively, sysmatic marker-gene-based cluster annotation can be also applied given that the clusters is well defined, which is 
# also a challenging problem (see Cellassign and Garnet)
##########################################
ms = readRDS(file = paste0(RdataDir, 'processed_6.5k.cells_scran.normalized.rds'))
source.my.script("scRNA_cluster_annotation_functions.R")

ms = reference.based.cluster.annotation(seurat.obj = ms, redefine.clusters = TRUE, predict.unassignedCells = FALSE)

########################################################
########################################################
# Section : here manual annotate BWM lineages using various information:
# 1) clusters (splitting and merging if necessay) 
# 2) predicted labels from seurat and scmap 
# 3) cell size info and estimated timing
# 3.5) marker genes
# 4) cluster connection by PAGA or VarID 
# 5) RNA velocity (not sure ...)
# 
########################################################
########################################################
rdsfile.saved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat.rds')
ms = readRDS(file = rdsfile.saved)

##########################################
# 1) overview of all given clusters and predicted labels
##########################################
## current 54 clusters were define using 3000 variable genes and resolution =3, 20 pcs and k = 10
p0 = DimPlot(ms, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
             na.value = "gray") + 
  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
  scale_colour_hue(drop = FALSE) + 
  NoLegend()

plot(p0)

## compare scmap and seurat reference-based annotation
source.my.script('scRNA_cluster_annotation_functions.R')
overview.and.compare.predicted.labels(seurat.obj = ms)

##########################################
# 2) focus short list of cell identities and manual annotate with other information
##########################################
source.my.script('scRNA_cluster_annotation_functions.R')

ms$manual.annot.ids = NA
ms = manual.annotation.for.BWM.clusters(seurat.obj = ms)


##########################################
# read annotated seurat.obj and save loom file for PAGA analysis
##########################################
# install scater https://bioconductor.org/packages/release/bioc/html/scater.html
library(scater)
# install loomR from GitHub using the remotes package remotes::install_github(repo =
# 'mojaveazure/loomR', ref = 'develop')
library(loomR)
library(Seurat)
library(patchwork)

nb.iteration = 36
RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_', 
                  nb.iteration, '.rds')
seurat.obj = readRDS(file = RDSsaved)

## select manually annotated bwm cells
ids.bwm = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.bwm))])
sub.obj = subset(seurat.obj, cells = cells.sels)

#saveRDS(sub.obj, file = paste0(RdataDir, 'manual_annotated_BWM_3.5k.cells.rds'))
source.my.script('scRNA_cluster_annotation_functions.R')
# compare.bwm.with.JMdata(sub.obj)

manual.annot.ids = sub.obj$manual.annot.ids
sub.obj@meta.data[, grep('pred|ids|scmap|previous.iteration.clusters|BWM.cells', colnames(sub.obj@meta.data))] = NULL
sub.obj$manual.ids = manual.annot.ids

sub.loom <- as.loom(sub.obj, filename = "seuratObj_BWM_manual_cellIds_iteration_36.loom", verbose = FALSE)
sub.loom
sub.loom$close_all()


########################################################
########################################################
# Section : regulator prediction using sc_TF_MARA.R
# ideally, at the end we will have a package easily to be install in Github
# inputs will be either mRNA or pre-mRNA matrix
# the data will be processed either cluster-wise (id-wise) or cell-wise
# imputation will be done for cell-wise analysis
# dynamic genes will be identified using gam package
# MARA is used for inferring TF actvities 
########################################################
########################################################
resDir = paste0("results/", version.analysis, '/annoted_BWM')
tabDir = paste0(resDir, "/tables/")

if(!dir.exists("results/")){dir.create("results/")}
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}

sub.obj = readRDS(file = paste0(RdataDir, 'manual_annotated_BWM_3.5k.cells.rds'))
sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa')] =
  'mixture_terminal_1'
sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap')] = 
  'mixture_terminal_2'

source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
require(tictoc)
tic()
test.umap.params.for.BWM.cells(sub.obj, 
                               pdfname = 'UMAP_param_TEST_BWM_all.pdf',
                               group.by = 'manual.annot.ids', with_legend = FALSE,
                               nfeatures.sampling = c(3000, 5000, 8000), nb.pcs.sampling = c(10, 20, 30, 50),
                               n.neighbors.sampling = c(5, 10, 30, 50), 
                               min.dist.sampling = c(0.01, 0.1)
)
toc()

##########################################
# run the sctf_MARA
##########################################
source.my.script('sctf_MARA.R')

res = predict.TF.MARA.for.scdata(sub.obj, mode = 'cluster.wise', id = 'manual.annot.ids')










