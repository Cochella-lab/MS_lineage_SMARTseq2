##########################################################################
##########################################################################
# Project:
# Script purpose:
# Usage example:
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Nov 19 16:30:34 2019
##########################################################################
##########################################################################
#Function for soursing functions
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
setwd(paste0("INSERT HERE WORKING DIRECTORY",user))

version.DATA = 'all_batches'
version.analysis =  paste0(version.DATA, '_20191115')
dataDir = paste0("../data/gene_counts/")
resDir = paste0("results/", version.analysis)
tabDir = paste0("results/", version.analysis, "/tables/")

RdataDir = paste0("results/", version.analysis, "/Rdata/")
if(!dir.exists("results/")){dir.create("results/")}
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

correct.cellCycle = FALSE

########################################################
########################################################
# Section : timingEst with cpm normalization and add it to the metadata
########################################################
########################################################
## import the R object from the previous step and double check the cells and genes from table
load(file=paste0("../data/R_processed_data/", version.DATA, '_QCed_cells_genes_filtered_SCE.Rdata'))

library(scater)
library(scran)
options(stringsAsFactors = FALSE)

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
# here estimate the timing with timer genes
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
  sce = estimate.timing.and.variance.with.timer.genes(sce)
  
  save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_timingEst_SCE.Rdata'))
  
}else{
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_timingEst_SCE.Rdata'))
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
load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_timingEst_SCE.Rdata'))

reducedDim(sce) <- NULL
endog_genes <- !rowData(sce)$is_feature_control

Normalization.Testing = FALSE
source.my.script("normalization_HVGs_cellCycle_batchCorrection_functions.R")

if(Normalization.Testing){
  pdfname = paste0(resDir, "/scRNAseq_filtered_normalization_testing.pdf")
  pdf(pdfname, width=14, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  compare.scran.seurat.sctransform(sce, using.HVGs = TRUE)
  
  dev.off()
}

## add some extra stat for sce (select normalization method: sctransform or scran ())
sce$library.size = apply(counts(sce), 2, sum)

library(Seurat)
library(ggplot2)

## convert sce to seurat object
ms = as.Seurat(sce, counts = 'counts', data = NULL, assay = "RNA")

nfeatures = 2000 # this can be variable
# new normalization from Seurat
# tried regress out the pct_counts_Mt but works less well
ms <- SCTransform(object = ms, variable.features.n = nfeatures) 
ms <- RunPCA(object = ms, features = VariableFeatures(ms), verbose = FALSE)
ElbowPlot(ms)

nb.pcs = 20; n.neighbors = 30; min.dist = 0.3;
ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(ms, reduction = "umap", group.by = 'request') + ggtitle('sctransform normalization')
DimPlot(ms, reduction = "umap", group.by = 'timingEst') + ggtitle('2000 HVGs')

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
# Not really important to call this function
##########################################
Correction.Batch.using.fastMNN = TRUE
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
    eval(parse(text= paste0('p', n, '=  DimPlot(ms, cols.highlight = "red", cells.highlight = as.list(which(ms$request =="', bcs[n], '"))) + NoLegend() + ggtitle("', bcs[n], '")')))
  
  CombinePlots(plots = list(p1, p2, p3, p4, p5, p6, p7, p8), ncol = 4)
  
  #order2correct = list(list(5, 8, 4), list(list(6, 7), list(2, 1, 3)))
  order2correct = list(2, 6, list(5, 8), 3, 4, 7, 1)
  # HVGs = find.HVGs(sce, Norm.Vars.per.batch = Norm.Vars.per.batch, method = "scran", ntop = 2000) # if use batch-specific HGVs
  # gene.chosen = match(HVGs, rownames(sce))
  # cat("nb of HGV : ", length(gene.chosen), "\n")
  
  ## note that after calling fastMNN, many features stored in Seurat object get lost
  msc <- RunFastMNN(object.list = SplitObject(ms, split.by = "request"), assay = "SCT", 
                        features = VariableFeatures(ms), reduction.name = 'mnn', 
                        cos.norm = TRUE, merge.order = order2correct, min.batch.skip = 0.6)
  metadata(msc@tools$RunFastMNN)$merge.info # check the mergint thresholds
  
  ms@reductions$mnn = Reductions(msc, slot = 'mnn')
  ms@tools$RunFastMNN = msc@tools$RunFastMNN
  
  nb.pcs = 20; n.neighbors = 30; min.dist = 0.3;
  ms <- RunUMAP(object = ms, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  ms <- RunUMAP(ms, reduction = "mnn", reduction.name = "umap_mnn", reduction.key = 'umap_mnn_',dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  
  p0 =DimPlot(ms, reduction = 'umap',  group.by = c("request"))
  p1 = DimPlot(ms, reduction = 'umap_mnn', group.by = c("request"))
  
  plot_grid(p0, p1)
  
  save(ms, file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_timingEst_Normed_bc_Seurat.Rdata'))

}




