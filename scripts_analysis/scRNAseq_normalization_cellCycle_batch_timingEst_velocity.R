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
setwd(paste0("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/",user))

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
#
########################################################
########################################################
## import the R object from the previous step and double check the cells and genes from table
load(file=paste0("../data/R_processed_data/", version.DATA, '_QCed_cells_genes_filtered_SCE.Rdata'))

library(scater)

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
  Test.Hashimshony_lineages = FALSE
  
  if(Test.Hashimshony_lineages){
    pdfname = paste0("../results/clustering_combining_variousInfos/test_timing_estimation_Hashimshony_lineages_test_with_improvedTimerGenes_v2.pdf")
    pdf(pdfname, width=10, height = 6)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    source('customizedClustering_functions.R')
    Test.timingEstimate.with.HashimshonyLineages(fastEstimate = TRUE, timerGenes.pval = 0.0001, lineageCorrs = 0.7,  loess.span = 0.5, lowFilter.threshold.target = 5,
                                                 PLOT.test = FALSE)
    
    dev.off()
    
  }
  
  ## Here we are sampling a range of parameters and timing estimation were done with each of them
  ## Whereby we assess the sensibility of our timingEst
  ## this will take some time to finish
  source('customizedClustering_functions.R')
  
  timingEst = c()
  for(pv in c(0.001, 0.0001, 0.00001))
  {
    for(cutoff.expr in c(4, 5, 6))
    {
      for(s in c(0.3, 0.5, 0.7))
      {
        cat('pv = ', pv, ' cutoff.expr = ', cutoff.expr, 's = ', s, "\n")
        sce.test = sc.estimateTiming.with.timer.genes(sce, fastEstimate = TRUE, timerGenes.pval = pv, lineageCorrs = 0.5, loess.span = s,
                                                      lowFilter.threshold.target = cutoff.expr)
        timingEst = rbind(timingEst, sce.test$timing.est)
      }
    }
  }
  
  #save(sce, timingEst, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrectedv2_facsInfos_timeEst_tmp.Rdata'))
  
  timingEst = t(timingEst)
  timingEst = as.matrix(timingEst)
  #colnames(timingEst) = c(as.vector(t(outer(c(0.001, 0.0001, 0.00001), c(4:6), paste, sep=""))))
  find.one.close.to.mean = function(x){
    # x = timingEst[1, ]
    difs = abs(x - mean(x))
    return(x[which(difs == min(difs))[1]])
  }
  sce$timingEst = apply(timingEst, 1, find.one.close.to.mean)
  sce$timingEst.sd = apply(timingEst, 1, sd)
  
  save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos_timingEst.Rdata'))
}else{
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos_timingEst.Rdata'))
}


par(mfrow = c(1, 3))
plot(sce$FSC_log2, sce$timingEst, type='p', cex = 0.5)
plot(sce$BSC_log2, sce$timingEst, type='p', cex = 0.5)
plot(sce$FSC_log2, sce$BSC_log2, type = 'p', cex = 0.5)

cdata = colData(sce)
cdata = data.frame(cdata[, c((ncol(cdata)-3): ncol(cdata))])
#cdata$timingEst = cdata$timingEst/50

cdata$timing.group = NA
cdata$timing.group[which(cdata$timingEst < 50)] = 1
cdata$timing.group[which(cdata$timingEst >= 450)] = 10
for(n in 2:9){cdata$timing.group[which(cdata$timingEst >= (n-1)*50 & cdata$timingEst < n*50)] = n}

cdata$timing.sd.group = 3
cdata$timing.sd.group[which(cdata$timingEst.sd<30)] = 1
cdata$timing.sd.group[which(cdata$timingEst.sd>=30 & cdata$timingEst.sd<60)] = 2
cdata$timing.sd.group = as.factor(cdata$timing.sd.group)

ggplot(cdata, aes(x=FSC_log2, y=BSC_log2, color=timingEst, shape = timing.sd.group)) +
  geom_point() +
  scale_color_gradientn(colours = rainbow(10))

sce$timingEst = as.factor(sce$timingEst)
sce$timingEst.group = as.factor(cdata$timing.group)
sce$timingEst.sd.group = as.factor(cdata$timing.sd.group)

save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos_timingEstGroups.Rdata'))

plotColData(sce,
            x = "FSC_log2",
            y = "BSC_log2",
            colour_by = "timingEst.group",
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
library(Seurat)
library(scRNA.seq.funcs)
library(scater)
library(scran)
options(stringsAsFactors = FALSE)

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

#source.my.script("scRNAseq_functions.R")
#gg.Mt = find.particular.geneSet("Mt")
#is.mito <- rownames(sce) %in% gg.Mt;

## convert sce to seurat object
ms = as.Seurat(sce, counts = 'counts', data = NULL, assay = "RNA")
nfeatures = 3000

# new normalization from Seurat
# tried regress out the pct_counts_Mt but works less well
ms <- SCTransform(object = ms, variable.features.n = nfeatures) 
ms <- RunPCA(object = ms, features = VariableFeatures(ms), verbose = FALSE)
ElbowPlot(ms)

nb.pcs = 20; n.neighbors = 30; min.dist = 0.4;
ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(ms, reduction = "umap", group.by = 'request') + ggtitle('sctransform')

#save(ms, file=paste0(RdataDir, version.DATA, '_QCleaned_sctransformNorm.Rdata'))
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
library(Seurat)
library(SeuratWrappers)

# ## set parameters for fastMNN 
# sce$timingEst = as.factor(sce$timingEst)
# sce$timingEst.group = as.factor(sce$timingEst.group)
# 
# # choose the batches (either plates or request)
# # here we choose the request as batch
# batches = sce$request
# bc.uniq = unique(batches)
# sce$batches <- batches
# 
# Use.fastMNN = TRUE
# Norm.Vars.per.batch = TRUE # HVGs for each batch or not
# Rescale.Batches = FALSE # scaling data in each batch or not
# k.mnn = 20
# cos.norm = TRUE
# nb.pcs = 50
# 
# batch.sequence.to.merge = c('R7130', 'R8729', 'R8612', 'R8526', 'R7926', # 3 plates for each request
#                             'R6875','R8613','R8348') # 1 plate for each request
# 
# order2correct = match(batch.sequence.to.merge, bc.uniq)
# #c(12, 13, # the same
# #                10, 9, 8, 5, 6, 7,  8, 11, 1, 2, 3, 4)
# #order2correct = c(15,14,13,12,11,10,9,8, 1, 2, 3, 4, 5,6,7 )
# #order2correct = c(3, 4, 1, 2)
# 
# ## double chekc  the mering order in the batch correction
# source('customizedClustering_functions.R')
# kk = match(sce$request, c('R7130', 'R8729', 'R8612', 'R8526', 'R7926'))
# kk = match(sce$request, c('R6875','R8613','R8348'))
# plotColData(sce[,which(!is.na(kk))],
#             x = "FSC_log2",
#             y = "BSC_log2",
#             colour_by = "request",
#             point_size = 1
# )
# 
# #cat("merging order for batch correction :\n", paste0(bc.uniq[order2correct], collapse = "\n"), "\n")
# for(n in 1:length(order2correct)){
#   
#   #n = 11
#   kk = order2correct[n]
#   
#   p = plotColData(sce[, which(sce$batches== bc.uniq[kk])],
#                   x = "FSC_log2",
#                   y = "BSC_log2",
#                   colour_by = "timingEst",
#                   point_size = 1
#                   
#   )
#   plot(p)
#   
#   cat('#', n, 'batch:',  bc.uniq[kk], ': ', length(which(sce$batches == bc.uniq[kk])), 'cells \n')
#   
# }
# 
# 
# HVGs = find.HVGs(sce, Norm.Vars.per.batch = Norm.Vars.per.batch, method = "scran", ntop = 2000)
# gene.chosen = match(HVGs, rownames(sce))
# cat("nb of HGV : ", length(gene.chosen), "\n")

pbmcsca <- RunFastMNN(object.list = SplitObject(ms, split.by = "request"), assay = "SCT", 
                      features = VariableFeatures(ms), reduction.name = 'mnn', 
                      cos.norm = TRUE, auto.merge = FALSE, min.batch.skip = 0.5)

nb.pcs = 20; n.neighbors = 30; min.dist = 0.25;
pbmcsca <- RunUMAP(pbmcsca, reduction = "mnn", dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)

p1 = DimPlot(pbmcsca, group.by = c("request"))
p0 =DimPlot(ms, reduction = 'umap',  group.by = c("request"))
plot_grid(p0, p1)

metadata(pbmcsca@tools$RunFastMNN)$merge.info
#library(harmony)
#pbmcsca <- RunHarmony(ms, group.by.vars = "request", assay.use="SCT")
#pbmcsca <- RunUMAP(pbmcsca, reduction = "harmony", dims = 1:20)
#pbmcsca <- FindNeighbors(pbmcsca, reduction = "harmony", dims = 1:30) %>% FindClusters()
#DimPlot(pbmcsca, reduction = 'harmony',  group.by = c("request"))
#DimPlot(pbmcsca, reduction = 'umap',  group.by = c("request"))


save(sce, file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_bcMNN.Rdata'))

##########################################
# Convert sce object to Seurat object
# check UMAP and tSNE
##########################################
require(Seurat)
pbmc = as.Seurat(sce)

#pbmc = Seurat::RunPCA(pbmc, pbmc, )
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = rownames(pbmc))
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#pbmc <- RunPCA(pbmc, features = HVGs)

pbmc <- Seurat::RunUMAP(pbmc, dims = 1:15, reduction = "MNN",
                        reduction.key = "umap", n.neighbors = 20, repulsion.strength = 1)

DimPlot(pbmc, reduction = "umap", group.by = 'timingEst')

pbmc <- Seurat::RunUMAP(pbmc, dims = 1:15, reduction = "pca",
                        reduction.key = "umap.pca", n.neighbors = 20, repulsion.strength = 1)

DimPlot(pbmc, reduction = "umap", group.by = 'timingEst')
