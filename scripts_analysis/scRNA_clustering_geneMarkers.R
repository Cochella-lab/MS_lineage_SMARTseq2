########################################################
########################################################
# Section : clustering and DE analysis or gene markers discovery
# feature selection and dimension reduction was done already in the batch correction part
# So we start with the PCs from fastMNN in scran, with which we define distance for clustering
# 1) different clustering methods will be tested  
# 2) special design matrix will be used for DE analysis 
# http://bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/de.html#2_blocking_on_uninteresting_factors_of_variation
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

library(Seurat)
library(SeuratWrappers)
library(cowplot)
library(ggplot2)
########################################################
########################################################
# Section : annotate scRNA clusters by mapping to reference
# 
########################################################
########################################################
dir.processed.data = "/Users/jiwang/workspace/imp/scRNAseq_MS_lineage_dev/results_aleks/results/all_batches_202005/Rdata/"
load(paste0(dir.processed.data, "all_batches_QCed_cells_genes_filtered_timingEst_Normed_bc_Seurat.Rdata"))

ms <- FindNeighbors(object = ms, reduction = "mnn", k.param = 20, dims = 1:20)
ms <- FindClusters(ms, resolution = 12, algorithm = 3)

DimPlot(object = ms, cells = colnames(ms), group.by = 'ident', label = TRUE, pt.size = 1)

all_ms.markers <- FindAllMarkers(ms, min.pct = 0.25, logfc.threshold = 0.25)
all_ms.top10 <- all_ms.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)


#length(unique(all_ms.top10$gene))

lin_sc_expr_190602 <- readRDS("../../../../data/JM_data/lin_sc_expr_190602.rds")
MS_names <- lin_sc_expr_190602[,grepl(pattern = "MS.", x = colnames(lin_sc_expr_190602))]
#lin_sc_expr_190602[,MS_names]

library(dplyr)

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

pdfname = paste0(paste0("../idents_new.pdf"))
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

saveRDS(ms_correlation_idents, file = './ms_correlation_idents.rds')


########################################################
# Section : Clustering section by integrating various informations: 
# gene expression, fac info, estimated timing and velocity 
########################################################
########################################################
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