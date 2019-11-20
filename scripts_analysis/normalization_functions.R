##########################################################################
##########################################################################
# Project:
# Script purpose: functions for scRNA-seq normalization
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Nov 20 12:15:18 2019
##########################################################################
##########################################################################
##########################################
# several common functions for normalizations 
# from Hemberg lab
# hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#normalization-theory
##########################################
calc_cpm <- function (expr_mat, spikes = NULL) 
{
  norm_factor <- colSums(expr_mat[-spikes, ])
  return(t(t(expr_mat)/norm_factor)) * 10^6
}

Down_Sample_Matrix <- function (expr_mat)
{
  min_lib_size <- min(colSums(expr_mat))
  down_sample <- function(x) {
    prob <- min_lib_size/sum(x)
    return(unlist(lapply(x, function(y) {
      rbinom(1, y, prob)
    })))
  }
  down_sampled_mat <- apply(expr_mat, 2, down_sample)
  return(down_sampled_mat)
}

cal_uq_Hemberg = function (expr_mat, spikes = NULL) 
{
  UQ <- function(x) {
    quantile(x[x > 0], 0.75)
  }
  if(!is.null(spikes)){
    uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
  }else{
    uq <- unlist(apply(expr_mat, 2, UQ))
  }
  
  norm_factor <- uq/median(uq)
  return(t(t(expr_mat)/norm_factor))
}

calculate.sizeFactors.DESeq2 = function(expr_mat)
{
  # expr_mat = counts(sce.qc)
  require('DESeq2')
  condition <- factor(rep("A", ncol(expr_mat)))
  dds <- DESeqDataSetFromMatrix(expr_mat, DataFrame(condition), design = ~ 1)
  dds <- estimateSizeFactors(dds)
  
  return(sizeFactors(dds))
  
}

test.normalization = function(sce, Methods.Normalization = c("cpm", "DESeq2", "scran", "seurat", "SCnorm", "sctransform"), min.size = 100)
{
  #Methods.Normalization = "DESeq2" 
  
  for(method in Methods.Normalization)
  {
    sce.qc = sce
    set.seed(1234567)
    
    method = 'cpm'
     
    cat('test normalization method -- ', method, "\n")
    main = paste0(method);
    
    if(method == "raw") { # raw log counts
      assay(sce.qc, "logcounts") <- log2(counts(sce.qc) + 1)
    }
    
    if(method == "cpm") { ### cpm
      #assay(sce.qc, "logcounts") <- log2(calculateCPM(sce.qc, use_size_factors = FALSE) + 1)
      sce.qc = logNormCounts(sce, size_factors = NULL, log = TRUE, pseudo_count =1)
    }
    if(method == "UQ"){
      logcounts(sce.qc) <- log2(cal_uq_Hemberg(counts(sce.qc)) + 1)
    }
    
    if(method == "DESeq2"){
      sizeFactors(sce.qc) = calculate.sizeFactors.DESeq2(counts(sce.qc))
      sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
    }
    if(method == "downsample") {
      assay(sce.qc, "logcounts") <- log2(Down_Sample_Matrix(counts(sce.qc)) + 1)
    }
    
    if(method == "scran"){
      ## scran normalization (not working here, because negative scaling factor found)
      qclust <- quickCluster(sce.qc, min.size = min.size,  method = 'igraph')
      sce.qc <- computeSumFactors(sce.qc, clusters = qclust)
      sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
    }
    
    if(method == "TMM"|method == "DESeq2"|method == "UQ"|method == "scran"){
      summary(sizeFactors(sce.qc))
      range(sizeFactors(sce.qc))
      
      plot(sce.qc$total_counts/1e6, sizeFactors(sce.qc), log="xy", main = paste0(method), 
           xlab="Library size (millions)", ylab="Size factor",
           pch=16)
      #legend("bottomright", col=c("black"), pch=16, cex=1.2, legend = "size factor from scran vs total library size")
    }
    
    if(method != 'Seurat' & method != "sctransform"){
      sce.qc = scater::runPCA(sce.qc)
      param.perplexity = 10; set.seed(100)
      sce.qc = scater::runTSNE(sce.qc,dimred = 'PCA', perplexity = param.perplexity)
      sce.qc = scater::runUMAP(sce.qc, dimred = 'PCA')
      #sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA", perplexity=5)
      #runUMAP(sce.zeisel, dimred="PCA")
      scater::plotReducedDim(sce.qc, dimred = "PCA", colour_by = 'request')
      scater::plotReducedDim(sce.qc, dimred = "UMAP", colour_by = 'request')
      scater::plotReducedDim(sce.qc, dimred = "TSNE", colour_by = 'request')
      
      # p1 = scater::plotPCA(
      #   sce.qc[endog_genes, ],
      #   run_args = list(exprs_values = "logcounts"), 
      #   size_by = "total_counts",
      #   #size_by = "total_features_by_counts",
      #   colour_by = "request"
      # ) + ggtitle(paste0("PCA -- ", main))
      # 
      # 
      # p2 = plotTSNE(
      #   sce.qc[endog_genes, ],
      #   run_args = list(exprs_values = "logcounts", perplexity = param.perplexity), 
      #   size_by = "total_counts",
      #   #size_by = "total_features_by_counts",
      #   colour_by = "seqInfos"  
      # ) + ggtitle(paste0("tSNE - perplexity = ", param.perplexity, "--", main))
      # 
      # p3 = plotUMAP(
      #   sce.qc[endog_genes, ],
      #   run_args = list(exprs_values = "logcounts"), 
      #   size_by = "total_counts",
      #   #size_by = "total_features_by_counts",
      #   colour_by = "seqInfos"
      # ) + ggtitle(paste0("UMAP -- ", main))
      # 
      # plot(p1); plot(p2); plot(p3)
    }
    
    
    if(method == 'seurat'| method == 'sctransform'){
      ms0 = as.Seurat(sce, counts = 'counts', data = NULL, assay = "RNA")
      
      if(method == 'seurat'){
        ms.logtransform <- NormalizeData(ms0, assay = "RNA" )
        ms.logtransform <- FindVariableFeatures(ms.logtransform, selection.method = "vst", nfeatures = 3000)
        all.genes <- rownames(ms.logtransform)
        ms.logtransform <- ScaleData(ms.logtransform, features = all.genes)
        
        ms.logtransform <- RunPCA(object = ms.logtransform, verbose = FALSE)
        
        #ms.logtransform <- FindNeighbors(object = ms.logtransform, dims = 1:20)
        #ms.logtransform <- FindClusters(object = ms.logtransform)
        
        ms.logtransform <- RunUMAP(object = ms.logtransform, reduction = 'pca', dims = 1:20, n.neighbors = 20)
        
        ms.logtransform = RunTSNE(ms.logtransform, reduction = 'pca', dims = 1:20, tsne.method = 'Rtsne')
        DimPlot(ms.logtransform, reduction = "umap", group.by = 'request')
        
        DimPlot(ms.logtransform, reduction = "tsne", group.by = 'request')
      }
      
      if(method == 'sctransform'){
        ms <- SCTransform(object = ms0) # new normalization from Seurat
        
        ms <- RunPCA(object = ms, verbose = FALSE)
        #ms <- FindNeighbors(object = ms, dims = 1:20)
        #ms <- FindClusters(object = ms)
        
        ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:20)
        DimPlot(ms, reduction = "umap", group.by = 'request')
        
        #ms <- RunUMAP(object = ms, reduction = 'MNN', dims = 1:20, n.neighbors = 30)
        #DimPlot(ms, reduction = "umap")
        
      }
    }
    
  }
}
