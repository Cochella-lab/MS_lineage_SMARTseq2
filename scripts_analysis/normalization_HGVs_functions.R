##########################################################################
##########################################################################
# Project:
# Script purpose: functions for scRNA-seq normalization
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Nov 20 12:15:18 2019
##########################################################################
##########################################################################
########################################################
########################################################
# Section : Compare normalization methods for scRNA-seqe 
# to choose the good one
########################################################
########################################################

# several common functions for normalizations 
# from Hemberg lab
# hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#normalization-theory
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

test.normalization = function(sce, Methods.Normalization = c("cpm", "DESeq2", "scran", "seurat", "sctransform"), use.HVGs = TRUE)
{
  #Methods.Normalization = "DESeq2"
  for(method in Methods.Normalization)
  {
    sce.qc = sce
    sce.qc$library.size = apply(counts(sce.qc), 2, sum)
    set.seed(1234567)
    
    method = 'cpm'
     
    cat('test normalization method -- ', method, "\n")
    main = paste0(method);
    
    if(method == "raw") { # raw log counts
      assay(sce.qc, "logcounts") <- log2(counts(sce.qc) + 1)
    }
    
    if(method == "cpm") { ### cpm
      #assay(sce.qc, "logcounts") <- log2(calculateCPM(sce.qc, use_size_factors = FALSE) + 1)
      sce.qc = logNormCounts(sce.qc, size_factors = NULL, log = TRUE, pseudo_count=1, center_size_factors = TRUE)
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
      #min.size = 100
      ## scran normalization (not working here, because negative scaling factor found)
      qclust <- quickCluster(sce.qc,  method = 'igraph')
      sce.qc <- computeSumFactors(sce.qc, clusters = qclust)
      sce.qc <- logNormCounts(sce.qc)
    }
    
    if(method == "TMM"|method == "DESeq2"|method == "UQ"|method == "scran"){
      summary(sizeFactors(sce.qc))
      range(sizeFactors(sce.qc))
      plot(sce.qc$library.size/1e6, sizeFactors(sce.qc), log="xy", main = paste0(method), 
           xlab="Library size (millions)", ylab="Size factor",
           pch=16)
      #legend("bottomright", col=c("black"), pch=16, cex=1.2, legend = "size factor from scran vs total library size")
    }
    
    if(method != 'Seurat' & method != "sctransform"){
      
      ## to have a fair comparison, here we also need use HVGs for PCAs, UMAP and t-SNE
      HVGs = find.HVGs(sce.qc, Norm.Vars.per.batch = FALSE, method = "scran", ntop = 2000)
      
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
      
      if(method == 'sctransform'){
        ms <- SCTransform(object = ms0) # new normalization from Seurat
        
        ms <- RunPCA(object = ms, verbose = FALSE)
        #ms <- FindNeighbors(object = ms, dims = 1:20)
        #ms <- FindClusters(object = ms)
        
        ElbowPlot(ms)
        
        ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:20, n.neighbors = 30, umap.method = "uwot",
                      metric = 'correlation', min.dist = 0.25)
        DimPlot(ms, reduction = "umap", group.by = 'request')
        
        
        ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:20, n.neighbors = 30, umap.method = "uwot",
                      metric = 'cosine', min.dist = 0.25)
        DimPlot(ms, reduction = "umap", group.by = 'request')
        
        
        ms <- RunTSNE(object = ms, reduction = 'pca', dims = 1:20)
        #ms <- RunUMAP(object = ms, reduction = 'MNN', dims = 1:20, n.neighbors = 30)
        DimPlot(ms, reduction = "tsne", group.by = 'request')
        
      }
      
      
      if(method == 'seurat'){
        #scale.factor = 10^4
        #scale.factor = mean(sce.qc$library.size)
        ms.logtransform <- NormalizeData(ms0, assay = "RNA", normalization.method = 'LogNormalize')
        ms.logtransform <- FindVariableFeatures(ms.logtransform, selection.method = "vst", nfeatures = 2000)
        ms.logtransform <- ScaleData(ms.logtransform, features = rownames(ms.logtransform))
        
        ms.logtransform <- RunPCA(object = ms.logtransform, verbose = FALSE)
        ElbowPlot(ms.logtransform)
        
        ms.logtransform <- RunUMAP(object = ms.logtransform, reduction = 'pca', dims = 1:20, n.neighbors = 30, min.dist = 0.25)
        DimPlot(ms.logtransform, reduction = "umap", group.by = 'request')
        
        ms.logtransform = RunTSNE(ms.logtransform, reduction = 'pca', dims = 1:10, tsne.method = 'Rtsne')
        DimPlot(ms.logtransform, reduction = "tsne", group.by = 'request')
      }
      
      ##########################################
      # Compare sctransform and standard normalization
      ##########################################
      Compare.two.norm.in.seurat.and.scran = FALSE
      if(Compare.two.norm.in.seurat){
        library(ggplot2)
        library(cowplot)
        p1 =DimPlot(ms.logtransform, reduction = "umap", group.by = 'request')
        p2 = DimPlot(ms, reduction = "umap", group.by = 'request')
        plot_grid(p1, p2)
      }
      
      Dissect.normalization.RunPCA.RunUMAP.in.Seurat.by.comparing.scater = FALSE
      if(Dissect.normalization.RunPCA.RunUMAP.in.Seurat.by.comparing.scater)
      {
        ## compare seurat normalizatioin vs 
        sce.qc = sce
        sce.qc$library.size = apply(counts(sce.qc), 2, sum)
        sce.qc = logNormCounts(sce.qc, size_factors = NULL, log = TRUE, pseudo_count=1, center_size_factors = FALSE)
        plot(sce.qc$library.size/1e6, sizeFactors(sce.qc), log="xy", main = paste0(method), 
             xlab="Library size (millions)", ylab="Size factor")
        
        HVGs = VariableFeatures(ms.logtransform)
        mm = match(VariableFeatures(ms.logtransform), rownames(sce.qc))
        sce.qc = scater::runPCA(sce.qc,  subset_row = mm, scale = TRUE)
        
        sce.qc = scater::runUMAP(sce.qc, dimred = 'PCA', n_dimred = 1:20,  n_neighbors = 20, metric = "cosine")
        scater::plotReducedDim(sce.qc, dimred = "UMAP", colour_by = 'request')
        param.perplexity = 50; set.seed(100)
        sce.qc = scater::runTSNE(sce.qc,dimred = 'PCA', perplexity = param.perplexity)
        scater::plotReducedDim(sce.qc, dimred = "TSNE", colour_by = 'request')
        
        kk = 1
        xx = log(counts(sce)[,kk]/sce.qc$library.size[kk]*10^4 +1)
        plot(xx, ms.logtransform@assays$RNA@data[, kk], log= 'xy')
        abline(0, 1, col='red')
        
        yy = log2((counts(sce)[,kk])/sizeFactors(sce.qc)[1] +1)
        plot(yy, logcounts(sce.qc)[, kk], log= 'xy')
        abline(0, 1, col='red')
        
        plot(logcounts(sce.qc)[, kk],  ms.logtransform@assays$RNA@data[, kk]/log(2))
        abline(0, 1, col='red')
        
        plot(reducedDim(sce.qc, "PCA")[, kk], Reductions(ms.logtransform, slot = 'pca')@cell.embeddings[,kk])
        
      }
      
    }
    
  }
}



########################################################
########################################################
# Section : HVGs selection
# 
########################################################
########################################################
find.HVGs = function(sce, Norm.Vars.per.batch = TRUE, method = "scran", ntop = NULL)
{
  if(method == "scran"){
    if(!Norm.Vars.per.batch){
      ## here we use batches as blockes, i.e. calculated mean and variances separately for each batches and then fit the trend
      ## In doing so, we implicitly assume that the trend is the same between batches
      ## https://bioconductor.org/packages/3.10/workflows/vignettes/simpleSingleCell/inst/doc/var.html#41_using_the_block=_argument
      fit <- trendVar(sce, block=sce$batches, parametric=TRUE, assay.type="logcounts", use.spikes=FALSE)
      dec <- decomposeVar(sce, fit)
      
      par(mfrow=c(1,2))
      cols = c(1:length(bc.uniq)) 
      nb.genes = min(table(sce$batches))
      
      matplot(fit$means[c(1:nb.genes), ], fit$vars[c(1:nb.genes), ], col=cols,
              xlab="Mean log-expression", ylab="Variance of log-expression")
      curve(fit$trend(x), add=TRUE, col="red")
      
      plot(dec$mean, dec$total, xlab="Mean log-expression", 
           ylab="Variance of log-expression", pch=16)
      curve(fit$trend(x), col="dodgerblue", add=TRUE)
      
      tmp.sce <- sce
      tmp.sce$log_size_factor <- log(sizeFactors(sce))
      plotColData(tmp.sce, x="batches", y="log_size_factor")
      #p2 <- plotColData(tmp.416B, x="Plate", y="log_size_factor_ERCC")
      #multiplot(p1, p2, cols=2)
      
      dec.sorted <- dec[order(dec$bio, decreasing=TRUE), ]
      head(dec.sorted)
      
      if(is.null(ntop)){
        # here HVGs selected with FDR < 0.01
        gene.chosen <- rownames(dec.sorted)[which(dec.sorted$FDR < 0.05)] 
      }else{
        gene.chosen = rownames(dec.sorted)[1:ntop] 
      }
      
    }else{
      ## here we search batch-specific features
      ## https://bioconductor.org/packages/3.10/workflows/vignettes/simpleSingleCell/inst/doc/var.html#42_fitting_batch-specific_trends
      sce.bc = sce
      ## recommendation in the function help: Centering should be performed by running multiBlockNorm before calling this function. 
      sce.bc <- multiBlockNorm(sce.bc, block = sce.bc$batches)
      
      par(mfrow=c(1,1))
      plot(sizeFactors(sce), sizeFactors(sce.bc), log='xy'); abline(0, 1, lwd =2, col = 'red') # did not change anything here, weired
      
      comb.out <- multiBlockVar(sce.bc, block=sce.bc$batches, assay.type="logcounts", trend.args=list(parametric=TRUE, use.spikes=FALSE))
      
      #comb.out = multiBlockVar(sce.bc, block = sce.bc$batches, trend.args = list(use.spikes = FALSE))
      head(comb.out[,1:6])
      
      par(mfrow=c(1, length(bc.uniq)))
      #is.spike <- isSpike(sce.416B.2)
      for (plate in unique(sce.bc$batches)) {
        cur.out <- comb.out$per.block[[plate]]
        plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
             ylab="Variance of log-expression", main=plate, ylim = c(0, 25), xlim = c(0, 15))
        curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
        #points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
      }
      
      dec.bc.sorted <- comb.out[order(comb.out$bio, decreasing=TRUE), ]
      head(dec.bc.sorted)
      
      if(is.null(ntop)){
        #length(which(dec.sorted$bio>0))
        # here HVGs selected with FDR<0.01
        #gene.chosen.bc <- which(dec.bc.sorted$p.value < 0.05)
        #gene.chosen.bc <- which(dec.bc.sorted$FDR < 0.1)
        gene.chosen.bc = rownames(dec.bc.sorted)[which(dec.bc.sorted$bio>0)]
        #length(which(dec.sorted$bio>0)) 
      }else{
        gene.chosen.bc = rownames(dec.bc.sorted)[1:ntop]
      }
      
      #length(gene.chosen.bc)
      
      ### compare the block and batch-specific HVGs selection
      #cat(length(gene.chosen), length(gene.chosen.bc), length(intersect(gene.chosen, gene.chosen.bc)), "\n")
      #library(VennDiagram)
      #venn.diagram(
      #  x = list(gene.chosen, gene.chosen.bc),
      #  category.names = c("batch-block" , "batch-specific"), filename = paste0(resDir, "/batch_block_specific.png"), output = TRUE)
      gene.chosen = gene.chosen.bc
      
    }
    cat("nb of HGVs : ", length(gene.chosen), "\n")
    
  }
  
  if(method == "Brenneck")
  {
    library(M3Drop)
    if(!Norm.Vars.per.batch){
      expr_matrix =  exp(logcounts(sce))
      Brennecke_HVG <- BrenneckeGetVariableGenes(
        expr_mat = expr_matrix,
        spikes = NA,
        fdr = 0.3,
        minBiolDisp = 0.
      )
      
      gene.chosen = Brennecke_HVG
    }
    
  }
  
  return(gene.chosen)
}
