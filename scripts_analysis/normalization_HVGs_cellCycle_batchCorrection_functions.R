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

compare.scran.seurat.sctransform = function(sce, using.HVGs = TRUE)
{
  sce.qc = sce
  sce.qc$library.size = apply(counts(sce.qc), 2, sum)
  set.seed(1234567)
  
  ms0 = as.Seurat(sce.qc, counts = 'counts', data = NULL, assay = "RNA")
  ss = sce.qc$library.size
  
  nfeatures = 3000
  n.neighbors = 30
  min.dist = 0.25
  nb.pcs = 20
  ### seurat normalization (cpm actually) 
  ms.logtransform <- NormalizeData(ms0, assay = "RNA", normalization.method = 'LogNormalize', scale.factor = mean(ss))
  ms.logtransform <- FindVariableFeatures(ms.logtransform, selection.method = "vst", nfeatures = nfeatures)
  
  ms.logtransform <- ScaleData(ms.logtransform, features = rownames(ms.logtransform))
  
  ms.logtransform <- RunPCA(object = ms.logtransform, verbose = FALSE)
  ElbowPlot(ms.logtransform)
  
  ms.logtransform <- RunUMAP(object = ms.logtransform, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  p0 = DimPlot(ms.logtransform, reduction = "umap", group.by = 'request') + ggtitle("seurat")
  
  ### cpm normalized in scater
  cpm = logNormCounts(sce.qc, size_factors = NULL, log = TRUE, pseudo_count=1, center_size_factors = FALSE)
  plot(cpm$library.size/1e6, sizeFactors(cpm), log="xy", main = 'cpm', xlab="Library size (millions)", ylab="Size factor")
  
  ms1 = as.Seurat(cpm, counts = 'counts', data = 'logcounts', assay = "RNA")
  ms1 = FindVariableFeatures(ms1, selection.method = 'vst', nfeatures = nfeatures)
  ms1 = ScaleData(ms1, features = rownames(ms1))
  ms1 = RunPCA(ms1, verbose = FALSE)
  ms1 = RunUMAP(ms1, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  p1 = DimPlot(ms1, reduction = "umap", group.by = 'request') + ggtitle('cpm')
  
  ### scran normalization
  qclust <- quickCluster(sce.qc)
  sce.norm <- computeSumFactors(sce.qc, clusters = qclust)
  sce.norm <- logNormCounts(sce.norm, log = TRUE, pseudo_count = 1)
  plot(sce.norm$library.size/1e6, sizeFactors(sce.norm), log="xy", xlab="Library size (millions)", ylab="Size factor")
  ms2 = as.Seurat(sce.norm, counts = 'counts', data = 'logcounts', assay = "RNA")
  ms2 = FindVariableFeatures(ms2, selection.method = 'vst', nfeatures = nfeatures)
  ms2 = ScaleData(ms2, features = rownames(ms1))
  ms2 = RunPCA(ms2, verbose = FALSE)
  ms2 = RunUMAP(ms2, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  p2 = DimPlot(ms2, reduction = "umap", group.by = 'request') + ggtitle('scran')
  
  
  ms3 <- SCTransform(object = ms0, variable.features.n =nfeatures) # new normalization from Seurat
  ms3 <- RunPCA(object = ms3, verbose = FALSE)
  #ElbowPlot(ms)
  ms3 <- RunUMAP(object = ms3, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  p3 = DimPlot(ms3, reduction = "umap", group.by = 'request') + ggtitle('sctransform')
  
  plot_grid(p0, p1, p2, p3, nrow = 2)
  
}

normalized.counts.using.scran.old = function(sce)
{
  set.seed(1000)
  clusters <- quickCluster(sce, min.size = 100, method="igraph")
  table(clusters)
  
  sce <- computeSumFactors(sce, clusters = clusters)
  
  ## quick check for size factors calculated by scran
  summary(sizeFactors(sce))
  plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
       xlab="Library size (millions)", ylab="Size factor", pch=16)
  
  sce <- normalize(sce, exprs_values = "counts", return_log = TRUE)
  
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

########################################################
########################################################
# Section : funcitons for cell cycle scoring and correction
# At the end, Seurat is used 
########################################################
########################################################
find.cc.markers.homologues = function()
{
  detach("package:Seurat", unload=TRUE)
  require(Seurat)
  #s.genes = c("cdk-4", "evl-18") # from GO:1901987 http://amigo1.geneontology.org/cgi-bin/amigo/term-assoc.cgi?term=GO:1902808&speciesdb=all&taxid=6239
  #g2m.genes = xx # GO:1902751
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  homologues = read.delim("/Volumes/groups/cochella/jiwang/annotations/cellCycle_genes_worm/BioMart_worm_human_homologe.txt", sep = "\t",
                          header = TRUE)
  #homologues = homologues[which(homologues$Human.orthology.confidence..0.low..1.high.==1), ]
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  s.genes = homologues$Gene.name[match(s.genes, homologues$Human.gene.name)]
  g2m.genes = homologues$Gene.name[match(g2m.genes, homologues$Human.gene.name)]
  s.genes = s.genes[which(!is.na(s.genes))]
  g2m.genes = g2m.genes[which(!is.na(g2m.genes))]
  
  c3.genes = list(s.genes = s.genes, g2m.genes = g2m.genes)
  return(c3.genes)
}

find.cc.markers.GO = function()
{
  s.genes = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_dnasynthesis.txt", 
                       sep = "\t", header = FALSE)
  gm1 = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_G2Mphase.txt", 
                   sep = "\t", header = FALSE)
  gm2 = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_Mphase.txt", 
                   sep = "\t", header = FALSE)
  gm3 = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_G2phase.txt", 
                   sep = "\t", header = FALSE)
  s.genes = unique(s.genes[, 3])
  g2m.genes = c(unique(as.character(gm1[,3])), unique(as.character(gm2[, 3])), unique(as.character(gm3[, 3])))
  
  c3.genes = list(s.genes = s.genes, g2m.genes = g2m.genes)
  return(c3.genes)
}

find.cellcycle.markers = function(list.sel = 'homologues')
{
  ##########################################
  # different strategies to find cell cycle markers of c. elegans for Seurat
  # 1st method): using homologue between human and c. elegans
  # 2nd method): find a list of 
  ##########################################
  if(list.sel == 'homologues') c3.genes = find.cc.markers.homologues() 
  if(list.sel == "curated") c3.genes = find.cc.markers.GO()
  s.genes = c3.genes$s.genes
  g2m.genes = c3.genes$g2m.genes
  
  # manually add genes from wormbook http://www.wormbook.org/chapters/www_cellcyclereguln/cellcyclereguln.html
  #s.genes = c(as.character(s.genes), c("cye-1")
  s.genes = unique(c(as.character(s.genes), c('cye-1', 'cya-1', 'evl-18')))
  g2m.genes = unique(c(as.character(g2m.genes), c('cdk-1', 'mat-1', 'mat-2', 'mat-3', 'emb-27', 'emb-30', 'mdf-1', 'san-1')))
  
  c3.genes = list(s.genes = s.genes, g2m.genes = g2m.genes)
  return(c3.genes)
  
  # test the code from https://github.com/hbc/macrae_ghazizadeh_zebrafish_heart_sc/blob/master/seurat_cluster/seurat_cluster_adapted_WT.Rmd
  # unfornately it does not work, because several pacakges can not be properly installed
  Test.query.cellCycle.markers = FALSE
  if(Test.query.cellCycle.markers){
    require(plotly)
    require(remotes)
    annot <- basejump::annotable("Danio rerio") %>% 
      dplyr::select(c(ensgene,symbol)) %>% 
      dplyr::mutate(symbol = toupper(symbol)) 
    cell_cycle_markers <- bcbioSingleCell::cellCycleMarkers[[camel("mus musculus")]] %>% 
      dplyr::mutate(symbol = toupper(symbol)) %>% dplyr::inner_join(annot,by = "symbol") %>% 
      dplyr::select(-c(ensgene.x)) %>% dplyr::rename(ensgene = ensgene.y)
    stopifnot(is.data.frame(cell_cycle_markers))
    markdownHeader("S phase markers", level = 3)
    s_genes <- cell_cycle_markers %>%
      filter(phase == "S") %>%
      pull("ensgene")
    print(s_genes)
    markdownHeader("G2/M phase markers", level = 3)
    g2m_genes <- cell_cycle_markers %>%
      filter(phase == "G2/M") %>%
      pull("ensgene")
    print(g2m_genes)
    saveData(cell_cycle_markers, s_genes, g2m_genes, dir = data_dir)
  }
  
}

cellCycle.correction = function(ms, method = "seurat")
{
  if(method == "seurat"){
    
    pdfname = paste0(resDir, "/scRNAseq_cellCycle_regression_Seurat.pdf")
    pdf(pdfname, width=12, height = 6)
    
    #library(scater)
    #library(scran)
    # install loomR from GitHub using the remotes package remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
    #library(loomR)
    #library(Seurat)
    # convert sce to seurat object (see https://satijalab.org/seurat/v3.0/conversion_vignette.html)
    #seurat = as.Seurat(sce, counts = "counts", data = "logcounts")
    #Idents(seurat) <- colnames(seurat) # quite important this assignment for cell identity
    #seurat <- FindVariableFeatures(seurat, selection.method = "vst")
    
    detach("package:scater", unload=TRUE)
    detach("package:scran", unload=TRUE)
    
    seurat = ms;
    
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(seurat), 25)
    plot1 <- VariableFeaturePlot(seurat)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    CombinePlots(plots = list(plot1, plot2))
    
    #seurat <- ScaleData(seurat, features = rownames(seurat), model.use = "linear") # standardize the data (x - mean(x))/sd(x)
    #seurat <- RunPCA(seurat, features = VariableFeatures(seurat), ndims.print = 6:10, nfeatures.print = 10)
    #DimPlot(seurat, reduction = "pca")
    # DimHeatmap(seurat, dims = c(1, 2))
    
    #source("scRNAseq_functions.R")
    c3.genes = find.cellcycle.markers(list.sel = "homologues")
    
    s.genes <- c3.genes$s.genes
    g2m.genes <- c3.genes$g2m.genes 
    s.genes = s.genes[which(!is.na(match(s.genes, rownames(seurat))))]
    g2m.genes = g2m.genes[which(!is.na(match(g2m.genes, rownames(seurat))))]
    
    seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    
    # view cell cycle scores and phase assignments
    head(seurat[[]])
    
    p00 = RidgePlot(seurat, features = c("cdk-1", "cdk-4", "cyd-1", "cye-1", "cya-1", "wee-1.3"), ncol = 2)
    plot(p00)
    
    #seurat <- RunPCA(seurat, features = VariableFeatures(seurat), verbose = TRUE)
    #DimPlot(seurat, reduction = 'pca')
    
    seurat <- RunPCA(seurat, features = c(as.character(s.genes), as.character(g2m.genes)))
    DimPlot(seurat, reduction = 'pca')
    
    # regress out the cell cycle
    seurat1 <- ScaleData(seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat))
    
    seurat1 <- RunPCA(seurat1, features = VariableFeatures(seurat1), nfeatures.print = 10)
    DimPlot(seurat1)
    
    # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
    seurat1 <- RunPCA(seurat1, features = c(s.genes, g2m.genes))
    p1 = DimPlot(seurat1)
    plot(p1)
    # regressing out the difference between the G2M and S phase scores
    seurat$CC.Difference <- seurat$S.Score - seurat$G2M.Score
    seurat2 <- ScaleData(seurat, vars.to.regress = "CC.Difference", features = rownames(seurat))
    
    # cell cycle effects strongly mitigated in PCA
    seurat2 <- RunPCA(seurat2, features = VariableFeatures(seurat2), nfeatures.print = 10)
    DimPlot(seurat2)
    
    # when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1
    # cells however, within actively proliferating cells, G2M and S phase cells group together
    seurat2 <- RunPCA(seurat2, features = c(s.genes, g2m.genes))
    p2 = DimPlot(seurat2)
    plot(p2)
    
    # save cell cycle scoring and corrected matrix
    library(scater)
    sce$S.Score = seurat$S.Score
    sce$G2M.Score = seurat$G2M.Score
    sce$Phase = seurat$Phase
    #sce$Phase.GO = seurat$old.ident
    sce$CC.Difference = seurat$CC.Difference
    
    #xx = as.data.frame(seurat@assays$RNA@scale.data); rownames(xx) = rownames(sce)
    #assay(sce, "logcounts_seurat") <- xx
    xx = as.data.frame(seurat1@assays$RNA@scale.data); rownames(xx) = rownames(sce)
    assay(sce, "logcounts_seurat_cellcycleCorrected") <- xx
    xx = as.data.frame(seurat2@assays$RNA@scale.data); rownames(xx) = rownames(sce)
    assay(sce, "logcounts_seurat_SG2MCorrected") <- xx
    
    save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2.Rdata'))
    
    dev.off()
    
    # return(sce)
    
    ##########################################
    # a peek of source code for cellCycle.scoring function from Seurat
    # https://github.com/ChristophH/in-lineage/blob/master/R/lib.R
    ##########################################
    Example.for.Seurat = FALSE
    if(Example.for.Seurat){
      library('Matrix')
      library('parallel')
      library('MASS')
      library('diffusionMap')
      library('FNN')
      library('igraph')
      library('princurve')
      library('ggplot2')
      library('inline')
      library('gplots')
      
      # for cell cycle score
      get.bg.lists <- function(goi, N, expr.bin) {
        res <- list()
        goi.bin.tab <- table(expr.bin[goi])
        for (i in 1:N) {
          res[[i]] <- unlist(lapply(names(goi.bin.tab), function(b) {
            sel <- which(expr.bin == as.numeric(b) & !(names(expr.bin) %in% goi))
            sample(names(expr.bin)[sel], goi.bin.tab[b])
          }))
        }
        return(res)
      }
      
      enr.score <- function(expr, goi, bg.lst) {
        goi.mean <- apply(expr[goi, ], 2, mean)
        bg.mean <- sapply(1:length(bg.lst), function(i) apply(expr[bg.lst[[i]], ], 2, mean))
        return((goi.mean - apply(bg.mean, 1, mean)) / apply(bg.mean, 1, sd))
      }
      
      get.cc.score <- function(cm, N=100, seed=42, 
                               s.gene.file='./annotation/s_genes.txt',
                               g2m.gene.file='./annotation/g2m_genes.txt')
      {
        set.seed(seed)
        cat('get.cc.score, ')
        cat('number of random background gene sets set to', N, '\n')
        
        min.cells <- 5
        
        cells.mols <- apply(cm, 2, sum)
        gene.cells <- apply(cm>0, 1, sum)
        cm <- cm[gene.cells >= min.cells, ]
        
        gene.mean <- apply(cm, 1, mean)
        
        breaks <- unique(quantile(log10(gene.mean), probs = seq(0,1, length.out = 50)))
        gene.bin <- cut(log10(gene.mean), breaks = breaks, labels = FALSE)
        names(gene.bin) <- rownames(cm)
        gene.bin[is.na(gene.bin)] <- 0
        
        regev.s.genes <- read.table(file=s.gene.file, header=FALSE, stringsAsFactors=FALSE)$V1
        regev.g2m.genes <- read.table(file=g2m.gene.file, header=FALSE, stringsAsFactors=FALSE)$V1
        
        goi.lst <- list('S'=rownames(cm)[!is.na(match(toupper(rownames(cm)), regev.s.genes))],
                        'G2M'=rownames(cm)[!is.na(match(toupper(rownames(cm)), regev.g2m.genes))])
        
        n <- min(40, min(sapply(goi.lst, length)))
        goi.lst <- lapply(goi.lst, function(x) x[order(gene.mean[x], decreasing = TRUE)[1:n]])
        
        bg.lst <- list('S'=get.bg.lists(goi.lst[['S']], N, gene.bin),
                       'G2M'=get.bg.lists(goi.lst[['G2M']], N, gene.bin))
        
        all.genes <- sort(unique(c(unlist(goi.lst, use.names=FALSE), unlist(bg.lst, use.names=FALSE))))
        
        expr <- log10(cm[all.genes, ]+1)
        
        s.score <- enr.score(expr, goi.lst[['S']], bg.lst[['S']])
        g2m.score <- enr.score(expr, goi.lst[['G2M']], bg.lst[['G2M']])
        
        phase <- as.numeric(g2m.score > 2 & s.score <= 2)
        phase[g2m.score <= 2 & s.score > 2] <- -1
        
        return(data.frame(score=s.score-g2m.score, s.score, g2m.score, phase))
      }
    }
    
  }
  
  ##########################################
  # scran method is based on trained classifier for mouse or human
  # so at the end it not usable
  ##########################################
  if(method == 'scran'){
    set.seed(100)
    library(scran)
    mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                    package="scran"))
    assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
    plot(assignments$score$G1, assignments$score$G2M, 
         xlab="G1 score", ylab="G2/M score", pch=16)
    sce$phases <- assignments$phases
    table(sce$phases)
  }
  
  if(method == "scLVM"){
    ##########################################
    # select the python verson to use for Rstudio
    # https://cran.r-project.org/web/packages/reticulate/vignettes/versions.html
    # still does not work at the end
    # we change the strategy: prepare the tables and run scLVM in the conda version
    ##########################################
    # system("python --version")
    # system("which python")
    # 
    # library(reticulate)
    # use_python("/Users/jiwang/anaconda3/envs/scLVM/bin/python")
    # #use_condaenv(condaenv = "scLVM", conda = "/Users/jiwang/anaconda3/condabin/conda", required = TRUE)
    # py_config()
    # 
    # system("python --version")
    # system("which python")
    # 
    # Sys.setenv(PATH = paste("/Users/jiwang/anaconda3/envs/scLVM/bin", Sys.getenv("PATH"),sep=":"))
    
    #install.packages("rPython", type = "source")
    #install.packages("/Users/jiwang/src_packages/scLVM_0.99.3.tar.gz", repos = NULL, type="source")
    
    ##########################################
    # example code from scLVM R tutorial
    # https://github.com/PMBio/scLVM/blob/master/R/tutorials/scLVM_vignette.Rmd
    ##########################################
    library(rPython)
    library(genefilter)
    library(statmod)
    require(ggplot2)
    library(gplots)
    require(DESeq2)
    library(scLVM)
    
    #limix_path = '/Users/jiwang/anaconda2/envs/scLVM/bin/python'
    #configLimix(limix_path)
    
    data(data_Tcells)
    help(data_Tcells)
    
    #dataMouse[ 1:5, 1:4 ]
    geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
      substr( rownames(dataMouse), 1, 4 ) ] )
    #2. calculate normalisation for counts
    countsMmus <- dataMouse[ which( geneTypes=="ENSM" ), ]
    countsERCC <- dataMouse[ which( geneTypes=="ERCC" ), ]
    lengthsMmus <- dataMouse[ which( geneTypes=="ENSM" ), 1 ]
    lengthsERCC <- dataMouse[ which( geneTypes=="ERCC" ), 1 ]
    sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
    sfMmus <- sfERCC #also use ERCC size factor for endogenous genes
    #normalise read counts
    nCountsERCC <- t( t(countsERCC) / sfERCC )
    nCountsMmus <- t( t(countsMmus) / sfERCC )
    
    countsMmus = counts(sce)
    sfERCC = estimateSizeFactorsForMatrix(countsMmus)
    sfMmus <- sfERCC
    nCountsMmus = t( t(countsMmus) / sfERCC )
    #use spike in to find tehcnical noise. 
    # If no spike-ins are available, we can also use the endogenous read counts for fitting the mean-CV2 relation using a log-linear fit in the log-space.
    # Alternatively, we can fit the mean-variance relationship in the log-space using local 2nd order polynomial regression (loess).
    #techNoise = fitTechnicalNoise(nCountsMmus,nCountsERCC=nCountsERCC, fit_type = 'counts')  
    techNoiseLogFit = fitTechnicalNoise(nCountsMmus, fit_type = 'log', use_ERCC = FALSE, plot=TRUE) 
    #techNoiseLogVarFit = fitTechnicalNoise(nCountsMmus, fit_type = 'logvar', use_ERCC = FALSE, plot=TRUE) 
    
    #call variable genes
    #is_het = getVariableGenes(nCountsMmus, techNoiseLogFit$fit, method = "fit", 
    #                          threshold = 0.1, fit_type="log",sfEndo=sfMmus, sfERCC=sfERCC)
    #table(is_het)
    
    #we an also do this for the other fits
    is_hetLog = getVariableGenes(nCountsMmus, techNoiseLogFit$fit, plot=TRUE)
    table(is_hetLog)
    #is_hetLogVar = getVariableGenes(nCountsMmus, techNoiseLogVarFit$fit, plot=TRUE)
    #table(is_hetLogVar)
    
    #get cell cycle genes from GO
    cc.genes = read.delim("../../../../annotations/cellCycle_genes_worm/GO_0007049_genes.txt", 
                          sep = "\t", header = FALSE)
    cc.genes = unique(cc.genes[,3])
    
    #rename a few variables
    Y = t(log10(nCountsMmus+1)) #normalised trandformed read counts
    genes_het_bool = as.vector(is_hetLog) #variable genes
    #genes_het_bool[]
    geneID = rownames(nCountsMmus) #gene IDs
    tech_noise = as.vector(techNoiseLogFit$techNoiseLog) #technical noise
    ens_ids_cc <- cc.genes
    
    index.cc = match(cc.genes, geneID)
    index.cc = index.cc[which(!is.na(index.cc))]
    ##########################################
    # can not proceed anymore and save tables for python in conda
    ##########################################
    #construct and initialize new scLVM object
    sclvm = new("scLVM")
    sclvm = init(sclvm,Y=Y,tech_noise = tech_noise)
    
    # CellCycleARD = fitFactor(sclvm,geneSet = ens_ids_cc, k=20,use_ard = TRUE)
    
    write.table(Y, file = paste0(tabDir, "gene_expression_matrx_4scLVM.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(tech_noise, file = paste0(tabDir, "tech_noise_4scLVM.txt"), sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE )
    write.table(geneID, file =paste0(tabDir, "geneNames_4scLVM.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE )
    write.table(which(genes_het_bool==TRUE), file =paste0(tabDir, "index_hetgenes_4scLVM.txt"), sep = "\t", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    write.table(index.cc, file =paste0(tabDir, "index_ccgenes_4scLVM.txt"), sep = "\t", 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  if(method == "ccRemover"){
    # see examplel from https://cran.r-project.org/web/packages/ccRemover/vignettes/ccRemover_tutorial.html
    # this method is pretty slow and probably risky (to test), because simply PCA was done and 
    # then compare cell cycle genes were compares with control genes and significantly different PCs were removed
    require(ccRemover)
    t.cell_data = logcounts(sce)
    head(t.cell_data[,1:5])
    
    summary(apply(t.cell_data,1, mean))
    mean_gene_exp <- rowMeans(t.cell_data)
    t_cell_data_cen <- t.cell_data - mean_gene_exp
    summary(apply(t_cell_data_cen,1,mean))
    
    gene_names <- rownames(t_cell_data_cen)
    # cell_cycle_gene_indices <- gene_indexer(gene_names, species = "mouse", 
    #                                         name_type = "symbols" )
    # length(cell_cycle_gene_indices)
    # if_cc <- rep(FALSE,nrow(t_cell_data_cen)) 
    # if_cc[cell_cycle_gene_indices] <- TRUE
    # summary(if_cc)
    
    cc.genes = read.delim("../../../../annotations/cellCycle_genes_worm/GO_0007049_genes.txt", header = FALSE)
    cc.genes = unique(cc.genes$V3)
    cc.index = match(cc.genes, gene_names)
    cc.index = cc.index[which(!is.na(cc.index))]
    
    if_cc <- rep(FALSE,nrow(t_cell_data_cen))
    if_cc[cc.index] <- TRUE
    summary(if_cc)
    
    dat <- list(x=t_cell_data_cen, if_cc=if_cc)
    xhat <- ccRemover(dat, bar=TRUE, max_it = 6, nboot = 100)
    
    xhat <- xhat + mean_gene_exp
    
    pca1 = prcomp(t(t.cell_data[if_cc,]), scale. = TRUE)
    pca2 = prcomp(t(xhat[if_cc,]), scale. = TRUE)
    par(mfrow = c(1, 2))
    plot(pca1$x[, c(2:3)])
    plot(pca2$x[, c(1:2)])
    
  }
  
}






