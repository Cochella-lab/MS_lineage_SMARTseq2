##########################################################################
##########################################################################
# Project: MS lineage embryogenesis in C elegans
# Script purpose: predicting TF regulation activity 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct  7 17:31:46 2020
##########################################################################
##########################################################################
predict.TF.MARA.for.scdata = function(sub.obj, mode = c('cluster.based', 'time.bin', 'cell.based'), 
                                      id = 'manual.annot.ids', process.motif.oc = FALSE, 
                                      Y_name = 'RNA')
{
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(Seurat)
  library(scater)
  library(SingleCellExperiment)
  library(scran)
  require(glmnet)
  source.my.script('scMARA_utility_functions.R')
  
  ##########################################
  # import processed tables : motif-tf mapping, motif.oc matrix,
  # known TF annotation
  ##########################################
  ll = readRDS(file = '../data/motifs_tfs/ce11_proteinCoding_genes_geneLength_transcriptLength.rds') # gene length and transript length
  #motif.oc = readRDS(file = '../data/motifs_tfs/motif_oc_all_proteinCodingGenes.rds')
  motif.oc = readRDS(file = '../data/motifs_tfs/motif_oc_scATACpeaks_all_proteinCodingGenes.rds')
  
  # modify this drosophila motif's name, because it is not for nhr-67 
  colnames(motif.oc)[which(colnames(motif.oc) == 'nhr-67.M141')] = 'M1471_1.02_M141' 
  motif.tf = readRDS(file = '../data/motifs_tfs/motif_tf_mapping.rds')
  tfs = readxl::read_xlsx('../data/motifs_tfs/Table-S2-wTF-3.0-Fuxman-Bass-Mol-Sys-Biol-2016.xlsx', sheet = 1)
  tf.mat = readRDS(file = paste0(RdataDir, 'TFs_expression_profiles_BWM.rds')) 
  
  ##########################################
  # normalized the single cell gene expression matrix with gene length 
  ##########################################
  ids = sub.obj$manual.annot.ids
  ids.uniq = unique(ids)
  ids.uniq = ids.uniq[order(ids.uniq)]
  ids.uniq = ids.uniq[order(nchar(ids.uniq))]
  
  ids.uniq = ids.uniq[grep('mixture_terminal_', ids.uniq, invert = TRUE)]
  
  # convert to SingleCellExperiment, recalculate scaling factor and normalized to fpkm
  sce = as.SingleCellExperiment(sub.obj)
  qclust <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters = qclust)
  sce <- logNormCounts(sce, log = TRUE, pseudo_count = 1)
  
  mm = match(rownames(sce), ll$gene.name)
  sce = sce[which(!is.na(mm)), ] # keep genes with correponding lengths
  transcript.length  = ll$transcript.length[mm[which(!is.na(mm) == TRUE)]]
  
  sce <- logNormCounts(sce, log = FALSE, size_factors = NULL)
  Y.fpkm <- log2(calculateFPKM(sce, lengths = transcript.length) + 1)
  
  remove(sce)
  
  ########################################################
  ########################################################
  # with manually annoted ids, there are seveval options to consider trajectories
  # 1) collapse cells in the same cell states to have pseudo-bulk
  # 2) pseudotime for trajectories, better resolution than pseud-bulk
  # 3) retain single cell 
  ########################################################
  ########################################################
  trajectory.resolution = 'pseudotime'
  
  if(trajectory.resolution == 'cluster.based'){
    ##########################################
    # average cells with the same ids (i.e. cluster-based motif activity )
    ##########################################
    Y.mat = matrix(NA, nrow = nrow(Y.fpkm), ncol = length(ids.uniq))
    colnames(Y.mat) = ids.uniq
    rownames(Y.mat) = rownames(Y.fpkm)
    for(n in 1:length(ids.uniq))
    {
      cat(ids.uniq[n], '\n')
      jj = which(ids == ids.uniq[n])
      if(length(jj) == 1) Y.mat[, n] = Y.fpkm[,jj]
      if(length(jj) > 1) Y.mat[, n] = apply(Y.fpkm[,jj], 1, mean)
    }
    
    # process.detected.tf.expression.profiles(Y.mat)
    remove(Y.fpkm) # free memory
    
    ##########################################
    # determine lineage-specific signatures 
    # i.e. selected gene sets (lineage-wide, specific, or restricted)
    # here we used markers from scRNA analysis
    ##########################################
    # gene.sels = define.modules.for.lineags(sub.obj, Y.fpkm, lineage = lineage)
    
    markers = readRDS(file = paste0(RdataDir,  'AllMarkers_MST_manual.annotation.rds'))
    markers.sels = markers[which(markers$p_val<10^-5 & markers$avg_logFC > 1), ]
    print(table(markers.sels$cluster))
    
    ids.groups = list(ids.uniq[which(nchar(ids.uniq) <=5)], 
                      ids.uniq[which(nchar(ids.uniq) == 6)],
                      ids.uniq[which(nchar(ids.uniq) == 7)], 
                      ids.uniq[which(nchar(ids.uniq) > 7)]) 
    
    ids.groups = as.list(ids.uniq)
    
    ids.groups = list(c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx'), 
                      c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp'))
    ids.groups = list(c('MSx', 'MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx', 
                        'MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp'))
    
    lineage = c('MSx', 'MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx', 
                'MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
    
    #lineage = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
    #lineage = c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
    # lineage = setdiff(ids.uniq, c("mixture_terminal_1", "mixture_terminal_2"))
    print(lineage)
    
    gene.sels = markers.sels[!is.na(match(markers.sels$cluster, lineage)), ]
    gene.sels = gene.sels[!is.na(match(gene.sels$gene, rownames(Y.mat))), ]
    #print(table(gene.sels$cluster))
    gene.sels = unique(gene.sels$gene)
    cat('nb of genes to use : ', length(gene.sels), '\n')
    
    #gene.sels = unique(markers$gene[which(!is.na(match(markers$cluster, lineage)) & markers$p_val_adj<10^-5 & markers$avg_logFC > 0.7)])
    #gene.sels = gene.sels[which(!is.na(match(gene.sels, rownames(Y.mat))))]
    ##########################################
    # prepare matrix A and reponse Y and run penalized.lm
    ##########################################
    index.sel = match(gene.sels, rownames(Y.mat))
    Y.sel = Y.mat[index.sel, match(lineage, colnames(Y.mat))]
    y = as.matrix(Y.sel)
    
    #mm = match(rownames(Y.sel), rownames(motif.oc))
    mm = match(gene.sels, rownames(motif.oc))
    y = y[!is.na(mm), ]
    x = as.matrix(motif.oc[mm[!is.na(mm)], ])
    x[which(is.na(x) == TRUE)] = 0
    
    source.my.script('scMARA_utility_functions.R')
    res = run.penelized.lm(x, y, alpha = 0, standardize = TRUE, intercept = TRUE, use.lambda.min = FALSE, 
                           Test = FALSE)
    
    print(res[grep('pha-4.mus|hnd-1..Tcf3|nhr-67.mus.M226|nhr-67.dm.M141|hlh-1.M175|unc-120.dm', rownames(res)),])
    
    
    #colnames(res) = ids.uniq
    print(res[grep('pha-4.mus|hnd-1..Tcf3|nhr-67.mus.M226|nhr-67.dm.M141|hlh-1.M175|unc-120.dm|tbx-|ceh-51', rownames(res)),])
    
    ss = apply(res, 1, function(x) length(which(abs(x)>1.2)))
    length(which(ss>=1))
    
    yy = res[which(ss>0), ] 
    yy[which(abs(yy)>2.5)] = 2.5
    
    
    pdfname = paste0(resDir, "/MARA_prediction_lineages_two_selected.pdf")
    pdf(pdfname, width=18, height = 16)
    par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    pheatmap(yy[grep('pha-4.mus|hnd-1..Tcf3|nhr-67.mus.M226|nhr-67.dm.M141|hlh-1.M175|unc-120.dm', rownames(yy)), ], 
             cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
             scale = 'none', cluster_cols=FALSE, main = paste0("MARA prediction for controls"), 
             na_col = "white", fontsize_col = 12) 
    
    pheatmap(yy, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
             scale = 'none', cluster_cols=FALSE, main = paste0("MARA prediction"), 
             na_col = "white", fontsize_col = 12) 
    
    dev.off()
    
  }
  
  if(trajectory.resolution == 'pseudotime'){
    
    library(tradeSeq)
    library(RColorBrewer)
    library(SingleCellExperiment)
    library(slingshot)
    
    # lineage or trajectory to consider
    lineage.list = list(c('MSx', 'MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx'),
                    c('MSx', 'MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
                    )
    
    pseudotime.method = 'diffusion.map'
    cat(length(lineage.list), ' lineages or trajectories to consider \n')
    
    # prepare input matrix for tradeSeq
    ids.sels = c()
    for(m in 1:length(lineage.list)) ids.sels = unique(c(ids.sels, lineage.list[[m]]))
    cells.sels = unique(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, ids.sels))])
    lineages.obj = subset(sub.obj, cells = cells.sels)
    
    counts.sel = as.matrix(lineages.obj@assays$RNA@counts)
    #crv <- newSlingshotDataSet(reducedDim = lineages.obj@reductions$umap@cell.embeddings[, c(1:2)], 
    #                           clusterLabels = lineages.obj$manual.annot.ids)
    
    pseudotime <- matrix(0, nrow = ncol(lineages.obj), ncol = length(lineage.list))
    rownames(pseudotime) = colnames(lineages.obj);
    colnames(pseudotime) = paste0('curve', c('MSxa', 'MSxp'))
    cellWeights <- pseudotime
    
    ##########################################
    # # infer pseudotime using slingshot or diffusion map
    # here the pseudotime were to be inferred and lineage-dependent genes will be identified with pseudotime
    # two options for pseudotime inferring: slingshot or diffusion map + princi_curve
    # diffusion map + princi_curve method is one trajectory one time
    ##########################################
    
    # the code were modified based on https://github.com/stevexniu/single-cell-ciona
    if(pseudotime.method == 'diffusion.map'){
      
      library(destiny)
      library(princurve)
      
      # save the pseudotime estimation 
      pdfname = paste0(resDir, "/pseudotime_estimation_v1.pdf")
      pdf(pdfname, width=12, height = 10)
      par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
      
      # loop over the trajectories
      for(n in 1:length(lineage.list)){
        
        #n = 2
        cat('lineage ', n, ': \n')
        lineage = lineage.list[[n]]
        print(lineage)
        
        ll.obj = subset(sub.obj, cells = unique(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, lineage))]))
        ll.obj = FindVariableFeatures(ll.obj, selection.method = "vst", nfeatures = 2000)
        
        ll.obj = ScaleData(ll.obj, features = rownames(ll.obj))
        ll.obj <- RunPCA(object = ll.obj, features = VariableFeatures(ll.obj), verbose = FALSE, weight.by.var = FALSE)
        ElbowPlot(ll.obj, ndims = 30)
        
        nb.pcs = 10; n.neighbors = 30; min.dist = 0.3;
        ll.obj <- RunUMAP(object = ll.obj, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
        
        p1 = DimPlot(ll.obj, reduction = 'pca', group.by = 'manual.annot.ids', label = TRUE) + NoLegend()
        p2 = DimPlot(ll.obj, reduction = 'umap', group.by = 'manual.annot.ids', label = TRUE) + NoLegend()
        p1 + p2
        
        ll.pca = ll.obj@reductions$pca@cell.embeddings[, c(1:5)]
        dm <- DiffusionMap(ll.pca, sigma = 'local', n_eigs = 5, k = 100, distance = 'euclidean')
        plot(dm)
        
        #plot(dm$DC1, dm$DC2)
        dcs = as.matrix(cbind(dm$DC1, dm$DC2))
        ll.obj[["DP"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(ll.obj))
        p1 = DimPlot(ll.obj, reduction = 'DP', group.by = 'manual.annot.ids')
        plot(p1)
        
        dcs = dcs[order(dcs[, 1]), ]
        princurve = principal_curve(dcs, start = dcs, smoother = 'smooth_spline', stretch = 2)
        
        plot(dcs)
        lines(princurve$s[order(princurve$lambda),], lty=1,lwd=4,col="purple",type = "l")
        whiskers(dcs, princurve$s)
        
        pseudotime.scaling = function(X) {
          return((X - min(X))/diff(range(X)))
        }
        
        pseudot = pseudotime.scaling(princurve$lambda)
        pseudot = pseudot[match(colnames(ll.obj), names(pseudot))] # match back with cell names
        
        #plot(pseudot, as.numeric(as.character(ll.obj$timingEst)),  cex = 0.5)
        Idents(ll.obj) = ll.obj$manual.annot.ids
        ll.obj$pseudotime = pseudot
        ll.obj$timingEst = as.numeric(as.character(ll.obj$timingEst))
        
        p2 = FeatureScatter(ll.obj, feature1 = "pseudotime", feature2 = "timingEst")
        plot(p2)
        # save the pseudotime in the initialized matrix
        mm = match(colnames(ll.obj), rownames(pseudotime))
        pseudotime[mm, n] = pseudot
        cellWeights[mm, n] = 1
        
      }
      
      dev.off()
      
    }
    
    # save(counts.sel, pseudotime, cellWeights, file = paste0(RdataDir, 'input_Matrix_for_tradeSeq.Rdata'))
    
    ##########################################
    # idnetify trajectory-associated genes using tradeSeq
    # the orignial code was from https://statomics.github.io/tradeSeq/articles/tradeSeq.html
    # updated version of analysis for multiple conditions were found 
    # https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
    ##########################################
    require(tictoc)
    load(file = paste0(RdataDir, 'input_Matrix_for_tradeSeq.Rdata'))
    
    palette(brewer.pal(8, "Dark2"))
    #data(countMatrix, package = "tradeSeq")
    #counts <- as.matrix(countMatrix)
    #rm(countMatrix)
    #data(crv, package = "tradeSeq")
    #data(celltype, package = "tradeSeq")
    set.seed(5)
    tic()
    icMat <- evaluateK(counts = counts.sel, k = 3:7,
                       pseudotime = pseudotime, cellWeights = cellWeights,
                       nGenes = 200, verbose = T)
    
    toc()
    
    
    # subsetting the genes of intest rather than all genes (too slow)
    nb.knots = 4;
    
    ss = apply(counts.sel, 1, function(x) length(which(x>10)))
    genes.sel = which(ss>50)
    length(genes.sel)
    
    BPPARAM <- BiocParallel::bpparam()
    BPPARAM # lists current options
    BPPARAM$workers <- 4 # use 2 cores
    
    tic()
    set.seed(7)
    #pseudotime <- slingPseudotime(crv, na = FALSE)
    #cellWeights <- slingCurveWeights(crv)
    sce <- fitGAM(counts = counts.sel, pseudotime = pseudotime, cellWeights = cellWeights,
                  nknots = nb.knots, verbose = TRUE, parallel=TRUE, BPPARAM = BPPARAM, genes = genes.sel)
    toc()
    
    save(sce, file = paste0(RdataDir, 'fitGAM_output_tradeSeq_v3.Rdata'))
    
    mean(rowData(sce)$tradeSeq$converged)
    
    # Assess DE along pseudotime, see more details in 
    # https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
    assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))
    
    Msxa.genes <-  rownames(assocRes)[
      which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.01)
      ]
    Msxp.genes <-  rownames(assocRes)[
      which(p.adjust(assocRes$pvalue_2, "fdr") <= 0.01)
      ]
    length(Msxa.genes)
    length(Msxp.genes)
    
    library(UpSetR)
    UpSetR::upset(fromList(list(Msxa = Msxa.genes, Msxp = Msxp.genes)))
    
    yhatSmooth <- predictSmooth(sce, gene = Msxa.genes, nPoints = 50, tidy = FALSE)
    heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           show_colnames = FALSE)
    
    shared.genes = intersect(Msxa.genes, Msxp.genes)
    Msxa.specific.genes = setdiff(Msxa.genes, shared.genes)
    Msxp.specific.genes = setdiff(Msxp.genes, shared.genes)
    
    gene.example = 'tbx-35'
    plotSmoothers(sce, assays(sce)$counts, gene = gene.example, alpha = 1, border = TRUE) + ggtitle(gene.example)
    
    pdfname = paste0(resDir, "/lineage_dependant_genes_MSxa_MSxp_shared.pdf")
    pdf(pdfname, width=16, height = 10)
    par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    for(g in shared.genes){
      cat(g, '\n')
      kk = which(tfs$`Public name` == g)
      gtitle = g
      if(length(kk) >0 ) gtitle = paste0(gtitle, '-- TF') 
      p1 = plotSmoothers(sce, assays(sce)$counts, gene = g, alpha = 1, border = TRUE,
                         nPoints = 100, size = 1) + ggtitle(gtitle)
      plot(p1)
    }
    
    dev.off()
    
    pdfname = paste0(resDir, "/lineage_dependant_genes_MSxa_specific.pdf")
    pdf(pdfname, width=16, height = 10)
    par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    for(g in Msxa.specific.genes){
      cat(g, '\n')
      kk = which(tfs$`Public name` == g)
      gtitle = g
      if(length(kk) >0 ) gtitle = paste0(gtitle, '-- TF') 
      p1 = plotSmoothers(sce, assays(sce)$counts, gene = g, alpha = 1, border = TRUE,
                         nPoints = 100, size = 1) + ggtitle(gtitle)
      plot(p1)
    }
    
    dev.off()
    
    pdfname = paste0(resDir, "/lineage_dependant_genes_MSxp_specific.pdf")
    pdf(pdfname, width=16, height = 10)
    par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    for(g in Msxp.specific.genes){
      cat(g, '\n')
      kk = which(tfs$`Public name` == g)
      gtitle = g
      if(length(kk) >0 ) gtitle = paste0(gtitle, '-- TF') 
      p1 = plotSmoothers(sce, assays(sce)$counts, gene = g, alpha = 1, border = TRUE,
                         nPoints = 100, size = 1) + ggtitle(gtitle)
      plot(p1)
    }
    
    dev.off()
    
    
  }
  
  
}
