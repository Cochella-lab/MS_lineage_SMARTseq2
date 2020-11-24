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
    
    ##########################################
    # here the pseudotime were to be inferred and lineage-dependent genes will be identified with pseudotime
    ##########################################
    library(tradeSeq)
    library(RColorBrewer)
    library(SingleCellExperiment)
    library(slingshot)
    
    
    # infer pseudotime using slingshot or diffusion map
    #lineage = c('MSx', 'MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx', 'MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
    lineage = c('MSx', 'MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
    
    cells.sels = unique(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, lineage))])
    ll.obj = subset(sub.obj, cells = cells.sels)
    ll.obj = FindVariableFeatures(ll.obj, selection.method = "vst", nfeatures = 2000)
    
    ll.obj = ScaleData(ll.obj, features = rownames(ll.obj))
    ll.obj <- RunPCA(object = ll.obj, features = VariableFeatures(ll.obj), verbose = FALSE, weight.by.var = FALSE)
    ElbowPlot(ll.obj, ndims = 50)
    
    DimPlot(ll.obj, reduction = 'pca', dims = c(1, 2),
            group.by = 'manual.annot.ids', label = TRUE) + NoLegend()
    
    nb.pcs = 10; n.neighbors = 30; min.dist = 0.3;
    ll.obj <- RunUMAP(object = ll.obj, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
    
    p1 = DimPlot(ll.obj, reduction = 'pca', group.by = 'manual.annot.ids', label = TRUE) + NoLegend()
    p2 = DimPlot(ll.obj, reduction = 'umap', group.by = 'manual.annot.ids', label = TRUE) + NoLegend()
    p1 + p2
    
    pseudotime.method = 'slingshot'
    if(pseudotime.method == 'slingshot'){
      pca <- prcomp(t(log1p(ll.obj@assays$RNA@data)[match(VariableFeatures(ll.obj), rownames(ll.obj)), ]), scale. = FALSE)
      rd1 <- pca$x[,1:2]
      plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
      #plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
      
      library(uwot)
      
      rd2 <- umap(t(log1p(ll.obj@assays$RNA@data)[match(VariableFeatures(ll.obj), rownames(ll.obj)), ]))
      colnames(rd2) <- c('UMAP1', 'UMAP2')
      
      plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
      DimPlot(ll.obj, reduction = 'umap', group.by = 'manual.annot.ids')
      
      
      
      
    }
    
    if(pseudotime.method == 'diffusion.map'){
      # the code were modified based on https://github.com/stevexniu/single-cell-ciona
      library(destiny)
      library(princurve)
      
      #ll.obj = RunPCA(ll.obj, features = VariableFeatures(ll.obj), npcs = 50, verbose = FALSE, weight.by.var = FALSE)
      
      ll.pca = ll.obj@reductions$pca@cell.embeddings[, c(1:5)]
      dm <- DiffusionMap(ll.pca, sigma = 'local', n_eigs = 5, k = 50, distance = 'euclidean')
      plot(dm)
      #plot(dm$DC1, dm$DC2)
      dcs = as.matrix(cbind(dm$DC1, dm$DC2))
      ll.obj[["DP"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(ll.obj))
      DimPlot(ll.obj, reduction = 'DP', group.by = 'manual.annot.ids')
      
      dcs = dcs[order(dcs[, 1]), ]
      princurve = principal_curve(dcs, start = dcs, smoother = 'smooth_spline', stretch = 2)
      
      plot(dcs)
      lines(princurve$s[order(princurve$lambda),], lty=1,lwd=4,col="purple",type = "l")
      
      whiskers(dcs, princurve$s)
      
      pseudot = princurve$lambda
      pseudot.scaled = (pseudot - min(pseudot))/diff(range(pseudot))
      #plot(order(-diffmap$X[,3]), ll.obj$timingEst,  cex = 0.5)  
      plot(pseudot, as.numeric(ll.obj$timingEst),  cex = 0.5)
      ll.obj$pseudotime = pseudot.scaled
      ll.obj$timingEst = as.numeric(as.character(ll.obj$timingEst))
      Idents(ll.obj) = ll.obj$manual.annot.ids
      FeatureScatter(ll.obj, feature1 = "pseudotime", feature2 = "timingEst")
      #mat.dist = as.matrix(dist(ll.pca, method = 'euclidean'))
      # Run diffusion map and return top 50 dimensions
      #set.seed(1)
      #diffmap = diffuse(mat.dist, maxdim=50)
      
      #diffmap.embedding = as.matrix(diffmap$X)
      #rownames(diffmap.embedding) = colnames(ll.obj)
      #colnames(diffmap.embedding) = paste0('diffumap_', c(1:ncol(diffmap.embedding)))
      
      #ll.obj[['diffmap']] = Seurat::CreateDimReducObject(embeddings=diffmap.embedding , key='diffmap_', assay='RNA')
      
      # Save first two diffusion map coordinators
      #shp@tsne.rot[1:2]=data.frame(shp.diff$X[,1:2],row.names = ll.obj@cell.names)
      #colnames(shp@tsne.rot)=c("tSNE_1","tSNE_2")
      # Visualize top two diffusion map components
      #tsne.pseudo(shp, do.label = F,label.cex.text = 1,name.y = "Diffusion Map Coordinator 2",name.x = "Diffusion Map Coordinator1",label.cols.use = c("green","yellow","orange1","orange4","orange4"),label.pt.size = 1.5,xlim=c(-0.1,0.05),ylim=c(-0.05,0.05))
      #legend("topleft",legend=c("12TVC","14STVC","16SHP","18SHP","20SHP"),col= c("green","yellow","orange1","orange4","orange4"),pch = 16,cex=0.5,pt.cex = 1)
      #plot(diffmap$X[, c(1:2)])
      # Fit the first two diffusion map components with principal curve
      
      
      df=data.frame(princurve$s[order(princurve$lambda),]);colnames(df) = c("x","y")
      
      ggplot(data=df,aes(x,y))+
        geom_line(size=1.5,colour="black")+
        geom_density2d(aes(colour=..level..),bins=6) + 
        scale_colour_gradient(low="darkgray",high="white",3) +
        xlim(-0.084,0.08) + ylim(-0.05,0.05) +
        geom_point(data=data.frame(shp@tsne.rot,color=shp@ident),aes(tSNE_1,tSNE_2),size=2,color=c(rep("green",table(shp@ident)[1]),
                                                                                                   rep("yellow",table(shp@ident)[2]),rep("orange1",table(shp@ident)[3]),rep("orange4",table(shp@ident)[4]),rep("orange4",table(shp@ident)[5]))) +
        theme_classic() +  
        theme(legend.position="none",axis.title=element_text(size = rel(1)),axis.text=element_blank(), axis.ticks = element_blank(),axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+ xlab(label = "Diffusion Component 1") + ylab(label = "Diffusion Component 2") +     geom_line(size=1.5,colour="black")
      
    }
    
    
    # For reproducibility
    #RNGversion("3.5.0")
    palette(brewer.pal(8, "Dark2"))
    data(countMatrix, package = "tradeSeq")
    counts <- as.matrix(countMatrix)
    rm(countMatrix)
    data(crv, package = "tradeSeq")
    data(celltype, package = "tradeSeq")
    
    set.seed(5)
    icMat <- evaluateK(counts = counts, sds = crv, k = 3:10, 
                       nGenes = 200, verbose = T)
    
    
    set.seed(7)
    pseudotime <- slingPseudotime(crv, na = FALSE)
    cellWeights <- slingCurveWeights(crv)
    sce <- fitGAM(counts = counts, pseudotime = pseudotime, cellWeights = cellWeights,
                  nknots = 6, verbose = FALSE)
    
    
  }
  
  
}
