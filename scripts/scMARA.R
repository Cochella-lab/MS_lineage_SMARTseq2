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
                                      id = 'manual.annot.ids', Y_name = 'RNA')
{
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(Seurat)
  library(scater)
  library(SingleCellExperiment)
  library(scran)
  source.my.script('scMARA_utility_functions.R')
  
  ll = readRDS(file = '../data/motifs_tfs/ce11_proteinCoding_genes_geneLength_transcriptLength.rds')
  motif.oc = readRDS(file = '../data/motifs_tfs/motif_oc_all_proteinCodingGenes.rds')
  # modify this drosophila motif's name, because it is not for nhr-67 
  colnames(motif.oc)[which(colnames(motif.oc) == 'nhr-67.M141')] = 'M1471_1.02_M141' 
  motif.tf = readRDS(file = '../data/motifs_tfs/motif_tf_mapping.rds')
  tfs = readxl::read_xlsx('../data/motifs_tfs/Table-S2-wTF-3.0-Fuxman-Bass-Mol-Sys-Biol-2016.xlsx', sheet = 1)
  tf.mat = readRDS(file = paste0(RdataDir, 'TFs_expression_profiles_BWM.rds')) 
  
  # mode = 'cluster.wise';
  ids = sub.obj$manual.annot.ids
  ids.uniq = unique(ids)
  ids.uniq = ids.uniq[order(nchar(ids.uniq))]
  
  # convert to SingleCellExperiment, recalculate scaling factor and normalized to fpkm
  sce = as.SingleCellExperiment(sub.obj)
  
  mm = match(rownames(sce), ll$gene.name)
  sce = sce[which(!is.na(mm)), ] # keep genes with correponding lengths
  ll = ll[mm[which(!is.na(mm))], ]
  
  sce <- logNormCounts(sce, log = FALSE, size_factors = NULL)
  Y.fpkm <- log2(calculateFPKM(sce, lengths = ll$transcript.length) + 1)
  
  remove(sce)
  
  if(mode == 'cluster.wise'){
    cat('-- averging the gene expression in clusters -- \n')
    
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
    remove(Y.fpkm)
    
    ##########################################
    # select dynamic genes for reponse Y
    # 1) with fano 
    # 2) ratio betwen daughter and mother, a lot of pairs to take care
    # 3) by lineage e.g. MSx, MSxa, MSxap, MSxapp, MSxappp, MSxapppp, MSxappppx (e.g. with gam)
    ##########################################
    select.dyn.genes.with.fano = FALSE
    if(select.dyn.genes.with.fano){
      ss = apply(Y.mat, 1, mean)
      fano = apply(Y.mat, 1, var)/ss
      plot(ss, fano, cex = 0.6);
      abline(h = c(0.5,  0.7, 1.0), col = 'blue', lwd=1.2)
      length(which(fano > 1.5))
      length(which(fano > 1.0))
      length(which(fano > 0.7))
      length(which(fano > 0.5))
      #length(which(fano > 0.3))
      
      Y.sel = Y.mat[which(fano > 1.5), ]
    }
    
    select.dyn.genes.with.pair.ratios = FALSE
    if(select.dyn.genes.with.pair.ratios){
      
      Y.mat = as.data.frame(Y.mat)
      
      rownames(Y.sel) = rownames(Y.mat)
      colnames(Y.sel) = c('MSxa', 'MSxp')
      
      hist(Y.sel, breaks = 100);abline(v = c(-1, 1))
      cutoff = 1;
      sels = apply(Y.sel, 1, function(x) sum(abs(x)> cutoff)>1)
      cat(sum(sels), ' gene were selected \n')
      Y.sel = Y.sel[sels, ]
       
    }
    
    select.dyn.genes.with.FindAllMarker.MST = FALSE
    if(select.dyn.genes.with.FindAllMarker.MST){
      run.FindAllMarkers = FALSE # it takes ~ 30 minutes
      if(run.FindAllMarkers){
        Idents(sub.obj) = sub.obj$manual.annot.ids
        markers.new <- FindAllMarkers(sub.obj, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
        saveRDS(markers.new, file = paste0(RdataDir,  'AllMarkers_MST_manual.annotation.rds'))
      }else{
        markers = readRDS(file = paste0(RdataDir,  'AllMarkers_MST_manual.annotation.rds'))
      }
      
      #top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
      #DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
      
    }
    
    pheatmap(Y.norm, cluster_rows=TRUE, show_rownames=FALSE, show_colnames = TRUE, breaks = NA,
             scale = 'row', cluster_cols=FALSE, main = paste0("dynamic genes"), 
             na_col = "white", fontsize_col = 10
    )
    
    ##########################################
    # prepare matrix A and reponse Y and run penalized.lm
    ##########################################
    require(glmnet)
    lineage = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
    lineage = c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
    
    gene.sel = markers$gene[which(!is.na(match(markers$cluster, lineage)) & markers$p_val_adj<0.001 & markers$avg_logFC>0.5)]
    gene.sel = gene.sel[which(!is.na(match(gene.sel, rownames(Y.mat))))]
    
    gene.sel_1 = gene.sel
    gene.sel_2 = gene.sel
    gene.shared = intersect(gene.sel_1, gene.sel_2)
    
    gene.sel_1 = setdiff(gene.sel_1, gene.shared)
    gene.sel_2 = setdiff(gene.sel_2, gene.shared)
    
    lineage = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
    gene.sel = gene.sel_1
    
    lineage = c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
    gene.sel = gene.sel_2
    index.sel = match(gene.sel, rownames(Y.mat))
    Y.sel = Y.mat[index.sel, match(lineage, colnames(Y.mat))]
    y = as.matrix(Y.sel)
    
    #mm = match(rownames(Y.sel), rownames(motif.oc))
    mm = match(gene.sel, rownames(motif.oc))
    y = y[!is.na(mm), ]
    x = as.matrix(motif.oc[mm[!is.na(mm)], ])
    x[which(is.na(x) == TRUE)] = 0
    
    source.my.script('scMARA_utility_functions.R')
    res = run.penelized.lm(x, y, alpha = 0, Test = TRUE)
    
  }
  
  if(mode == 'time.bin'){
    library(destiny)
    library(princurve)
        
    lineage = c('MSx', 'MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
    lineage = c('MSx', 'MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
    
    cells.sels = unique(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, lineage))])
    ll.obj = subset(sub.obj, cells = cells.sels)
    ll.obj = FindVariableFeatures(ll.obj, selection.method = "vst", nfeatures = 1000)
    ll.obj = RunPCA(ll.obj, features = VariableFeatures(ll.obj), npcs = 50, verbose = FALSE, weight.by.var = TRUE)
    
    ll.pca = ll.obj@reductions$pca@cell.embeddings[, c(1:50)]
    dm <- DiffusionMap(ll.pca, sigma = 'local', n_eigs = 5)
    #plot(dm)
    #plot(dm$DC1, dm$DC2)
    dcs = as.matrix(cbind(dm$DC1, dm$DC2))
    ll.obj[["DP"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(ll.obj))
    DimPlot(ll.obj, reduction = 'DP', group.by = 'manual.annot.ids')
    
    princurve = principal_curve(dcs, start = dcs, smoother = 'lowess', stretch = 2)
    
    plot(dcs)
    lines(princurve$s[order(princurve$lambda),], lty=1,lwd=4,col="purple",type = "l")
    
    whiskers(dcs, princurve$s)
    
    pseudot = princurve$lambda
    #plot(order(-diffmap$X[,3]), ll.obj$timingEst,  cex = 0.5)  
    plot(pseudot, ll.obj$timingEst,  cex = 0.5)  
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
      geom_point(data=data.frame(shp@tsne.rot,color=shp@ident),aes(tSNE_1,tSNE_2),size=2,color=c(rep("green",table(shp@ident)[1]),rep("yellow",table(shp@ident)[2]),rep("orange1",table(shp@ident)[3]),rep("orange4",table(shp@ident)[4]),rep("orange4",table(shp@ident)[5]))) +
      theme_classic() +  
      theme(legend.position="none",axis.title=element_text(size = rel(1)),axis.text=element_blank(), axis.ticks = element_blank(),axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+ xlab(label = "Diffusion Component 1") + ylab(label = "Diffusion Component 2") +     geom_line(size=1.5,colour="black")
    
    
  }
}

