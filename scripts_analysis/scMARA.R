##########################################################################
##########################################################################
# Project: MS lineage embryogenesis in C elegans
# Script purpose: predicting TF regulation activity 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct  7 17:31:46 2020
##########################################################################
##########################################################################
predict.TF.MARA.for.scdata = function(sub.obj, mode = 'cluster.wise', id = 'manual.annot.ids', Y_name = 'RNA')
{
  library("pheatmap")
  library("RColorBrewer")
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
    
        
  }else{
    Y.mat = Y.fpkm
  }
  
}



########################################################
########################################################
# Section : package installation
# 
########################################################
########################################################
install.magic = FALSE
if(install.magic){
  install.packages("Rmagic")
  system('which python')
  system('python --version')
  system('pip install --user magic-impute') 
  # at the same time the magic-impute was installed in default python 
  # /usr/local/bin/python
  
}

##########################################
# test MAGIC and it works for the example
# original code from https://github.com/KrishnaswamyLab/MAGIC
##########################################
Test.Rmagic = FALSE
if(Test.Rmagic){
  library(Rmagic)
  library(ggplot2)
  data(magic_testdata)
  
  ss = apply(magic_testdata, 2, sum)
  magic_testdata = magic_testdata[, ss>0]
  MAGIC_data <- Rmagic::magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"), verbose = 1)
  
  ggplot(MAGIC_data) +
    geom_point(aes(x=VIM, y=CDH1, color=ZEB1))
  
}

Test.phateR = FALSE
if(Test.phateR){
  ##########################################
  # test PHATE from https://github.com/KrishnaswamyLab/phateR
  # FAQ
  # 
  # Should genes (features) by rows or columns?
  #   
  #   To be consistent with common dimensionality reductions such as PCA (stats::prcomp) and t-SNE (Rtsne::Rtsne), we require that cells (observations) be rows and genes (features) be columns of your input data.
  # 
  # Can I run PHATE with Seurat?
  #   
  #   PHATE was removed from Seurat in version 3. You can install a version of Seurat with RunPHATE included by following the instructions at https://github.com/satijalab/seurat/pull/1172#issuecomment-564782167.
  # 
  # I have installed PHATE in Python, but phateR says it is not installed!
  #   
  #   Check your reticulate::py_discover_config("phate") and compare it to the version of Python in which you installed PHATE (run which python and which pip in a terminal.) Chances are reticulate can’t find the right version of Python; you can fix this by adding the following line to your ~/.Renviron:
  #   
  #   PATH=/path/to/my/python
  # 
  # You can read more about Renviron at https://CRAN.R-project.org/package=startup/vignettes/startup-intro.html.
  # Help
  # 
  # Please let us know of any issues at the GitHub repository. If you have any questions or require assistance using PHATE, please read the documentation at https://CRAN.R-project.org/package=phateR/phateR.pdf or by running help(phateR::phate) or contact us at https://krishnaswamylab.org/get-help.
  ##########################################
  library(phateR)
  #> Loading required package: Matrix
  data(tree.data)
  plot(prcomp(tree.data$data)$x, col=tree.data$branches)
  # runs phate
  tree.phate <- phate(tree.data$data)
  summary(tree.phate)
  #> PHATE embedding
  #> k = 5, alpha = 40, t = auto
  #> Data: (3000, 100)
  #> Embedding: (3000, 2)
  
  # plot embedding
  palette(rainbow(10))
  plot(tree.phate, col = tree.data$branches)
  
}
