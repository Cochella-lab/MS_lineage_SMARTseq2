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
  
  ll = readRDS(file = '../data/motifs_tfs/ce11_proteinCoding_genes_geneLength_transcriptLength.rds')
  motif.oc = readRDS(file = '../data/motifs_tfs/motif_oc_all_proteinCodingGenes.rds')
  motif.tf = readRDS(file = '../data/motifs_tfs/motif_tf_mapping.rds')
  tfs = readxl::read_xlsx('../data/motifs_tfs/Table-S2-wTF-3.0-Fuxman-Bass-Mol-Sys-Biol-2016.xlsx', sheet = 1)
  
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
  
  if(mode == 'cluster.wise'){
    cat('-- averging the gene expression in clusters -- \n')
    
    Y.mat = matrix(NA, nrow = nrow(Y.fpkm), ncol = length(ids.uniq))
    colnames(Y.mat) = ids.uniq
    rownames(Y.mat) = rownames(Y.fpkm)
    for(n in 1:length(ids.uniq))
    {
      jj = which(ids == ids.uniq[n])
      if(length(jj) == 1) Y.mat[, n] = Y.fpkm[,jj]
      if(length(jj) > 1) Y.mat[, n] = apply(Y.fpkm[,jj], 1, mean)
    }
    
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
      
      Y.sel = cbind(Y.mat$MSxa - Y.mat$MSx,
                    Y.mat$MSxp - Y.mat$MSx
                    #Y.mat$MSxap - Y.mat$MSxa, 
                    #Y.mat$MSxapp - Y.mat$MSxap, 
                    #Y.mat$MSxappp - Y.mat$MSxapp
                    )
      rownames(Y.sel) = rownames(Y.mat)
      colnames(Y.sel) = c('MSxa', 'MSxp')
      
      hist(Y.sel, breaks = 100);abline(v = c(-1, 1))
      cutoff = 1;
      sels = apply(Y.sel, 1, function(x) sum(abs(x)> cutoff)>1)
      cat(sum(sels), ' gene were selected \n')
      Y.sel = Y.sel[sels, ]
       
    }
    cal_z_score <- function(x){ (x - mean(x)) / sd(x)}
    Y.norm <- t(apply(Y.sel, 1, cal_z_score))
    #cols = c(colorRampPalette((brewer.pal(n = 7, name="RdYlBu")))(100))
    pheatmap(Y.norm, cluster_rows=TRUE, 
             show_rownames=FALSE, show_colnames = TRUE, breaks = NA,
             scale = 'none',
             cluster_cols=FALSE, 
             main = paste0("gene dynamics"), 
             na_col = "white",
             #color = cols, 
             #annotation_col = my_sample_col,
             #gaps_row = c(1:nrow(map)-1),
             fontsize_col = 10,
             height = 8,
             width = 30
    )
    
    ##########################################
    # prepare matrix A and reponse Y and run elastic-net
    ##########################################
    require(glmnet)
    y = as.matrix(Y.sel)
    #y = y[]
    mm = match(rownames(Y.sel), rownames(motif.oc))
    y = y[!is.na(mm), ]
    x = as.matrix(motif.oc[mm[!is.na(mm)], ])
    x[which(is.na(x) == TRUE)] = 0
    
    cat('nb of motifs which were not found among all considered genes ', length(which(apply(x, 2, sum) == 0)), '\n')
    cat('nb of genes without motifs ', length(which(apply(x, 1, sum) == 0)), '\n')
    
    #y = Y.sel[, c(1, 2)] 
    #y = Y.sel[, c(1:10)]
    
    alpha = 0.5
    #binary = 1;  binary.matrix = binary== 0
    intercept=FALSE
    ### standardize matrix of motif occurrence makes more sense because the absolute number of motif occurrence is not precise.
    standardize=TRUE
    standardize.response=FALSE
    
    #if(binary.matrix){x = x >0; standardize=FALSE}
    ### use Cross-validation to select tuning paprameter
    cv.fit=cv.glmnet(x, y, family='mgaussian', grouped=FALSE, 
                     alpha=alpha, nlambda=100, standardize=standardize, 
                     standardize.response=standardize.response, intercept=intercept, relax = TRUE)
    plot(cv.fit)
    #cv.fit$lambda
    
    #optimal = which(cv.fit$lambda==cv.fit$lambda.min)
    s.optimal = cv.fit$lambda.1se
    fit=glmnet(x,y,alpha=alpha, lambda=cv.fit$lambda,family='mgaussian', 
               standardize=standardize, standardize.response=standardize.response, intercept=intercept)
    #plot(fit, xvar = "lambda", label = TRUE, type.coef = "coef")
    ## collect result from the elastic-net
    #colnames(x)[which(fit$beta[[1]][,optimal]!=0)]
    #colnames(x)[which(fit$beta[[2]][,optimal]!=0)]
    keep = as.data.frame(coef.glmnet(fit, s = s.optimal))
    keep = keep[-1, ] # remove intecept
    colnames(keep) = names(fit$beta)
    ss = apply(keep, 1, function(x) all(x==0))
    keep = keep[!ss, ] 
    
    #keep[which(abs(keep)<0.05)] = 0
    #keep = data.frame(keep, stringsAsFactors = FALSE)
    head(rownames(keep)[order(-abs(keep$MSxp))], 10)
    
  }else{
    Y.mat = Y.fpkm
  }
  
}


########################################################
########################################################
# Section : utility functions for sctf_MARA
# 
########################################################
########################################################
##########################################
# process worm gene tss to have unique tss for each gene. 
# to do this, we just pick the tss furthest from the gene start so that the promoter can cover as much as regulatory elements
# in addition, save the gene length, transcript length for the scRNA-seq length normalization
##########################################
process.worm.gene.tss = function()
{
    
  rm(list=ls())
  setwd('/Volumes/groups/cochella/jiwang/annotations')
  
  tss = read.table('ce11_tss.bed', header = FALSE, sep = '\t')
  load('BioMart_WBcel235.Rdata')
  annot = annot[which(annot$Gene.type == 'protein_coding'), ]
  
  # filter non-protein-coding genes
  mm = match(tss$V4, annot$Gene.stable.ID)
  tss = tss[which(!is.na(mm)), ]
  
  gg.uniq = unique(tss$V4)
  keep = rep(NA, length(gg.uniq))
  names(keep) = gg.uniq
  gg.counts = table(tss$V4)
  
  gg.with.singleTss = names(gg.counts)[which(gg.counts == 1)]
  #gg.with.multiTss = names(gg.counts)[which(gg.counts > 1)]
  keep[match(gg.with.singleTss, names(keep))] = match(gg.with.singleTss, tss$V4) 
  
  nn = which(is.na(keep))
  for(n in nn)
  {
     kk = which(tss$V4 == names(keep)[n])
     if(length(kk) == 1){
        cat('Error \n')
     }else{
       cat(n, '--', as.character(names(keep)[n]), '\n')
       #jj = which(annot$)
       if(unique(tss$V6[kk]) == '+'){
         keep[n] = kk[which(tss$V2[kk] == max(tss$V2[kk]))][1]
       }else{
         keep[n] = kk[which(tss$V2[kk] == min(tss$V2[kk]))][1]
       }
     }
  }
  
  tss = tss[keep, ]
  #write.table(tss, file = 'ce11_tss_curated_singleTss_perGene_proteinCoding.bed', sep = '\t', col.names = FALSE, row.names = FALSE,
  #            quote = FALSE)
  
  
  ## save the gene length and transcript length (averge of isoform lengths)
  aa = data.frame(annot$Gene.stable.ID, annot$Gene.Start..bp., annot$Gene.End..bp., annot$Transcript.length..including.UTRs.and.CDS., 
                  annot$Gene.name, annot$Gene.type, stringsAsFactors = FALSE)
  
  aa$gene.length = abs(aa$annot.Gene.End..bp. - aa$annot.Gene.Start..bp.)
  
  gg.uniq = unique(aa$annot.Gene.stable.ID)
  keep = aa[match(gg.uniq, aa$annot.Gene.stable.ID), ]
  colnames(keep) = c('wormbase.id', 'gene.start', 'gene.end', 'transcript.length', 'gene.name', 'gene.type', 'gene.length')
  
  for(n in 1:nrow(keep))
  {
    # n = 1
    jj = which(aa$annot.Gene.stable.ID == keep$wormbase.id[n])
    if(length(jj) > 1){
      cat(n, '--', as.character(keep$wormbase.id[n]), '\n')
      keep$transcript.length[n] = as.integer(median(aa$annot.Transcript.length..including.UTRs.and.CDS.[jj]))
    }
  }
  #saveRDS(keep, file = 'ce11_proteinCoding_genes_geneLength_transcriptLength.rds')
}

##########################################
# after running FIMO, make motif occurrency matrix  
##########################################
make.motif.oc.matrix.from.fimo.output = function()
{
  library(data.table)
  fimo.out = '../data/motifs_tfs/fimo.tsv'
  fimo = fread(fimo.out, header = TRUE)
  motif.oc = table(fimo$motif_id, fimo$sequence_name, useNA = 'ifany')
  motif.oc = t(motif.oc)
  
  ##########################################
  # import gene annotation and convert ensgene 
  ##########################################
  load(file = '../data/Hashimsholy_et_al/annotMapping_ensID_Wormbase_GeneName.Rdata')
  mm = match(rownames(motif.oc), geneMapping$Wormbase)
  rownames(motif.oc) = geneMapping$Gene.name[mm]
  #kk = match(rownames(motif.oc, ))
    
  ## modify the motif names with associated TFs
  motif.tf = read.table('../data/motifs_tfs/motifs_tfs_mapping.txt', header = FALSE, sep = ' ')
  motif.tf = motif.tf[, c(2:3)]
  colnames(motif.tf) = c('motifs', 'tfs')
  
  # manually modify the motif names
  motif.tf = data.frame(motif.tf, stringsAsFactors = FALSE)
  motif.tf$motifs.new = motif.tf$motifs
  motif.tf$tfs.new = motif.tf$tfs
  
  xx = motif.tf
  #xx$motifs.new = gsub('1.02', '', xx$motifs.new)
  #xx$motifs.new = gsub('1.02', '', xx$motifs.new)
  #xx$tfs.new = paste0(xx$tfs.new, '_', xx$tfs)
  #xx$tfs.new = gsub('-', '', xx$tfs.new)
  xx$tfs.new = gsub(':', '.', xx$tfs.new)
  xx$tfs.new = gsub('/', '.', xx$tfs.new)
  xx$tfs.new = gsub("\\(","", xx$tfs.new)
  xx$tfs.new = gsub("\\)","", xx$tfs.new)
  xx$tfs.new = gsub("_Homo_sapiens_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Caenorhabditis_briggsae_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Drosophila_melanogaster_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Mus_musculus_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Brugia_pahangi_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Wuchereria_bancrofti_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_PBM_CONSTRUCTS_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Tetraodon_nigroviridis_DBD*.*","", xx$tfs.new)
  
  xx$motifs.new = paste0(xx$motifs.new, '_', xx$tfs.new)
  
  motif.tf = xx
  
  saveRDS(motif.tf, file = '../data/motifs_tfs/motif_tf_mapping.rds')
  
  mm = match(colnames(motif.oc), motif.tf$motifs)
  colnames(motif.oc) = motif.tf$motifs.new[mm]
  
  saveRDS(motif.oc, file = '../data/motifs_tfs/motif_oc_all_proteinCodingGenes.rds')
  
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
