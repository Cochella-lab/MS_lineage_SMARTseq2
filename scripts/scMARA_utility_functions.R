##########################################################################
##########################################################################
# Project: C elegans embryogenesis for MS lineage
# Script purpose: utility functions for scMARA
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Oct 23 18:27:14 2020
##########################################################################
##########################################################################
run.penelized.lm = function(x, y, alpha = 0, intercept=TRUE, standardize=FALSE,  standardize.response=FALSE, 
                            Test = FALSE )
{
  # check the x and y
  cat(length(which(apply(x, 2, sum) == 0)), '  motifs which were not found among all considered genes \n')
  cat(length(which(apply(x, 1, sum) == 0)), ' genes without motifs \n')
  
  cat(ncol(x), ' motifs \n')
  cat(nrow(x), ' gene selected \n')
  
  ### standardize matrix of motif occurrence makes more sense because the absolute number of motif occurrence is not precise.
  if(is.null(ncol(y))){
    family = 'gaussian'
    cat('responses nb : 1  \n')
  }else{
    family = 'mgaussian'
    cat('responses nb :',  ncol(y), '\n')
  }
  
  cat('-- start the penalized linear regression -- \n')
  ### use Cross-validation to select tuning paprameter
  cv.fit=cv.glmnet(x, y, family=family, grouped=FALSE, 
                   alpha=alpha, nlambda=100, standardize=standardize, 
                   standardize.response=standardize.response, intercept=intercept, relax = FALSE)
  plot(cv.fit)
  #s.optimal = cv.fit$lambda.1se
  s.optimal = cv.fit$lambda.min
  
  fit=glmnet(x,y,alpha=alpha, lambda=s.optimal, family=family, 
             standardize=standardize, standardize.response=standardize.response, intercept=intercept, 
             relax = FALSE)
  #plot(fit, xvar = "lambda", label = TRUE, type.coef = "coef")
  ## collect result from the elastic-net
  #colnames(x)[which(fit$beta[[1]][,optimal]!=0)]
  #colnames(x)[which(fit$beta[[2]][,optimal]!=0)]
  
  # extract fitting results for either multiple response or single response; 
  # in particular, we decided to use only ridge penalized linear regression here
  if(family == 'mgaussian'){
    keep = as.data.frame(coef.glmnet(fit, s = s.optimal))
    keep = keep[-1, ] # remove intecept
    colnames(keep) = names(fit$beta)
    keep = apply(keep, 2, scale)
    rownames(keep) = rownames(fit$beta[[2]])
    #ss = apply(keep, 1, function(x) !all(x==0))
    #keep = keep[ss, ]
    #head(rownames(keep)[order(-abs(keep$MSxp))], 10)
    #head(rownames(keep)[order(-abs(keep$MSxa))], 10)
    res = keep
    if(Test){
      ss = apply(keep, 1, function(x) length(which(abs(x)>2.0)))
      keep = keep[which(ss>0), ]
      print(keep)
    }
    
  }else{
    keep = coef.glmnet(fit, s = s.optimal)
    motif.names = keep@Dimnames[[1]]
    motif.names = motif.names[-1]
    keep = keep[-1]
    names(keep) = motif.names
    
    keep = scale(keep)
    o1 = order(-abs(keep))
    keep = keep[o1]
    motif.names = motif.names[o1]
    res = data.frame(motif = motif.names, scores = keep, stringsAsFactors = FALSE)
    if(Test) print(res[which(abs(res$scores) > 2.0), ])
    
    # if(alpha > 0.0) {
    #   keep = keep[which(keep != 0)]
    #   print(names(keep))
    # }else{
     
  }
  
  #for(n in 1:ncol(keep)) keep[,n] = scale(keep[,n], center = TRUE, scale = TRUE)
  #keep[which(abs(keep)<0.05)] = 0
  #keep = data.frame(keep, stringsAsFactors = FALSE)
  return(res)
  
}

########################################################
########################################################
# Section : utility functions for sctf_MARA
# part-1 : worm gene annotation convertion
# part-2 : motif occurrency matrix preparation from fimo output
# part-3 : motif-tf mapping table
# part-4 : pwm clustering based on sequence similarity (here the binding site overlapping is not considered as before) 
# and modify motif-tf mapping 
# part-5: prepare the detecte tf expression profiles
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

display.gene.expression.MS.lineage = function(tf.mat)
{
  library("treeio")
  library("ggtree")
  
  bwm_tree = readRDS(file = paste0(dataDir, 'BWM_tree_for_visualization.rds'))
  
  source.my.script('make_lineage_ggtree_Viscello.R')
  ids.names = colnames(tf.mat)
  ids.names[which(ids.names == "MSxppapp/MSxpappp")] = 'MSxpappp'
  
  pdfname = paste0(resDir, "/BWM_lineage_expressed_TFs_profiles_v2.pdf")
  pdf(pdfname, width=10, height = 8)
  par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  options(warn=-1)
  for(n in 1:nrow(tf.mat))
  #for(n in 1:10)
  {
    cat(n, ' -- ', rownames(tf.mat)[n], '\n')
    
    tf.expr = tf.mat[n, ]
    bwm_tree$value = tf.expr[match(bwm_tree$lineage, ids.names)]
    out.tree = make_lineage_ggtree(bwm_tree, root = 'MS', color.annot = "value") + 
      ggtitle(rownames(tf.mat)[n])
    plot(out.tree)
    
  }
  
  options(warn=0)
  dev.off()
  
  #nwk <- system.file("extdata", "sample.nwk", package="treeio")
  #tree <- read.tree(nwk)
  # ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
  # ggtree(tree, color="firebrick", size=2, linetype="dotted")
  # ggtree(tree, ladderize=TRUE)
  # ggtree(tree, branch.length="none")
  # 
  # 
  # beast_file <- system.file("examples/MCC_FluA_H3.tree", 
  #                           package="ggtree")
  # beast_tree <- read.beast(beast_file)
  # ggtree(beast_tree, mrsd="2013-01-01") + theme_tree2()
  
  #ggtree(tree_tbl, )
  
}

##########################################
# after running FIMO, make motif occurrency matrix  
##########################################
make.motif.oc.matrix.from.fimo.output = function()
{
  library(data.table)
  motif.tf = readRDS( '../data/motifs_tfs/motif_tf_mapping.rds')
  
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
  
  ss1 = apply(motif.oc, 1, sum)
  cat(length(which(ss1 == 0)), 'genes without scanned motifs \n')
  ss2 = apply(motif.oc, 2, sum)
  
  ## merge occurrence for motifs in the same cluster
  #mm = match(colnames(motif.oc), motif.tf$motifs)
  #colnames(motif.oc) = motif.tf$motifs.new[mm]
  names = unique(motif.tf$names[match(colnames(motif.oc), motif.tf$motifs)])
  xx = matrix(0, ncol = length(names), nrow = nrow(motif.oc))
  rownames(xx) = rownames(motif.oc)
  colnames(xx) = names
  
  for(n in 1:ncol(xx))
  {
    # n = 2
    cat(n, '\n')
    mtf = motif.tf$motifs[which(motif.tf$names == colnames(xx)[n])]
    kk = match(mtf, colnames(motif.oc))
    kk = kk[!is.na(kk)]
    
    if(length(kk) == 0){
      cat('Error : no motif found \n')
    }else{
      if(length(kk) == 1){
        xx[,n] = motif.oc[, kk]
      }else{
        xx[,n] = ceiling(apply(motif.oc[,kk], 1, median))
      }
    }
  }
  
  motif.oc = xx;
  saveRDS(motif.oc, file = '../data/motifs_tfs/motif_oc_all_proteinCodingGenes.rds')
  
}

## modify the motif names with associated TFs
process.motif.tf.mapping = function()
{
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
  
}

cluster.pwm.based.similarity = function()
{
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra) # clustering visualization
  require(graphics)
  
  pwm.corr = read.table(file = '../data/motifs_tfs/pwm_similarity_correction_PCC.txt', header = TRUE, 
                        row.names = 1)
  motif.tf$motifs.new = paste0(motif.tf$tfs.new, '_', motif.tf$motifs)
  
  newName = motif.tf$motifs.new[match(rownames(pwm.corr), motif.tf$motifs)]
  rownames(pwm.corr) = newName
  colnames(pwm.corr) = newName
  
  comparisons <- 1 - pwm.corr
  dd <- as.dist(comparisons)
  
  # Hierarchical clustering using Complete Linkage
  hc <- hclust(dd, method = "ward.D2" )
  
  # Plot the obtained dendrogram
  #plot(hc, cex = 0.6, hang = -1)
  #sub_grp <- cutree(hc, h = 0.1)
  pdfname = paste0(resDir, "/pwm_celegans_similarity_clustering.pdf")
  pdf(pdfname, width=20, height = 30)
  par(cex =0.5, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  hc.cutoff = 0.1
  #plot(hc, cex = 0.5, hang = -1)
  plot(as.dendrogram(hc), cex=0.5, horiz=TRUE)
  abline(v = c(0.1, 0.15, 0.2), col = 'red')
  #rect.hclust(hc, h = hc.cutoff, border="darkred")
  #groups <- 
  length(unique(cutree(hc, h = 0.1)))
  length(unique(cutree(hc, h = 0.15)))
  length(unique(cutree(hc, h = 0.2)))
  length(unique(cutree(hc, h = 0.25)))
  
  dev.off()
  
  change.pwm.logo.names = FALSE
  if(change.pwm.logo.names){
    setwd('../data/motifs_tfs/pwm_logo')
    logo.file = list.files(path = '.', pattern = '*.pdf', full.names = FALSE)
    for(n in 1:nrow(motif.tf))
    {
      # n = 1 
      cmd = paste0('mv ', logo.file[grep(motif.tf$motifs[n], logo.file)], ' ',  motif.tf$motifs.new[n], '.pdf')
      system(cmd)
    }
    
  }
  #fviz_nbclust(diss = comparisons, FUN = hcut, method = "wss")
  #fviz_nbclust(df, FUN = hcut, method = "silhouette")
  
  ##########################################
  # merge motifs using height = 0.1 and change motif names
  ##########################################
  groups <- cutree(hc, h = 0.1)
  motif.tf = data.frame(motif.tf, group = groups, stringsAsFactors = FALSE)
  motif.tf$names = NA
  for(nn in unique(motif.tf$group))
  {
    # nn = 5
    kk = which(motif.tf$group == nn)
    motif.tf$names[kk] = paste0(paste0(unique(motif.tf$tfs.new[kk]), collapse = '_'), '.M', nn)
    
  }
  
  saveRDS(motif.tf, file = '../data/motifs_tfs/motif_tf_mapping.rds') # motif-to-tf mapping for non-redundant motifs (to some extent)
  
}

process.detected.tf.expression.profiles = function(Y.mat)
{
  # subset Y.mat for TFs
  jj = match(tfs$`Public name`, rownames(Y.mat))
  jj = jj[!is.na(jj)]
  tf.mat = Y.mat[jj, ]
  cutoff.tf = 1;
  ss = apply(tf.mat, 1, function(x) !all(x<cutoff.tf))
  tf.mat = tf.mat[ss, ]
  
  save.tf.profiles.across.lineage = FALSE
  if(save.tf.profiles.across.lineage){
    pdfname = paste0(resDir, "/TFs_standardized_fpkm_in_BWM.pdf")
    pdf(pdfname, width=12, height = 50)
    par(cex =0.3, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    pheatmap(tf.mat, cluster_rows=TRUE, 
             show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
             scale = 'row',
             cluster_cols=FALSE, 
             main = paste0("standardized fpkm of TFs "), 
             na_col = "white",
             #color = cols, 
             #annotation_col = my_sample_col,
             #gaps_row = c(1:nrow(map)-1),
             fontsize_col = 10
    )
    dev.off()
    
    
    write.csv(tf.mat, file = paste0(tabDir, 'detected_TFs_in_BWM.csv'), row.names = TRUE)
  }
  
  saveRDS(tf.mat, file = paste0(RdataDir, 'TFs_expression_profiles_BWM.rds')) 
  
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
