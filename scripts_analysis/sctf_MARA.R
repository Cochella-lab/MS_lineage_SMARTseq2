##########################################################################
##########################################################################
# Project: MS lineage embryogenesis in C elegans
# Script purpose: predicting TF regulation activity 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct  7 17:31:46 2020
##########################################################################
##########################################################################
library("pheatmap")
library("RColorBrewer")
library(grid)


predict.TF.MARA.for.scdata = function(sub.obj, mode = 'cluster.wise', id = 'manual.annot.ids', Y_name = 'RNA')
{
  # mode = 'cluster.wise';
  ids = sub.obj$manual.annot.ids
  Y.data = sub.obj@assays$RNA@data
  ids.uniq = unique(ids)
  
  if(mode == 'cluster.wise'){
    cat('-- averging the gene expression in clusters -- \n')
    
    Y.mat = matrix(NA, nrow = nrow(Y.data), ncol = length(ids.uniq))
    colnames(Y.mat) = ids.uniq
    rownames(Y.mat) = rownames(Y.data)
    
    for(n in 1:length(ids.uniq)) 
    {
      jj = which(ids == ids.uniq[n])
      if(length(jj) == 1) Y.mat[, n] = Y.data[,jj]
      if(length(jj) > 1) Y.mat[, n] = apply(Y.data[,jj], 1, mean)
    }
    
    cal_z_score <- function(x){ (x - mean(x)) / sd(x)}
    Y.norm <- t(apply(Y.mat, 1, cal_z_score))
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
    
  }else{
    Y.mat = Y.data
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
