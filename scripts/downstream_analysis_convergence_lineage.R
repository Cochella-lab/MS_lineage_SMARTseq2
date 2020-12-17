##########################################################################
##########################################################################
# Project: MS lineage 
# Script purpose: analyze convergence lineages and identify potential regulators
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Dec 16 11:55:32 2020
##########################################################################
##########################################################################
library(pheatmap)
library(RColorBrewer)
library(grid)
library(Seurat)
library(scater)
library(SingleCellExperiment)
library(scran)

aggregate.cells.across.ids = function(seurat.obj)
{
  ids.current = names(table(seurat.obj$manual.annot.ids))
  cat(length(ids.current), ' annotated ids in total \n')
  
  ids.to.drop = ids.current[grep('doublets|like', ids.current)]
  if(length(ids.to.drop) >0) {
    ids.sel = setdiff(ids.current, ids.to.drop)
    cat('ids to drop : \n')
    print(ids.to.drop)
  }else{
    cat('no ids dropped \n')
    ids.sel = ids.current
  }
  
  cat(length(ids.sel), ' annotated ids were selected to aggregate \n')
  
  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sel))])
  
  cat(length(cells.sels), ' cells selected to annotate \n')
  
  sub.obj = subset(seurat.obj, cells = cells.sels)
  
  # convert seurat object to SingleCellExperiment 
  sce = as.SingleCellExperiment(sub.obj)
  
  summed <- aggregateAcrossCells(sce, 
                                 id=colData(sce)[, c("manual.annot.ids")])
  summed
  
  # Creating up a DGEList object for use in edgeR:
  library(edgeR)
  y <- DGEList(counts(summed), samples=colData(summed))
  #y
  
  # normalize pseudo-bulk 
  y <- calcNormFactors(y)
  y$samples
  
  par(mfrow=c(2,2))
  for (i in seq_len(ncol(y))) {
    plotMD(y, column=i)
  }
  
  #plotMDS(cpm(y, log=TRUE), 
  #        col=ifelse(y$samples$tomato, "red", "blue"))
  
  # output is DGEList object from edgeR 
  return(y) 
  
}

compare.convergence.lineages.with.others = function(y)
{
  library("pheatmap")
  library("RColorBrewer")
  
  # here the input is the DGEList object from edgeR
  cpm = edgeR::cpm(y, log = TRUE, prior.count = 1)
  
  sampleDists <- dist(t(cpm))
  
  sampleDistMatrix <- as.matrix(sampleDists)
  #rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
  #colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  
  
  # Grouped Bar Plot
  counts <- table(mtcars$vs, mtcars$gear)
  barplot(counts, main="Car Distribution by Gears and VS",
          xlab="Number of Gears", col=c("darkblue","red"),
          legend = rownames(counts), beside=TRUE)
  ids.convergence = c('MSx', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
  list2.ids = c('MSxp', 'MSxpp', 'MSxppp', 'MSxapa', 'MSxpppp', 'MSxapap', 'MSxapapp', 
                'MSxppppp', "MSpaaappp/MSxapappa")
  
  compares =  matrix(NA, nrow = length(list2.ids), ncol = length(list1.ids))
  colnames(compares) = list1.ids
  rownames(compares) = list2.ids
  ii1 = match(list1.ids, colnames(sampleDistMatrix))
  ii2 = match(list2.ids, rownames(sampleDistMatrix))
  compares = sampleDistMatrix[ii2, ii1]
  
  pdfname = paste0(resDir, "/convergence_lineage_distance_to_others.pdf")
  pdf(pdfname, width=12, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  barplot(compares, beside = TRUE, col = c(1:nrow(compares)), 
          legend.text = rownames(compares), args.legend = c(x = 'topleft', bty = 'n'), 
          ylim = c(0, 800))
          
  dev.off()
  
  
}



