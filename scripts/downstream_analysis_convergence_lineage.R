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
  #counts <- table(mtcars$vs, mtcars$gear)
  #barplot(counts, main="Car Distribution by Gears and VS",
  #        xlab="Number of Gears", col=c("darkblue","red"),
  #        legend = rownames(counts), beside=TRUE)
  ids.convergence = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
  ids.bwm = c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
  ids.phrx = c('MSxapa', 'MSxapap', 'MSxapapp', "MSpaaappp/MSxapappa")
  
  pdfname = paste0(resDir, "/convergence_lineage_compared_with_one.BWM.lineage_one.Pharynx.lineage.pdf")
  pdf(pdfname, width=12, height = 8)
  par(cex =0.7, mar = c(5,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  cc = matrix(NA, ncol = length(ids.convergence), nrow = length(c(ids.bwm, ids.phrx)))
  colnames(cc) = ids.convergence
  rownames(cc) = c(ids.bwm, ids.phrx)
  
  #ii1 = match(ids.convergence[-c(1:2)], colnames(sampleDistMatrix))
  #ii2 = match(ids.phrx, rownames(sampleDistMatrix))
  cc.vec = c()
  cc.vec.names = c()
  for(n in 1:ncol(cc))
  {
    # n = 1
    ii1 = match(colnames(cc)[n], colnames(sampleDistMatrix))
    ids2compare = c(ids.bwm[which(nchar(ids.bwm) <= (nchar(colnames(cc)[n]) + 1))], 
                    ids.phrx[which(nchar(ids.phrx) == nchar(colnames(cc)[n]))])
    ii2 = match(ids2compare, rownames(sampleDistMatrix))
    
    cc[match(ids2compare, rownames(cc)), n] = sampleDistMatrix[ii2, ii1]
    cc.vec = c(cc.vec, sampleDistMatrix[ii2, ii1])
    cc.dist = sampleDistMatrix[ii2, ii1]
    cc.vec.names = c(cc.vec.names, paste0(colnames(cc)[n], '_vs_', ids2compare))
    
    barplot(cc.dist, beside = TRUE, col=c(1:length(cc.dist)), main = colnames(cc)[n], 
            #legend.text = rownames(cc.vec), args.legend = c(x = 'topleft', bty = 'n'), 
            names.arg = names(cc.dist), las = 2,
            ylim = c(0, 600))
    
  }
  names(cc.vec) = cc.vec.names
  
  dev.off()
  
  ids.mothers = c('MSx', 'MSxa', 'MSxap')
   
  pdfname = paste0(resDir, "/convergence_lineage_compared_mothers_sisters.pdf")
  pdf(pdfname, width=10, height = 6)
  par(cex =0.7, mar = c(4,10,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(n in 1:length(ids.mothers))
  {
    # n = 1
    ids2compare = c(ids.mothers[n], paste0(ids.mothers[n], c('p', 'a')))
    #ii1 = match(list1.ids, colnames(sampleDistMatrix))
    #ii2 = match(list2.ids, rownames(sampleDistMatrix))
    compares = sampleDistMatrix[match(ids2compare, rownames(sampleDistMatrix)), 
                                match(ids2compare, colnames(sampleDistMatrix))]
    compares[lower.tri(compares, diag = TRUE)] <- NA
    cc.names = expand.grid(colnames(compares), '_vs_', rownames(compares))
    cc.names = paste0(cc.names[, 3], cc.names[,2], cc.names[,1]) 
    cc = as.vector(compares)
    names(cc) = cc.names
    cc = cc[!is.na(cc)]
    
    barplot(cc, beside = TRUE, col = c(1:length(cc)), horiz = TRUE,
            names.arg = names(cc),
           las = 2, xlim = c(0, 500))
    
  }
  
  dev.off()
  
}



