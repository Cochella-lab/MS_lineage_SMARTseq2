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

compare.convergence.lineages.with.others = function(y, method = c('euclidean', 'correlation', 'jsd'))
{
  library("pheatmap")
  library("RColorBrewer")
  #library(philentropy)
  
  # here the input is the DGEList object from edgeR
  cpm = edgeR::cpm(y, log = TRUE, prior.count = 1)
  
  method = 'jsd'
  #colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  if(method == 'euclidean'){ sampleDists <- dist(t(cpm), method = 'euclidean'); xlim = c(0, 500)}
  
  if(method == 'correlation'){ sampleDists = cor(cpm, method = 'pearson'); xlim = c(0, 1)}
  
  if(method == 'jsd') {
    library(cummeRbund)
    cpm2 = log2(2^cpm + 1)
    probs <- cummeRbund::makeprobs(cpm2)
    sampleDists <-cummeRbund::JSdist(probs)
    #sampleDists = JSD(t(as.matrix(cpm)), unit = 'log2')
  }
  
  sampleDistMatrix <- as.matrix(sampleDists)
  
  ids.convergence = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
  ids.bwm = c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
  ids.phrx = c('MSxapa', 'MSxapap', 'MSxapapp', "MSpaaappp/MSxapappa")
  
  pdfname = paste0(resDir, "/convergence_lineage_compared_with_one.BWM.lineage_one.Pharynx.lineage_", method, ".pdf")
  pdf(pdfname, width=12, height = 8)
  par(cex =0.7, mar = c(5,8,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
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
            names.arg = names(cc.dist), las = 2, horiz = TRUE,
            xlim = xlim)
    
  }
  
  names(cc.vec) = cc.vec.names
  
  dev.off()
  
  ids.mothers = c('MSx', 'MSxa', 'MSxap')
   
  pdfname = paste0(resDir, "/convergence_lineage_compared_mothers_sisters_", method, ".pdf")
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
           las = 2, xlim = xlim)
    
  }
  
  dev.off()
  
}

find.regulators.for.convergence.lineage = function(y)
{
  library("pheatmap")
  library("RColorBrewer")
    
  # here the input is the DGEList object from edgeR
  cpm = edgeR::cpm(y, log = TRUE, prior.count = 1)
  
  ids.convergence = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
  ids.bwm = c('MSx', 'MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppx', 'MSxppppp')
  ids.phrx = c('MSxapa', 'MSxapap', 'MSxapapp', "MSpaaappp/MSxapappa")
  #ids.phrx2 = c('MSxaa')
  
  ids.sel = c(ids.bwm, ids.convergence, ids.phrx)
  #indexs = c(c(1:7), c(2:7), c(3:6))
  #cols = c('black', rep('darkblue', 6), rep('red', 6), rep('orange', 4))
  
  tfs = readxl::read_xlsx('../data/motifs_tfs/Table-S2-wTF-3.0-Fuxman-Bass-Mol-Sys-Biol-2016.xlsx', sheet = 1)
  gene2plot = unique(c('pha-4', 'hnd-1', 'hlh-1', 'unc-120', 'nhr-67', tfs$`Public name`))
  gene2plot = gene2plot[!is.na(match(gene2plot, rownames(cpm)))]
  
  cpmxx = cpm[match(gene2plot, rownames(cpm)), match(ids.sel, colnames(cpm))]
  kk = apply(as.matrix(cpmxx), 1, function(x) !all(as.numeric(x)<2))
  cpmxx = cpmxx[kk, ]
  
  pdfname = paste0(resDir, "/convergence_lineage_regulators_profiles_v2.pdf")
  pdf(pdfname, width=10, height = 6)
  par(cex =1.0, mar = c(4,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(n in 1:nrow(cpmxx))
  {
    # n = 2
    if(!all(is.na(cpmxx[n]))){
      cat(n, ' -- ', rownames(cpmxx)[n], '\n')
      ylims = range(cpmxx[n, ]) + c(-1, 1)
      plot(c(1:7), cpmxx[n, c(1:7)], ylim = ylims, col = 'darkblue', type='b', main = rownames(cpmxx)[n], pch = 15, lwd = 1.5,
           ylab = 'log2(cpm)', xlab = 'cell ids')
      points(c(1:7), cpmxx[n, c(1, 8:13)], col = 'red', type = 'b', lwd = 3.0, pch = 16)
      points(c(3:7), cpmxx[n, c(9, 14:17)], col = 'black', type = 'b', pch =17)
      abline(v = c(4, 5), lwd =2.0, lty = 'dashed', col = 'gray')
      legend('topright', c('bwm', 'convergence', 'pharynx'),lwd = c(1, 2, 1), pch = c(15, 16, 17),
             col = c('darkblue', 'red', 'black'),  adj = c(0, 0.6), bty = 'n')
      # add text
      text(c(1:7), cpmxx[n, c(1:7)], ids.bwm, pos = 3, cex = 0.8)
      text(c(2:7), cpmxx[n, c(8:13)], ids.convergence, pos = 3, cex = 0.8)
      text(c(4:7), cpmxx[n, c(14:17)], ids.phrx, pos = 3, cex = 0.8)
    }
   
  }
  dev.off()
  
}

