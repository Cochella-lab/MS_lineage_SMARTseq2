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
  markers.sels = markers[which(markers$p_val<10^-3 & markers$avg_logFC > 0.5), ]
  print(table(markers.sels$cluster))
  
  ids.groups = list(ids.uniq[which(nchar(ids.uniq) <=5)], 
                    ids.uniq[which(nchar(ids.uniq) == 6)],
                    ids.uniq[which(nchar(ids.uniq) == 7)], 
                    ids.uniq[which(nchar(ids.uniq) > 7)]) 
  
  ids.groups = as.list(ids.uniq)
  
  for(n in 1:length(ids.groups)){
    
    # n = 1;
    lineage = ids.groups[[n]]
    cat(n, ' --- ')
    print(lineage)
    # lineage = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
    # lineage = c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
    # lineage = setdiff(ids.uniq, c("mixture_terminal_1", "mixture_terminal_2"))
    
    gene.sels = markers.sels[!is.na(match(markers.sels$cluster, lineage)), ]
    gene.sels = gene.sels[!is.na(match(gene.sels$gene, rownames(Y.mat))), ]
    #print(table(gene.sels$cluster))
    gene.sels = unique(gene.sels$gene)
    
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
    res = run.penelized.lm(x, y, alpha = 0, standardize = TRUE, intercept = TRUE, use.lambda.min = TRUE, 
                           Test = FALSE)
    
    print(res[grep('pha-4|hnd-1..Tcf|nhr-67.homo.M227|hlh-1.M175|unc-120.dm', rownames(res)),])
    
    if(n == 1) {
      keep = res[match(colnames(x), rownames(res)), ]
    }else{
      keep = cbind(keep, res[match(colnames(x), rownames(res)), ])
    }
        
  }
  
  print(keep[grep('pha-4|hnd-1..Tcf|nhr-67.homo.M227|hlh-1.M175|unc-120.dm', rownames(keep)),])
  
  ss = apply(keep, 1, function(x) length(which(abs(x)>1.5)))
  length(which(ss>=1))
  
  yy = keep[which(ss>0), ] 
  yy[which(abs(yy)>2.5)] = 2.5
  
  pdfname = paste0(resDir, "/MARA_prediction_all_lineages_v1.pdf")
  pdf(pdfname, width=18, height = 16)
  par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  pheatmap(yy[grep('pha-4|hnd-1|nhr-67|hlh-1|unc-120.dm', rownames(yy)), ], 
           cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
           scale = 'none', cluster_cols=FALSE, main = paste0("MARA prediction for controls"), 
           na_col = "white", fontsize_col = 12) 
  
  pheatmap(yy, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
           scale = 'none', cluster_cols=FALSE, main = paste0("MARA prediction"), 
           na_col = "white", fontsize_col = 12) 
  
  
  dev.off()
  
}
