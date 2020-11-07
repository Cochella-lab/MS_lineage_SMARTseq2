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
  ids.uniq = ids.uniq[order(nchar(ids.uniq))]
  
  # convert to SingleCellExperiment, recalculate scaling factor and normalized to fpkm
  sce = as.SingleCellExperiment(sub.obj)
  
  mm = match(rownames(sce), ll$gene.name)
  sce = sce[which(!is.na(mm)), ] # keep genes with correponding lengths
  ll = ll[mm[which(!is.na(mm))], ]
  
  sce <- logNormCounts(sce, log = FALSE, size_factors = NULL)
  Y.fpkm <- log2(calculateFPKM(sce, lengths = ll$transcript.length) + 1)
  
  remove(sce)
  
  ##########################################
  # determine lineage-specific signatures 
  # i.e. selected gene sets (lineage-wide, specific, or restricted)
  ##########################################
  define.modules.for.lineags(Y.fpkm)
  
  
  ##########################################
  # prepare matrix A and reponse Y and run penalized.lm
  ##########################################
  lineage = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
  
  lineage = c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
  lineage = setdiff(ids.uniq, c("mixture_terminal_1", "mixture_terminal_2"))
  
  gene.sel = unique(markers$gene[which(!is.na(match(markers$cluster, lineage)) & markers$p_val_adj<10^-5 & markers$avg_logFC > 0.7)])
  gene.sel = gene.sel[which(!is.na(match(gene.sel, rownames(Y.mat))))]
  
  remove.shared.genes = FALSE
  if(remove.shared.genes){
    gene.sel_1 = gene.sel
    gene.sel_2 = gene.sel
    gene.shared = intersect(gene.sel_1, gene.sel_2)
    
    gene.sel_1 = setdiff(gene.sel_1, gene.shared)
    gene.sel_2 = setdiff(gene.sel_2, gene.shared)
    
    lineage = c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
    gene.sel = gene.sel_1
    
    lineage = c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
    gene.sel = gene.sel_2
    
  }
  
  index.sel = match(gene.sel, rownames(Y.mat))
  Y.sel = Y.mat[index.sel, match(lineage, colnames(Y.mat))]
  y = as.matrix(Y.sel)
  
  #mm = match(rownames(Y.sel), rownames(motif.oc))
  mm = match(gene.sel, rownames(motif.oc))
  y = y[!is.na(mm), ]
  x = as.matrix(motif.oc[mm[!is.na(mm)], ])
  x[which(is.na(x) == TRUE)] = 0
  
  source.my.script('scMARA_utility_functions.R')
  res = run.penelized.lm(x, y, alpha = 0.1, standardize = FALSE, intercept = TRUE, use.lambda.min = FALSE, 
                         Test = TRUE)
  
  print(res[grep('pha-4.mus|hnd-1..Tcf|nhr-67.homo.M227|hlh-1.M175|unc-120.dm', rownames(res)),])
  
  #ss = apply(res, 1, function(x) length(which(abs(x)>1.3)))
  pheatmap(res[which(ss>0), ], cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
           scale = 'none', cluster_cols=FALSE, main = paste0("MARA prediction"), 
           na_col = "white", fontsize_col = 10) 
  
  
  
  return(res)
  
}
