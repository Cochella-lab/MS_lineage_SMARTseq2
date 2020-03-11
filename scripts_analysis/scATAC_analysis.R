##########################################################################
##########################################################################
# Project: Aleks' single-cell MS project 
# Script purpose: analyze scATAC-seq with cisTopic
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Mon Mar  2 14:43:04 2020
##########################################################################
##########################################################################
source.my.script <- function(name.of.function){
  tryCatch(path <- rstudioapi::getSourceEditorContext()$path,
           error = function(e){
             install.packages("rstudioapi")
             path <-  rstudioapi::getSourceEditorContext()$path})
  source.path <- sub(basename(path), "", path)
  source(paste0(source.path,name.of.function))
}

## set up the paths for the data and results
tryCatch(path <- rstudioapi::getSourceEditorContext()$path, 
         error = function(e){
           install.packages("rstudioapi")
           path <-  rstudioapi::getSourceEditorContext()$path})
source.path <- sub(basename(path), "", path)


user <- "results_jiwang/"
setwd(paste0("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/",user))

version.DATA = 'scATAC_earlyEmbryo'
version.analysis =  paste0(version.DATA, '_20200302')
dataDir = paste0("../data/gene_counts/")
resDir = paste0("results/", version.analysis)
tabDir = paste0("results/", version.analysis, "/tables/")

RdataDir = paste0("results/", version.analysis, "/Rdata/")
if(!dir.exists("results/")){dir.create("results/")}
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

source.my.script('scATAC_functions.R')
if(!require("DropletUtils")) BiocManager::install('DropletUtils')
library(DropletUtils)
library(data.table)
library(Matrix)
library(tictoc)

scATACpro_output_path = '/Volumes/groups/cochella/jiwang/Projects/Aleks/R8898_scATAC/output/'
raw_mtx_file = paste0(scATACpro_output_path, '/raw_matrix/MACS2/matrix.mtx')
filtered_mtx_dir = paste0(scATACpro_output_path, 'filtered_matrix_peaks_barcodes')
raw_mtx_dir = dirname(raw_mtx_file)

if(!dir.exists(filtered_mtx_dir)){dir.create(filtered_mtx_dir)}

########################################################
########################################################
# Section : prepare annotation files
# 
########################################################
########################################################
# # args = commandArgs(T)
Generate.enhancer.annotation = FALSE
if(Generate.enhancer.annotation){
  enhancers = read.table('/Volumes/groups/cochella/jiwang/annotations/promoter_enhancer_annotation-elife-37344-fig2-data1-v2.txt', header = TRUE)
  enhancers = enhancers[ which(enhancers$annot == 'putative_enhancer'), ]
  enhancers = data.frame(enhancers$chrom_ce11, enhancers$start_ce11, enhancers$end_ce11, enhancers$associated_gene_id, stringsAsFactors = FALSE)
  enhancers$scores = '.'
  enhancers$strand = '.'
  write.table(enhancers, file = '/Volumes/groups/cochella/jiwang/annotations/ce11_enhancers_Jaenes_et_al_elife_2018.bed',
              row.names = FALSE, col.names = FALSE, sep = '\t', quote = FALSE)
}

########################################################
########################################################
# Section : Cell filtering using DropletUtils
#  prepare filtered matrix, barcodes and peaks
########################################################
########################################################
mat = readMM(raw_mtx_file)
features = fread(paste0(raw_mtx_dir, '/features.txt'), header = F)
barcodes = fread(paste0(raw_mtx_dir, '/barcodes.txt'), header = F)
rownames(mat) = features$V1
colnames(mat) = barcodes$V1

# metrics = read.csv('res_counts/outs/singlecell.csv', header = TRUE)

set.seed(2019)

tic('emptyDrops running :')
cell.out <- emptyDrops(mat, BPPARAM = SerialParam())
toc()

PLOT.barcode.rank = FALSE
if(PLOT.barcode.rank){
  
  pdfname = paste0(resDir, "/EmptyDrop_barcode_filtering.pdf")
  pdf(pdfname, width=14, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #my.counts <- DropletUtils:::simCounts()
  br.out <- barcodeRanks(mat)

  # Making a plot.
  plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
  o <- order(br.out$rank)
  lines(br.out$rank[o], br.out$fitted[o], col="red")
  abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
  abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
  legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
         legend=c("knee", "inflection"))
  
  
  #set.seed(100)
  #cell.out <- cell.out
  cell.out
  
  is.cell <- cell.out$PValue < 0.0001
  sum(is.cell, na.rm=TRUE)
  
  table(Limited=cell.out$Limited, Significant=is.cell)
  plot(cell.out$Total, -cell.out$LogProb, col=ifelse(is.cell, "red", "black"),
       xlab="Total UMI count", ylab="-Log Probability", log='x')
  
  dev.off()
  
}


filter.out <- cell.out[complete.cases(cell.out), ]

saveRDS(filter.out, file = paste0(filtered_mtx_dir, '/EmptyDrop_obj.rds'))

fdr = 0.01
#filter.out = filter.out[filter.out$FDR <= fdr, ]
kk = filter.out$Total>10^3
sum(kk)

filter.out = filter.out[kk, ]
select.cells = rownames(filter.out)
#length(select.cells)

out_mat = mat[, colnames(mat) %in% select.cells]
barcodes = colnames(out_mat)
peaks = rownames(out_mat)

peaks = t(sapply(peaks, FUN = function(x){
  ccord = unlist(strsplit(as.character(x), "-"))
  if(!startsWith(ccord[1], 'chr')){
    return(c(paste0('chr', ccord[1]), ccord[-1]))
  }else{
    return(ccord)
  }  
}))

#metrics$barcode = sapply(metrics$barcode, function(x) gsub('-1', '', x))
#mm = match(barcodes, metrics$barcode)
#metrics = metrics[]

## save filtered matrix, peaks and barcodes
system(paste('mkdir -p', filtered_mtx_dir))
writeMM(out_mat, file = paste0(filtered_mtx_dir, '/matrix.mtx'))  

write.table(barcodes, file = paste0(filtered_mtx_dir, '/barcodes.tsv'), sep = '\t', 
            row.names = F, quote = F, col.names = F)
write.table(peaks, file = paste0(filtered_mtx_dir, '/peaks.bed'), sep = '\t',
            row.names = F, quote = F, col.names = F)


########################################################
########################################################
# Section :  lsi and lsi_log transform
# 
########################################################
########################################################
source.my.script('scATAC_functions.R')
tenx.bmat = load_tenx_atac(paste0(filtered_mtx_dir, '/matrix.mtx'), 
                                 paste0(filtered_mtx_dir, '/peaks.bed'), 
                                 paste0(filtered_mtx_dir, '/barcodes.tsv'))

# select features showing in > 50 cells
tenx.bmat = filter_features(tenx.bmat, cells=50)

# Binarize the matrix for consistency
tenx.bmat@x[tenx.bmat@x > 1] = 1

tenx.seurat.lsi = lsi_workflow(tenx.bmat, dims=2:50, log_scale_tf=FALSE, reduction='pca', resolution=1.5)
tenx.seurat.lsi_log = lsi_workflow(tenx.bmat, dims=2:50, log_scale_tf=TRUE, reduction='pca.l2', resolution=1.5)

plot_clustering_comparison(tenx.seurat.lsi,
                           tenx.seurat.lsi_log,
                           reduction='umap',
                           description1='LSI',
                           description2='LSI logTF',
                           cluster_column1='RNA_snn_res.1.5',
                           cluster_column2='RNA_snn_res.1.5')


nb.pcs = 50; n.neighbors = 20; min.dist = 0.2;
tenx.seurat.lsi_log <- RunUMAP(object = tenx.seurat.lsi_log, reduction = 'pca.l2', dims = 2:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(object = tenx.seurat.lsi_log, label = TRUE, reduction = 'umap') + NoLegend()


##########################################
# compare LAD from cistopic with lsi-log
##########################################
source.my.script('scATAC_functions.R')

Run.cistopic_workflow = FALSE
if(Run.cistopic_workflow){
  tic('run model for cisTopic')
  tenx.cistopic = cistopic_workflow(tenx.bmat, topic = seq(20, 100, by=5))
  saveRDS(tenx.cistopic, file =  paste0(RdataDir, 'atac_cisTopic_runWrapLDAModel_localRunning.rds'))
  
  toc 
}

tenx.cistopic = readRDS(file =  paste0(RdataDir, 'atac_cisTopic_runWrapLDAModel_localRunning.rds'))
#tenx.cistopic = readRDS('data_downloads/atac_v1_adult_brain_fresh_5k.cistopic.rds')

cisTopicObject = tenx.cistopic
par(mfrow=c(3,3))
cisTopicObject <- selectModel(cisTopicObject, type='maximum')
cisTopicObject <- selectModel(cisTopicObject, type='perplexity')
#cisTopicObject <- selectModel(cisTopicObject, type='derivative')

seurat.cistopic = seurat_workflow_on_cistopic(tenx.cistopic, resolution=0.8, method='Z-score')
# plot_clustering_comparison(tenx.seurat.lsi_log,
#                            seurat.cistopic,
#                            reduction='umap',
#                            description1='LSI logTF',
#                            description2='cisTopic',
#                            cluster_column1='RNA_snn_res.1.5',
#                            cluster_column2='RNA_snn_res.1.5')

nb.pcs = 35; n.neighbors = 30; min.dist = 0.25;
seurat.cistopic <- RunUMAP(object = seurat.cistopic, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(object = seurat.cistopic, label = TRUE, reduction = 'umap') + NoLegend()

saveRDS(seurat.cistopic, file =  paste0(RdataDir, 'atac_LDA_seurat_object.rds'))

########################################################
########################################################
# Section : cluster annotations using cell-type specific features (motifs) or cis-regulatory elements 
# or using scRNA-seq data
########################################################
########################################################
seurat.cistopic = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object.rds'))

p1 <- DimPlot(seurat.cistopic, label = TRUE, pt.size = 0.1, label.size = 5) + NoLegend()
p1


##########################################
# generate bigwigs for each cluster for visualization
##########################################
barcode.cluster = data.frame(Barcode = colnames(seurat.cistopic),
                             Cluster = seurat.cistopic$seurat_clusters, 
                             stringsAsFactors = FALSE)
write.table(barcode.cluster, file = paste0(filtered_mtx_dir, '/barcodes_and_clusters.txt'), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

##########################################
# compute gene activity score to annotate clusters
##########################################
peaks.chrM = grep('chrM', rownames(seurat.cistopic))
seurat.cistopic = seurat.cistopic[-peaks.chrM, ]

source.my.script('scATAC_functions.R')
fragment.file = '/Volumes/groups/cochella/jiwang/Projects/Aleks/R8898_scATAC/cellranger_atac_wbcel235/outs/fragments.tsv.gz'

seurat.cistopic = compute.gene.acitivity.scores(seurat.cistopic, fragment.file = fragment.file)

##########################################
# motif enrichment analysis  
##########################################
seurat.cistopic = compute.motif.enrichment(seurat.cistopic)


# tbx-38
p1 = FeaturePlot(
  object = seurat.cistopic,
  features = rownames(seurat.cistopic)[201],
  #cols = c('red', 'blue'),
  min.cutoff = 0,
  max.cutoff = 'q95',
  pt.size = 0.1
)

# med-2 
p2 = FeaturePlot(
  object = seurat.cistopic,
  features = rownames(seurat.cistopic)[106],
  #cols = c('red', 'blue'),
  min.cutoff = 0,
  max.cutoff = 'q95',
  #blend = TRUE,
  pt.size = 0.1
)

# elt-3
p3 = FeaturePlot(
  object = seurat.cistopic,
  features = rownames(seurat.cistopic)[280],
  #cols = c('red', 'blue'),
  min.cutoff = 0,
  max.cutoff = 'q95',
  #blend = TRUE,
  pt.size = 0.1
)

# pal-1 
p4 = FeaturePlot(
  object = seurat.cistopic,
  features = rownames(seurat.cistopic)[292],
  #cols = c('green', 'blue'),
  min.cutoff = 0,
  max.cutoff = 'q95',
  #blend = TRUE,
  pt.size = 0.1
)

CombinePlots(list(p1, p2, p3, p4), ncol = 2)


