##########################################################################
##########################################################################
# Project: single-cell analysis for scATAC-seq 
# Script purpose: main function and annotate clusters
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri May  8 15:24:55 2020
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
setwd(paste0("../", user))

version.DATA = 'scATAC_earlyEmbryo'
version.analysis =  paste0(version.DATA, '_20200302')
resDir = paste0("results/", version.analysis)
tabDir = paste0("results/", version.analysis, "/tables/")
RdataDir = paste0("results/", version.analysis, "/Rdata/")
if(!dir.exists("results/")){dir.create("results/")}
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

source.my.script('scATAC_functions.R')
source.my.script('scATAC_cluster_annotation_functions.R')

########################################################
########################################################
# Section : main function
# 
########################################################
########################################################
source.my.script('scATAC_cluster_annotation_functions.R')

## import peak, RNA (gene activity matrix) as Assays
# gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter.geneBody_2000bp.rds')
# gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.500bpDownstream.rds')
# gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.2000bpDownstream.rds')
gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_nonOverlapped.withGenes_promoter_2000bpUpstream.500bpDownstream.rds')
seurat.cistopic = readRDS(file = gamat)

# add motif activities into seurat.cistopic
xx = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object_motifClassChromVar.rds'))
seurat.cistopic[['chromvar']] = CreateAssayObject(data = xx@assays$chromvar@data)

# remove cluster0 cells because they are not real cells
seurat.cistopic = subset(seurat.cistopic, cells = which(seurat.cistopic$peaks_snn_res.0.8 != 0) )

## normalize the gene activity matrix
DefaultAssay(seurat.cistopic) = 'RNA'
seurat.cistopic <- NormalizeData(seurat.cistopic, 
                                 #normalization.method = 'CLR',
                                 scale.factor = median(Matrix::colSums(seurat.cistopic@assays$RNA@counts))
)
seurat.cistopic <- ScaleData(seurat.cistopic) 

gene.check = 'hlh-2'
range(seurat.cistopic@assays$RNA@data[which(rownames(seurat.cistopic) == gene.check),])
range(seurat.cistopic@assays$RNA@scale.data[which(rownames(seurat.cistopic) == gene.check),])

seurat.cistopic = Manual.update.cluster.annotation(seurat.cistopic);

# test.umap.param(seurat.cistopic)
#nb.pcs = ncol(seurat.cistopic[['pca']]);
nb.pcs = 80;
n.neighbors = 50; min.dist = 0.1;

seurat.cistopic <- RunUMAP(object = seurat.cistopic, reduction = 'pca', 
                           dims = 1:nb.pcs, 
                           n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(seurat.cistopic, label = TRUE, pt.size = 2, label.size = 6,  repel = TRUE) + 
  NoLegend()

source.my.script('scATAC_functions.R')

xx = detect_doubletCell(seurat.cistopic)

#DimPlot(seurat.cistopic, label = TRUE, group.by = 'peaks_snn_res.0.8', pt.size = 2, label.size =8,  repel = TRUE) +  NoLegend()

##########################################
# 1) first identify cluster-specific marker genes
##########################################
Run.findMarkers = FALSE
if(Run.findMarkers){
  DefaultAssay(seurat.cistopic) = 'RNA'
  Idents(seurat.cistopic) = seurat.cistopic$seurat_clusters
  
  cms <- FindMarkers(seurat.cistopic, ident.1 = 6, only.pos = TRUE,
                     min.pct = 0.1, logfc.threshold = 0.1)
  
  #head(cms, n = 20)
  
  #FeaturePlot(seurat.cistopic, features = rownames(cms)[1:4])
  #FeaturePlot(seurat.cistopic, features = rownames(cms)[c(nrow(cms):(nrow(cms)-4))])
  #cms[grep('end-1', rownames(cms)), ]
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  ee.markers <- FindAllMarkers(seurat.cistopic, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
  saveRDS(ee.markers, file = paste0(RdataDir, 'scATAC_seq_geneActivityMatrix.inPromoter.resolution.0.8_markerGenes.rds'))
}

##########################################
#  ## heatmap for all marker genes
##########################################
make.heatmap.for.markerGenes = FALSE
if(make.heatmap.for.markerGenes)
{
  ee.markers = readRDS(file = paste0(RdataDir, 'scATAC_seq_geneActivityMatrix.inPromoter.resolution.0.8_markerGenes.rds'))
  
  # ee.markers$cluster = lineages[(1+as.integer(ee.markers$cluster))]
  top.markers <- ee.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
  
  #Idents(seurat.cistopic) = $lineage
  DoHeatmap(seurat.cistopic, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend() 
  
  head(as.data.frame(top.markers[which(top.markers$cluster==11), ]), 50)
}

##########################################
# heatmap for TF marker genes 
##########################################
make.heatmap.for.tf.markerGenes = FALSE
if(make.heatmap.for.tf.markerGenes)
{
  ee.markers = readRDS(file = paste0(RdataDir, 'scATAC_seq_geneActivityMatrix.inPromoter.resolution.0.8_markerGenes.rds'))
  library(openxlsx)
  tfs = openxlsx::read.xlsx('data/wTF_3.0.xlsx', sheet = 1)
  
  
  tf.markers = ee.markers[!is.na(match(ee.markers$gene, tfs$Public.name)), ]
  #tf.markers = tf.markers[tf.markers$p_val_adj<0.01, ]
  tf.markers <- tf.markers %>% group_by(cluster) %>% top_n(n = 5 , wt = avg_logFC)
  
  cat(length(unique(tf.markers$gene)), ' TFs in marker genes \n')
  library("pheatmap")
  library("RColorBrewer")
  cols = c(colorRampPalette((brewer.pal(n = 7, name="RdYlBu")))(100))
  
  DefaultAssay(seurat.cistopic) = 'RNA'
  
  seurat.cistopic <- ScaleData(seurat.cistopic, do.center = FALSE, scale.max = 2)
  gene.check = 'hlh-2'
  range(seurat.cistopic@assays$RNA@data[which(rownames(seurat.cistopic) == gene.check),])
  range(seurat.cistopic@assays$RNA@scale.data[which(rownames(seurat.cistopic) == gene.check),])
  
  groups.orders = levels(seurat.cistopic)
  new.orders = unique(c(grep('7.P0', groups.orders), 
                        grep('21.germline', groups.orders),
                        grep('intestine', groups.orders),
                        grep('hypodermis', groups.orders),
                        grep('BWM', groups.orders), 
                        grep('pharynx', groups.orders),
                        grep('9.neurons', groups.orders),
                        grep('neurons', groups.orders)))
  
  levels(seurat.cistopic) <- groups.orders[new.orders]
  
  
  pdfname = paste0(resDir, "/cluster_annotations/marker_genes_TFs_top5_log2normalized_v5.pdf")
  pdf(pdfname, width=12, height = 16)
  par(cex =1.0, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  # 
  DoHeatmap(seurat.cistopic, features = tf.markers$gene, 
            slot = 'data',
            size = 5, hjust = 0, label = TRUE, 
            disp.max = 1., 
            lines.width = 50) + scale_fill_gradientn(colors = rev(cols))
  # scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", 
  #                       high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
  #                       midpoint = 0, guide = "colourbar", aesthetics = "fill")
  dev.off()
  
  DefaultAssay(seurat.cistopic) = 'peaks'
  VlnPlot(seurat.cistopic, features = 'nFeature_peaks', idents = groups.orders, log = TRUE, pt.size = 0.3) + 
    NoLegend()
  
  VlnPlot(seurat.cistopic, features = 'nFeature_peaks', idents = groups.orders, log = TRUE, pt.size = 0.3) + 
    NoLegend()
  
  VlnPlot(seurat.cistopic, features = 'nFeature_peaks', idents = groups.orders[c(7, grep('intestine', groups.orders))], 
          sort = TRUE, log = TRUE, pt.size = 0.3) + 
    NoLegend()
  
  VlnPlot(seurat.cistopic, features = 'nFeature_peaks', idents = groups.orders[c(7, grep('hypodermis', groups.orders))], 
          sort = TRUE, log = TRUE, pt.size = 0.3) + 
    NoLegend()
  VlnPlot(seurat.cistopic, features = 'nFeature_peaks', idents = groups.orders[c(7, grep('BWM', groups.orders))], 
          sort = TRUE, log = TRUE, pt.size = 0.3) + 
    NoLegend()
  
  VlnPlot(seurat.cistopic, features = 'nFeature_peaks', idents = groups.orders[c(7, grep('pharynx', groups.orders))], 
          sort = TRUE, log = TRUE, pt.size = 0.3) + 
    NoLegend()
  
  VlnPlot(seurat.cistopic, features = 'nFeature_peaks', idents = groups.orders[c(7, grep('neurons', groups.orders))], 
          sort = TRUE, log = TRUE, pt.size = 0.3) + 
    NoLegend()
  
}

##########################################
# check lineage-specifc genes or tf or tf motif activities 
##########################################
Check.lineage.tissue.specific.markerGenes = FALSE
if(Check.lineage.tissue.specific.markerGenes){
  # import processed tintori from raw data
  tintori = import.tintori.et.al.dataset(nfeatures = 3000, npcs = 30)
  # import Aleks' scRNA-seq data
  aleks = readRDS(file = paste0(RdataDir, 'EE_Aleks_rawCounts_processed_sizefactorNormalization.rds'))
  
  DefaultAssay(seurat.cistopic) <- 'peaks'
  FeaturePlot(
    object = seurat.cistopic,
    features = "chrV:10647260-10647485", # lsy-6 downstream peaks specific for ABaxx
    pt.size = 1,
    #max.cutoff = ,
    #min.cutoff = 'q20',
    ncol = 1
  )
  
  options(warn=-1)
  
  # marker genes for intestine
  feature.example.intestine = c('sys-1','med-1', 'med-2', 'end-3', 'end-1', 'elt-2', 'elt-7', 'elt-4',
                                'nit-1', 'sdz-31', 'sdz-1',
                                'sdz-38', 'tbx-35',
                                'sdz-26', 'sdz-23',
                                'ges-1' 
  )
  
  # marker genes (main tfs) for body wall muscle
  feature.example.bwm = c('myo-3', 'atn-1', 'unc-89', 'pat-6', 'pal-1', 'hlh-8', 'hnd-1', 'hlh-1', 'unc-120')
  
  # marker genes for epidermis 
  feature.example.epidermis = c('elt-1', 'lin-26', 'elt-3', 'nhr-25', 'nhr-23', 'ajm-1', 'dlg-1', 'che-14',
                                'pal-1', 'tbx-8', 'tbx-9', 'vab-7',
                                'die-1','che-16', 
                                'egl-18', # elt-5
                                'elt-6', 'rnd-1')
  # marker genes for pharynx
  feature.example.pharynx = c('skn-1', 'med-1', 'med-2', 'tbx-35', 
                              'pha-4', 'tbx-2', 'tbx-37', 'tbx-38', 
                              
                              'glp-1', # maternal specific to AB lineage
                              'lag-1', # tf regulator of pha-4 in EMS
                              'ref-1', # 'hlh-25', 'hlh-26', 'hlh-27', 'hlh-28', 'hlh-29', # target of lag-1
                              'sem-2', #'sox-1' # med-1/2 target
                              'pha-1', 'ceh-24',
                              'ceh-22','peb-1','daf-12', 'myo-2'
  )
  feature.example.neurons = c('lin-22', 'lin-26', 'cpb-1', 'hlh-14',  'glp-1', # genes determining neuron or non-neurons
                              'unc-86', 'nhr-67', 'egl-5', 'lin-32', 'cnd-1', 'hlh-2', 'pag-3', # genes required for the type of neuronal lineage program
                              'daf-19', # ciliated sensory neurons specification
                              'cnd-1' # vnc motorneuron specification
                              #'unc-25', 'unc-47', 'unc-30', 'lim-6', # GABAergic neurons
                              #'ceh-10', 'ttx-3', 'che-1', 'ceh-36', 'unc-30', 'ceh-14', 'unc-42', 'mex-3',
                              #'unc-119', 'unc-33', 'hen-1', 'unc-17', 'sra-11'
  )
  
  feature.example.germline = c('pgl-1', 'glh-1', 'glh-4', 'nos-1', 'htp-3', 'rec-8',
                               'nos-2', # maternal mRNA maintained and translated
                               'par-1', 'par-2',
                               'mex-5', 'mex-6'
  )
  
  #feature.example = feature.example.intestine
  feature.example = c('end-1', 'end-3', 'elt-2', 'elt-7')
  feature.example = c('skn-1', 'med-1', 'med-2', 'tbx-35', 'che-51', 'end-1', 'end-3', 'elt-2', 'elt-7')
  #
  #feature.example = c('pal-1', 'hnd-1', 'hlh-1', 'unc-120')
  
  # feature.example = feature.example.pharynx
  
  feature.example = c('mec-3', 'ceh-43', 'hlh-14', 'unc-30', 'ceh-36')
  DefaultAssay(seurat.cistopic) <- 'RNA'
  FeaturePlot(seurat.cistopic, features = feature.example, ncol = 3, reduction = 'umap')
  
  DefaultAssay(seurat.cistopic) <- 'chromvar'
  
  motifs2check = rownames(seurat.cistopic)[grep(paste0(feature.example, collapse = '|'), rownames(seurat.cistopic))]
  
  if(length(motifs2check)>0){
    FeaturePlot(
      object = seurat.cistopic,
      features = motifs2check,
      pt.size = 0.1,
      #cols = c('red', 'orange', 'blue'),
      #max.cutoff = 1,
      min.cutoff = 0,
      ncol = 4
      #cols =  cols
    )
  }else{
    cat('no motif found !\n')
  }
  
  
}



########################################################
########################################################
# OLD Section : scATAC-seq cluster annotation
# 
########################################################
########################################################

##########################################
# using marker genes
##########################################
gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.500bpDownstream.rds')
seurat.cistopic = readRDS(file = gamat)
#seurat.cistopic = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object_geneBody.promoter.activityscores.rds'))

#source.my.script('scATAC_functions.R')
#seurat.cistopic = annotate.scATAC.clusters.with.marker.genes(seurat.cistopic)
  
##########################################
# label transferring from scRNA-seq data using liger and seurat 
# see details in https://satijalab.org/signac/articles/mouse_brain_vignette.html#integrating-with-scrna-seq-data
##########################################

#seurat.cistopic = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object_promoterOnly.activityscores.rds'))
#seurat.cistopic = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object_geneBody.promoter.activityscores.rds'))
gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.500bpDownstream.rds')
seurat.cistopic = readRDS(file =  gamat)

source.my.script('scATAC_functions.R')
transfer.labels.from.scRNA.to.scATAC(seurat.cistopic, tintori, method = 'liger')






