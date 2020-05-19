##########################################################################
##########################################################################
# Project: single-cell analysis for scATAC-seq 
# Script purpose: annotate clusters
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri May  8 15:24:55 2020
##########################################################################
##########################################################################
Manual.update.cluster.annotation = function(seurat.cistopic)
{
  DefaultAssay(seurat.cistopic) <- 'peaks'
  Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
  new.cluster.ids <- c("1.neurons", 
                       "2.neurons", 
                       "3.neurons", 
                       "4.hypodermis", 
                       "5.intestine",
                       '6.intestine', 
                       '7.P0/AB/P1/ABa.p/EMS/P2',
                       '8.hypodermis', 
                       '9.neurons/BWM', 
                       '10.pharynx/neurons', 
                       '11.BWM', 
                       '12.hypodermis', 
                       '13.pharynx/neurons', 
                       '14.pharynx/neurons', 
                       '15.BWM',
                       '16.pharynx/neurons', 
                       '17.neurons', 
                       '18.neurons',
                       '19.pharynx/neurons', 
                       '20.hypodermis', 
                       '21.germline.P4')
  
  names(new.cluster.ids) <- levels(seurat.cistopic)
  seurat.cistopic <- RenameIdents(seurat.cistopic, new.cluster.ids)
  
  #DimPlot(seurat.cistopic, reduction = "umap", label = TRUE, pt.size = 2,  label.size = 5, repel = TRUE) + NoLegend()
  
  return(seurat.cistopic)
  
}


##########################################
# integrate scRNA-seq data from published dataset in https://github.com/qinzhu/VisCello.celegans
# Packer, J. S J. I. Murray (2019). A lineage-resolved molecular atlas of C. elegans embryogenesis at single-cell resolution. Science: eaax1971.
# here we are using Seurat label transferring method
##########################################
##########################################
compare.gene.activity.matrix.computing = function()
{
  ##########################################
  # there are different ideas to compute the gene activity matrix for scATAC-seq data
  # here we compare 1) the native promoter(-2kb, +0.5/2kb), promoter.genebody (-2kb) method,
  # 2) promoter, promoter.genebody with only integenetic regions
  # 3) probably also cicero method
  ##########################################
  gam.file = list.files(path = 'results/scATAC_earlyEmbryo_20200302/Rdata', pattern = '*.rds', full.names = TRUE)
  gam.file = gam.file[grep('atac_LDA_seurat_object_geneActivityMatrix', gam.file)]
  gam.file = gam.file[grep('promoter.geneBody_2000bp|promoter_2000bpUpstream.2000bpDownstream|promoter_2000bpUpstream.500bpDownstream', 
                           gam.file)]
  
  feature.example = unique(c('pha-4', 'hnd-1', 'skn-1', 'hlh-1', 'tbx-8', 'unc-120', 'med-2', 'med-1', 'sdz-38', 'pal-1',
                      'end-1', 'elt-2', 'end-1', 'end-3', 'elt-7',
                      'hnd-1', 'nrh-25', 
                      'pie-1', 'mex-3', 'lin-53', 'gld-1','ceh-32', 'myo-2', 'cnd-1', 'ceh-32', 'tbx-37', 'tbx-38', 'tbx-8',
                      'erm-1', # AB enriched
                      'cpg-2',
                      'med-2', # marker for EMS
                      'era-1', #'W02F12.3',
                      'tbx-35',
                      'end-3',
                      'mex-5',
                      'pigv-1', #'T09B4.1'
                      'ceh-51',
                      'end-1',
                      'pal-1',
                      'hlh-27',
                      'ref-1',
                      'tbx-38',
                      'sdz-38',
                      'pha-4', 'hnd-1',
                      'lsy-6', # marker for ABaxx
                      'tbx-37','tbx-38', 'ceh-27',
                      'ceh-36', 'elt-1',
                      'dmd-4', 
                      'elt-7', 'end-1', 'nrh-79', 
                       "med-1","med-2", 'pal-1', 'pha-4', 
                      'pes-1', 'nos-1', 'nos-2'))

  for(gamat in gam.file){
    
    gamat = gam.file[grep('nonOverlapped.withGenes_promoter_2000bpUpstream.500bpDownstream', gam.file)]
    seurat.cistopic = readRDS(file = gamat)
    
    cat('gene activity matrix -- ', gamat, '\n')
    Manual.modify.cluster.annotation = TRUE
    if(Manual.modify.cluster.annotation){
      Manual.modify.cluster.annotation(seurat.cistopic)
    }
    
    # remove cluster0 cells because they are not real cells
    seurat.cistopic = subset(seurat.cistopic, cells = which(seurat.cistopic$peaks_snn_res.0.8 != 0) )
    
    # pdfname = paste0(resDir, "/cluster_annotations/markerGenes_examples_", gsub('.rds', '', basename(gamat)), ".pdf")
    # pdf(pdfname, width=10, height = 8)
    # par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    # 
    DefaultAssay(seurat.cistopic) <- 'peaks'
    # test.umap.param(seurat.cistopic)
    #nb.pcs = ncol(seurat.cistopic[['pca']]);
    nb.pcs = 80;
    n.neighbors = 30; min.dist = 0.1;
    
    seurat.cistopic <- RunUMAP(object = seurat.cistopic, reduction = 'pca', 
                               dims = 1:nb.pcs, 
                               n.neighbors = n.neighbors, min.dist = min.dist)
    p1 = DimPlot(seurat.cistopic, label = TRUE, pt.size = 1, label.size = 6) + 
      NoLegend()
    plot(p1)
    
    Normalize.Gene.Activity.Matrix = TRUE
    if(Normalize.Gene.Activity.Matrix){
      DefaultAssay(seurat.cistopic) = 'RNA'
      seurat.cistopic <- NormalizeData(seurat.cistopic, 
                                       #normalization.method = 'CLR',
                                       scale.factor = median(Matrix::colSums(seurat.cistopic@assays$RNA@counts))
      )
      seurat.cistopic <- ScaleData(seurat.cistopic)
    }
    
    f = c('mes-2', 'mes-3', 'mes-6', 'med-4')
    f = c('mep-1', 'let-418', 'hda-1', 'apx-1')
    f = c('ref-1', 'hlh-25', 'pop-1', 'glp-1', 'lag-2', 'pha-4', 'lin-12', 'tbx-37', 'tbx-38')
    f = c('che-1', 'unc-3', 'che-14', 'unc-86', 'unc-42')
    FeaturePlot(seurat.cistopic, features = f, ncol = 2, pt.size = 0.7)
    
    for(f in feature.example){
      cat(f, '\n')
      if(length(which(rownames(seurat.cistopic) == f))>0)
        plot(FeaturePlot(seurat.cistopic, features = f, ncol = 1, pt.size = 1))
    }
    
    #dev.off()
    
  }
  
  ##########################################
  #  the following cod is to test the influence of nb of features into the UMAP in gene activity matrix 
  ##########################################
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneBody.promoter.activityscores.rds')
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.2000bpDownstream.rds')
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_promoterOnly.activityscores.rds')
  
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.200bpDownstream.rds')
  
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.2000bpDownstream.rds')
  
  gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.500bpDownstream.rds')
  seurat.cistopic = readRDS(file =  gamat)
  #seurat.cistopic = readRDS(file =  paste0(RdataDir, 'atac_LDA_seurat_object_geneBody.promoter.activityscores.rds'))
  
  DefaultAssay(seurat.cistopic) <- 'RNA'
  
  hist(log2(as.numeric(GetAssayData(seurat.cistopic, slot = 'counts'))+1), breaks = 10000)
  
  pdfname = paste0(resDir, "/UMAP_geneActivityMatrix_promoter.genebody_test_nbHVGs.pdf")
  pdf(pdfname, width=16, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(nb.variableFeatures in c(seq(1000, 5000, by = 1000), 8000, 10000, 12000, 15000, nrow(seurat.cistopic)))
  {
    # nb.variableFeatures = 20000
    cat('nb of HVGs -- ', nb.variableFeatures, '\n')
    seurat.cistopic <- FindVariableFeatures(seurat.cistopic, nfeatures = nb.variableFeatures)
    seurat.cistopic <- NormalizeData(seurat.cistopic, 
                                     #normalization.method = 'CLR',
                                     scale.factor = median(Matrix::colSums(seurat.cistopic@assays$RNA@counts)))
    seurat.cistopic <- ScaleData(seurat.cistopic)
    
    labels = Idents(seurat.cistopic)
    Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
    
    seurat.cistopic <- RunPCA(seurat.cistopic, npcs = 50, verbose = FALSE, reduction.name = 'pca.ga')
    
    nb.pcs = 50; n.neighbors = 20; min.dist = 0.1;
    seurat.cistopic <- RunUMAP(seurat.cistopic, reduction = "pca.ga", dims = 1:nb.pcs,
                               n.neighbors = n.neighbors, min.dist = min.dist, reduction.name = "umap.ga")
    p1 = DimPlot(seurat.cistopic, reduction = "umap.ga", label = TRUE, pt.size = 2, label.size = 5, repel = FALSE) + 
      ggtitle(paste0('nb of HVGs - ', nb.variableFeatures))
    plot(p1)
  }
  
  
  
  dev.off()
}


transfer.labels.from.scRNA.to.scATAC = function(seurat.cistopic, tintori, method = c('seurat', 'liger'))
{
  if(method == 'seurat'){
    
    #gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.500bpDownstream.rds')
    gamat = paste0(RdataDir, 'atac_LDA_seurat_object_geneActivityMatrix_seurat_promoter_2000bpUpstream.2000bpDownstream.rds')
    seurat.cistopic = readRDS(file =  gamat)
    
    # drop cells in cluster 0
    seurat.cistopic = subset(seurat.cistopic, cells = which(seurat.cistopic$peaks_snn_res.0.8 != 0) )
    
    # process.tintori.et.al()
    # processed raw data
    tintori = readRDS(file = paste0(RdataDir, 'Tintori.et.al_rawCounts_processed_sizefactorNormalization.rds')) 
    #xx <- readRDS(paste0(RdataDir, 'Tintori_et_al_highQualtiyCells.rds')) # processed data from cello 
    tintori = tintori[!is.na(match(rownames(tintori), rownames(seurat.cistopic)))] # select only protein-coding genes and miRNAs
    
    tintori <- FindVariableFeatures(
      object = tintori,
      nfeatures = 3000
    )
    
    tintori <- ScaleData(object = tintori)
    tintori <- RunPCA(object = tintori, npcs = 50, verbose = FALSE)
    
    Idents(tintori) = tintori$lineage
    tintori <- RunUMAP(object = tintori, dims = 1:30, n.neighbors = 10, min.dist = 0.2)
    DimPlot(tintori, reduction = "umap", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
    
    resolution = 0.8;
    
    pdfname = paste0(resDir, "/label_transferring_seurat_resolution_rawTintoriLineage_", resolution,  ".pdf")
    pdf(pdfname, width=12, height = 8)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    DefaultAssay(seurat.cistopic) <- 'peaks'
    seurat.cistopic = FindClusters(object = seurat.cistopic, 
                                   n.start=20, resolution=resolution,
                                   algorithm = 3)
    #seurat.cistopic <- RunPCA(seurat.cistopic, npcs = 50, verbose = FALSE, reduction.name = 'pca.ga',  features = tops$gene)
    
    Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
    #nb.pcs = ncol(seurat.cistopic[['pca']]); n.neighbors = 20; min.dist = 0.2;
    nb.pcs = ncol(seurat.cistopic[['pca']]); n.neighbors = 30; min.dist = 0.3;
    seurat.cistopic <- RunUMAP(seurat.cistopic, reduction = "pca", dims = 1:nb.pcs,
                               n.neighbors = n.neighbors, min.dist = min.dist, reduction.name = "umap")
    DimPlot(seurat.cistopic, reduction = "umap", label = TRUE, pt.size = 2,  label.size = 8, repel = FALSE) + NoLegend()
    
    # Here, we process the gene activity matrix 
    # in order to find anchors between cells in the scATAC-seq dataset 
    # and the scRNA-seq dataset.
    DefaultAssay(seurat.cistopic) <- 'RNA'
    nb.variableFeatures = 3000
    
    seurat.cistopic <- FindVariableFeatures(seurat.cistopic, nfeatures = nb.variableFeatures)
    seurat.cistopic <- NormalizeData(seurat.cistopic)
    seurat.cistopic <- ScaleData(seurat.cistopic)
    
    transfer.anchors <- FindTransferAnchors(
      reference = seurat.cistopic,
      query = tintori,
      features = unique(c(VariableFeatures(object = tintori), VariableFeatures(seurat.cistopic))),
      #features = features.to.use,
      reference.assay = 'RNA',
      query.assay = 'RNA',
      reduction = 'cca',
      k.anchor = 30, # k.anchor is neighborhood size for MNN big k.anchor, the bigger, the more anchors found
      k.filter = 150, # retain the anchor (cell from one dataset to annother) if within k.filter neighbors, the bigger, the more retained  
      max.features = 5000, # max nb of features used for anchor filtering
      k.score = 30, 
      npcs = 50, 
      dims = 1:50
    )
    
    cat('nb of cells in scATAC and in tintori as anchors : ', 
        length(unique(transfer.anchors@anchors[, 1])), '--',  length(unique(transfer.anchors@anchors[, 2])), '\n')
    
    Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
    predicted.labels <- TransferData(
      anchorset = transfer.anchors,
      #refdata = Idents(tintori),
      refdata = Idents(seurat.cistopic),
      #refdata = as.vector(Idents(seurat.cistopic)),
      #weight.reduction = seurat.cistopic[['pca']],
      weight.reduction = tintori[['pca']],
      dims = 1:30,
      k.weight = 50
    )
    
    tintori <- AddMetaData(object = tintori, metadata = predicted.labels)
    
    p1 = DimPlot(tintori, group.by = "predicted.id", reduction = 'umap', pt.size = 2,  label.size = 6,
                 label = TRUE, repel = TRUE) + ggtitle("scRNA-seq cells") + 
      scale_colour_hue(drop = FALSE)
    p2 = DimPlot(tintori,reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,  label.size = 6) + 
      ggtitle("scRNA-seq cells") + 
      scale_colour_hue(drop = FALSE)
    
    p1 + p2
    
    res = data.frame(lineage = Idents(tintori), predicted.labels, stringsAsFactors = FALSE)
    
    map = res[, match(paste0('prediction.score.', levels(seurat.cistopic)), colnames(res))]
    rownames(map) = sapply(rownames(map),function(x) gsub('-cell_r*', '',x)) 
    map = t(map)
    
    # process the rowname and colnames
    rownames(map) = sapply(rownames(map), function(x) gsub('prediction.score.', '', x))
    colnames(map) = sapply(colnames(map), function(x) paste0(unlist(strsplit(x, '_'))[c(2,1)], collapse = "_"))
    
    celltypes = sapply(colnames(map), function(x) unlist(strsplit(x, '_'))[1])

    # manually sort the cell types
    #cellrefs = levels(tintori)
    cellrefs = c("P0",  "AB", "P1", "ABa",  "ABp" , "EMS",   "P2",
                  "ABal",  "ABar","ABpl", "ABpr","MS",  "E" , 
                 "ABalx", "ABarx", "ABplx", "ABprx",  "MSx1", "MSx2",  "Cx1" ,"Cx2", "C",  "P3", "D", "P4", "Ea", "Ep")   
    
    kk = c()
    for(cc in cellrefs)
    {
      kk = c(kk, which(celltypes==cc))
    }
    
    map = map[, kk]
    my_sample_col <- data.frame(sample = celltypes[kk])
    rownames(my_sample_col) <- colnames(map)
    #colnames(map) = celltypes[kk]
    
    library("pheatmap")
    library("RColorBrewer")
    library(grid)
    
    pdfname = paste0(resDir, "/cluster_annotations/label_transferring_seurat_resolution_rawTintoriLineage_", 
                     resolution,  "_dropCluster0_v3.pdf")
    pdf(pdfname, width=24, height = 10)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    cols = c(colorRampPalette((brewer.pal(n = 7, name="Reds")))(10))
    pheatmap(map, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, breaks = seq(0, 1, by = 0.1),
             cluster_cols=FALSE, main = paste0("resolution -- ",  resolution), na_col = "white",
             color = cols, 
             annotation_col = my_sample_col,
             gaps_col = match(unique(my_sample_col$sample), my_sample_col$sample)[-1] -1,
             gaps_row = c(1:nrow(map)-1),
             fontsize_col = 10,
             height = 8,
             width = 30
    )
    
    #grid.text(levels(my_sample_col), x=c(12,0.6),y=c(12,0.89), gp=gpar(fontsize=10))
    dev.off()
    
    map.fitered = map
    map.fitered[which(map.fitered<0.1)] = 0
    pheatmap(map.fitered, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
             cluster_cols=FALSE, main = paste0("fitlered < 0.1 with resolution -- ",  resolution), na_col = "white",
             color = cols)
    
    dev.off()
    
    # Manual.modify.cluster.annotation(seurat.cistopic)
    
    # xx = as.matrix(map)
    # pheatmap(xx, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
    #          cluster_cols=FALSE, main = paste0("alpha = ",  " -- "), na_col = "white",
    #          color = cols)
    #table(seurat.cistopic$prediction.score.max > 0.5)
    
    # 
    # seurat.cistopic <- AddMetaData(object = seurat.cistopic, metadata = predicted.labels)
    # DimPlot(seurat.cistopic, group.by = "predicted.id",reduction = 'umap', label = TRUE, repel = FALSE) + ggtitle("scATAC-seq cells") + 
    #   scale_colour_hue(drop = FALSE)
    # 
    # table(seurat.cistopic$prediction.score.max > 0.5)
    
    hist(seurat.cistopic$prediction.score.max)
    abline(v = 0.5, col = "red")
    
    seurat.cistopic.filtered <- subset(seurat.cistopic, subset = prediction.score.max > 0.5)
    # to make the colors match
    seurat.cistopic.filtered$predicted.id <- factor(seurat.cistopic.filtered$predicted.id, levels = levels(tintori))  
    DimPlot(seurat.cistopic.filtered, group.by = "predicted.id", label = TRUE, repel = TRUE) + ggtitle("scATAC-seq cells") + 
      NoLegend() + scale_colour_hue(drop = FALSE)
    
    
  }
  
  if(method == 'liger'){
    transfer.labels.from.scRNA.to.scATAC.with.liger()
  }
  
}

transfer.labels.from.Murray.scRNA.to.scATAC.seurat. = function(seurat.cistopic)
{
  # process Parker et al. data
  library(VisCello.celegans)
  eset = readRDS(file = paste0(RdataDir, 'Parker_et_al_dataSet_afterFiltering_89701cell.rds'))
  pmeda = data.frame(pData(eset))
  
  
  ## select the cells for early embryo
  kk = which(pmeda$embryo.time < 270 & pmeda$cell.type != 'unannotated' & pmeda$cell.type != 'XXX')
  cat('nb of cell in reference -- ', length(kk), '\n')
  
  ee = CreateSeuratObject(counts = eset@assayData$exprs[,kk], assay = 'RNA', meta.data = pmeda[kk, ])
  ee@assays$RNA@data = eset@assayData$norm_exprs[,kk]
  
  ee <- FindVariableFeatures(
    object = ee,
    nfeatures = 3000
  )
  
  ee <- ScaleData(object = ee)
  
  Idents(ee) = ee$cell.type
  
  # Here, we process the gene activity matrix 
  # in order to find anchors between cells in the scATAC-seq dataset 
  # and the scRNA-seq dataset.
  Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
  DimPlot(seurat.cistopic, reduction = "umap", label = TRUE, pt.size = 2,  label.size = 8, repel = FALSE) + NoLegend()
  
  DefaultAssay(seurat.cistopic) <- 'RNA'
  nb.variableFeatures = 5000
  seurat.cistopic <- FindVariableFeatures(seurat.cistopic, nfeatures = nb.variableFeatures)
  seurat.cistopic <- NormalizeData(seurat.cistopic)
  seurat.cistopic <- ScaleData(seurat.cistopic)
  
  seurat.cistopic <- RunPCA(seurat.cistopic, npcs = 50, verbose = FALSE, reduction.name = 'pca.ga')
  
  transfer.anchors <- FindTransferAnchors(
    reference = ee,
    query = seurat.cistopic,
    features = unique(c(VariableFeatures(object = ee), VariableFeatures(seurat.cistopic))),
    #features = features.to.use,
    reference.assay = 'RNA',
    query.assay = 'RNA',
    reduction = 'cca',
    k.anchor = 10, # k.anchor is neighborhood size for MNN big k.anchor, the bigger, the more anchors found
    k.filter = 200, # retain the anchor (cell from one dataset to annother) if within k.filter neighbors, the bigger, the more retained  
    max.features = 200, # max nb of features used for anchor filtering
    k.score = 30, 
    npcs = 50, 
    dims = 1:50
  )
  
  cat('nb of cells in scATAC and in tintori as anchors : ', 
      length(unique(transfer.anchors@anchors[, 1])), '--',  length(unique(transfer.anchors@anchors[, 2])), '\n')
  
  #Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
  
  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    #refdata = Idents(tintori),
    refdata = Idents(ee),
    #refdata = as.vector(Idents(seurat.cistopic)),
    #weight.reduction = seurat.cistopic[['pca']],
    weight.reduction = seurat.cistopic[['pca.ga']],
    dims = 1:30,
    k.weight = 50
  )
  
  seurat.cistopic <- AddMetaData(object = seurat.cistopic, metadata = predicted.labels)
  DimPlot(seurat.cistopic, group.by = "predicted.id",reduction = 'umap', label = TRUE, repel = FALSE) + ggtitle("scATAC-seq cells") +
    scale_colour_hue(drop = FALSE)
  
  table(seurat.cistopic$prediction.score.max > 0.5)
  
  hist(seurat.cistopic$prediction.score.max)
  abline(v = 0.5, col = "red")
  
  seurat.cistopic.filtered <- subset(seurat.cistopic, subset = prediction.score.max > 0.5)
  
  # to make the colors match
  seurat.cistopic.filtered$predicted.id <- factor(seurat.cistopic.filtered$predicted.id, levels = levels(tintori))  
  DimPlot(seurat.cistopic.filtered, group.by = "predicted.id", reduction = 'umap',  label = TRUE, repel = TRUE) + 
    ggtitle("scATAC-seq cells") + 
    scale_colour_hue(drop = FALSE)
  
  ##########################################
  # make summary of cluster-to-predicted label mapping
  ##########################################
  res = data.frame(clusters = Idents(seurat.cistopic), predicted.labels, stringsAsFactors = FALSE)
  
  labels.pred =  unique(res$predicted.id)
  clusters = unique(res$clusters)
  res.map = matrix(0, nrow = length(labels.pred), ncol = length(clusters))
  rownames(res.map) = labels.pred
  colnames(res.map) = clusters
  res.map = res.map[, sort(colnames(res.map))]
  
  for(n in 1:ncol(res.map))
  {
    # n = 1
    predicted.ids = res$predicted.id[which(res$clusters == colnames(res.map)[n])]
    stats.ids = table(predicted.ids)
    stats.ids = stats.ids[match(rownames(res.map), names(stats.ids))]
    stats.ids[is.na(stats.ids)] = 0
    res.map[,n] = stats.ids/sum(stats.ids)
  }
  
  library("pheatmap")
  library("RColorBrewer")
  library(grid)
  
  cols = c(colorRampPalette((brewer.pal(n = 7, name="Reds")))(10))
  pheatmap(res.map, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, breaks = seq(0, 1, by = 0.1),
           cluster_cols=TRUE, main = paste0("resolution -- 0.8"), na_col = "white",
           color = cols, 
           #annotation_col = my_sample_col,
           #gaps_row = c(1:nrow(map)-1),
           fontsize_col = 10,
           height = 8,
           width = 30
  )
  
  cols =  c(colorRampPalette((brewer.pal(n = 7, name="Set3")))(ncol(res.map)))
  barplot(res.map, horiz = TRUE, legend = rownames(res.map), 
          col = cols )
  
  map.fitered = map
  map.fitered[which(map.fitered<0.1)] = 0
  pheatmap(map.fitered, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
           cluster_cols=FALSE, main = paste0("fitlered < 0.1 with resolution -- ",  resolution), na_col = "white",
           color = cols)
  
}


process.tintori.et.al = function()
{
  # DefaultAssay(seurat.cistopic) <- 'peaks'
  # seurat.cistopic <- RunTFIDF(seurat.cistopic)
  # seurat.cistopic <- FindTopFeatures(seurat.cistopic, min.cutoff = 'q25')
  # seurat.cistopic <- RunSVD(
  #   object = seurat.cistopic,
  #   assay = 'peaks',
  #   reduction.key = 'LSI_',
  #   reduction.name = 'lsi'
  # )
  
  #VariableFeatures(seurat.cistopic) <- names(which(Matrix::rowSums(seurat.cistopic) > 100))
  #seurat.cistopic <- RunLSI(seurat.cistopic, n = 50, scale.max = NULL)
  #seurat.cistopic <- RunUMAP(seurat.cistopic, reduction = "lsi", dims = 1:50)
  
  # Here, we process the gene activity matrix 
  # in order to find anchors between cells in the scATAC-seq dataset 
  # and the scRNA-seq dataset.
  DefaultAssay(seurat.cistopic) <- 'RNA'
  nb.variableFeatures = nrow(seurat.cistopic)
  seurat.cistopic <- FindVariableFeatures(seurat.cistopic, nfeatures = nb.variableFeatures)
  seurat.cistopic <- NormalizeData(seurat.cistopic)
  seurat.cistopic <- ScaleData(seurat.cistopic)
  
  Refine.HGV.Gene.Activity = FALSE
  if(Refine.HGV.Gene.Activity){
    
    labels = Idents(seurat.cistopic)
    Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
    
    seurat.cistopic <- RunPCA(seurat.cistopic, npcs = 50, verbose = FALSE, reduction.name = 'pca.ga')
    
    nb.pcs = 30; n.neighbors = 20; min.dist = 0.3;
    seurat.cistopic <- RunUMAP(seurat.cistopic, reduction = "pca.ga", dims = 1:nb.pcs,
                               n.neighbors = n.neighbors, min.dist = min.dist, reduction.name = "umap.ga")
    DimPlot(seurat.cistopic, reduction = "umap.ga", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
    
    
    # find all markers of cluster 1
    markers.ga <- FindAllMarkers(seurat.cistopic, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    tops <- markers.ga %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
    #length(intersect(tops$gene, VariableFeatures(seurat.cistopic)))
    #VariableFeatures(seurat.cistopic) = tops$gene
    seurat.cistopic <- RunPCA(seurat.cistopic, npcs = 50, verbose = FALSE, reduction.name = 'pca.ga',  features = tops$gene)
    
    nb.pcs = 30; n.neighbors = 20; min.dist = 0.3;
    seurat.cistopic <- RunUMAP(seurat.cistopic, reduction = "pca.ga", dims = 1:nb.pcs,
                               n.neighbors = n.neighbors, min.dist = min.dist, reduction.name = "umap.ga")
    DimPlot(seurat.cistopic, reduction = "umap.ga", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
    
  }
  
  # Load the pre-processed scRNA-seq data 
  tintori <- readRDS(paste0(RdataDir, 'Tintori_et_al_highQualtiyCells.rds'))
  tintori = tintori[!is.na(match(rownames(tintori), rownames(seurat.cistopic)))] # select only protein-coding genes and miRNAs
  
  tintori <- FindVariableFeatures(
    object = tintori,
    nfeatures = 5000
  )
  
  DimPlot(tintori, reduction = "umap", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
  
  Merge.cell.types.tintori = FALSE
  if(Merge.cell.types.tintori){
    #Idents(seurat.cistopic) = seurat.cistopic$peaks_snn_res.0.8
    new.cluster.ids <- c("P0/AB/P1", "P0/AB/P1", "P0/AB/P1","ABa/p/EMS", "ABa/p/EMS", "ABa/p/EMS",
                         "P2/P3/C","ABa/pr/l",  "ABa/pr/l",  "ABa/pr/l", "ABa/pr/l",
                         "P2/P3/C", "MS/E","P2/P3/C", "MS/E",
                         "ABa/pr/lx", "ABa/pr/lx", "ABa/pr/lx", "ABa/pr/lx",
                         "MSx1/2", "MSx1/2",  "Cx1/2",   
                         "Ea/p",    "Ea/p",    "Cx1/2",   "P4/D",    "P4/D" )   
    names(new.cluster.ids) <- levels(tintori)
    tintori <- RenameIdents(tintori, new.cluster.ids)
    
    DimPlot(tintori, reduction = "umap", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
  }
  
  DimPlot(tintori, reduction = "umap", label = TRUE, pt.size = 2, label.size = 5, repel = TRUE)
  
}

##########################################
# test liger to transfer lables 
# see details in https://macoskolab.github.io/liger/walkthrough_rna_atac.html
##########################################
transfer.labels.from.scRNA.to.scATAC.with.liger = function(seurat.cistopic, tintori)
{
  
  # rna_clusts = readRDS("../liger-rna-atac-vignette/rna_cluster_assignments.RDS")
  # atac_clusts = readRDS("../liger-rna-atac-vignette/atac_cluster_assignments.RDS")
  # pbmc.atac <- readRDS('../liger-rna-atac-vignette/pbmc.atac.expression.mat.RDS')
  # pbmc.rna <- readRDS('../liger-rna-atac-vignette/pbmc.rna.expression.mat.RDS')
  
  # sels = sample(c(1:length(rna_clusts)), 200, replace = FALSE)
  # rna_clusts = rna_clusts[sels]
  #pbmc.rna = pbmc.rna[, sels]
  
  atac_clusts = seurat.cistopic$peaks_snn_res.0.8
  rna_clusts = Idents(tintori)
  pbmc.atac = GetAssayData(seurat.cistopic, slot = 'counts', assay = 'RNA')
  pbmc.rna = GetAssayData(tintori, slot = 'counts', assay = 'RNA')
  
  # data preprocessing
  library(liger)
  
  convertSeurt_toLiger = FALSE
  if(convertSeurt_toLiger){
    ligerex2 <- seuratToLiger(list(seurat.cistopic, tintori), assays.use = c('RNA', 'RNA'), num.hvg.info = 2500, renormalize = TRUE)
    
    cat('nb of HVGs : ', length(ligerex2@var.genes), '\n')
    int.pbmc = ligerex2
    
  }else{
    ggplot2::theme_set(theme_cowplot())
    #xx = pbmc.atac[,names(atac_clusts)]
    pbmc.data = list(atac=pbmc.atac[,names(atac_clusts)], rna=pbmc.rna[,names(rna_clusts)])
    int.pbmc <- createLiger(pbmc.data)
    
    int.pbmc <- liger::normalize(int.pbmc)
    
    par(mfrow=c(1,2))
    int.pbmc <- selectGenes(int.pbmc, datasets.use = c(1, 2), 
                            num.genes = c(20000, 5000), do.plot = TRUE)
    cat('nb of HVGs : ', length(unique(int.pbmc@var.genes)), '\n')
    
  }
  
  int.keep = int.pbmc
  
  pdfname = paste0(resDir, "/Test_liger_lambda.pdf")
  pdf(pdfname, width=16, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(lambda in c(1:4, seq(5, 100, by = 5)))
  {
    lambda = 50
    int.pbmc <- scaleNotCenter(int.keep)
    
    Suggest.k.lambda = FALSE
    if(Suggest.k.lambda){
      suggestK(int.pbmc) # plot entropy metric to find an elbow that can be used to select the number of factors
      suggestLambda(int.pbmc, k = 20, num.cores = 5) # plot alignment metric to find an elbow that can be used to select the value of lambda
    }
    
    # Factorization and Quantile Normalization
    int.pbmc <- optimizeALS(int.pbmc, k=30, lambda = lambda, nrep = 1)
    
    #int.pbmc <- liger::runUMAP(int.pbmc, use.raw = T)
    int.pbmc <- liger::quantile_norm(int.pbmc)
    
    # visualization
    int.pbmc <- liger::runUMAP(int.pbmc, use.raw = FALSE,  dims.use = 1:30, distance = 'euclidean', n_neighbors = 20, min_dist = 0.2)
    
    plots1 <- plotByDatasetAndCluster(int.pbmc, return.plots = T, clusters=rna_clusts) 
    #print(plots1[[1]])
    p1 =  print(plots1[[2]])
    plots2 <- plotByDatasetAndCluster(int.pbmc, return.plots = T, clusters=atac_clusts) 
    #print(plots2[[1]])
    p2 = print(plots2[[2]])
    
    plot(p1 + p2) + ggtitle(paste0('lambda -', lambda))
    
  }
  
  dev.off()
  
  calcAlignment(int.pbmc)
  
  xx = ligerToSeurat(int.pbmc, use.liger.genes = T)
  
  # NKG7 = plotGene(int.pbmc,gene="NKG7",return.plots=T)
  # MS4A1 = plotGene(int.pbmc,gene="MS4A1",return.plots=T)
  # plot_grid(NKG7[[2]],MS4A1[[2]],NKG7[[1]],MS4A1[[1]],ncol=2)
  
  #seurat_liger = ligerToSeurat(int.pbmc, use.liger.genes = T)
  
  
}

########################################################
########################################################
# Section : functions for scATAC-seq cluster annotation
# 
########################################################
########################################################
annotate.scATAC.clusters.with.marker.genes = function(seurat.cistopic)
{
  
}

import.tintori.et.al.dataset = function(nfeatures = 5000, npcs = 30, n.neighbors = 10, min.dist = 0.2)
{
  tintori = readRDS(file = paste0(RdataDir, 'Tintori.et.al_rawCounts_processed_sizefactorNormalization.rds')) 
  tintori = tintori[!is.na(match(rownames(tintori), rownames(seurat.cistopic)))] # select only protein-coding genes and miRNAs
  
  tintori <- FindVariableFeatures(
    object = tintori,
    nfeatures = 5000
  )
  
  tintori <- ScaleData(object = tintori)
  tintori <- RunPCA(object = tintori, npcs = 50, verbose = FALSE)
  
  Idents(tintori) = tintori$lineage
  tintori <- RunUMAP(object = tintori, dims = 1:30, n.neighbors = 10, min.dist = 0.2)
  
  #DimPlot(tintori, reduction = "umap", label = TRUE, pt.size = 2, label.size = 5, repel = FALSE)
  
  return(tintori)
}

plot.feature.example.scRNA.scATAC.tf.motif = function(seurat.cistopic, feature.example, tintori, aleks)
{
  plots = list()
  if(length(which(rownames(tintori) == feature.example))>0){
    p01 = FeaturePlot(tintori, features = feature.example, ncol = NULL, by.col = FALSE,
                      label = TRUE, label.size = 4, pt.size = 2, 
                      reduction = 'umap', repel = TRUE) + 
      ggtitle(paste0('tintori - ' , feature.example))
    plots[['tintori']] = p01
  }  
  if(length(which(rownames(aleks) == feature.example))>0){
    p02 = FeaturePlot(aleks, features = feature.example, ncol = NULL, by.col = FALSE) +
      ggtitle(paste0('aleks - ', feature.example))
    plots[['aleks']] = p02
  }
  DefaultAssay(seurat.cistopic) <- 'RNA'
  p1 = FeaturePlot(
    object = seurat.cistopic,
    features = feature.example,
    pt.size = 0.5,
    #max.cutoff = 'q90',
    #min.cutoff = 'q5',
    ncol = 1
  )
  
  plots[['scatac']] = p1
  DefaultAssay(seurat.cistopic) <- 'chromvar'
  
  motifs2check = rownames(seurat.cistopic)[grep(feature.example, rownames(seurat.cistopic))]
  if(length(motifs2check)>0){
    p2 = FeaturePlot(
      object = seurat.cistopic,
      features = motifs2check[1],
      pt.size = 0.1,
      #cols = c('red', 'blue'),
      #max.cutoff = 1,
      min.cutoff = 0.2,
      ncol = 1
    )
    plots[['tf.motif']] = p2 
  }
  
  CombinePlots(plots, ncol = 2)
}

test.umap.params = function(seurat.obj)
{
  pdfname = paste0(resDir, "/cluster_annotations/test_umap_params_for_seurat.cistopic.pdf")
  pdf(pdfname, width=12, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #nb.pcs = ncol(seurat.obj[['pca']]);
  #n.neighbors = 50; min.dist = 0.01;
  
  for(nb.pcs in c(20, 30, 50, 70, 80, ncol(seurat.obj[['pca']])))
  {
    if(nb.pcs <= ncol(seurat.obj[['pca']])){
     for(n.neighbors in seq(10, 50, by = 10))
     {
       for(min.dist in c(0.01, 0.05, seq(0.1, 0.6, by = 0.1)))
       {
         cat('nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, ', min.dist - ', min.dist, '\n')
         seurat.obj <- RunUMAP(object = seurat.obj, reduction = 'pca', 
                                    dims = 1:nb.pcs, 
                                    n.neighbors = n.neighbors, min.dist = min.dist)
        p1 =  DimPlot(seurat.obj, label = TRUE, pt.size = 1, label.size = 6) + 
           NoLegend() + 
           ggtitle(paste0('nb.pcs - ', nb.pcs, '; n.neighbors - ', n.neighbors, ', min.dist - ', min.dist))
        plot(p1)
       }
     }
    }
  }
  dev.off()
  
}

cluster.annotation.using_marker.genes_tf.motif_peaks = function(seurat.cistopic)
{
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
  seurat.cistopic <- ScaleData(seurat.cistopic, do.center = FALSE, scale.max = 2)
  
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
                          grep('11.BWM|15.BWM', groups.orders), 
                          grep('pharynx/neurons', groups.orders),
                          grep('9.neurons/BWM', groups.orders),
                          grep('neurons', groups.orders)))
    
    levels(seurat.cistopic) <- groups.orders[new.orders]
    
    
    pdfname = paste0(resDir, "/cluster_annotations/marker_genes_TFs_top5_log2normalized_v4.pdf")
    pdf(pdfname, width=12, height = 16)
    par(cex =1.0, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    # 
    DoHeatmap(seurat.cistopic, features = tf.markers$gene, 
              slot = 'scale.data',
              size = 5, hjust = 0, label = TRUE, 
              disp.max = 2, 
              lines.width = 50) + scale_fill_gradientn(colors = rev(cols))
      # scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", 
      #                       high = rev(c('#b2182b','#ef8a62','#fddbc7')), 
      #                       midpoint = 0, guide = "colourbar", aesthetics = "fill")
    dev.off()
    
  }
  
  
  ##########################################
  # annotat those clusters using known lineage-specifc genes or tf or tf motif activities 
  ##########################################
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
  
  
  feature.example = feature.example.neurons
  
  DefaultAssay(seurat.cistopic) <- 'RNA'
  FeaturePlot(seurat.cistopic, features = feature.example, ncol = 4, reduction = 'umap')
  
  DefaultAssay(seurat.cistopic) <- 'chromvar'
  motifs2check = rownames(seurat.cistopic)[grep(paste0(feature.example, collapse = '|'), rownames(seurat.cistopic))]
  if(length(motifs2check)>0){
    FeaturePlot(
      object = seurat.cistopic,
      features = motifs2check,
      pt.size = 0.1,
      #cols = c('red', 'orange', 'blue'),
      #max.cutoff = 1,
      #min.cutoff = 0,
      ncol = 2
      #cols =  cols
    )
  }
  
  
}

wormmine.database.download = function(seurat.cistopic)
{
  library(InterMineR)
  listMines()
  
  # load HumaMine
  im <- initInterMine(mine=listMines()["WormMine"])
  im
  
  # Get template (collection of pre-defined queries)
  template = getTemplates(im)
  head(template)
  
  # Get gene-related templates
  template[grep("gene", template$name, ignore.case=TRUE),]
  
  # Query for gene pathways
  queryGenePath = getTemplateQuery(
    im = im, 
    name = "anatomical_expression_patterns_for_genes"
    #name = 'genes_in_expression_cluster'
  )
  
  
  queryGenePath
  
  DefaultAssay(seurat.cistopic) = 'RNA'
  
  
  wmm.data = c()
  for(n in 15884:nrow(seurat.cistopic))
  {
    # n = 7300
    g = rownames(seurat.cistopic)[n];
    
    cat(g)
    queryGenePath$where[[1]][["value"]] <- g
    
    resGenePath <- runQuery(im, queryGenePath)
    if(!is.null(resGenePath)) {
      wmm.data = rbind(wmm.data, resGenePath[, c(3,5,6)])
      cat('\n')
    }else{
      cat(' not found\n')
    }
  }
  
  saveRDS(wmm.data, file = paste0(RdataDir, 'wormbase_anatomical_expression_patterns_part.rds'))
  anatomy = rbind(anatomy, wmm.data)
  
  saveRDS(anatomy, file = paste0(RdataDir, 'wormbase_anatomical_expression_patterns_all.rds'))
  
  length(unique(anatomy$Gene.symbol[which(anatomy$Gene.expressionPatterns.anatomyTerms.name == 'Aba')]))
  
  
  ##########################################
  # start to edit the marker-gene-lineage table 
  ##########################################
  cells.0.28 = c('P0', # 1-cell  
               'AB', 'P1', # 2-cell
               'ABa','ABp', 'EMS', 'P2', # 4-cell
               'ABal', 'ABar', 'ABpl', 'ABpr', 'MS', 'E', 'C', 'P3', # 8-cell
               # 16-cell
               'ABala', 'ABalp', 'ABara', 'ABarp', 
               'ABpla', 'ABplp', 'ABpra', 'ABprp', 
               'MSa', 'MSp', 'Ea', 'Ep', 'Ca', 'Cp', 'D', 'P4', 
               
               # 28-cell
               'ABalaa', 'ABalap', 'ABalpa', 'ABalpp', 'ABaraa', 'ABarap','ABarpa', 'ABarpp',
               'ABplaa', 'ABplap', 'ABplpa', 'ABplpp', 'ABpraa', 'ABprap','ABprpa', 'ABprpp',
               'MSaa', 'MSap', 'MSpa', 'MSpp',
               'Ea', 'Ep',
               'Caa', 'Cap', 'Cpa', 'Cpp', 
               'D', 'P4'
  )
  cells.51 = c(# 51-cell
    'ABalaaa', 'ABalapa', 'ABalpaa', 'ABalppa', 'ABaraaa', 'ABarapa','ABarpaa', 'ABarppa',
    'ABalaap', 'ABalapp', 'ABalpap', 'ABalppp', 'ABaraap', 'ABarapp','ABarpap', 'ABarppp',
    'ABplaaa', 'ABplapa', 'ABplpaa', 'ABplppa', 'ABpraaa', 'ABprapa','ABprpaa', 'ABprppa',
    'ABplaap', 'ABplapp', 'ABplpap', 'ABplppp', 'ABpraap', 'ABprapp','ABprpap', 'ABprppp',
    'MSaaa', 'MSapa', 'MSpaa', 'MSppa', 
    'MSaap', 'MSapp', 'MSpap', 'MSppp',
    'Eaa', 'Eap', 'Epa', 'Epp',
    'Caa', 'Cap', 'Cpa', 'Cpp', 
    'Da', 'Dp', 
    'P4')
  
  lineages = unique(cells.0.28)
  lineages = lineages[grep('P1|P2|P3|P4|P0|EMS', lineages, invert = TRUE)]
  lineages = lineages[which(lineages != 'AB')]
  
  anatomy = readRDS(file = paste0(RdataDir, 'wormbase_anatomical_expression_patterns_all.rds'))
  
  marker.annatom = matrix(0, ncol = length(lineages), nrow = length(unique(anatomy$Gene.symbol)))
  colnames(marker.annatom) = lineages
  rownames(marker.annatom) = unique(anatomy$Gene.symbol)
  
  for(n in 1:ncol(marker.annatom))
  {
    #n = 1
    ll = colnames(marker.annatom)[n]
    cat(n, ' -- ', ll, '\n')
    ggs = unique(anatomy$Gene.symbol[which(anatomy$Gene.expressionPatterns.anatomyTerms.name == ll)])
    
    marker.annatom[match(ggs, rownames(marker.annatom)) ,n]  = 1 
  }
  
  ss = apply(marker.annatom, 1, sum)
  marker.annatom = marker.annatom[ss>0, ]
  
  #DefaultAssay(seurat.cistopic) = 'RNA'
  counts = GetAssayData(seurat.cistopic, slot = 'counts', assay = 'RNA')
  counts = counts > 0
  rownames(counts) = rownames(seurat.cistopic)
  counts = counts[match(rownames(marker.annatom), rownames(seurat.cistopic)), ]
  
  counts = as.data.frame(counts)
  
  ee.markers = readRDS(file = paste0(RdataDir, 'scATAC_seq_geneActivityMatrix.inPromoter.resolution.0.8_markerGenes.rds'))
  mm = match(rownames(counts), ee.markers$gene)
  length(which(!is.na(mm)))
  mm = !is.na(mm)
  
  counts = counts[mm, ]
  marker.annatom = marker.annatom[mm, ]
  
  library("pheatmap")
  library("RColorBrewer")
  cols = c(colorRampPalette((brewer.pal(n = 7, name="Reds")))(10))
  pheatmap(t(marker.annatom), cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = seq(0, 1, by = 0.2),
           cluster_cols=TRUE, 
           color = cols, 
           fontsize_col = 10,
           height = 8,
           width = 30
  )
  
  corefs = cor(counts, marker.annatom)
  
  components = unlist(apply(corefs, 1, function(x) colnames(corefs)[which.max(x)]))
  
  seurat.cistopic$anatamy.inferred = components
  
  #DimPlot(seurat.cistopic, label = TRUE, pt.size = 1, label.size = 6) 
  DimPlot(seurat.cistopic, label = FALSE, group.by = 'anatamy.inferred', pt.size = 0.5, label.size = 6) 
  
  nn = c(1:ncol(seurat.cistopic))
  
  nn = which(seurat.cistopic$peaks_snn_res.0.8 == 7)[c(1:100)]
  components = c()
  for(n in  nn)
  {
    #n = 56
    cors = cor(counts[match(rownames(marker.annatom), rownames(seurat.cistopic)),n], marker.annatom, method = 'spearman')
    names(cors) = colnames(marker.annatom)
    cat(names(cors)[which.max(cors)], '\n')
    components = c(components, names(cors)[which.max(cors)])
    
  }
  
  barplot(table(components), las = 2)

}


