Save.and.Clean.annotation.from.bwm = function(seurat.obj)
{
  # import the most updated annotation of bwm that is 37th iteration
  nb.iteration = 37
  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_',
                    nb.iteration, '.rds')
  seurat.obj = readRDS(file = RDSsaved)


  # 'manual.annot.ids.2' have the annotated non-BWM cell ids saved
  DimPlot(seurat.obj, group.by = "manual.annot.ids.2", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_annotedIDs")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  seurat.obj$manual.annot.ids.backupBWM = seurat.obj$manual.annot.ids


  # put back non-BWM cluster annotations; but one can select bwm with metatdata 'BWM.cells'
  kk = which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$manual.annot.ids.2))

  cat(length(kk), ' cells annotated as non-BWM cell ids; \nhowever due to early iteration some of them have BWM annotation  \n')

  ids.1 = names(table(seurat.obj$manual.annot.ids)) # lastest BWM annotations
  ids.2 = names(table(seurat.obj$manual.annot.ids.2[kk])) # cells with non-BWM annotations
  print(ids.1);
  print(ids.2)

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()


  ##########################################
  # here is what we wanna do:
  # if ids.2 have BWM annotations, it means that the annotations were outdated and put back to NA
  # if ids.2 are not annotated as BWM, they were put back to clean annotations
  ##########################################
  ids2clean = intersect(ids.2, ids.1)
  ids2add = setdiff(ids.2, ids.1)

  jj = which(!is.na(match(seurat.obj$manual.annot.ids.2, ids2add)))
  cells2add = unique(colnames(seurat.obj)[jj])
  cat(length(cells2add), 'cells will be added non-BWM annotatation\n')

  seurat.obj$manual.annot.ids[jj] = seurat.obj$manual.annot.ids.2[jj]

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()


  idntt = c("MSapaapp")
  cells_to_show <- list(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, idntt))])
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
          cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1,
          pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)

  idntt = c("MSppaapp")
  cells_to_show <- list(colnames(seurat.obj)[!is.na(match(seurat.obj$predicted.ids.scmap, idntt))])
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
          cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1,
          pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)

  saveRDS(seurat.obj, file = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                 '_cleanedBWM_and_Pharynx_iteration_0.rds'))


}

################################################################################################
################################################################################################
# some parts from BWM annotation could be reused later
################################################################################################
################################################################################################
highligh.cells.dissociate.close.cell.ids.with = FALSE
 if(highligh.cells.check.markers){
   idntt = c('MSxpapap')
   cells_to_show <- list(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, idntt))])
   p1 = DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
                cells.highlight = cells_to_show, cols.highlight = 'blue', sizes.highlight = 1,
                pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
   #idntt = 'MSxpappp'
   cells_to_show <-  list(colnames(sub.obj)[which(sub.obj$predicted.ids.seurat.keep == idntt)])
   p2 = DimPlot(sub.obj, group.by = "pred.ids", reduction = 'umap', label = TRUE, repel = TRUE,
                cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1.5,
                pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
   p1 + p2

   idntt = c('likely_MSxpappp')
   cells_to_show <- list(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, idntt))])
   p3 = DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
                cells.highlight = cells_to_show, cols.highlight = 'blue', sizes.highlight = 1,
                pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
   idntt = 'MSxpappp'
   cells_to_show <-  list(colnames(sub.obj)[which(sub.obj$predicted.ids.seurat.keep == idntt)])
   p4 = DimPlot(sub.obj, group.by = "pred.ids", reduction = 'umap', label = TRUE, repel = TRUE,
                cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1.5,
                pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
   p3 + p4


   # MSx
   extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSx', group_2 = 'MSxp', ntop = 10)
   features.sels = c( 'sdz-1', 'hnd-1', 'pha-4', 'fbxb-70'
                      # 'clec-91', 'clec-87',  'clec-87','arrd-2','clec-88', 'otub-3'
   )
   FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

   # MSxpa
   FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-39', 'R11A5.3', 'egl-43'))
   FeaturePlot(sub.obj, reduction = 'umap', features = c('zig-8','ham-1', 'F19C6.4',  'shc-2', 'F07C6.4',
                                                         'unc-120', 'abts-1' # MSxpaaa
   ))
   # MSxpappa
   #source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
   extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpappp', group_2 = 'MSxppapp', ntop = 20)
   features.sels = c('stn-2', 'Y9C2UA.1', 'hsp-12.1', 'unc-62',
                     'ceh-34', 'zig-6',
                     'fbxb-88', 'fbxc-24', 'E01G4.5',  'igcm-4'
   )
   FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

   extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpappp', group_2 = 'MSxppapp', ntop = 20)
   features.sels = c('maph-1.3', 'shc-2', 'twk-31', 'mec-2',   # MSxppapp
                     'D1086.12', 'tbx-2', 'K09G1.1', 'zig-6', 'stn-2', 'T25G12.11', 'unc-5') # MSxpappp
   FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

   # MSxpapap (almost sure)
   features.sels = c('maph-1.2', 'maph-1.3', 'K09G1.1', 'T25G12.11', 'stn-2', 'hnd-1', 'Y37E3.30', 'zig-6', 'Y9C2UA.1')
   FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)


   DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) + NoLegend()

   p1 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2) + NoLegend() +
     ggtitle('manual ids')

   p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat.keep', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2) +
     NoLegend() + ggtitle('seurat.pred.ids.keep')

   p3 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2) +
     NoLegend() + ggtitle('seurat.pred.ids.bwm.all')

   p1 /(p2 + p3)

   FeaturePlot(sub.obj, reduction = 'umap', features = 'maph-1.2')

 }
 ##########################################
   # rerun the seurat for label transferring
   ##########################################
   RErun.seurat.transferring.labels = FALSE
   if(RErun.seurat.transferring.labels){
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

     ids.refs = c(terminals,
                  'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxpapp', 'MSxpapa', 'MSxpaaa', 'MSxappp', 'MSpappa')
     #ids.refs = terminals
     #sub.obj = find.reference.mapped.ids.for.terminalCells.scmap(sub.obj, nfeatures = 2000, terminals = terminals)
     sub.obj = seurat.transfer.labels.from.Murray.scRNA.to.scRNA.terminalCells(sub.obj, nfeatures = 3000, npcs = 30,
                                                                               reduction = 'cca',
                                                                               k.anchor = 5, k.filter = 200, max.features = 200,
                                                                               terminals = ids.refs)
     sub.obj$predicted.ids = sub.obj$predicted.ids.seurat.terminal
     sub.obj$predicted.ids.prob = sub.obj$predicted.ids.seurat.terminal.prob

     hist(sub.obj$predicted.ids.prob)
     abline(v = 0.5, col = 'red', lwd=2.0)

     seurat.obj$pred.ids.terminals.mothers.seurat = NA
     seurat.obj$pred.ids.terminals.mothers.seurat.prob = NA
     mm = match(colnames(sub.obj), colnames(seurat.obj))
     seurat.obj$pred.ids.terminals.mothers.seurat[mm] = sub.obj$predicted.ids.seurat.terminal
     seurat.obj$pred.ids.terminals.mothers.seurat.prob[mm] = sub.obj$predicted.ids.seurat.terminal.prob


     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     sub.obj = reference.based.cell.projection.rf.svm(sub.obj, nfeatures = 3000, cost = 0.5, ntree = 50, scale.cell = TRUE,
                                                      terminals = ids.refs)

     sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.seurat.terminal
     sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.5] = NA
     sub.obj$pred.ids.svm.filter = sub.obj$pred.ids.svm
     sub.obj$pred.ids.svm.filter[sub.obj$pred.ids.svm.prob <0.5] = NA

     p1 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
       ggtitle('manual ids')
     #p1 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6)
     p2 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
       ggtitle('predicted ids. seurat')
     p3 = DimPlot(sub.obj, group.by = 'pred.ids.svm', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
       ggtitle('predicted ids svm')
     p4 = DimPlot(sub.obj, group.by = 'predicted.ids.fitered', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
       ggtitle('seurat.pred.ids seurat. filtered')
     p5 = DimPlot(sub.obj, group.by = 'pred.ids.svm.filter', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
       ggtitle('predicted ids svm filtered')

     (p2 + p3) / (p4 + p5)

     p2 + p3
     p1 + p2
     p1 + p3

     DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
     DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE) + NoLegend()
     #DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 4, repel = TRUE)

     # save the seurat prediction for all BWM cells
     seurat.obj$pred.ids.terminals.mothers.seurat = NA
     seurat.obj$pred.ids.terminals.mothers.svm = NA
     #sub.obj$predicted.ids.seurat.keep.bwm.all = NA
     #sub.obj$pred.ids.terminals.mothers.seurat = sub.obj$predicted.ids.seurat.terminal
     #seurat.obj$pred.ids.seurat.terminals.mothers = NA
     mm = match(colnames(sub.obj), colnames(seurat.obj))
     seurat.obj$pred.ids.terminals.mothers.seurat[mm] = sub.obj$predicted.ids.seurat.terminal
     seurat.obj$pred.ids.terminals.mothers.svm[mm] = sub.obj$pred.ids.svm

   }

   p0 = DimPlot(sub.obj, group.by = 'timingEst', reduction = 'umap', label = FALSE, label.size = 5)
   p1 = DimPlot(sub.obj, group.by = 'request', reduction = 'umap', label = FALSE, label.size = 5)
   p0 + p1

   #FeaturePlot(sub.obj, reduction = 'umap', features = c('unc-120', 'pha-4', 'hnd-1'))
   #FeaturePlot(seurat.obj, reduction = 'umap', features = c('unc-120', 'pha-4', 'hnd-1'))
   for(idntt in unique(sub.obj$manual.annot.ids))
     {
       cells_to_show <- list(c(WhichCells(sub.obj, idents = idntt)))
       p1  = DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
                     cells.highlight = cells_to_show, cols.highlight = 'blue', sizes.highlight = 1,
                     pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)
       plot(p1)
     }

     idntt = c('MSxpapap'); cells_to_show <- list(c(WhichCells(sub.obj, idents = idntt)))
     DimPlot(sub.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
             cells.highlight = cells_to_show, cols.highlight = 'red', sizes.highlight = 1,
             pt.size = 2, label.size = 5) + NoLegend() + ggtitle(idntt)

     # check info in JM data for specific lineage
     #' source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     #' extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppppp', group_2 = 'MSxpppaa', ntop = 10)
     #' features.sels = c(#'hnd-1',  'pha-4',
     #'   #'fbxb-70',
     #'   #'ceh-37', 'C45G7.4', # MSxap
     #'   #'tbx-8', 'hlh-1', 'hnd-1', 'unc-39' # MSxpap
     #'   #'ham-1',  'nhr-67', # MSxapp
     #'   #'unc-39', 'irx-1', 'egl-43'
     #'   #'unc-120', 'tab-1', 'rpm-1', 'F55C5.10', 'hil-3' # MSpappa
     #'   #'tbx-8', 'asic-2', 'skpo-1', 'pxn-2', 'D1086.12','tab-1'  # MSxappp
     #'   #'mnr-1', 'strl-1', 'zig-10', 'gana-1', 'fbxb-88', # MSxapppa
     #'   #'lin-39', 'F54D5.5'  # MSxapppp
     #'   #'gsnl-1', 'abts-1', 'stn-2', # MSxapppax
     #'   #'lin-39',  'tbx-2', 'B0379.1'  #MSxappppx
     #'   #'egl-43', 'R11A5.3' #MSxpa
     #'   #'zip-7', 'hlh-1', 'col-111' # MSxpap
     #'   #'tbx-7', 'unc-120', 'ref-2', 'hlh-16', 'unc-39', 'ttr-50', 'clec-266' # MSxpaa
     #'   #'F41D9.2', 'hnd-1', 'unc-120', 'tbx-7', 'abts-1', 'Y66D12A.13' # MSxpaaa
     #'   #'let-381', 'ins-2', 'F40H3.3', 'ZK183.5', # MSxpapa
     #'   #'unc-120', 'ceh-51', 'lag-2', 'fbxb-22', 'T02G6.11', 'C06A8.3', 'T05D4.2',   # MSxpp
     #'   #'zip-7', 'hnd-1', 'tbx-8', 'fkh-2', 'tbx-11', # MSxppp
     #'   #'hlh-16', 'rgs-7', 'tbx-8', 'unc-120', 'T24C2.2', 'T24C2.3'  #MSxppa
     #'   #'col-118', 'let-381', 'C03B1.1',  'T11B7.2' #MSxppaa
     #'   #'F40H3.3', 'Y42H9B.3', 'unc-120', 'ZK180.5', 'sfrp-1', 'ceh-34', 'irx-1', 'Y116A8C.3',
     #'   #'let-381', 'ccg-1' # MSxpapa
     #'   #'zig-8', 'pat-9', 'F19C6.4', 'tbx-8', 'her-1' # MSxpapp
     #'   #'sul-2', 'camt-1', 'irx-1', 'ctg-1',  # MSxpppa
     #'   #'tbx-2', 'D1086.12' #MSxpppp
     #'   #'ceh-34', 'ham-1', 'F19C6.4', 'eya-1', 'shc-2', 'F07C6.4', 'ten-1', 'hil-7' # MSxppap
     #'   #'maph-1.3', 'fkh-2', 'ccg-1', 'zig-6',
     #'   #'ham-1', 'shc-2', 'lam-3', 'T25G12.11','F57F4.4', 'cnt-2', 'Y17G7B.23', 'Y37E3.30'  # MSxpaaap
     #' )
     #'
     # MSxppppp
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppppp', group_2 = 'MSxpppaa', ntop = 10)
     features.sels = c('clec-264', 'ceh-13', 'gana-1','D1086.12', 'tbx-2', 'F37H8.5', 'B0379.1'
     )
     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

     # MSxpppap
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpppap', group_2 = 'MSxppppp', ntop = 10)
     features.sels = c('ceh-13', 'fbxc-24', 'hot-1', 'pde-6',  'hmg-1.1' )
     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

     # MSxpppaa
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpppaa', group_2 = 'MSxppppa', ntop = 10)
     features.sels = c('bgal-1', 'sul-1', 'zig-7')
     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

     # MSxppppa
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppppa', group_2 = 'MSxpppaa', ntop = 10)
     features.sels = c('lntl-1', 'D1086.12', 'tbx-2', 'ZK512.1')
     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

     # MSxpppp and MSxpppa
     features.sels = c( 'sul-2', 'camt-1', 'irx-1', 'ctg-1',  # MSxpppa
                        'tbx-2', 'D1086.12' #MSxpppp
     )
     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

     #FeaturePlot(sub.obj, reduction = 'umap', features = c('zig-6', 'eya-1'))

     # MSxpaaap
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpaaap', group_2 = 'MSxpapap', ntop = 20, test.use = 'roc')

     features.sels = c('maph-1.3', 'maph-1.2', 'ZC449.5', 'shc-2',
                       'T25G12.11', 'zig-8', 'ceh-34', 'ham-1', 'lam-3', 'twk-31', 'stn-2')
     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

     # MSxppapp
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxppapp', group_2 = 'MSxpappa', ntop = 10, test.use = 'roc')

     features.sels = c('K09G1.1', 'tbx-2', 'D1086.12', 'maph-1.3', 'shc-2', 'twk-31',
                       'zig-6', 'mec-2')
     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

     # MSxpappp
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpappp', group_2 = 'MSxpapap', ntop = 10)

     features.sels = c('zig-8', 'D1086.12', 'K09G1.1', 'tbx-2', 'zig-6', 'stn-2',
                       'lipl-7')
     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

     # MSxpappa
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpappa', group_2 = 'MSxpppaa', ntop = 10)

     features.sels = c('Y9C2UA.1', 'hsp-12.1', 'fbxc-24',
                       'unc-62', 'cpn-3', 'lipl-7', 'fbxb-88', 'tnt-3', 'lbp-1', 'F37H8.5', 'E01G4.5',
                       'ceh-34', 'zig-6'

     )

     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

     # MSxpapap
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpapap', group_2 = 'MSxpappp', ntop = 10)

     features.sels = c('maph-1.2', 'maph-1.3', 'K09G1.1', 'T25G12.11', 'stn-2', 'hnd-1', 'Y37E3.30', 'zig-6', 'Y9C2UA.1')
     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)


     # MSxpapa and MSxpapp
     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     extrack.markers.from.JM(markers = markers, eet = eet, group_1 = 'MSxpapa', group_2 = 'MSxpapp', ntop = 10)

     features.sels = c('zig-8', 'tbx-8', 'hlh-1', 'pat-9','tag-196', 'fkh-2', 'spp-15', 'cnt-2',
                       'her-1','ceh-34', 'sfrp-1','let-381', 'irx-1', 'ins-2'
     )
     FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

     

 ##############################################################################################################################


# iteration of manual annotation for pharynx
library(ggplot2)
 library(patchwork)
 library("pheatmap")
 library("RColorBrewer")
 library(grid)
 #library(RaceID) # refer to the vignett https://cran.r-project.org/web/packages/RaceID/vignettes/RaceID.html
 library(Matrix)
 #library(lsa)
 library(dplyr)
 library(openxlsx)
