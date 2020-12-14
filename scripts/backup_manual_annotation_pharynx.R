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
##############################################################################################################################
# iteration of manual annotation for pharynx
##############################################################################################################################
##############################################################################################################################
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


 ##########################################
   # Main aim:  to find MSxaa, while MSxa and MSxap were found already
   # Notes:
   # adding whole MSxap actually does not help
   ##########################################
   nb.iteration = 1
   Refine.annotated.ids = FALSE;

   resDir = paste0("results/", version.analysis, '/annoted_pharynx')
   if(!dir.exists(resDir)){dir.create(resDir)}

   if(Refine.annotated.ids){by.group = 'manual.annot.ids';
   }else{by.group = 'seurat_clusters'}

   RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                  '_cleanedBWM_and_Pharynx_iteration_', nb.iteration -1, '.rds')

   RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                      '_cleanedBWM_and_Pharynx_iteration_', nb.iteration, '.rds')

   seurat.obj = readRDS(file = RDSsaved)

   pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

   #seurat.obj$predicted.ids.scmap = seurat.obj$scmap.pred.id.500
   #seurat.obj$predicted.ids.seurat = seurat.obj$seurat.pred.id



   pdf(pdfname, width=18, height = 10)
   par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
   ##########################################
   # select subset of cells to annotate
   ##########################################
   # cluster.index = '25'
   # table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
   #
   # table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
   #
   # table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]
   #
   # xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
   # xx[which(xx > 0)]
   #
   # ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
   # xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
   # table(seurat.obj$manual.annot.ids[ii1])
   # #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
   #FeaturePlot(sub.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))

   # select cells with cluster index
   ##########################################
   cluster.sels = c('29', '32', '35', '40', '42')
   #cluster.sels = c('25', '36', '8', '39', '2', '19', '27', '13', '1', '11', '33', '48', '18', '46', '15', '26')

   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
   #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
   #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) & is.na(seurat.obj$manual.annot.ids)]
   #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))]
   #                                   & is.na(seurat.obj$manual.annot.ids))]

   # select BWM terminal and middle time points cells
   #' ##########################################
   #' cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
   #'                  '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
   #'                  '25', # possible transition clusters
   #'                  '24', # also transitions cells and many of them are not annotated
   #'                  #'44', '31', '52', '28', '50', # cluster '44', '31', '52', '28', '50' were not included here
   #'                  '3', '5', '16', '30', '22', '4' # all middle time points
   #' )
   #'
   #cluster.sels = c('29', '32', '35', '40', '42')
   #sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
   # cells.sels = unique(colnames(seurat.obj)[seurat.obj$BWM.cells == 'BWM' &
   #   (!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
   #    !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
   #     )])
   # cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
   #                                               !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
   #                                            ])
   #seurat.obj$manual.annot.ids.6 = seurat.obj$manual.annot.ids
   #seurat.obj$BWM.cells[which(seurat.obj$manual.annot.ids == 'unknown_MSxpppaa_MSxppppa_later')] = NA

   #kk = which(seurat.obj$pred.ids.terminals.mothers.seurat == 'MSpappax')
   #ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
   #ids.sels = ids.current[which(nchar(ids.current)>6)]
   #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')

   # ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)],
   #                    c('MSxppaa'))
   #
   # ids.left = setdiff(ids.current, ids.sels)
   # print(ids.left)
   # nchar(ids.left)

   #cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

   ##########################################
   # subset seurat object with selected cells
   ##########################################
   sub.obj = subset(seurat.obj, cells = cells.sels)

   #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
   sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
   sub.obj$pred.ids = sub.obj$predicted.ids.seurat.keep
   xx = table(sub.obj$predicted.ids.seurat.keep)
   xx[xx>10]
   sub.obj$pred.ids.filtered = sub.obj$pred.ids
   sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

   DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

   barplot(table(sub.obj$seurat_clusters)/table(seurat.obj$seurat_clusters), ylim = c(0, 1), cex.names = 1.0, las=2)

   # sub.obj$manual.annot.ids = sub.obj$predicted.ids.seurat.keep
   ##########################################
   # check potential ids for selected clusters
   ##########################################
   # #DimPlot(sub.obj, reduction = 'umap', group.by = 'scmap.pred.id.500')
   # threshold = 0.7
   # predicted.ids = sub.obj$scmap.pred.id.500
   # #predicted.ids[which(sub.obj$scmap.corr.500 < threshold)] = 'unassigned'
   #
   # if(Refine.annotated.ids){
   #   counts = table(predicted.ids, sub.obj$manual.annot.ids)
   #   counts.seurat = table(as.character(sub.obj$seurat.pred.id), sub.obj$manual.annot.ids)
   # }else{
   #   counts = table(predicted.ids, as.character(sub.obj$seurat_clusters))
   #   counts.seurat = table(as.character(sub.obj$seurat.pred.id), as.character(sub.obj$seurat_clusters))
   # }
   counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters))
   barplot(counts, main="cluster compositions by scmap ",
           xlab=NULL, col=c(1:nrow(counts)), las = 2,
           legend = rownames(counts))
   #
   # barplot(counts.seurat, main="cluster compositions by seurat ",
   #         xlab=NULL, col=c(1:nrow(counts)), las = 2,
   #         legend = rownames(counts))

   #counts[, match(c('31', '28', '52'), colnames(counts))]
   #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
   ##########################################
   # find new set of variable genes and redo pca and umap
   ##########################################
   Explore.umap.parameters.for.BWMcells = FALSE
   if(Explore.umap.parameters.for.BWMcells){

     source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
     require(tictoc)
     tic()
     test.umap.params.for.BWM.cells(sub.obj,
                                    pdfname = 'UMAP_param_TEST_BWM_searching_for_MSx.pdf',
                                    group.by = 'manual.annot.ids', with_legend = TRUE,
                                    nfeatures.sampling = c(500, 1000), nb.pcs.sampling = c(10, 20, 30),
                                    n.neighbors.sampling = c(5, 10, 30, 50),
                                    min.dist.sampling = c(0.01, 0.1)
     )
     toc()

   }

   nfeatures = 1000;
   sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
   #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
   sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
   sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
   ElbowPlot(sub.obj, ndims = 50)

   nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
   n.neighbors = 30;
   min.dist = 0.01; spread = 1
   sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                      spread = spread, n.neighbors = n.neighbors,
                      min.dist = min.dist, verbose = TRUE)
   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()

   #VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'manual.annot.ids') + NoLegend()
   #xx = table(sub.obj$predicted.ids.seurat.keep)
   #xx[xx>10]
   #sub.obj$pred.ids.filtered = sub.obj$pred.ids
   #sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA
   #jj2 = which(!is.na(match(sub.obj$predicted.ids.seurat.keep, c('MSxppppx', 'MSxpppax'))) == TRUE)
   p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()
   #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
   p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
   p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)

   p0 + p1

   p0 + p2
   p1 + p2


   ##########################################
   # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
   ##########################################
   FindClusters_subclusters = function(sub.obj, resolution = 0.4)
   {
     sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
     return(sub.obj$seurat_clusters)
   }
   sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
   sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.5)
   DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

   p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
     ggtitle('manual.annot.ids') + NoLegend()
   p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
                label.size = 6, na.value = "gray", combine = TRUE)
   p1 + p2

   p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
                group.by = 'seurat_clusters_split')

   p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                group.by = 'seurat_clusters_split') + NoLegend()

   VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
           group.by = 'manual.annot.ids') + NoLegend()

   p1 + p4
   p2 + p4
   plot(p3)


   dev.off()

   ##########################################
   # check the counts of predicted ids for newly split clusters
   ##########################################
   sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
   sub.obj$predicted.ids.prob = sub.obj$predicted.scores
   sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
   sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

   Idents(sub.obj) = sub.obj$seurat_clusters_split
   counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
   counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
   #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
   counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)

   p1 = DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
   p2 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()
   p3 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()
   p4 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()


   (p1 + p2)/(p3 + p4)


   Idents(sub.obj) = sub.obj$seurat_clusters_split
   idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
   idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

   idents.sel = c('4', '8', '10')

   ## chcek the reference-mapped ids for the rest of clusters
   counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
   counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
   counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
   counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
   counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
   counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

   features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                     'ceh-36', # MSxaap
                     'tbx-2', 'ceh-27') # MSxaaa

   VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)

   FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

   # check info in JM data for specific lineage
   # ee = process.import.Murray.scRNA()
   # murray.ids = unique(ee$lineage)
   #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
   markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
   #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
   #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))

   ids.sel = c('MSxaa')
   source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
   find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)

   top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
   DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

   # to find new marker genes
   top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
   top.markers[top.markers$cluster == '13',]


   FeaturePlot(sub.obj, reduction = 'umap',
               features = c('zig-6', 'frpr-8','C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-4', 'maph-1.2', 'irx-1')) # MSxpapa

   ##########################################
   # update the manual annotation if good marker genes or mapped ids were found
   ##########################################
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
   # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
   #mm = match(colnames(sub.obj), colnames(seurat.obj))
   #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

   # sub.obj$ids.backup = sub.obj$manual.annot.ids

   cluster.assingment = list(#c('0', 'MSxp'),
     c('4', 'MSxaa'),
     c('8', 'MSxaap'),
     c('10', 'MSxaaa')
     #c('3', 'mixture_BWM_terminal_2'),
     #c('4', 'mixture_BWM_terminal_2')
   )

   for(n in 1:length(cluster.assingment)){
     cluster.index = cluster.assingment[[n]][1];
     id2assign =  cluster.assingment[[n]][2];
     cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
     cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
     sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
     seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
     # if(is.na(id2assign)) {
     #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
     # }else{
     #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
     # }
   }

   #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
           pt.size = 2) + NoLegend()

   DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
     scale_colour_hue(drop = FALSE) + NoLegend()


   saveRDS(seurat.obj, file = RDS2save)


   ##########################################
    # Main aim: because MSxaa, MSxaap, MSxaaa found, this iteration will devoted to MSxapa and to explore the whole MSxa lineages
    # in this iteration, the pharynx-related lineages were reclustered, so that the cell selection will be done with clusters_pharynx
    # Notes:
    #
    ##########################################
    nb.iteration = 2
    Refine.annotated.ids = FALSE;

    resDir = paste0("results/", version.analysis, '/annoted_pharynx')
    if(!dir.exists(resDir)){dir.create(resDir)}

    if(Refine.annotated.ids){by.group = 'manual.annot.ids';
    }else{by.group = 'seurat_clusters'}

    RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                   '_cleanedBWM_and_Pharynx_iteration_', nb.iteration -1, '.rds')

    RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                       '_cleanedBWM_and_Pharynx_iteration_', nb.iteration, '.rds')

    seurat.obj = readRDS(file = RDSsaved)

    pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

    DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
            na.value = "gray") +
      ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
      scale_colour_hue(drop = FALSE) +
      NoLegend()

    cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')

    pdf(pdfname, width=18, height = 10)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    ##########################################
    # select subset of cells to annotate
    ##########################################
    # cluster.index = '25'
    # table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
    #
    # table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
    #
    # table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]
    #
    # xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
    # xx[which(xx > 0)]
    #
    # ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
    # xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
    # table(seurat.obj$manual.annot.ids[ii1])
    # #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
    #FeaturePlot(sub.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))

    # select cells with cluster index
    ##########################################
    #cluster.sels = c('29', '32', '35', '40', '42')
    #ids.sel = c()

    cluster.sels = c('6', '14', '0', '10', '38', '9', '12', '7', '20', '49', '47')
    ids.sel = c('MSxaa', 'MSxaap', 'MSxaaa', 'MSxap')
    ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')

    if(length(ids.sel)>0){
      cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
                                               !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
                                                 is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
    }else{
      cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
    }

    # select BWM terminal and middle time points cells
    #' ##########################################
    #' cluster.sels = c('36', '8', '39', '2', '19', '27', # BWM_terminal_1 without transition
    #'                  '13', '1', '11', '33', '48', '18', '46', '15', '26', # BWM_terminal_2
    #'                  '25', # possible transition clusters
    #'                  '24', # also transitions cells and many of them are not annotated
    #'                  #'44', '31', '52', '28', '50', # cluster '44', '31', '52', '28', '50' were not included here
    #'                  '3', '5', '16', '30', '22', '4' # all middle time points
    #' )
    #'
    #cluster.sels = c('29', '32', '35', '40', '42')
    #sub.obj = subset(seurat.obj, cells = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
    # cells.sels = unique(colnames(seurat.obj)[seurat.obj$BWM.cells == 'BWM' &
    #   (!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
    #    !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
    #     )])
    # cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels)) |
    #                                               !is.na(match(seurat.obj$manual.annot.ids, ids.sels))
    #                                            ])
    #seurat.obj$manual.annot.ids.6 = seurat.obj$manual.annot.ids
    #seurat.obj$BWM.cells[which(seurat.obj$manual.annot.ids == 'unknown_MSxpppaa_MSxppppa_later')] = NA

    #kk = which(seurat.obj$pred.ids.terminals.mothers.seurat == 'MSpappax')
    #ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
    #ids.sels = ids.current[which(nchar(ids.current)>6)]
    #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')

    # ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)],
    #                    c('MSxppaa'))
    #
    # ids.left = setdiff(ids.current, ids.sels)
    # print(ids.left)
    # nchar(ids.left)

    #cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

    ##########################################
    # subset seurat object with selected cells
    ##########################################
    cat(length(cells.sels), ' cells selected to annotate \n')
    sub.obj = subset(seurat.obj, cells = cells.sels)

    #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
    sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
    sub.obj$pred.ids = sub.obj$predicted.ids.seurat
    xx = table(sub.obj$predicted.ids.seurat.keep)
    xx[xx>10]
    sub.obj$pred.ids.filtered = sub.obj$pred.ids
    sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

    DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()


    ##########################################
    # check potential ids for selected clusters
    ##########################################
    counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters))
    barplot(counts, main="cluster compositions by scmap ",
            xlab=NULL, col=c(1:nrow(counts)), las = 2,
            legend = rownames(counts))

    #counts[, match(c('31', '28', '52'), colnames(counts))]
    #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
    ##########################################
    # find new set of variable genes and redo pca and umap
    ##########################################
    Explore.umap.parameters.for.BWMcells = FALSE
    if(Explore.umap.parameters.for.BWMcells){

      source.my.script('scRNA_functions.R')
      require(tictoc)
      tic()
      test.umap.params.for.BWM.cells(sub.obj,
                                     pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                     group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                     nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                     n.neighbors.sampling = c(10, 30, 50),
                                     min.dist.sampling = c(0.01, 0.1)
      )
      toc()

    }

    nfeatures = 1000;
    sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
    #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
    sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
    sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
    ElbowPlot(sub.obj, ndims = 50)

    nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
    n.neighbors = 30;
    min.dist = 0.1; spread = 1
    sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                       spread = spread, n.neighbors = n.neighbors,
                       min.dist = min.dist, verbose = TRUE)
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()


    p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend()
    #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
    p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend() + ggtitle('scmap prediction')

    p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
      NoLegend() + ggtitle('scran prediction')

    p0 + p2

    p1 + p2


    ##########################################
    # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
    ##########################################
    FindClusters_subclusters = function(sub.obj, resolution = 0.4)
    {
      sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
      return(sub.obj$seurat_clusters)
    }
    sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:20, compute.SNN = TRUE)
    sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 2.0)
    DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

    p1  = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,
                  pt.size = 2) +
      ggtitle('manual.annot.ids') + NoLegend()
    p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
                 label.size = 8, na.value = "gray", combine = TRUE)

    p1 + p2 + ggsave(paste0(resDir, '/UMAP_pharynx_seurat_prediction_reclustering_newBase.pdf'), width = 22, height = 10)


    p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
                 group.by = 'seurat_clusters_split')

    p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                 group.by = 'seurat_clusters_split') + NoLegend()

    VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
            group.by = 'manual.annot.ids') + NoLegend()

    p1 + p4
    p2 + p4
    plot(p3)

    seurat.obj$seurat_clusters_pharynx = NA
    seurat.obj$seurat_clusters_pharynx[match(colnames(sub.obj), colnames(seurat.obj))] =
      as.integer(as.character(sub.obj$seurat_clusters_split))


    dev.off()

    ##########################################
    # check the counts of predicted ids for newly split clusters
    ##########################################
    sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
    sub.obj$predicted.ids.prob = sub.obj$predicted.scores
    sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
    sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

    Idents(sub.obj) = sub.obj$seurat_clusters_split
    counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
    counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
    #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
    counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)

    p1 = DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
    p2 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()
    p3 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()
    p4 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()


    (p1 + p2)/(p3 + p4)


    Idents(sub.obj) = sub.obj$seurat_clusters_split
    idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
    idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

    idents.sel = c('4', '8', '10')

    ## chcek the reference-mapped ids for the rest of clusters
    counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
    counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
    counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
    counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
    counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
    counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

    features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                      'ceh-36', # MSxaap
                      'tbx-2', 'ceh-27') # MSxaaa

    features.sels = c('dmd-4', 'cnd-1', 'swt-3')

    FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

    VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)
    # check info in JM data for specific lineage
    # ee = process.import.Murray.scRNA()
    # murray.ids = unique(ee$lineage)
    #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
    markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
    #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
    #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))

    ids.sel = c('MSxaa')
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
    find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)

    top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

    # to find new marker genes
    top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
    top.markers[top.markers$cluster == '13',]


    FeaturePlot(sub.obj, reduction = 'umap',
                features = c('zig-6', 'frpr-8','C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-4', 'maph-1.2', 'irx-1')) # MSxpapa

    ##########################################
    # update the manual annotation if good marker genes or mapped ids were found
    ##########################################
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
    # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
    #mm = match(colnames(sub.obj), colnames(seurat.obj))
    #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

    # sub.obj$ids.backup = sub.obj$manual.annot.ids

    cluster.assingment = list(#c('0', 'MSxp'),
      c('4', 'MSxaa'),
      c('8', 'MSxaap'),
      c('10', 'MSxaaa')
      #c('3', 'mixture_BWM_terminal_2'),
      #c('4', 'mixture_BWM_terminal_2')
    )

    for(n in 1:length(cluster.assingment)){
      cluster.index = cluster.assingment[[n]][1];
      id2assign =  cluster.assingment[[n]][2];
      cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
      cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
      sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
      seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
      # if(is.na(id2assign)) {
      #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
      # }else{
      #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
      # }
    }

    #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
            pt.size = 2) + NoLegend()

    DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
            na.value = "gray") +
      ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
      scale_colour_hue(drop = FALSE) + NoLegend()


    saveRDS(seurat.obj, file = RDS2save)



    ##########################################
  # Main aim:
  # there are three main lineages: MSxapa, MSxaap, MSxaaa; I start the most separated lineages, MSxapa
  # Notes:
  # we have for the first time identified MSxapa, MSxapap, MSxapaa and MSxapaaa/MSxapaap
  # seurat_clusters_pharynx 22 were not annotated due to its mixed predicted id
  ##########################################
  nb.iteration = 3
  Refine.annotated.ids = FALSE;

  resDir = paste0("results/", version.analysis, '/annoted_pharynx')
  if(!dir.exists(resDir)){dir.create(resDir)}

  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}

  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                 '_cleanedBWM_and_Pharynx_iteration_', nb.iteration -1, '.rds')

  RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                     '_cleanedBWM_and_Pharynx_iteration_', nb.iteration, '.rds')

  seurat.obj = readRDS(file = RDSsaved)

  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()


  cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
  cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  ##########################################
  # select subset of cells to annotate
  ##########################################
  # cluster.index = '25'
  # table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
  #
  # table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
  #
  # table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]
  #
  # xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
  # xx[which(xx > 0)]
  #
  # ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
  # xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
  # table(seurat.obj$manual.annot.ids[ii1])
  # #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))

  # select cells with cluster index
  ##########################################
  #cluster.sels = c('29', '32', '35', '40', '42')
  #ids.sel = c()

  cluster.sels = c('22', '23', '3', '0', '21')
  ids.sel = c('MSxaa', 'MSxaap', 'MSxaaa', 'MSxap')
  ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')

  # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]

  if(length(ids.sel)>0){
    cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
                                             !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
                                               is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
  }else{
    cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  }

  #ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
  #ids.sels = ids.current[which(nchar(ids.current)>6)]
  #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')

  # ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)],
  #                    c('MSxppaa'))
  #
  # ids.left = setdiff(ids.current, ids.sels)
  # print(ids.left)
  # nchar(ids.left)

  #cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

  ##########################################
  # subset seurat object with selected cells
  ##########################################
  cat(length(cells.sels), ' cells selected to annotate \n')
  sub.obj = subset(seurat.obj, cells = cells.sels)

  sub.obj$seurat_clusters = sub.obj$seurat_clusters_pharynx
  #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  sub.obj$pred.ids = sub.obj$predicted.ids.seurat
  xx = table(sub.obj$predicted.ids.seurat.keep)
  xx[xx>10]
  sub.obj$pred.ids.filtered = sub.obj$pred.ids
  sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()


  ##########################################
  # check potential ids for selected clusters
  ##########################################
  counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters))
  barplot(counts, main="cluster compositions by scmap ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))

  #counts[, match(c('31', '28', '52'), colnames(counts))]
  #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  Explore.umap.parameters.for.BWMcells = FALSE
  if(Explore.umap.parameters.for.BWMcells){

    source.my.script('scRNA_functions.R')
    require(tictoc)
    tic()
    test.umap.params.for.BWM.cells(sub.obj,
                                   pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                   group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                   nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                   n.neighbors.sampling = c(10, 30, 50),
                                   min.dist.sampling = c(0.01, 0.1)
    )
    toc()

  }

  nfeatures = 1000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 10;
  min.dist = 0.1; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist, verbose = TRUE)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
    NoLegend()


  p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend()

  #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
  p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend() + ggtitle('scamp prediction')

  p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
    NoLegend() + ggtitle('scran prediction')

  p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend()

  p0 + p2

  p1 + p2

  p2 + p3

  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p1  = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    ggtitle('manual.annot.ids') + NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 6, na.value = "gray", combine = TRUE)

  p1 + p2  #+ ggsave(paste0(resDir, '/UMAP_pharynx_seurat_prediction_reclustering_newBase.pdf'), width = 22, height = 10)


  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
               group.by = 'seurat_clusters_split')

  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
               group.by = 'seurat_clusters_split') + NoLegend()

  VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
          group.by = 'manual.annot.ids') + NoLegend()

  p1 + p4
  p2 + p4
  plot(p3)

  seurat.obj$seurat_clusters_pharynx = NA
  seurat.obj$seurat_clusters_pharynx[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$seurat_clusters_split

  dev.off()

  ##########################################
  # check the counts of predicted ids for newly split clusters
  ##########################################
  sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
  sub.obj$predicted.ids.prob = sub.obj$predicted.scores
  sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
  sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

  Idents(sub.obj) = sub.obj$seurat_clusters_split
  counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
  counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
  #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)

  p1 = DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
  p2 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()
  p3 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()
  p4 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()


  (p1 + p2)/(p3 + p4)


  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
  idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

  idents.sel = c('8', '1', '0', '2', '5')

  ## chcek the reference-mapped ids for the rest of clusters
  counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
  counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
  counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
  counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
  counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
  counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

  features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                    'ceh-36', # MSxaap
                    'tbx-2', 'ceh-27') # MSxaaa

  features.sels = c(#'dmd-4', 'cnd-1', 'swt-3', # MSxapa
                    #'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                    #'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                    'gst-20', 'sdz-21', 'ceh-13', #MSxapaap
                    'pha-2', 'ceh-22' # MSxapaaa
  )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)
  # check info in JM data for specific lineage
  # ee = process.import.Murray.scRNA()
  # murray.ids = unique(ee$lineage)
  #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
  #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))

  ids.sel = c('MSxapaap')
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')
  find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)

  top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

  # to find new marker genes
  top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
  top.markers[top.markers$cluster == '13',]


  FeaturePlot(sub.obj, reduction = 'umap',
              features = c('zig-6', 'frpr-8','C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-4', 'maph-1.2', 'irx-1')) # MSxpapa

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
  # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
  #mm = match(colnames(sub.obj), colnames(seurat.obj))
  #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

  # sub.obj$ids.backup = sub.obj$manual.annot.ids

  cluster.assingment = list(#c('0', 'MSxp'),
    c('0', 'MSxapaa'),
    c('1', 'MSxapa'),
    c('5', 'MSxapap'),
    c('2', 'MSxapaap/MSxapaaa')
    #c('3', 'mixture_BWM_terminal_2'),
    #c('4', 'mixture_BWM_terminal_2')
  )

  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    # if(is.na(id2assign)) {
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
    # }else{
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
    # }
  }

  #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
          pt.size = 2) + NoLegend()

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE) + NoLegend()


  saveRDS(seurat.obj, file = RDS2save)


  ##########################################
   # Main aim:
   # try to annotate the middle stage of MSxaap
   # one daughter MSxaapa is symetric; but the other daughter MSxaapp is asymetric, because the cell death of MSpaapp;
   # So there are MSaaapp and following single granddaughter MSaaappp, because the other one is also dead
   #
   # Notes:
   # MSxaapa is well annotated for sure; but cluster_pharyx 8 may have some other ids;
   # we have a mixture of MSxaapap/MSpaapaa/MSaaapaa to tease down later
   ##########################################
   nb.iteration = 4
   Refine.annotated.ids = FALSE;

   resDir = paste0("results/", version.analysis, '/annoted_pharynx')
   if(!dir.exists(resDir)){dir.create(resDir)}

   if(Refine.annotated.ids){by.group = 'manual.annot.ids';
   }else{by.group = 'seurat_clusters'}

   RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                  '_cleanedBWM_and_Pharynx_iteration_', nb.iteration -1, '.rds')

   RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                      '_cleanedBWM_and_Pharynx_iteration_', nb.iteration, '.rds')

   seurat.obj = readRDS(file = RDSsaved)

   pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

   DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
     scale_colour_hue(drop = FALSE) +
     NoLegend()

   DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
     scale_colour_hue(drop = FALSE) +
     NoLegend()


   cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
   cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

   pdf(pdfname, width=18, height = 10)
   par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
   ##########################################
   # select subset of cells to annotate
   ##########################################
   # cluster.index = '25'
   # table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
   #
   # table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
   #
   # table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]
   #
   # xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
   # xx[which(xx > 0)]
   #
   # ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
   # xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
   # table(seurat.obj$manual.annot.ids[ii1])
   # #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
   #FeaturePlot(sub.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))

   # select cells with cluster index
   ##########################################
   #cluster.sels = c('29', '32', '35', '40', '42')

   cluster.sels = c('8', '4')
   ids.sel = c('MSxaa', 'MSxap','MSxaap','MSxaaa')
   ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')

   # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]

   if(length(ids.sel)>0){
     cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
                                              !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
                                                is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
   }else{
     cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
   }

   #ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
   #ids.sels = ids.current[which(nchar(ids.current)>6)]
   #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')

   # ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)],
   #                    c('MSxppaa'))
   #
   # ids.left = setdiff(ids.current, ids.sels)
   # print(ids.left)
   # nchar(ids.left)

   #cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

   ##########################################
   # subset seurat object with selected cells
   ##########################################
   cat(length(cells.sels), ' cells selected to annotate \n')
   sub.obj = subset(seurat.obj, cells = cells.sels)

   sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
   #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
   sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
   sub.obj$pred.ids = sub.obj$predicted.ids.seurat
   xx = table(sub.obj$predicted.ids.seurat.keep)
   xx[xx>10]
   sub.obj$pred.ids.filtered = sub.obj$pred.ids
   sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

   DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()


   ##########################################
   # check potential ids for selected clusters
   ##########################################
   counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters))
   barplot(counts, main="cluster compositions by scmap ",
           xlab=NULL, col=c(1:nrow(counts)), las = 2,
           legend = rownames(counts))

   #counts[, match(c('31', '28', '52'), colnames(counts))]
   #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
   ##########################################
   # find new set of variable genes and redo pca and umap
   ##########################################
   Explore.umap.parameters.for.BWMcells = FALSE
   if(Explore.umap.parameters.for.BWMcells){

     source.my.script('scRNA_functions.R')
     require(tictoc)
     tic()
     test.umap.params.for.BWM.cells(sub.obj,
                                    pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                    group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                    nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                    n.neighbors.sampling = c(10, 30, 50),
                                    min.dist.sampling = c(0.01, 0.1)
     )
     toc()

   }

   nfeatures = 1000;
   sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
   #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
   sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
   sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
   ElbowPlot(sub.obj, ndims = 50)

   nb.pcs = 5 # nb of pcs depends on the considered clusters or ids
   n.neighbors = 10;
   min.dist = 0.1; spread = 1
   sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                      spread = spread, n.neighbors = n.neighbors,
                      min.dist = min.dist, verbose = TRUE)
   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()


   p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()

   #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
   p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend() + ggtitle('scamp prediction')

   p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
     NoLegend() + ggtitle('seurat prediction')

   p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()

   p0 + p2
   p0 + p1
   p1 + p2

   p2 + p3

   ##########################################
   # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
   ##########################################
   FindClusters_subclusters = function(sub.obj, resolution = 0.4)
   {
     sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
     return(sub.obj$seurat_clusters)
   }
   sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
   sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
   DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

   p1  = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
     ggtitle('seurat.pred.ids') + NoLegend()
   p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
                label.size = 6, na.value = "gray", combine = TRUE)

   p1 + p2  #+ ggsave(paste0(resDir, '/UMAP_pharynx_seurat_prediction_reclustering_newBase.pdf'), width = 22, height = 10)


   p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
                group.by = 'seurat_clusters_split')

   p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                group.by = 'seurat_clusters_split') + NoLegend()

   VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
           group.by = 'manual.annot.ids') + NoLegend()

   p1 + p4
   p2 + p4
   plot(p3)

   seurat.obj$seurat_clusters_pharynx = NA
   seurat.obj$seurat_clusters_pharynx[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$seurat_clusters_split

   dev.off()

   ##########################################
   # check the counts of predicted ids for newly split clusters
   ##########################################
   sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
   sub.obj$predicted.ids.prob = sub.obj$predicted.scores
   sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
   sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

   Idents(sub.obj) = sub.obj$seurat_clusters_split
   counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
   counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
   #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
   counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)

   p1 = DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
   p2 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()
   p3 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()
   p4 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()


   (p1 + p2)/(p3 + p4)


   Idents(sub.obj) = sub.obj$seurat_clusters_split
   idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
   idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

   idents.sel = c('8', '1', '0', '2', '5')

   ## chcek the reference-mapped ids for the rest of clusters
   counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
   counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
   counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
   counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
   counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
   counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

   features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                     'ceh-36', # MSxaap
                     'tbx-2', 'ceh-27') # MSxaaa
   features.sels = c(#'dmd-4', 'cnd-1', 'swt-3', # MSxapa
                     #'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                     #'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                     'gst-20', 'sdz-21', 'ceh-13', #MSxapaap
                     'pha-2', 'ceh-22' # MSxapaaa
   )
   features.sels = c('lin-12', 'ceh-36', 'ceh-22','tbx-2', 'tbx-7', 'ceh-34', # MSxaapa
                     'fos-1', 'ceh-36', 'ceh-34', # MSaaapp
                     'lin-12', 'pax-1', # MSxaapap
                     'ngn-1', 'ces-1', # MSpaapaa
                     'sptf-1', 'ceh-22' #MSaaapaa

   )

   FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

   VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)
   # check info in JM data for specific lineage
   # ee = process.import.Murray.scRNA()
   # murray.ids = unique(ee$lineage)
   #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
   markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
   #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
   #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
   source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

   ids.sel = c('MSpaapaa')
   find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)

   top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
   DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()

   # to find new marker genes
   top.markers <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
   top.markers[top.markers$cluster == '13',]


   FeaturePlot(sub.obj, reduction = 'umap',
               features = c('zig-6', 'frpr-8','C14B4.2', 'C06A1.2', 'F48C5.2', 'gst-4', 'maph-1.2', 'irx-1')) # MSxpapa

   ##########################################
   # update the manual annotation if good marker genes or mapped ids were found
   ##########################################
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
   # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
   #mm = match(colnames(sub.obj), colnames(seurat.obj))
   #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

   # sub.obj$ids.backup = sub.obj$manual.annot.ids

   cluster.assingment = list(#c('0', 'MSxp'),
     #c('1', 'MSaaap'),
     c('1', 'MSaaapp'),
     c('6', 'MSxaapap/MSpaapaa/MSaaapaa'),
     c('2', 'MSxaapa')
     #c('3', 'mixture_BWM_terminal_2'),
     #c('4', 'mixture_BWM_terminal_2')
   )

   for(n in 1:length(cluster.assingment)){
     cluster.index = cluster.assingment[[n]][1];
     id2assign =  cluster.assingment[[n]][2];
     cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
     cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
     sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
     seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
     # if(is.na(id2assign)) {
     #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
     # }else{
     #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
     # }
   }

   #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
           pt.size = 2) + NoLegend()

   DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
     scale_colour_hue(drop = FALSE) + NoLegend()


   saveRDS(seurat.obj, file = RDS2save)



   ##########################################
  # Main aim:
  # iteration 5 is to tackle lineage MSxaaa that is totally asymmetric
  #
  # Notes:
  # cluster 21 is considered for lineage MSxapa, so it is not included in this iteration
  # four daugthers of MSxaaa were annotated together with one grand daguther (likely)
  ##########################################
  nb.iteration = 5
  Refine.annotated.ids = FALSE;

  resDir = paste0("results/", version.analysis, '/annoted_pharynx')
  if(!dir.exists(resDir)){dir.create(resDir)}

  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}

  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                 '_cleanedBWM_and_Pharynx_iteration_', nb.iteration -1, '.rds')

  RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                     '_cleanedBWM_and_Pharynx_iteration_', nb.iteration, '.rds')

  seurat.obj = readRDS(file = RDSsaved)

  pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()


  cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
  cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

  pdf(pdfname, width=18, height = 10)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  ##########################################
  # select subset of cells to annotate
  ##########################################
  # cluster.index = '25'
  # table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
  #
  # table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
  #
  # table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]
  #
  # xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
  # xx[which(xx > 0)]
  #
  # ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
  # xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
  # table(seurat.obj$manual.annot.ids[ii1])
  # #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
  #FeaturePlot(sub.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))

  # select cells with cluster index
  ##########################################
  #cluster.sels = c('29', '32', '35', '40', '42')

  cluster.sels = c('15', '2', '5', '23')
  ids.sel = c('MSxaaa')
  ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')

  # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]

  if(length(ids.sel)>0){
    cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
                                             !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
                                               is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
  }else{
    cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  }

  #ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
  #ids.sels = ids.current[which(nchar(ids.current)>6)]
  #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')

  # ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)],
  #                    c('MSxppaa'))
  #
  # ids.left = setdiff(ids.current, ids.sels)
  # print(ids.left)
  # nchar(ids.left)

  #cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

  ##########################################
  # subset seurat object with selected cells
  ##########################################
  cat(length(cells.sels), ' cells selected to annotate \n')
  sub.obj = subset(seurat.obj, cells = cells.sels)

  sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
  #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  sub.obj$pred.ids = sub.obj$predicted.ids.seurat
  xx = table(sub.obj$predicted.ids.seurat.keep)
  xx[xx>10]
  sub.obj$pred.ids.filtered = sub.obj$pred.ids
  sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()


  ##########################################
  # check potential ids for selected clusters
  ##########################################
  counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters))
  barplot(counts, main="cluster compositions by scmap ",
          xlab=NULL, col=c(1:nrow(counts)), las = 2,
          legend = rownames(counts))

  #counts[, match(c('31', '28', '52'), colnames(counts))]
  #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  Explore.umap.parameters.for.BWMcells = FALSE
  if(Explore.umap.parameters.for.BWMcells){

    source.my.script('scRNA_functions.R')
    require(tictoc)
    tic()
    test.umap.params.for.BWM.cells(sub.obj,
                                   pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                   group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                   nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                   n.neighbors.sampling = c(10, 30, 50),
                                   min.dist.sampling = c(0.01, 0.1)
    )
    toc()

  }

  nfeatures = 1000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 5 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 20;
  min.dist = 0.1; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist, verbose = TRUE)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
    NoLegend()


  p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend()

  #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
  p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend() + ggtitle('scamp prediction')

  p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
    NoLegend() + ggtitle('seurat prediction')

  p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend()

  p0 + p2
  p0 + p1
  p1 + p2

  p2 + p3

  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 0.4)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:5, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p1  = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
    ggtitle('seurat.pred.ids') + NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
               label.size = 6, na.value = "gray", combine = TRUE)

  p1 + p2  #+ ggsave(paste0(resDir, '/UMAP_pharynx_seurat_prediction_reclustering_newBase.pdf'), width = 22, height = 10)


  p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
               group.by = 'seurat_clusters_split')

  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
               group.by = 'seurat_clusters_split') + NoLegend()

  VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
          group.by = 'manual.annot.ids') + NoLegend()

  p1 + p4
  p2 + p4
  plot(p3)

  seurat.obj$seurat_clusters_pharynx = NA
  seurat.obj$seurat_clusters_pharynx[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$seurat_clusters_split

  #dev.off()

  ##########################################
  # check the counts of predicted ids for newly split clusters
  ##########################################
  sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
  sub.obj$predicted.ids.prob = sub.obj$predicted.scores
  sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
  sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

  Idents(sub.obj) = sub.obj$seurat_clusters_split
  counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
  counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
  #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)

  p1 = DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label =TRUE, label.size = 6, repel = TRUE)
  p2 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label =TRUE, label.size = 4, repel = TRUE) + NoLegend()
  p3 = DimPlot(sub.obj, group.by = 'predicted.ids', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()
  p4 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label =TRUE, label.size = 5, repel = TRUE) + NoLegend()


  (p1 + p2)/(p3 + p4)


  Idents(sub.obj) = sub.obj$seurat_clusters_split
  idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
  idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

  idents.sel = c('8', '1', '0', '2', '5')

  ## chcek the reference-mapped ids for the rest of clusters
  counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
  counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
  counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
  counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
  counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
  counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

  features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                    'ceh-36', # MSxaap
                    'tbx-2', 'ceh-27') # MSxaaa
  features.sels = c(#'dmd-4', 'cnd-1', 'swt-3', # MSxapa
                    #'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                    #'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                    'gst-20', 'sdz-21', 'ceh-13', #MSxapaap
                    'pha-2', 'ceh-22' # MSxapaaa
  )
  features.sels = c('lin-12', 'ceh-36', 'ceh-22','tbx-2', 'tbx-7', 'ceh-34', # MSxaapa
                    'fos-1', 'ceh-36', 'ceh-34', # MSaaapp
                    'lin-12', 'pax-1', # MSxaapap
                    'ngn-1', 'ces-1', # MSpaapaa
                    'sptf-1', 'ceh-22' #MSaaapaa

  )

  features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                    'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                    #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                    'ceh-13', # MSaaaap not MSpaaap
                    'hlh-3', 'fem-1'

                    )

  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  #VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)

  # update the annotation with marker genes
  cluster.assingment = list(#c('0', 'MSxp'),
    c('2', 'MSxaaa'),
    c('1', 'MSaaaaa'),
    c('0', 'MSpaaaa'),
    c('5', 'MSxaaa.like.or.others'),
    c('3', 'MSpaaap'),
    c('7', 'MSaaaap'),
    c('4', 'MSaaaaap.others')
    #c('3', 'mixture_BWM_terminal_2'),
    #c('4', 'mixture_BWM_terminal_2')
  )

  # check info in JM data for specific lineage
  # ee = process.import.Murray.scRNA()
  # murray.ids = unique(ee$lineage)
  #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
  #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

  ids.sel = c('MSpaaaa'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)



  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
  # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
  #mm = match(colnames(sub.obj), colnames(seurat.obj))
  #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

  # sub.obj$ids.backup = sub.obj$manual.annot.ids

  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    # if(is.na(id2assign)) {
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
    # }else{
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
    # }
  }

  #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
          pt.size = 2) + NoLegend()

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE) + NoLegend()


  saveRDS(seurat.obj, file = RDS2save)



  ##########################################
    # Main aim: try to annotate cluster_pharynx 6, 20, 24, 7
    #
    # Notes:
    # MSpaaappp.MSxapappa and MSxaapapa.ABalpappapp.or.others.not.sure to verify later,
    # because the marker genes were not working greatly to comfirm the predicted annotations
    ##########################################
    nb.iteration = 6
    Refine.annotated.ids = FALSE;

    resDir = paste0("results/", version.analysis, '/annoted_pharynx')
    if(!dir.exists(resDir)){dir.create(resDir)}

    if(Refine.annotated.ids){by.group = 'manual.annot.ids';
    }else{by.group = 'seurat_clusters'}

    RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                   '_cleanedBWM_and_Pharynx_iteration_', nb.iteration -1, '.rds')

    RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                       '_cleanedBWM_and_Pharynx_iteration_', nb.iteration, '.rds')

    seurat.obj = readRDS(file = RDSsaved)

    pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

    DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
            na.value = "gray") +
      ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
      scale_colour_hue(drop = FALSE) +
      NoLegend()

    DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
            na.value = "gray") +
      ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
      scale_colour_hue(drop = FALSE) +
      NoLegend()

    cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
    cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

    #pdf(pdfname, width=18, height = 10)
    #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    ##########################################
    # select subset of cells to annotate
    ##########################################
    # cluster.index = '25'
    # table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
    #
    # table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
    #
    # table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]
    #
    # xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
    # xx[which(xx > 0)]
    #
    # ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
    # xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
    # table(seurat.obj$manual.annot.ids[ii1])
    # #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
    #FeaturePlot(sub.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))

    # select cells with cluster index
    ##########################################
    #cluster.sels = c('29', '32', '35', '40', '42')
    cluster.sels = c('6', '24', '20', '7')
    #ids.sel = c('MSxaaa')
    #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')

    cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
    #
    # if(length(ids.sel)>0){
    #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
    #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
    #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
    # }else{
    #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
    # }

    #ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
    #ids.sels = ids.current[which(nchar(ids.current)>6)]
    #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')

    # ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)],
    #                    c('MSxppaa'))
    #
    # ids.left = setdiff(ids.current, ids.sels)
    # print(ids.left)
    # nchar(ids.left)

    #cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

    ##########################################
    # subset seurat object with selected cells
    ##########################################
    cat(length(cells.sels), ' cells selected to annotate \n')
    sub.obj = subset(seurat.obj, cells = cells.sels)

    sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
    #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
    sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
    sub.obj$pred.ids = sub.obj$predicted.ids.seurat
    xx = table(sub.obj$predicted.ids.seurat.keep)
    xx[xx>10]
    sub.obj$pred.ids.filtered = sub.obj$pred.ids
    sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

    DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()


    ##########################################
    # check potential ids for selected clusters
    ##########################################
    counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters))
    barplot(counts, main="cluster compositions by scmap ",
            xlab=NULL, col=c(1:nrow(counts)), las = 2,
            legend = rownames(counts))

    #counts[, match(c('31', '28', '52'), colnames(counts))]
    #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
    ##########################################
    # find new set of variable genes and redo pca and umap
    ##########################################
    Explore.umap.parameters.for.BWMcells = FALSE
    if(Explore.umap.parameters.for.BWMcells){

      source.my.script('scRNA_functions.R')
      require(tictoc)
      tic()
      test.umap.params.for.BWM.cells(sub.obj,
                                     pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                     group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                     nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                     n.neighbors.sampling = c(10, 30, 50),
                                     min.dist.sampling = c(0.01, 0.1)
      )
      toc()

    }

    nfeatures = 1000;
    sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
    #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
    sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
    sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
    ElbowPlot(sub.obj, ndims = 50)

    nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
    n.neighbors = 10;
    min.dist = 0.01; spread = 1
    sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                       spread = spread, n.neighbors = n.neighbors,
                       min.dist = min.dist, verbose = TRUE)
    DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()


    p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend()

    #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
    p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend() + ggtitle('scamp prediction')

    p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
      NoLegend() + ggtitle('seurat prediction')

    p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend()

    p0 + p2
    p0 + p1
    p1 + p2

    p2 + p3

    ##########################################
    # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
    ##########################################
    FindClusters_subclusters = function(sub.obj, resolution = 0.4)
    {
      sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
      return(sub.obj$seurat_clusters)
    }
    sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
    sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
    DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

    p1  = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
      ggtitle('seurat.pred.ids') + NoLegend()
    p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
                 label.size = 8, na.value = "gray", combine = TRUE)

    p1 + p2  #+ ggsave(paste0(resDir, '/UMAP_pharynx_seurat_prediction_reclustering_newBase.pdf'), width = 22, height = 10)


    p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
                 group.by = 'seurat_clusters_split')

    p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                 group.by = 'seurat_clusters_split') + NoLegend()

    VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
            group.by = 'manual.annot.ids') + NoLegend()

    p2 + p4 + ggsave(paste0(resDir, '/splitcluster_timingEstimation_iteration_6.pdf'), width = 18, height = 10)

    p1 + p4


    plot(p3)

    #seurat.obj$seurat_clusters_pharynx = NA
    #seurat.obj$seurat_clusters_pharynx[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$seurat_clusters_split

    #dev.off()

    ##########################################
    # check the counts of predicted ids for newly split clusters
    ##########################################
    sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
    sub.obj$predicted.ids.prob = sub.obj$predicted.scores
    sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
    sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

    Idents(sub.obj) = sub.obj$seurat_clusters_split
    counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
    counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
    #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
    counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)


    Idents(sub.obj) = sub.obj$seurat_clusters_split
    idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
    idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

    idents.sel = c('8', '1', '0', '2', '5')

    ## chcek the reference-mapped ids for the rest of clusters
    counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
    counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
    counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
    counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
    counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
    counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

    features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                      'ceh-36', # MSxaap
                      'tbx-2', 'ceh-27') # MSxaaa
    features.sels = c(#'dmd-4', 'cnd-1', 'swt-3', # MSxapa
                      #'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                      #'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                      'gst-20', 'sdz-21', 'ceh-13', #MSxapaap
                      'pha-2', 'ceh-22' # MSxapaaa
    )
    features.sels = c('lin-12', 'ceh-36', 'ceh-22','tbx-2', 'tbx-7', 'ceh-34', # MSxaapa
                      'fos-1', 'ceh-36', 'ceh-34', # MSaaapp
                      'lin-12', 'pax-1', # MSxaapap
                      'ngn-1', 'ces-1', # MSpaapaa
                      'sptf-1', 'ceh-22' #MSaaapaa

    )

    features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                      'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                      #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                      'ceh-13', # MSaaaap not MSpaaap
                      'hlh-3', 'fem-1'

                      )

    features.sels = c(
      'nfki-1', 'unc-62', 'ser-2', 'tnc-2',
      'pax-1', 'ttx-1', 'agr-1', 'irx-1',
       'aff-1',  'ceh-27', 'unc-129', 'ceh-22'
    )
    FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

    #VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)
    # update the annotation with marker genes
    cluster.assingment = list(
      c('3', 'MSaaapaa'),
      c('0', 'MSxapaaa'),
      c('7', 'MSaaaappp/MSxapaapp'),
      c('1', 'MSpaaaap'),
      c('5', 'MSxapapp'),
      c('2', 'MSpaaappp.MSxapappa'),
      c('6', 'MSpaaappp.MSxapappa'),
      c('4', 'MSxaapapa.ABalpappapp.or.others.not.sure')

    )

    # check info in JM data for specific lineage
    # ee = process.import.Murray.scRNA()
    # murray.ids = unique(ee$lineage)
    #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
    markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
    #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
    #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

    ids.sel = c('MSpaaappp.MSxapappa'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


    ##########################################
    # update the manual annotation if good marker genes or mapped ids were found
    ##########################################
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
    # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
    #mm = match(colnames(sub.obj), colnames(seurat.obj))
    #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

    # sub.obj$ids.backup = sub.obj$manual.annot.ids

    for(n in 1:length(cluster.assingment)){
      cluster.index = cluster.assingment[[n]][1];
      id2assign =  cluster.assingment[[n]][2];
      cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
      cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
      sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
      seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
      # if(is.na(id2assign)) {
      #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
      # }else{
      #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
      # }
    }

    #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
            pt.size = 2) + NoLegend()

    DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
            na.value = "gray") +
      ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
      scale_colour_hue(drop = FALSE) + NoLegend()


    saveRDS(seurat.obj, file = RDS2save)




    ##########################################
      # Main aim: try to annotate cluster_pharynx 13, 17, 18, 11, 25, 10
      # those well split clusters for terminal cells
      # Notes:
      # second MSpaaappp.MSxapappa and MSaaaappp.MSxapaapp were found more sure
      # those two ids annotated in the previous iteration will be revised later
      ##########################################
      nb.iteration = 7
      Refine.annotated.ids = FALSE;

      resDir = paste0("results/", version.analysis, '/annoted_pharynx')
      if(!dir.exists(resDir)){dir.create(resDir)}

      if(Refine.annotated.ids){by.group = 'manual.annot.ids';
      }else{by.group = 'seurat_clusters'}

      RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                     '_cleanedBWM_and_Pharynx_iteration_', nb.iteration -1, '.rds')

      RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                         '_cleanedBWM_and_Pharynx_iteration_', nb.iteration, '.rds')

      seurat.obj = readRDS(file = RDSsaved)

      pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

      DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
              na.value = "gray") +
        ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
        scale_colour_hue(drop = FALSE) +
        NoLegend()

      DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
              na.value = "gray") +
        ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
        scale_colour_hue(drop = FALSE) +
        NoLegend()

      cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
      cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

      #pdf(pdfname, width=18, height = 10)
      #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
      ##########################################
      # select subset of cells to annotate
      ##########################################
      # cluster.index = '25'
      # table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
      #
      # table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
      #
      # table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]
      #
      # xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
      # xx[which(xx > 0)]
      #
      # ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
      # xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
      # table(seurat.obj$manual.annot.ids[ii1])
      # #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
      #FeaturePlot(sub.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))

      # select cells with cluster index
      ##########################################
      #cluster.sels = c('29', '32', '35', '40', '42')
      #cluster.sels = c('6', '24', '20', '7')
      cluster.sels = c('13', '17', '18', '11', '25', '10')
      #ids.sel = c('MSxaaa')
      #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')

      cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
      #
      # if(length(ids.sel)>0){
      #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
      #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
      #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
      # }else{
      #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
      # }

      #ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
      #ids.sels = ids.current[which(nchar(ids.current)>6)]
      #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')

      # ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)],
      #                    c('MSxppaa'))
      #
      # ids.left = setdiff(ids.current, ids.sels)
      # print(ids.left)
      # nchar(ids.left)

      #cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

      ##########################################
      # subset seurat object with selected cells
      ##########################################
      cat(length(cells.sels), ' cells selected to annotate \n')
      sub.obj = subset(seurat.obj, cells = cells.sels)

      sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
      #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
      sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
      sub.obj$pred.ids = sub.obj$predicted.ids.seurat
      xx = table(sub.obj$predicted.ids.seurat.keep)
      xx[xx>10]
      sub.obj$pred.ids.filtered = sub.obj$pred.ids
      sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

      DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()


      ##########################################
      # check potential ids for selected clusters
      ##########################################
      counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters))
      barplot(counts, main="cluster compositions by scmap ",
              xlab=NULL, col=c(1:nrow(counts)), las = 2,
              legend = rownames(counts))

      #counts[, match(c('31', '28', '52'), colnames(counts))]
      #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
      ##########################################
      # find new set of variable genes and redo pca and umap
      ##########################################
      Explore.umap.parameters.for.BWMcells = FALSE
      if(Explore.umap.parameters.for.BWMcells){

        source.my.script('scRNA_functions.R')
        require(tictoc)
        tic()
        test.umap.params.for.BWM.cells(sub.obj,
                                       pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                       group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                       nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                       n.neighbors.sampling = c(10, 30, 50),
                                       min.dist.sampling = c(0.01, 0.1)
        )
        toc()

      }

      nfeatures = 1000;
      sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
      #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
      sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
      sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
      ElbowPlot(sub.obj, ndims = 50)

      nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
      n.neighbors = 10;
      min.dist = 0.01; spread = 1
      sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                         spread = spread, n.neighbors = n.neighbors,
                         min.dist = min.dist, verbose = TRUE)
      DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
        NoLegend()


      p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
        NoLegend()

      #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
      p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
        NoLegend() + ggtitle('scamp prediction')

      p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
        NoLegend() + ggtitle('seurat prediction')

      p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
        NoLegend()

      p0 + p2
      p0 + p1
      p1 + p2

      p2 + p3

      ##########################################
      # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
      ##########################################
      FindClusters_subclusters = function(sub.obj, resolution = 1.0)
      {
        sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
        return(sub.obj$seurat_clusters)
      }
      sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
      sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.5)
      DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

      p1  = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
        ggtitle('seurat.pred.ids') + NoLegend()
      p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
                   label.size = 8, na.value = "gray", combine = TRUE)

      p1 + p2  #+ ggsave(paste0(resDir, '/UMAP_pharynx_seurat_prediction_reclustering_newBase.pdf'), width = 22, height = 10)


      p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
                   group.by = 'seurat_clusters_split')

      p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                   group.by = 'seurat_clusters_split') + NoLegend()

      VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
              group.by = 'manual.annot.ids') + NoLegend()

      p2 + p4 + ggsave(paste0(resDir, '/splitcluster_timingEstimation_iteration_', nb.iteration, '.pdf'), width = 18, height = 10)

      p1 + p4


      plot(p3)

      #seurat.obj$seurat_clusters_pharynx = NA
      #seurat.obj$seurat_clusters_pharynx[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$seurat_clusters_split

      #dev.off()

      ##########################################
      # check the counts of predicted ids for newly split clusters
      ##########################################
      sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
      sub.obj$predicted.ids.prob = sub.obj$predicted.scores
      sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
      sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

      Idents(sub.obj) = sub.obj$seurat_clusters_split
      counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
      counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
      #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
      counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)

      Idents(sub.obj) = sub.obj$seurat_clusters_split
      idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
      idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

      idents.sel = c('8', '1', '0', '2', '5')

      ## chcek the reference-mapped ids for the rest of clusters
      counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
      counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
      counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
      counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
      counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
      counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

      features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                        'ceh-36', # MSxaap
                        'tbx-2', 'ceh-27') # MSxaaa
      features.sels = c(#'dmd-4', 'cnd-1', 'swt-3', # MSxapa
                        #'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                        #'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                        'gst-20', 'sdz-21', 'ceh-13', #MSxapaap
                        'pha-2', 'ceh-22' # MSxapaaa
      )
      features.sels = c('lin-12', 'ceh-36', 'ceh-22','tbx-2', 'tbx-7', 'ceh-34', # MSxaapa
                        'fos-1', 'ceh-36', 'ceh-34', # MSaaapp
                        'lin-12', 'pax-1', # MSxaapap
                        'ngn-1', 'ces-1', # MSpaapaa
                        'sptf-1', 'ceh-22' #MSaaapaa

      )

      features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                        'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                        #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                        'ceh-13', # MSaaaap not MSpaaap
                        'hlh-3', 'fem-1'
                        )
      features.sels = c(
        'nfki-1', 'unc-62', 'ser-2', 'tnc-2',
        'pax-1', 'ttx-1', 'agr-1', 'irx-1',
         'aff-1',  'ceh-27', 'unc-129', 'ceh-22'
      )

      features.sels = c('ceh-32', 'ceh-34', 'tbx-7', 'ngn-1', 'ces-1')

      features.sels = c('nhr-67', 'pha-4', 'cwn-2',
                        'ttx-1', 'pax-1', 'agr-1', 'irx-1', 'lin-12')
      features.sels = c('hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', 'ceh-53', 'dmd-4', 'asp-4', 'ceh-6', 'C39E9.8')

      features.sels = c('nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1',
                        'ceh-22', 'spp-7')
      FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

      #VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)
      # update the annotation with marker genes
      cluster.assingment = list(
        c('0', 'MSpaapaa'),
        c('1', 'MSxaapapa.ABalpappapp'),
        c('3', 'MSaappa'),
        c('7', 'MSaappa'),
        c('8', 'MSxapaapa'),
        c('2', 'MSaaaaapa.and.MSxapaapa'),
        c('5', 'MSxapapa.and.MSxapapaa'),
        c('4', 'MSaaaappp.MSxapaapp.sure'),
        c('6', 'MSpaaappp.MSxapappa.sure')

      )

      # check info in JM data for specific lineage
      # ee = process.import.Murray.scRNA()
      # murray.ids = unique(ee$lineage)
      #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
      markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
      #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
      #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
      source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

      ids.sel = c('MSaaaaapa'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


      ##########################################
      # update the manual annotation if good marker genes or mapped ids were found
      ##########################################
      # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
      # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
      # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
      # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

      # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
      # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
      #mm = match(colnames(sub.obj), colnames(seurat.obj))
      #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

      # sub.obj$ids.backup = sub.obj$manual.annot.ids

      for(n in 1:length(cluster.assingment)){
        cluster.index = cluster.assingment[[n]][1];
        id2assign =  cluster.assingment[[n]][2];
        cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
        cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
        sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
        seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
        # if(is.na(id2assign)) {
        #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
        # }else{
        #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
        # }
      }

      #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

      DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
              pt.size = 2) + NoLegend()

      DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
              na.value = "gray") +
        ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
        scale_colour_hue(drop = FALSE) + NoLegend()


      saveRDS(seurat.obj, file = RDS2save)



      ##########################################
        # Main aim:
        # try to annotate pharynx clusters 1, 12, 9, 19
        # Notes:
        #
        # many unknown mixture appear and also previously predicted ids show up in this iteration again while some ids never found
        # this iteration is the last one for the first round of rough manual annotation
        ##########################################
        nb.iteration = 8
        Refine.annotated.ids = FALSE;

        resDir = paste0("results/", version.analysis, '/annoted_pharynx')
        if(!dir.exists(resDir)){dir.create(resDir)}

        if(Refine.annotated.ids){by.group = 'manual.annot.ids';
        }else{by.group = 'seurat_clusters'}

        RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                       '_cleanedBWM_and_Pharynx_iteration_', nb.iteration -1, '.rds')

        RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                           '_cleanedBWM_and_Pharynx_iteration_', nb.iteration, '.rds')

        seurat.obj = readRDS(file = RDSsaved)

        #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

        DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                na.value = "gray") +
          ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
          scale_colour_hue(drop = FALSE) +
          NoLegend()

        DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                na.value = "gray") +
          ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
          scale_colour_hue(drop = FALSE) +
          NoLegend()

        cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
        cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

        #pdf(pdfname, width=18, height = 10)
        #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
        ##########################################
        # select subset of cells to annotate
        ##########################################
        # cluster.index = '25'
        # table(seurat.obj$manual.annot.ids[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
        #
        # table(seurat.obj$predicted.ids.seurat[seurat.obj$seurat_clusters == cluster.index], useNA = 'ifany')
        #
        # table(seurat.obj$manual.annot.ids)[grep('MSppaap', names(table(seurat.obj$manual.annot.ids)))]
        #
        # xx = table(seurat.obj$seurat_clusters[which(seurat.obj$manual.annot.ids == 'MSxpa')])
        # xx[which(xx > 0)]
        #
        # ii1 = which(seurat.obj$predicted.ids.seurat == 'MSxpa')
        # xx = table(seurat.obj$seurat_clusters[ii1]); xx[which(xx > 0)]
        # table(seurat.obj$manual.annot.ids[ii1])
        # #FeaturePlot(seurat.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))
        #FeaturePlot(sub.obj, reduction = 'umap', features = c('lin-39', 'clec-264', 'zig-6', 'ceh-34'))

        # select cells with cluster index
        ##########################################
        #cluster.sels = c('29', '32', '35', '40', '42')
        #cluster.sels = c('6', '24', '20', '7')
        #cluster.sels = c('13', '17', '18', '11', '25', '10')
        cluster.sels = c('1', '9', '12', '19')
        #ids.sel = c('MSxaaa')
        #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')

        cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
        #
        # if(length(ids.sel)>0){
        #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
        #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
        #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
        # }else{
        #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
        # }

        #ids.current = names(table(seurat.obj$manual.annot.ids[!is.na(seurat.obj$BWM.cells)], useNA = 'ifany'))
        #ids.sels = ids.current[which(nchar(ids.current)>6)]
        #ids.sels = c('MSx', 'MSxp', 'MSxa', 'MSxpp', 'MSxpa', 'MSxap')

        # ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)],
        #                    c('MSxppaa'))
        #
        # ids.left = setdiff(ids.current, ids.sels)
        # print(ids.left)
        # nchar(ids.left)

        #cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

        ##########################################
        # subset seurat object with selected cells
        ##########################################
        cat(length(cells.sels), ' cells selected to annotate \n')
        sub.obj = subset(seurat.obj, cells = cells.sels)

        sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
        #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
        sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
        sub.obj$pred.ids = sub.obj$predicted.ids.seurat
        xx = table(sub.obj$predicted.ids.seurat.keep)
        xx[xx>10]
        sub.obj$pred.ids.filtered = sub.obj$pred.ids
        sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

        DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()


        ##########################################
        # check potential ids for selected clusters
        ##########################################
        counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters))
        barplot(counts, main="cluster compositions by scmap ",
                xlab=NULL, col=c(1:nrow(counts)), las = 2,
                legend = rownames(counts))

        #counts[, match(c('31', '28', '52'), colnames(counts))]
        #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
        ##########################################
        # find new set of variable genes and redo pca and umap
        ##########################################
        Explore.umap.parameters.for.BWMcells = FALSE
        if(Explore.umap.parameters.for.BWMcells){

          source.my.script('scRNA_functions.R')
          require(tictoc)
          tic()
          test.umap.params.for.BWM.cells(sub.obj,
                                         pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                         group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                         nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                         n.neighbors.sampling = c(10, 30, 50),
                                         min.dist.sampling = c(0.01, 0.1)
          )
          toc()

        }

        nfeatures = 1000;
        sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
        #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
        sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
        sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
        ElbowPlot(sub.obj, ndims = 50)

        nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
        n.neighbors = 20;
        min.dist = 0.01; spread = 1
        sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                           spread = spread, n.neighbors = n.neighbors,
                           min.dist = min.dist, verbose = TRUE)
        DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
          NoLegend()


        p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
          NoLegend()

        #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
        p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
          NoLegend() + ggtitle('scamp prediction')

        p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
          NoLegend() + ggtitle('seurat prediction')

        p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
          NoLegend()

        p0 + p2
        p0 + p1
        p1 + p2

        p2 + p3

        ##########################################
        # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
        ##########################################
        FindClusters_subclusters = function(sub.obj, resolution = 1.0)
        {
          sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
          return(sub.obj$seurat_clusters)
        }
        sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
        sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
        DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

        p1  = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
          ggtitle('seurat.pred.ids') + NoLegend()
        p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
                     label.size = 8, na.value = "gray", combine = TRUE)

        p1 + p2  #+ ggsave(paste0(resDir, '/UMAP_pharynx_seurat_prediction_reclustering_newBase.pdf'), width = 22, height = 10)


        p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
                     group.by = 'seurat_clusters_split')

        p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                     group.by = 'seurat_clusters_split') + NoLegend()

        VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                group.by = 'manual.annot.ids') + NoLegend()

        p2 + p4 + ggsave(paste0(resDir, '/splitcluster_timingEstimation_iteration_', nb.iteration, '.pdf'), width = 18, height = 10)

        p1 + p4


        plot(p3)

        #seurat.obj$seurat_clusters_pharynx = NA
        #seurat.obj$seurat_clusters_pharynx[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$seurat_clusters_split

        #dev.off()

        ##########################################
        # check the counts of predicted ids for newly split clusters
        ##########################################
        sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
        sub.obj$predicted.ids.prob = sub.obj$predicted.scores
        sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
        sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

        Idents(sub.obj) = sub.obj$seurat_clusters_split
        counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
        counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
        #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
        counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)

        Idents(sub.obj) = sub.obj$seurat_clusters_split
        idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
        idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

        idents.sel = c('8', '1', '0', '2', '5')

        ## chcek the reference-mapped ids for the rest of clusters
        counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
        counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
        counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
        counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
        counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
        counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

        features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                          'ceh-36', # MSxaap
                          'tbx-2', 'ceh-27') # MSxaaa
        features.sels = c(#'dmd-4', 'cnd-1', 'swt-3', # MSxapa
                          #'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                          #'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                          'gst-20', 'sdz-21', 'ceh-13', #MSxapaap
                          'pha-2', 'ceh-22' # MSxapaaa
        )
        features.sels = c('lin-12', 'ceh-36', 'ceh-22','tbx-2', 'tbx-7', 'ceh-34', # MSxaapa
                          'fos-1', 'ceh-36', 'ceh-34', # MSaaapp
                          'lin-12', 'pax-1', # MSxaapap
                          'ngn-1', 'ces-1', # MSpaapaa
                          'sptf-1', 'ceh-22' #MSaaapaa

        )

        features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                          'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                          #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                          'ceh-13', # MSaaaap not MSpaaap
                          'hlh-3', 'fem-1'
                          )
        features.sels = c(
          'nfki-1', 'unc-62', 'ser-2', 'tnc-2',
          'pax-1', 'ttx-1', 'agr-1', 'irx-1',
           'aff-1',  'ceh-27', 'unc-129', 'ceh-22'
        )

        features.sels = c('ceh-32', 'ceh-34', 'tbx-7', 'ngn-1', 'ces-1')
        features.sels = c('nhr-67', 'pha-4', 'cwn-2',
                          'ttx-1', 'pax-1', 'agr-1', 'irx-1', 'lin-12')
        features.sels = c('hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', 'ceh-53', 'dmd-4', 'asp-4', 'ceh-6', 'C39E9.8')
        features.sels = c('nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1',
                          'ceh-22', 'spp-7')
        features.sels = c('lim-7', 'ces-1', 'ceh-27', 'dod-6') # MSpaaapa

        features.sels = c('ceh-32', 'ceh-34', 'ces-1', 'tbx-7', 'ngn-1',
                          'fem-1', 'C45G7.4', 'K04G2.12')
        features.sels = c('ceh-34', 'ceh-27', 'ceh-32', 'K04G2.12', 'C45G7.4', 'fem-1')

        features.sels = c('srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'ceh-27', 'ceh-22')
        FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

        #VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)
        # update the annotation with marker genes
        cluster.assingment = list(
          c('5', 'MSpaaapa'),
          c('6', 'MSpaapaa'), # for second time
          c('0', 'MSpaaaap'), # second time
          c('3', 'MSpaaaaa.to.confirm'),
          c('4', 'MSpaaaaa.to.confirm'),
          c('2', 'mixture.unknown'),
          c('1', 'mixture.unknown')

        )

        # check info in JM data for specific lineage
        # ee = process.import.Murray.scRNA()
        # murray.ids = unique(ee$lineage)
        #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
        markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
        #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
        #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
        source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

        ids.sel = c('MSpaaapp'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


        ##########################################
        # update the manual annotation if good marker genes or mapped ids were found
        ##########################################
        # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
        # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
        # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
        # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

        # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
        # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
        #mm = match(colnames(sub.obj), colnames(seurat.obj))
        #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

        # sub.obj$ids.backup = sub.obj$manual.annot.ids

        for(n in 1:length(cluster.assingment)){
          cluster.index = cluster.assingment[[n]][1];
          id2assign =  cluster.assingment[[n]][2];
          cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
          cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
          sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
          seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
          # if(is.na(id2assign)) {
          #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
          # }else{
          #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
          # }
        }

        #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

        DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                pt.size = 2) + NoLegend()

        DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                na.value = "gray") +
          ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
          scale_colour_hue(drop = FALSE) + NoLegend()


        saveRDS(seurat.obj, file = RDS2save)


        ##########################################
          # Main aim:
          # iteration R1, beginning of second round of annotations
          # start with the branch of MSxapa
          # this iteration focuses on the beginning of MSxapa with ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
          # Notes:
          # 1) newly annotated MSaaappp not quite sure due to the lacking of confident markers
          # 2) MSxapa and MSxapaa are still a little confusing
          ##########################################
          Rnb.iteration = 1
          Refine.annotated.ids = TRUE;

          resDir = paste0("results/", version.analysis, '/annoted_pharynx')
          if(!dir.exists(resDir)){dir.create(resDir)}

          if(Refine.annotated.ids){by.group = 'manual.annot.ids';
          }else{by.group = 'seurat_clusters'}

          RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                         '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration -1, '.rds')

          RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                             '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration, '.rds')

          seurat.obj = readRDS(file = RDSsaved)

          #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

          DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                  na.value = "gray") +
            ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
            scale_colour_hue(drop = FALSE) +
            NoLegend()

          DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                  na.value = "gray") +
            ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
            scale_colour_hue(drop = FALSE) +
            NoLegend()

          cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
          cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

          #pdf(pdfname, width=18, height = 10)
          #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

          ##########################################
          # select subset of cells
          ##########################################

          # select cells with cluster index
          ##########################################
          #cluster.sels = c('29', '32', '35', '40', '42')
          #cluster.sels = c('6', '24', '20', '7')
          #cluster.sels = c('13', '17', '18', '11', '25', '10')
          cluster.sels = c('1', '9', '12', '19')
          #ids.sel = c('MSxaaa')
          #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')

          cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
          #
          # if(length(ids.sel)>0){
          #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
          #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
          #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
          # }else{
          #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
          # }

          # select cells with ids
          ##########################################
          index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
          table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

          ids.current = names(table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany'))
          ids.current = ids.current[grep('MSxapa', ids.current)]
          ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
          #ids.sels = ids.current
          # ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)],
          #                    c('MSxppaa'))
          #
          ids.left = setdiff(ids.current, ids.sels)
          # print(ids.left)
          # nchar(ids.left)

          cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels)) |
                                                     (seurat.obj$seurat_clusters_pharynx == '22')])

          ##########################################
          # subset seurat object with selected cells
          ##########################################
          cat(length(cells.sels), ' cells selected to annotate \n')
          sub.obj = subset(seurat.obj, cells = cells.sels)

          sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
          #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
          sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
          sub.obj$pred.ids = sub.obj$predicted.ids.seurat
          xx = table(sub.obj$predicted.ids.seurat.keep)
          xx[xx>10]
          sub.obj$pred.ids.filtered = sub.obj$pred.ids
          sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

          DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()


          ##########################################
          # check potential ids for selected clusters
          ##########################################
          counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters))
          barplot(counts, main="cluster compositions by scmap ",
                  xlab=NULL, col=c(1:nrow(counts)), las = 2,
                  legend = rownames(counts))

          #counts[, match(c('31', '28', '52'), colnames(counts))]
          #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
          ##########################################
          # find new set of variable genes and redo pca and umap
          ##########################################
          Explore.umap.parameters.for.BWMcells = FALSE
          if(Explore.umap.parameters.for.BWMcells){

            source.my.script('scRNA_functions.R')
            require(tictoc)
            tic()
            test.umap.params.for.BWM.cells(sub.obj,
                                           pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                           group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                           nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                           n.neighbors.sampling = c(10, 30, 50),
                                           min.dist.sampling = c(0.01, 0.1)
            )
            toc()

          }

          nfeatures = 1000;
          sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
          #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
          sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
          sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
          ElbowPlot(sub.obj, ndims = 50)

          nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
          n.neighbors = 10;
          min.dist = 0.01; spread = 1
          sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                             spread = spread, n.neighbors = n.neighbors,
                             min.dist = min.dist, verbose = TRUE)
          DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
            NoLegend()

          DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
            NoLegend()

          p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
            NoLegend()

          #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
          p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
            NoLegend() + ggtitle('scamp prediction')

          p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
            NoLegend() + ggtitle('seurat prediction')

          p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
            NoLegend()

          p0 + p2
          p0 + p1
          p1 + p2

          p2 + p3

          ##########################################
          # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
          ##########################################
          FindClusters_subclusters = function(sub.obj, resolution = 1.0)
          {
            sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
            return(sub.obj$seurat_clusters)
          }
          sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:10, compute.SNN = TRUE)
          sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
          DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

          p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
            ggtitle('manual.ids') + NoLegend()
          p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
                       label.size = 8, na.value = "gray", combine = TRUE)

          p2 + p1 #+ ggsave(paste0(resDir, '/splitcluster_manual_IDs_iteration_', Rnb.iteration, '.pdf'), width = 18, height = 10)


          p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
                       group.by = 'seurat_clusters_split')

          p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                       group.by = 'seurat_clusters_split') + NoLegend()

          VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                  group.by = 'manual.annot.ids') + NoLegend()

          (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_R', Rnb.iteration, '.pdf'),
                                   width = 18, height = 16)

          p1 + p4


          plot(p3)

          #seurat.obj$seurat_clusters_pharynx = NA
          #seurat.obj$seurat_clusters_pharynx[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$seurat_clusters_split

          #dev.off()

          ##########################################
          # check the counts of predicted ids for newly split clusters
          ##########################################
          sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
          sub.obj$predicted.ids.prob = sub.obj$predicted.scores
          sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
          sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

          Idents(sub.obj) = sub.obj$seurat_clusters_split
          counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
          counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
          #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
          counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)

          Idents(sub.obj) = sub.obj$seurat_clusters_split
          idents.sel = as.character(levels(sub.obj$seurat_clusters_split))
          idents.sel = setdiff(idents.sel, c('0', '2', '5', '4', '7', '8', '3', '12', '13', '14', '11'))

          idents.sel = c('8', '1', '0', '2', '5')

          ## chcek the reference-mapped ids for the rest of clusters
          counts.sel = counts[, !is.na(match(colnames(counts), idents.sel))]
          counts.sel = counts.sel[apply(as.matrix(counts.sel), 1, sum) >0, ]
          counts.seurat.sel = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
          counts.seurat.sel = counts.seurat.sel[apply(as.matrix(counts.seurat.sel), 1, sum)>0, ]
          counts.annot.sel = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
          counts.annot.sel = counts.annot.sel[apply(as.matrix(counts.annot.sel), 1, sum) >0, ]

          features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                            'ceh-36', # MSxaap
                            'tbx-2', 'ceh-27') # MSxaaa
          features.sels = c('lin-12', 'ceh-36', 'ceh-22','tbx-2', 'tbx-7', 'ceh-34', # MSxaapa
                            'fos-1', 'ceh-36', 'ceh-34', # MSaaapp
                            'lin-12', 'pax-1', # MSxaapap
                            'ngn-1', 'ces-1', # MSpaapaa
                            'sptf-1', 'ceh-22' #MSaaapaa

          )

          features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                            'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                            #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                            'ceh-13', # MSaaaap not MSpaaap
                            'hlh-3', 'fem-1'
                            )

          features.sels = c('dmd-4', 'cnd-1', 'swt-3', # MSxapa
            'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
            'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
            'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
            'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1' # MSxapaaa
          )

          features.sels = c('ceh-36', 'ceh-32', 'tbx-2', 'ceh-27', 'hphd-1', 'let-381',
                            'fos-1', 'jun-1', 'ttx-1'
                            )
          FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

          #VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)
          # update the annotation with marker genes


          cluster.assingment = list(
            c('2', 'MSxapa'),
            c('8', 'MSxapa'),
            c('1', 'MSxapap'),
            c('3', 'MSxapaa'),
            c('4', 'MSxapaa'),

            c('0', 'MSxapaap'),
            c('6', 'MSxapaap'),
            c('5', 'MSxapaaa'),
            c('7', 'MSaaappp')
          )


          # check info in JM data for specific lineage
          # ee = process.import.Murray.scRNA()
          # murray.ids = unique(ee$lineage)
          #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
          markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
          #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
          #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
          source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

          ids.sel = c('MSaaapp'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


          ##########################################
          # update the manual annotation if good marker genes or mapped ids were found
          ##########################################
          # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
          # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
          # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
          # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

          # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
          # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
          #mm = match(colnames(sub.obj), colnames(seurat.obj))
          #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

          # sub.obj$ids.backup = sub.obj$manual.annot.ids

          for(n in 1:length(cluster.assingment)){
            cluster.index = cluster.assingment[[n]][1];
            id2assign =  cluster.assingment[[n]][2];
            cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
            cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
            sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
            seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
            # if(is.na(id2assign)) {
            #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
            # }else{
            #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
            # }
          }

          #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

          DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                  pt.size = 2) + NoLegend()

          DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                  na.value = "gray") +
            ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
            scale_colour_hue(drop = FALSE) + NoLegend()


          saveRDS(seurat.obj, file = RDS2save)




          ##########################################
            # Main aim:
            # iteration R2,
            # still with the branch of MSxapa
            # this iteration focuses on MSxapa with ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
            # Notes:
            # MSxapa is almost done;
            # in particular, MSxapapa and daughter cell can be well separated with ces-1 and C39E9.8, but now they were just kept together.
            ##########################################
            Rnb.iteration = 2
            Refine.annotated.ids = TRUE;

            resDir = paste0("results/", version.analysis, '/annoted_pharynx')
            if(!dir.exists(resDir)){dir.create(resDir)}

            if(Refine.annotated.ids){by.group = 'manual.annot.ids';
            }else{by.group = 'seurat_clusters'}

            RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                           '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration -1, '.rds')

            RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                               '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration, '.rds')

            seurat.obj = readRDS(file = RDSsaved)

            #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

            DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                    na.value = "gray") +
              ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
              scale_colour_hue(drop = FALSE) +
              NoLegend()

            DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                    na.value = "gray") +
              ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
              scale_colour_hue(drop = FALSE) +
              NoLegend()

            cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
            cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

            #pdf(pdfname, width=18, height = 10)
            #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

            ##########################################
            # select subset of cells
            ##########################################

            # select cells with cluster index
            ##########################################
            #cluster.sels = c('29', '32', '35', '40', '42')
            #cluster.sels = c('6', '24', '20', '7')
            #cluster.sels = c('13', '17', '18', '11', '25', '10')
            # cluster.sels = c('1', '9', '12', '19')
            # #ids.sel = c('MSxaaa')
            # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
            #
            # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
            #
            # if(length(ids.sel)>0){
            #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
            #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
            #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
            # }else{
            #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
            # }

            # select cells with ids
            ##########################################
            index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
            table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

            ids.current = names(table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany'))
            ids.current = ids.current[grep('MSxapa', ids.current)]

            #ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
            ids.sels = ids.current
            # ids.sels = setdiff(ids.current[which(nchar(ids.current)>5)],
            #                    c('MSxppaa'))
            #
            ids.left = setdiff(ids.current, ids.sels)
            # print(ids.left)
            # nchar(ids.left)

            cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

            ##########################################
            # subset seurat object with selected cells
            ##########################################
            cat(length(cells.sels), ' cells selected to annotate \n')
            sub.obj = subset(seurat.obj, cells = cells.sels)

            sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
            #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
            sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
            sub.obj$pred.ids = sub.obj$predicted.ids.seurat
            xx = table(sub.obj$predicted.ids.seurat.keep)
            xx[xx>10]
            sub.obj$pred.ids.filtered = sub.obj$pred.ids
            sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

            DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()


            ##########################################
            # check potential ids for selected clusters
            ##########################################
            counts = table(sub.obj$manual.annot.ids, as.character(sub.obj$seurat_clusters))
            barplot(counts, main="cluster compositions by scmap ",
                    xlab=NULL, col=c(1:nrow(counts)), las = 2,
                    legend = rownames(counts))

            #counts[, match(c('31', '28', '52'), colnames(counts))]
            #counts.seurat[, match(c('31', '28', '52'), colnames(counts.seurat))]
            ##########################################
            # find new set of variable genes and redo pca and umap
            ##########################################
            Explore.umap.parameters.for.BWMcells = FALSE
            if(Explore.umap.parameters.for.BWMcells){

              source.my.script('scRNA_functions.R')
              require(tictoc)
              tic()
              test.umap.params.for.BWM.cells(sub.obj,
                                             pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                             group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                             nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                             n.neighbors.sampling = c(10, 30, 50),
                                             min.dist.sampling = c(0.01, 0.1)
              )
              toc()

            }

            nfeatures = 1000;
            sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
            #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
            sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
            sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
            ElbowPlot(sub.obj, ndims = 50)

            nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
            n.neighbors = 20;
            min.dist = 0.1; spread = 1
            sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                               spread = spread, n.neighbors = n.neighbors,
                               min.dist = min.dist, verbose = TRUE)
            DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
              NoLegend()

            DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
              NoLegend()

            p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
              NoLegend()

            #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
            p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
              NoLegend() + ggtitle('scamp prediction')

            p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
              NoLegend() + ggtitle('seurat prediction')

            p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
              NoLegend()

            p0 + p2
            p0 + p1
            p1 + p2

            p2 + p3

            ##########################################
            # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
            ##########################################
            FindClusters_subclusters = function(sub.obj, resolution = 1.0)
            {
              sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
              return(sub.obj$seurat_clusters)
            }
            sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:10, compute.SNN = TRUE)
            sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
            DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

            p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, repel = TRUE,  pt.size = 2) +
              ggtitle('manual.ids') + NoLegend()
            p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2,
                         label.size = 8, na.value = "gray", combine = TRUE)

            p2 + p1 #+ ggsave(paste0(resDir, '/splitcluster_manual_IDs_iteration_', Rnb.iteration, '.pdf'), width = 18, height = 10)


            p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
                         group.by = 'seurat_clusters_split')

            p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                         group.by = 'seurat_clusters_split') + NoLegend()

            VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                    group.by = 'manual.annot.ids') + NoLegend()

            (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_R', Rnb.iteration, '.pdf'),
                                     width = 18, height = 16)

            p1 + p4

            plot(p3)

            #seurat.obj$seurat_clusters_pharynx = NA
            #seurat.obj$seurat_clusters_pharynx[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$seurat_clusters_split

            #dev.off()

            ##########################################
            # check the counts of predicted ids for newly split clusters
            ##########################################
            sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
            sub.obj$predicted.ids.prob = sub.obj$predicted.scores
            sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
            sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

            Idents(sub.obj) = sub.obj$seurat_clusters_split
            counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
            counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
            #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
            counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)

            check.specific.cluster = TRUE
            if(check.specific.cluster){
              idents.sel = c('14', '11')
              counts = counts[, !is.na(match(colnames(counts), idents.sel))]
              counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
              counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
              counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
              counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
              counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
            }

            features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                              'ceh-36', # MSxaap
                              'tbx-2', 'ceh-27') # MSxaaa
            features.sels = c('lin-12', 'ceh-36', 'ceh-22','tbx-2', 'tbx-7', 'ceh-34', # MSxaapa
                              'fos-1', 'ceh-36', 'ceh-34', # MSaaapp
                              'lin-12', 'pax-1', # MSxaapap
                              'ngn-1', 'ces-1', # MSpaapaa
                              'sptf-1', 'ceh-22' #MSaaapaa

            )

            features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                              'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                              #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                              'ceh-13', # MSaaaap not MSpaaap
                              'hlh-3', 'fem-1'
                              )

            # features of lineage MSxapa
            features.sels = c('dmd-4', 'cnd-1', 'swt-3', # MSxapa
              'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
              'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
              'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
              'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
              'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
              'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
              'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

            )

            features.sels = c('Y51H7C.10', 'hlh-6', 'ces-1', 'dmd-4', 'ceh-6', 'asp-4',
                              'C39E9.8'
                              )
            FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

            #VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)
            # update the annotation with marker genes


            cluster.assingment = list(
              c('12', 'MSaaaappp.MSxapaapp'),
              c('9', 'MSaaaappp.MSxapaapp'),
              c('6', 'MSxapaaa'),
              c('13', 'MSxapaaa'),
              c('10', "MSpaaappp.MSxapappa"),
              c('3', 'MSpaaapp'),
              c('7', 'MSxapap'),
              c('8', 'MSxapap'),
              c('5', 'MSxapaapa'),
              c('15', 'MSxapaapa'),
              c('11', 'MSxapapa.MSxapapaa'),
              c('14', 'MSxapapp')
            )


            # check info in JM data for specific lineage
            # ee = process.import.Murray.scRNA()
            # murray.ids = unique(ee$lineage)
            #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
            markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
            #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
            #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
            source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

            ids.sel = c('MSxapapp'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


            ##########################################
            # update the manual annotation if good marker genes or mapped ids were found
            ##########################################
            # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
            # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
            # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
            # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

            # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
            # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
            #mm = match(colnames(sub.obj), colnames(seurat.obj))
            #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

            # sub.obj$ids.backup = sub.obj$manual.annot.ids

            for(n in 1:length(cluster.assingment)){
              cluster.index = cluster.assingment[[n]][1];
              id2assign =  cluster.assingment[[n]][2];
              cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
              cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
              sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
              seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
              # if(is.na(id2assign)) {
              #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
              # }else{
              #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
              # }
            }

            #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

            DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                    pt.size = 2) + NoLegend()

            DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                    na.value = "gray") +
              ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
              scale_colour_hue(drop = FALSE) + NoLegend()


            saveRDS(seurat.obj, file = RDS2save)




            ##########################################
              # Main aim:
              # iteration R3,
              # just started branch of MSxaap
              #
              # Notes:
              # MSaaapp and MSaaappp were annotated with ceh-36 and absence of lin-12, but without id-specific markers
              #
              ##########################################
              Rnb.iteration = 3
              Refine.annotated.ids = TRUE;

              resDir = paste0("results/", version.analysis, '/annoted_pharynx')
              if(!dir.exists(resDir)){dir.create(resDir)}

              if(Refine.annotated.ids){by.group = 'manual.annot.ids';
              }else{by.group = 'seurat_clusters'}

              RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                             '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration -1, '.rds')

              RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                 '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration, '.rds')

              seurat.obj = readRDS(file = RDSsaved)

              #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

              DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                      na.value = "gray") +
                ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
                scale_colour_hue(drop = FALSE) +
                NoLegend()

              DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                      na.value = "gray") +
                ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
                scale_colour_hue(drop = FALSE) +
                NoLegend()

              cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
              cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

              #pdf(pdfname, width=18, height = 10)
              #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

              ##########################################
              # select subset of cells
              ##########################################

              # select cells with cluster index
              ##########################################
              #cluster.sels = c('29', '32', '35', '40', '42')
              #cluster.sels = c('6', '24', '20', '7')
              #cluster.sels = c('13', '17', '18', '11', '25', '10')
              # cluster.sels = c('1', '9', '12', '19')
              # #ids.sel = c('MSxaaa')
              # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
              #
              # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
              #
              # if(length(ids.sel)>0){
              #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
              #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
              #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
              # }else{
              #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
              # }

              # select cells with ids
              ##########################################
              index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
              table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

              ids.current = names(table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany'))

              ids.current = ids.current[grep('MSxaap|MSaaapp', ids.current)]
              #ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
              ids.sels = ids.current
              #ids.sels = setdiff(ids.current, c('MSxaap'))
              #
              ids.left = setdiff(ids.current, ids.sels)
              # print(ids.left)
              # nchar(ids.left)

              cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

              ##########################################
              # subset seurat object with selected cells
              ##########################################
              cat(length(cells.sels), ' cells selected to annotate \n')
              sub.obj = subset(seurat.obj, cells = cells.sels)

              sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
              #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
              sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
              sub.obj$pred.ids = sub.obj$predicted.ids.seurat
              xx = table(sub.obj$predicted.ids.seurat.keep)
              xx[xx>10]
              sub.obj$pred.ids.filtered = sub.obj$pred.ids
              sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

              DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

              ##########################################
              # find new set of variable genes and redo pca and umap
              ##########################################
              Explore.umap.parameters.for.BWMcells = FALSE
              if(Explore.umap.parameters.for.BWMcells){

                source.my.script('scRNA_functions.R')
                require(tictoc)
                tic()
                test.umap.params.for.BWM.cells(sub.obj,
                                               pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                               group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                               nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                               n.neighbors.sampling = c(10, 30, 50),
                                               min.dist.sampling = c(0.01, 0.1)
                )
                toc()

              }

              nfeatures = 1000;
              sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
              #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
              sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
              sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
              ElbowPlot(sub.obj, ndims = 50)

              nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
              n.neighbors = 5;
              min.dist = 0.01; spread = 1
              sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                                 spread = spread, n.neighbors = n.neighbors,
                                 min.dist = min.dist, verbose = TRUE)
              if(!Refine.annotated.ids){
                DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
                  NoLegend()
              }else{
                DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                  NoLegend()
              }

              p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                NoLegend()

              #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
              p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                NoLegend() + ggtitle('scamp prediction')

              p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
                NoLegend() + ggtitle('seurat prediction')

              p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                NoLegend()

              p0 + p2
              p0 + p1
              p1 + p2

              p2 + p3

              ##########################################
              # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
              ##########################################
              FindClusters_subclusters = function(sub.obj, resolution = 1.0)
              {
                sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
                return(sub.obj$seurat_clusters)
              }
              sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
              sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1)
              DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

              p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
                ggtitle('manual.ids') + NoLegend()
              p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
                           label.size = 8, na.value = "gray", combine = TRUE)

              #p2 + p1 #+ ggsave(paste0(resDir, '/splitcluster_manual_IDs_iteration_', Rnb.iteration, '.pdf'), width = 18, height = 10)

              p3 = VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
                           group.by = 'seurat_clusters_split')
              p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1,
                           group.by = 'seurat_clusters_split') + NoLegend()

              #VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'manual.annot.ids') + NoLegend()

              (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_R', Rnb.iteration, '.pdf'),
                                       width = 18, height = 16)

              #p1 + p4

              #plot(p3)

              #seurat.obj$seurat_clusters_pharynx = NA
              #seurat.obj$seurat_clusters_pharynx[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$seurat_clusters_split

              #dev.off()

              ##########################################
              # check the counts of predicted ids for newly split clusters
              ##########################################

              idents.sel = c('3', '4', '7')
              Idents(sub.obj) = sub.obj$seurat_clusters_split
              sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
              sub.obj$predicted.ids.prob = sub.obj$predicted.scores
              sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
              sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

              counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
              counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
              #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
              counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
              if(exists('idents.sel')){
                if(length(setdiff(colnames(counts), idents.sel)) > 0){
                  counts = counts[, !is.na(match(colnames(counts), idents.sel))]
                  counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
                  counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
                  counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
                  counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
                  counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
                }
              }


              features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                                'ceh-36', # MSxaap
                                'tbx-2', 'ceh-27') # MSxaaa
              # features of lineage MSxapa
              features.sels = c('dmd-4', 'cnd-1', 'swt-3', # MSxapa
                                'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                                'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                                'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                                'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                                'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                                'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                                'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

              )

              features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                                'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                                #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                                'ceh-13', # MSaaaap not MSpaaap
                                'hlh-3', 'fem-1'
                                )
              # features of MSxaap lineage
              features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                                'sptf-1', 'ceh-22',  #MSaaapaa
                                'lin-12', 'pax-1', # MSxaapap
                                'ngn-1', 'ces-1', # MSpaapaa
                                'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                                'jun-1', 'ttx-1', 'irx-1' # MSaaappp

              )


              features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258'
                                )
              FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)
              #VlnPlot(sub.obj, features = features.sels,  group.by = 'seurat_clusters_split', idents = idents.sel)

              # update the annotation with marker genes
              cluster.assingment = list(
                c('5', 'MSxaap'),
                c('3', 'MSaaapp'),
                c('4', 'MSaaapp'),
                c('7', 'MSaaappp')

              )


              # check info in JM data for specific lineage
              # ee = process.import.Murray.scRNA()
              # murray.ids = unique(ee$lineage)
              #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
              markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
              #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
              #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
              source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

              ids.sel = c('MSaaapp'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


              ##########################################
              # update the manual annotation if good marker genes or mapped ids were found
              ##########################################
              # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
              # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
              # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
              # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

              # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
              # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
              #mm = match(colnames(sub.obj), colnames(seurat.obj))
              #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

              # sub.obj$ids.backup = sub.obj$manual.annot.ids

              for(n in 1:length(cluster.assingment)){
                cluster.index = cluster.assingment[[n]][1];
                id2assign =  cluster.assingment[[n]][2];
                cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
                cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
                sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
                seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
                # if(is.na(id2assign)) {
                #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
                # }else{
                #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
                # }
              }

              #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

              DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                      pt.size = 2) + NoLegend()

              DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                      na.value = "gray") +
                ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
                scale_colour_hue(drop = FALSE) + NoLegend()


              saveRDS(seurat.obj, file = RDS2save)




              ##########################################
                # Main aim:
                # revise the lineage MSxaap
                #
                # Notes:
                # it is still tricky to dissect MSxaapap, MSaaapapp and MSxaapapa
                # and also the dissection between MSpaapaa and MSaaapaa
                ##########################################
                Rnb.iteration = 4
                Refine.annotated.ids = TRUE;

                resDir = paste0("results/", version.analysis, '/annoted_pharynx')
                if(!dir.exists(resDir)){dir.create(resDir)}

                if(Refine.annotated.ids){by.group = 'manual.annot.ids';
                }else{by.group = 'seurat_clusters'}

                RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                               '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration -1, '.rds')

                RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                   '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration, '.rds')

                seurat.obj = readRDS(file = RDSsaved)

                #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

                DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                        na.value = "gray") +
                  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
                  scale_colour_hue(drop = FALSE) +
                  NoLegend()

                DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                        na.value = "gray") +
                  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
                  scale_colour_hue(drop = FALSE) +
                  NoLegend()

                cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
                cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

                #pdf(pdfname, width=18, height = 10)
                #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

                ##########################################
                # select subset of cells
                ##########################################

                # select cells with cluster index
                ##########################################
                #cluster.sels = c('29', '32', '35', '40', '42')
                #cluster.sels = c('6', '24', '20', '7')
                #cluster.sels = c('13', '17', '18', '11', '25', '10')
                # cluster.sels = c('1', '9', '12', '19')
                # #ids.sel = c('MSxaaa')
                # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
                #
                # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
                #
                # if(length(ids.sel)>0){
                #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
                #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
                #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
                # }else{
                #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
                # }

                # select cells with ids
                ##########################################
                index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
                table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

                ids.current = names(table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany'))

                ids.current = ids.current[grep('MSxaap|MSaaapp', ids.current)]
                #ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
                #ids.sels = ids.current
                ids.sels = setdiff(ids.current, c('MSxaap', 'MSaaapp', 'MSaaappp'))

                ids.left = setdiff(ids.current, ids.sels)
                # print(ids.left)
                # nchar(ids.left)

                cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

                ##########################################
                # subset seurat object with selected cells
                ##########################################
                cat(length(cells.sels), ' cells selected to annotate \n')
                sub.obj = subset(seurat.obj, cells = cells.sels)

                sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
                #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
                sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
                sub.obj$pred.ids = sub.obj$predicted.ids.seurat
                xx = table(sub.obj$predicted.ids.seurat.keep)
                xx[xx>10]
                sub.obj$pred.ids.filtered = sub.obj$pred.ids
                sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

                DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

                ##########################################
                # find new set of variable genes and redo pca and umap
                ##########################################
                Explore.umap.parameters.for.BWMcells = FALSE
                if(Explore.umap.parameters.for.BWMcells){

                  source.my.script('scRNA_functions.R')
                  require(tictoc)
                  tic()
                  test.umap.params.for.BWM.cells(sub.obj,
                                                 pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                                 group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                                 nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                                 n.neighbors.sampling = c(10, 30, 50),
                                                 min.dist.sampling = c(0.01, 0.1)
                  )
                  toc()

                }

                nfeatures = 1000;
                sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
                #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
                sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
                sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
                ElbowPlot(sub.obj, ndims = 50)

                nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
                n.neighbors = 10;
                min.dist = 0.01; spread = 1
                sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                                   spread = spread, n.neighbors = n.neighbors,
                                   min.dist = min.dist, verbose = TRUE)
                if(!Refine.annotated.ids){
                  DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
                    NoLegend()
                }else{
                  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                    NoLegend()
                }

                p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                  NoLegend()

                #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
                p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                  NoLegend() + ggtitle('scamp prediction')

                p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
                  NoLegend() + ggtitle('seurat prediction')

                p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                  NoLegend()

                p0 + p2
                p0 + p1
                p1 + p2

                p2 + p3

                ##########################################
                # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
                ##########################################
                FindClusters_subclusters = function(sub.obj, resolution = 1.0)
                {
                  sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
                  return(sub.obj$seurat_clusters)
                }
                sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:10, compute.SNN = TRUE)
                sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1)
                DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

                p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
                  ggtitle('manual.ids') + NoLegend()
                p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
                             label.size = 8, na.value = "gray", combine = TRUE)
                p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

                (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_R', Rnb.iteration, '.pdf'),
                                         width = 18, height = 16)

                #dev.off()
                ##########################################
                # check the counts of predicted ids for newly split clusters
                ##########################################

                #idents.sel = c('3', '4', '7')
                Idents(sub.obj) = sub.obj$seurat_clusters_split
                sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
                sub.obj$predicted.ids.prob = sub.obj$predicted.scores
                sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
                sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

                counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
                counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
                #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
                counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
                if(exists('idents.sel')){
                  if(length(setdiff(colnames(counts), idents.sel)) > 0){
                    counts = counts[, !is.na(match(colnames(counts), idents.sel))]
                    counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
                    counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
                    counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
                    counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
                    counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
                  }
                }


                features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                                  'ceh-36', # MSxaap
                                  'tbx-2', 'ceh-27') # MSxaaa
                # features of lineage MSxapa
                features.sels = c('dmd-4', 'cnd-1', 'swt-3', # MSxapa
                                  'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                                  'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                                  'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                                  'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                                  'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                                  'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                                  'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

                )

                features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                                  'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                                  #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                                  'ceh-13', # MSaaaap not MSpaaap
                                  'hlh-3', 'fem-1'
                                  )
                # features of MSxaap lineage
                features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                                  'sptf-1', 'ceh-22',  #MSaaapaa
                                  'lin-12', 'pax-1', # MSxaapap
                                  'ngn-1', 'ces-1', # MSpaapaa
                                  'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                                  'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                                  'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                                  'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
                )


                features.sels = c('lin-12', 'ceh-22', 'ngn-1', 'ces-1', 'ceh-34', 'sptf-1'
                                  )
                FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

                # update the annotation with marker genes
                cluster.assingment = list(
                  c('6', 'MSxaapap'),
                  c('5', 'MSxaapapa/ABalpappapp'),
                  c('1', 'MSxaapapa/ABalpappapp'),
                  c('4', 'MSaaapapp'),
                  c('2', 'MSaaapaa'),
                  c('0', 'MSaaapaa'),
                  c('3', 'MSpaapaa')

                )


                # check info in JM data for specific lineage
                # ee = process.import.Murray.scRNA()
                # murray.ids = unique(ee$lineage)
                #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
                markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
                #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
                #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
                source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

                ids.sel = c('MSaaapapp'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


                ##########################################
                # update the manual annotation if good marker genes or mapped ids were found
                ##########################################
                # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
                # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
                # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
                # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

                # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
                # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
                #mm = match(colnames(sub.obj), colnames(seurat.obj))
                #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

                # sub.obj$ids.backup = sub.obj$manual.annot.ids

                for(n in 1:length(cluster.assingment)){
                  cluster.index = cluster.assingment[[n]][1];
                  id2assign =  cluster.assingment[[n]][2];
                  cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
                  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
                  sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
                  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
                  # if(is.na(id2assign)) {
                  #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
                  # }else{
                  #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
                  # }
                }

                #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

                DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                        pt.size = 2) + NoLegend()

                DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                        na.value = "gray") +
                  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
                  scale_colour_hue(drop = FALSE) + NoLegend()


                saveRDS(seurat.obj, file = RDS2save)



                ##########################################
                  # Main aim:
                  # revise the MSxaaa lineage
                  #
                  # Notes:
                  # in this iteration we deal with MSxaaa and some well-split clusters
                  #
                  ##########################################
                  Rnb.iteration = 5
                  Refine.annotated.ids = TRUE;

                  resDir = paste0("results/", version.analysis, '/annoted_pharynx')
                  if(!dir.exists(resDir)){dir.create(resDir)}

                  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
                  }else{by.group = 'seurat_clusters'}

                  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                                 '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration -1, '.rds')

                  RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                     '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration, '.rds')

                  seurat.obj = readRDS(file = RDSsaved)

                  #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

                  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                          na.value = "gray") +
                    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
                    scale_colour_hue(drop = FALSE) +
                    NoLegend()

                  DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                          na.value = "gray") +
                    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
                    scale_colour_hue(drop = FALSE) +
                    NoLegend()

                  cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
                  cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

                  #pdf(pdfname, width=18, height = 10)
                  #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

                  ##########################################
                  # select subset of cells
                  ##########################################

                  # select cells with cluster index
                  ##########################################
                  #cluster.sels = c('29', '32', '35', '40', '42')
                  #cluster.sels = c('6', '24', '20', '7')
                  #cluster.sels = c('13', '17', '18', '11', '25', '10')
                  # cluster.sels = c('1', '9', '12', '19')
                  # #ids.sel = c('MSxaaa')
                  # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
                  #
                  # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
                  #
                  # if(length(ids.sel)>0){
                  #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
                  #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
                  #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
                  # }else{
                  #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
                  # }

                  # select cells with ids
                  ##########################################
                  index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
                  table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

                  ids.current = names(table(seurat.obj$manual.annot.ids[index.pharynx.cells]))
                  ids.current = ids.current[grep('MSxaap|MSaaapp|MSxapa', ids.current, invert = TRUE)]
                  #ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
                  #ids.sels = ids.current
                  ids.sels = setdiff(ids.current,
                                     c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp'))

                  ids.left = setdiff(ids.current, ids.sels)
                  # print(ids.left)
                  # nchar(ids.left)

                  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

                  ##########################################
                  # subset seurat object with selected cells
                  ##########################################
                  cat(length(cells.sels), ' cells selected to annotate \n')
                  sub.obj = subset(seurat.obj, cells = cells.sels)

                  sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
                  #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
                  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
                  sub.obj$pred.ids = sub.obj$predicted.ids.seurat
                  xx = table(sub.obj$predicted.ids.seurat.keep)
                  xx[xx>10]
                  sub.obj$pred.ids.filtered = sub.obj$pred.ids
                  sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

                  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

                  ##########################################
                  # find new set of variable genes and redo pca and umap
                  ##########################################
                  Explore.umap.parameters.for.BWMcells = FALSE
                  if(Explore.umap.parameters.for.BWMcells){

                    source.my.script('scRNA_functions.R')
                    require(tictoc)
                    tic()
                    test.umap.params.for.BWM.cells(sub.obj,
                                                   pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                                   group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                                   nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                                   n.neighbors.sampling = c(10, 30, 50),
                                                   min.dist.sampling = c(0.01, 0.1)
                    )
                    toc()

                  }

                  nfeatures = 1000;
                  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
                  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
                  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
                  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
                  ElbowPlot(sub.obj, ndims = 50)

                  nb.pcs = 20 # nb of pcs depends on the considered clusters or ids
                  n.neighbors = 10;
                  min.dist = 0.01; spread = 1
                  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                                     spread = spread, n.neighbors = n.neighbors,
                                     min.dist = min.dist, verbose = TRUE)
                  if(!Refine.annotated.ids){
                    DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
                      NoLegend()
                  }else{
                    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                      NoLegend()
                  }

                  p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                    NoLegend()

                  #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
                  p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                    NoLegend() + ggtitle('scamp prediction')

                  p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
                    NoLegend() + ggtitle('seurat prediction')

                  p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                    NoLegend()

                  p0 + p2
                  p0 + p1
                  p1 + p2

                  p2 + p3

                  ##########################################
                  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
                  ##########################################
                  FindClusters_subclusters = function(sub.obj, resolution = 1.0)
                  {
                    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
                    return(sub.obj$seurat_clusters)
                  }
                  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
                  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
                  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

                  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
                    ggtitle('manual.ids') + NoLegend()
                  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
                               label.size = 8, na.value = "gray", combine = TRUE)
                  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

                  (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_R', Rnb.iteration, '.pdf'),
                                           width = 18, height = 16)

                  #dev.off()
                  ##########################################
                  # check the counts of predicted ids for newly split clusters
                  ##########################################

                  idents.sel = c('6', '7', '8', '11')
                  Idents(sub.obj) = sub.obj$seurat_clusters_split
                  sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
                  sub.obj$predicted.ids.prob = sub.obj$predicted.scores
                  sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
                  sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

                  counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
                  counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
                  #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
                  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
                  if(exists('idents.sel')){
                    if(length(setdiff(colnames(counts), idents.sel)) > 0){
                      counts = counts[, !is.na(match(colnames(counts), idents.sel))]
                      counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
                      counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
                      counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
                      counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
                      counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
                    }
                  }

                  # features of early stage
                  features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                                    'ceh-36', # MSxaap
                                    'tbx-2', 'ceh-27') # MSxaaa
                  # features of lineage MSxapa
                  features.sels = c('dmd-4', 'cnd-1', 'swt-3', # MSxapa
                                    'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                                    'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                                    'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                                    'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                                    'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                                    'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                                    'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

                  )

                  # features of MSxaap lineage
                  features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                                    'sptf-1', 'ceh-22',  #MSaaapaa
                                    'lin-12', 'pax-1', # MSxaapap
                                    'ngn-1', 'ces-1', # MSpaapaa
                                    'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                                    'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                                    'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                                    'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
                  )
                  # features of MSxaaa
                  features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                                    'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                                    #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                                    'ceh-13', # MSaaaap not MSpaaap
                                    'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                                    'lim-7', 'dod-6', 'ces-1' # MSpaaapa
                                    )

                  features.sels = c('cwn-2', 'ttx-1'

                                    )
                  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

                  # update the annotation with marker genes
                  cluster.assingment = list(
                    c('5', 'MSxaaa'),
                    c('11', 'MSpaaapa'),
                    c('6', 'MSpaaapp'),
                    c('8', 'MSpaaaap.sure'),
                    c('7', 'MSaappa')
                  )


                  # check info in JM data for specific lineage
                  # ee = process.import.Murray.scRNA()
                  # murray.ids = unique(ee$lineage)
                  #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
                  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
                  #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
                  #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
                  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

                  ids.sel = c('MSpaaapp'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


                  ##########################################
                  # update the manual annotation if good marker genes or mapped ids were found
                  ##########################################
                  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
                  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
                  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
                  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

                  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
                  # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
                  #mm = match(colnames(sub.obj), colnames(seurat.obj))
                  #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

                  # sub.obj$ids.backup = sub.obj$manual.annot.ids

                  for(n in 1:length(cluster.assingment)){
                    cluster.index = cluster.assingment[[n]][1];
                    id2assign =  cluster.assingment[[n]][2];
                    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
                    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
                    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
                    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
                    # if(is.na(id2assign)) {
                    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
                    # }else{
                    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
                    # }
                  }

                  #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

                  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                          pt.size = 2) + NoLegend()

                  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                          na.value = "gray") +
                    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
                    scale_colour_hue(drop = FALSE) + NoLegend()


                  saveRDS(seurat.obj, file = RDS2save)




                  ##########################################
                  # Main aim:
                  # revise close clusters in the MSxaaa lineage
                  #
                  # Notes:
                  # here we annotated again the four daughters of MSxaaa
                  #
                  ##########################################
                  Rnb.iteration = 6
                  Refine.annotated.ids = TRUE;

                  resDir = paste0("results/", version.analysis, '/annoted_pharynx')
                  if(!dir.exists(resDir)){dir.create(resDir)}

                  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
                  }else{by.group = 'seurat_clusters'}

                  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                                 '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration -1, '.rds')

                  RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                     '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration, '.rds')

                  seurat.obj = readRDS(file = RDSsaved)

                  #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

                  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                          na.value = "gray") +
                    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
                    scale_colour_hue(drop = FALSE) +
                    NoLegend()

                  DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                          na.value = "gray") +
                    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
                    scale_colour_hue(drop = FALSE) +
                    NoLegend()

                  cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
                  cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

                  #pdf(pdfname, width=18, height = 10)
                  #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

                  ##########################################
                  # select subset of cells
                  ##########################################

                  # select cells with cluster index
                  ##########################################
                  #cluster.sels = c('29', '32', '35', '40', '42')
                  #cluster.sels = c('6', '24', '20', '7')
                  #cluster.sels = c('13', '17', '18', '11', '25', '10')
                  # cluster.sels = c('1', '9', '12', '19')
                  # #ids.sel = c('MSxaaa')
                  # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
                  #
                  # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
                  #
                  # if(length(ids.sel)>0){
                  #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
                  #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
                  #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
                  # }else{
                  #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
                  # }

                  # select cells with ids
                  ##########################################
                  index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
                  table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

                  ids.current = names(table(seurat.obj$manual.annot.ids[index.pharynx.cells]))
                  ids.current = ids.current[grep('MSxaap|MSaaapp|MSxapa', ids.current, invert = TRUE)]
                  #ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
                  #ids.sels = ids.current
                  ids.sels = setdiff(ids.current,
                                     c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp',
                                       'MSxaaa', 'MSpaaapp', 'MSpaaapa', 'MSaappa', 'MSpaaaap.sure'))

                  ids.sels = c("MSaaaaa", "MSaaaap", "MSpaaaa", "MSpaaap", "MSxaaa.like.or.others")
                  ids.left = setdiff(ids.current, ids.sels)
                  # print(ids.left)
                  # nchar(ids.left)

                  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

                  ##########################################
                  # subset seurat object with selected cells
                  ##########################################
                  cat(length(cells.sels), ' cells selected to annotate \n')
                  sub.obj = subset(seurat.obj, cells = cells.sels)

                  sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
                  #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
                  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
                  sub.obj$pred.ids = sub.obj$predicted.ids.seurat
                  xx = table(sub.obj$predicted.ids.seurat.keep)
                  xx[xx>10]
                  sub.obj$pred.ids.filtered = sub.obj$pred.ids
                  sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

                  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

                  ##########################################
                  # find new set of variable genes and redo pca and umap
                  ##########################################
                  Explore.umap.parameters.for.BWMcells = FALSE
                  if(Explore.umap.parameters.for.BWMcells){

                    source.my.script('scRNA_functions.R')
                    require(tictoc)
                    tic()
                    test.umap.params.for.BWM.cells(sub.obj,
                                                   pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                                   group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                                   nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                                   n.neighbors.sampling = c(10, 30, 50),
                                                   min.dist.sampling = c(0.01, 0.1)
                    )
                    toc()

                  }

                  nfeatures = 2000;
                  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
                  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
                  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
                  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
                  ElbowPlot(sub.obj, ndims = 50)

                  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
                  n.neighbors = 10;
                  min.dist = 0.01; spread = 1
                  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                                     spread = spread, n.neighbors = n.neighbors,
                                     min.dist = min.dist, verbose = TRUE)
                  if(!Refine.annotated.ids){
                    DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
                      NoLegend()
                  }else{
                    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
                      NoLegend()
                  }

                  p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                    NoLegend()

                  #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
                  p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                    NoLegend() + ggtitle('scamp prediction')

                  p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
                    NoLegend() + ggtitle('seurat prediction')

                  p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
                    NoLegend()

                  p0 + p2
                  p0 + p1
                  p1 + p2

                  p2 + p3

                  ##########################################
                  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
                  ##########################################
                  FindClusters_subclusters = function(sub.obj, resolution = 1.0)
                  {
                    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
                    return(sub.obj$seurat_clusters)
                  }
                  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
                  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.)
                  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

                  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
                    ggtitle('manual.ids') + NoLegend()
                  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
                               label.size = 8, na.value = "gray", combine = TRUE)
                  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

                  (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_R', Rnb.iteration, '.pdf'),
                                           width = 18, height = 16)

                  #dev.off()
                  ##########################################
                  # check the counts of predicted ids for newly split clusters
                  ##########################################

                  idents.sel = c('2', '4', '3')
                  Idents(sub.obj) = sub.obj$seurat_clusters_split
                  sub.obj$predicted.ids = sub.obj$predicted.ids.scmap
                  sub.obj$predicted.ids.prob = sub.obj$predicted.scores
                  sub.obj$predicted.ids.fitered = sub.obj$predicted.ids.scmap
                  sub.obj$predicted.ids.fitered[sub.obj$predicted.ids.prob < 0.7] = NA

                  counts = table(sub.obj$predicted.ids, sub.obj$seurat_clusters_split)
                  counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
                  #counts.seurat.filter = table(sub.obj$predicted.ids.fitered, sub.obj$seurat_clusters_split)
                  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
                  if(exists('idents.sel')){
                    if(length(setdiff(colnames(counts), idents.sel)) > 0){
                      counts = counts[, !is.na(match(colnames(counts), idents.sel))]
                      counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
                      counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
                      counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
                      counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
                      counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
                    }
                  }

                  # features of early stage
                  features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                                    'ceh-36', # MSxaap
                                    'tbx-2', 'ceh-27') # MSxaaa
                  # features of lineage MSxapa
                  features.sels = c('dmd-4', 'cnd-1', 'swt-3', # MSxapa
                                    'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                                    'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                                    'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                                    'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                                    'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                                    'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                                    'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

                  )

                  # features of MSxaap lineage
                  features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                                    'sptf-1', 'ceh-22',  #MSaaapaa
                                    'lin-12', 'pax-1', # MSxaapap
                                    'ngn-1', 'ces-1', # MSpaapaa
                                    'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                                    'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                                    'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                                    'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
                  )
                  # features of MSxaaa
                  features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                                    'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                                    #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                                    'ceh-13', # MSaaaap not MSpaaap
                                    'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                                    'lim-7', 'dod-6', 'ces-1' # MSpaaapa
                                    )

                  features.sels = c('fem-1', 'ceh-32', 'ceh-27', 'let-381',
                                    'fkh-2', 'ceh-53', 'hlh-3', 'spi-1', 'F54E2.2',
                                    'igcm-4', 'K04G2.12', 'F26B1.1',
                                    'ceh-34', 'C45G7.4', 'unc-30'

                                    )
                  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

                  # update the annotation with marker genes
                  cluster.assingment = list(
                    c('0', 'MSpaaaa'),
                    c('3', 'MSaaaaa.to.confirm'),
                    c('5', 'MSaaaap'),
                    c('2', 'MSpaaap'),
                    c('4', 'MSpaaap'),
                    c('1', 'MSaaaaaa.to.confirm')

                  )


                  # check info in JM data for specific lineage
                  # ee = process.import.Murray.scRNA()
                  # murray.ids = unique(ee$lineage)
                  #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
                  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
                  #markers = markers[!is.na(match(markers$Lineage, bwms)), ]
                  #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
                  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

                  ids.sel = c('MSaaaaaa'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


                  ##########################################
                  # update the manual annotation if good marker genes or mapped ids were found
                  ##########################################
                  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
                  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
                  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
                  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

                  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == '')] = 'mixture_terminal_mothers'
                  # seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids
                  #mm = match(colnames(sub.obj), colnames(seurat.obj))
                  #seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

                  # sub.obj$ids.backup = sub.obj$manual.annot.ids

                  for(n in 1:length(cluster.assingment)){
                    cluster.index = cluster.assingment[[n]][1];
                    id2assign =  cluster.assingment[[n]][2];
                    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
                    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
                    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
                    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
                    # if(is.na(id2assign)) {
                    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
                    # }else{
                    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
                    # }
                  }

                  #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

                  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
                          pt.size = 2) + NoLegend()

                  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
                          na.value = "gray") +
                    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
                    scale_colour_hue(drop = FALSE) + NoLegend()

                  saveRDS(seurat.obj, file = RDS2save)

                  ##########################################
  # Main aim:
  # revise terminal cells of MSxaaa
  #
  # Notes:
  # annotated "MSxpaaaa/MSxpaaaax" changed to "MSxpaaaa" ring gang cells, because MSxpaaaax is together with MSaaaaaax,
  # i.e. "MSaaaaaax/MSxpaaaax"
  # all terminal cells of MSxaaa were kind of annotated, even though many of them were annotated based on a small set of markers,
  # while other markers were controversial and confusing.
  # further comfirmed will be good by more reliable markers
  # after this iteration, the first annotation revision is finished.
  ##########################################
  Rnb.iteration = 7
  Refine.annotated.ids = TRUE;

  resDir = paste0("results/", version.analysis, '/annoted_pharynx')
  if(!dir.exists(resDir)){dir.create(resDir)}

  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}

  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                 '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration -1, '.rds')

  RDS2save =  paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                     '_cleanedBWM_and_Pharynx_iteration_R', Rnb.iteration, '.rds')

  seurat.obj = readRDS(file = RDSsaved)

  #jj = which(seurat.obj$manual.annot.ids == "MSxpaaaa/MSxpaaaax")
  #seurat.obj$manual.annot.ids[jj] = 'MSxpaaaa'

  #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  DimPlot(seurat.obj, group.by = "seurat_clusters_pharynx", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
  cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

  #pdf(pdfname, width=18, height = 10)
  #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

  ##########################################
  # select subset of cells
  ##########################################

  # select cells with cluster index
  ##########################################
  #cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('6', '24', '20', '7')
  #cluster.sels = c('13', '17', '18', '11', '25', '10')
  # cluster.sels = c('1', '9', '12', '19')
  # #ids.sel = c('MSxaaa')
  # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
  #
  # cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
  #
  # if(length(ids.sel)>0){
  #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
  #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
  #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
  # }else{
  #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  # }

  # select cells with ids
  ##########################################
  index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
  table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

  ids.current = names(table(seurat.obj$manual.annot.ids[index.pharynx.cells]))
  ids.current = ids.current[grep('MSxaap|MSaaapp|MSxapa', ids.current, invert = TRUE)]
  #ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
  #ids.sels = ids.current
  ids.sels = setdiff(ids.current,
                     c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp',
                       'MSxaaa', 'MSpaaapp', 'MSpaaapa', 'MSaappa', 'MSpaaaap.sure',
                        "MSaaaaa",  "MSaaaap", "MSpaaaa", "MSpaaap"))

  #ids.sels = c("MSaaaaa", "MSaaaap", "MSpaaaa", "MSpaaap", "MSxaaa.like.or.others")
  ids.left = setdiff(ids.current, ids.sels)
  # print(ids.left)
  # nchar(ids.left)

  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

  ##########################################
  # subset seurat object with selected cells
  ##########################################
  cat(length(cells.sels), ' cells selected to annotate \n')
  sub.obj = subset(seurat.obj, cells = cells.sels)

  sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
  #sub.obj$predicted.ids.fitered[is.na(sub.obj$predicted.ids.fitered)] = 'unassigned'
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))
  sub.obj$pred.ids = sub.obj$predicted.ids.seurat
  xx = table(sub.obj$predicted.ids.seurat.keep)
  xx[xx>10]
  sub.obj$pred.ids.filtered = sub.obj$pred.ids
  sub.obj$pred.ids.filtered[is.na(match(sub.obj$pred.ids, names(xx[xx>10])))] = NA

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  Explore.umap.parameters.for.BWMcells = FALSE
  if(Explore.umap.parameters.for.BWMcells){

    source.my.script('scRNA_functions.R')
    require(tictoc)
    tic()
    test.umap.params.for.BWM.cells(sub.obj,
                                   pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                   group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                   nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                   n.neighbors.sampling = c(10, 30, 50),
                                   min.dist.sampling = c(0.01, 0.1)
    )
    toc()

  }

  nfeatures = 1000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 10;
  min.dist = 0.01; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist, verbose = TRUE)
  if(!Refine.annotated.ids){
    DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()
  }else{
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()
  }

  p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend()

  #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
  p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend() + ggtitle('scamp prediction')

  p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
    NoLegend() + ggtitle('seurat prediction')

  p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend()

  p0 + p2
  p0 + p1
  p1 + p2

  p2 + p3

  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 1.0)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:10, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
    ggtitle('manual.ids') + NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
               label.size = 8, na.value = "gray", combine = TRUE)
  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

  (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_R', Rnb.iteration, '.pdf'),
                           width = 18, height = 16)

  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'seurat_clusters_split')
  #dev.off()
  ##########################################
  # check the counts of predicted ids for newly split clusters
  ##########################################

  idents.sel = c('0', '2')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
  counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
  if(exists('idents.sel')){
    if(length(setdiff(colnames(counts), idents.sel)) > 0){
      counts = counts[, !is.na(match(colnames(counts), idents.sel))]
      counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
      counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
      counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
      counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
      counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
    }
  }

  # features of early stage
  features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                    'ceh-36', # MSxaap
                    'tbx-2', 'ceh-27') # MSxaaa
  # features of lineage MSxapa
  features.sels = c('dmd-4', 'cnd-1', 'swt-3', # MSxapa
                    'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                    'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                    'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                    'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                    'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                    'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                    'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

  )
  # features of MSxaap lineage
  features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                    'sptf-1', 'ceh-22',  #MSaaapaa
                    'lin-12', 'pax-1', # MSxaapap
                    'ngn-1', 'ces-1', # MSpaapaa
                    'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                    'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                    'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                    'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
  )
  # features of MSxaaa
  features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                    'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                    #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                    'ceh-13', # MSaaaap not MSpaaap
                    'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                    'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                    'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                    'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                    'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                    'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                    'igcm-4', 'F54E2.2', 'K04G2.12',
                    'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                    'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6' # MSaaaaapa
                    )

  features.sels = c('hlh-3', 'hlh-6'
                    )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  # update the annotation with marker genes
  cluster.assingment = list(
    c('4', 'MSaaaaaa'),
    c('6', 'MSaaaaaaa'),
    c('5', 'MSaaaaap'),
    c('3', 'MSaaaapa'),
    c('7', 'MSaaaapa'),
    c('1', 'MSaaaapp'),
    c('0', 'MSpaaaaa'),
    c('2', 'MSaaaaapa')
  )

  cat(length(cluster.assingment), 'clusters assigned \n')
  cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
  nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
  if(length(nb.unassigned) > 1){
    cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
  }

  # check info in JM data for specific lineage
  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
  #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

  ids.sel = c('MSaaaaapa'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    # if(is.na(id2assign)) {
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
    # }else{
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
    # }
  }

  #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
          pt.size = 2) + NoLegend()

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE) + NoLegend()

  saveRDS(seurat.obj, file = RDS2save)












################################################################################
## start the global revision (GR) iteration,
## try to consolidate the annotations and to find some missing terminal cells for pharynx
################################################################################
##########################################
 # Main aim:
 # revise the pharynx cells but with annotation missing
 #
 # Notes:
 # here we just resplit the cell into tiny clusters first, and then check assigned ids of tiny clusters and assign them to the cells
 # with assigned ids.
 ##########################################
 GR.iteration = 1 # RG (revison global)
 Refine.annotated.ids = TRUE;

 resDir = paste0("results/", version.analysis, '/annoted_pharynx')
 if(!dir.exists(resDir)){dir.create(resDir)}

 if(Refine.annotated.ids){by.group = 'manual.annot.ids';
 }else{by.group = 'seurat_clusters'}

 RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
 RDS2save =  paste0(RdataDir,
                    'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                    '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

 seurat.obj = readRDS(file = RDSsaved)

 #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")

 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
   scale_colour_hue(drop = FALSE) +
   NoLegend()

 DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_BWM_manual.annoted.IDs")) +
   scale_colour_hue(drop = FALSE) +
   NoLegend()

 cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' cells not annotated \n')
 cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))), ' pharynx cells left to annotate \n')

 #pdf(pdfname, width=18, height = 10)
 #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)

 ##########################################
 # select subset of cells
 ##########################################
 jj = which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))
 table(seurat.obj$seurat_clusters_pharynx[jj])

 # select cells with cluster index
 ##########################################
 #cluster.sels = c('29', '32', '35', '40', '42')
 #cluster.sels = c('6', '24', '20', '7')
 #cluster.sels = c('13', '17', '18', '11', '25', '10')
 # cluster.sels = c('1', '9', '12', '19')
 # #ids.sel = c('MSxaaa')
 # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
 cluster.sels = c('3', '4', '8', '22', '0', '21', '23')

 cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
 #
 # if(length(ids.sel)>0){
 #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
 #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
 #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
 # }else{
 #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
 # }

 # select cells with ids
 ##########################################
 index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
 table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

 ids.current = names(table(seurat.obj$manual.annot.ids[index.pharynx.cells]))
 ids.current = ids.current[grep('MSxaap|MSaaapp|MSxapa', ids.current, invert = TRUE)]
 #ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
 #ids.sels = ids.current
 ids.sels = setdiff(ids.current,
                    c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp',
                      'MSxaaa', 'MSpaaapp', 'MSpaaapa', 'MSaappa', 'MSpaaaap.sure',
                       "MSaaaaa",  "MSaaaap", "MSpaaaa", "MSpaaap"))

 #ids.sels = c("MSaaaaa", "MSaaaap", "MSpaaaa", "MSpaaap", "MSxaaa.like.or.others")
 ids.left = setdiff(ids.current, ids.sels)
 # print(ids.left)
 # nchar(ids.left)

 cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

 ##########################################
 # subset seurat object with selected cells
 ##########################################
 cat(length(cells.sels), ' cells selected to annotate \n')
 sub.obj = subset(seurat.obj, cells = cells.sels)

 sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
 sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

 DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

 ##########################################
 # find new set of variable genes and redo pca and umap
 ##########################################
 Explore.umap.parameters.for.BWMcells = FALSE
 if(Explore.umap.parameters.for.BWMcells){

   source.my.script('scRNA_functions.R')
   require(tictoc)
   tic()
   test.umap.params.for.BWM.cells(sub.obj,
                                  pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                  group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                  nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                  n.neighbors.sampling = c(10, 30, 50),
                                  min.dist.sampling = c(0.01, 0.1)
   )
   toc()

 }

 nfeatures = 1000;
 sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
 #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
 sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
 sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
 ElbowPlot(sub.obj, ndims = 50)

 nb.pcs = 1 # nb of pcs depends on the considered clusters or ids
 n.neighbors = 10;
 min.dist = 0.01; spread = 1
 sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                    spread = spread, n.neighbors = n.neighbors,
                    min.dist = min.dist, verbose = TRUE)
 if(!Refine.annotated.ids){
   DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()
 }else{
   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()
 }

 p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
   NoLegend()

 #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
 p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
   NoLegend() + ggtitle('scamp prediction')

 p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
   NoLegend() + ggtitle('seurat prediction')

 p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
   NoLegend()

 p0 + p2
 p0 + p1
 p1 + p2

 p2 + p3

 ##########################################
 # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
 ##########################################
 FindClusters_subclusters = function(sub.obj, resolution = 1.0)
 {
   sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
   return(sub.obj$seurat_clusters)
 }
 sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 5, dims = 1:10, compute.SNN = TRUE)
 sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 3.0)
 DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

 p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
   ggtitle('manual.ids') + NoLegend()
 p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
              label.size = 8, na.value = "gray", combine = TRUE)
 p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

 (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration, '.pdf'),
                          width = 18, height = 16)

 VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
         group.by = 'seurat_clusters_split')
 #dev.off()
 ##########################################
 # check the counts of predicted ids for newly split clusters
 ##########################################

 idents.sel = c('0', '2')
 Idents(sub.obj) = sub.obj$seurat_clusters_split
 counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
 counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
 counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
 if(exists('idents.sel')){
   if(length(setdiff(colnames(counts), idents.sel)) > 0){
     counts = counts[, !is.na(match(colnames(counts), idents.sel))]
     counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
     counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
     counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
     counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
     counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
   }
 }

 # features of early stage
 features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                   'ceh-36', # MSxaap
                   'tbx-2', 'ceh-27') # MSxaaa
 # features of lineage MSxapa
 features.sels = c('dmd-4', 'cnd-1', 'swt-3', # MSxapa
                   'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                   'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                   'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                   'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                   'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                   'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                   'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

 )
 # features of MSxaap lineage
 features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                   'sptf-1', 'ceh-22',  #MSaaapaa
                   'lin-12', 'pax-1', # MSxaapap
                   'ngn-1', 'ces-1', # MSpaapaa
                   'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                   'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                   'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                   'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
 )
 # features of MSxaaa
 features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                   'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                   #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                   'ceh-13', # MSaaaap not MSpaaap
                   'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                   'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                   'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                   'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                   'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                   'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                   'igcm-4', 'F54E2.2', 'K04G2.12',
                   'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                   'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6' # MSaaaaapa
                   )

 features.sels = c('hlh-3', 'hlh-6'
                   )
 FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

 # update the annotation with marker genes
 cluster.assingment = list(
   c('4', 'MSaaaaaa'),
   c('6', 'MSaaaaaaa'),
   c('5', 'MSaaaaap'),
   c('3', 'MSaaaapa'),
   c('7', 'MSaaaapa'),
   c('1', 'MSaaaapp'),
   c('0', 'MSpaaaaa'),
   c('2', 'MSaaaaapa')
 )

 cat(length(cluster.assingment), 'clusters assigned \n')
 cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
 nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
 if(length(nb.unassigned) > 1){
   cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
 }

 # check info in JM data for specific lineage
 markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
 #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
 #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
 source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

 ids.sel = c('MSaaaaapa'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


 ##########################################
 # update the manual annotation if good marker genes or mapped ids were found
 ##########################################
 add.missing.annotation.for.pharynx.left.over = FALSE
 if(add.missing.annotation.for.pharynx.left.over){
   jj = which(is.na(sub.obj$manual.annot.ids))

   for(j in jj){
     #j = 1;
     cat('cell : ',  j, '\n')
     index.cluster = sub.obj$seurat_clusters_split[j]
     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
   }

   mm = match(colnames(sub.obj), colnames(seurat.obj))
   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids

 }

 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

 for(n in 1:length(cluster.assingment)){
   cluster.index = cluster.assingment[[n]][1];
   id2assign =  cluster.assingment[[n]][2];
   cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
   cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
   sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
   seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
   # if(is.na(id2assign)) {
   #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
   # }else{
   #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
   # }
 }

 #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

 DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
         pt.size = 2) + NoLegend()

 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
   scale_colour_hue(drop = FALSE) + NoLegend()

 saveRDS(seurat.obj, file = RDS2save)




 ##########################################
# Main aim:
# I realized that there are many cells in seurat_clusters 42 and 44 without annotation and I will look into them
#
# Notes:
# 17 cells in cluster 42 were annotated as MSxapa
# cluster 44 were not annotated at all and failed to be annotated again, because no clear prediction were found in cluster 44
# not much progress done for cluster 28, 52 and 31 except a likely MSxapappp
##########################################
GR.iteration = 2 # RG (revison global)
Refine.annotated.ids = FALSE;

resDir = paste0("results/", version.analysis, '/annoted_pharynx')
if(!dir.exists(resDir)){dir.create(resDir)}
if(Refine.annotated.ids){by.group = 'manual.annot.ids';
}else{by.group = 'seurat_clusters'}


RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                               '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
RDS2save =  paste0(RdataDir,
                   'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                   '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

seurat.obj = readRDS(file = RDSsaved)

#pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
        na.value = "gray") +
  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
  scale_colour_hue(drop = FALSE) +
  NoLegend()

DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
        na.value = "gray") +
  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
  scale_colour_hue(drop = FALSE) +
  NoLegend()

#cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
#    ' pharynx cells left to annotate \n')

cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

#pdf(pdfname, width=18, height = 10)
#par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
##########################################
# select subset of cells
##########################################
jj = which(is.na(seurat.obj$manual.annot.ids))
print(table(seurat.obj$seurat_clusters[jj]))


# select cells with cluster index
##########################################
#cluster.sels = c('29', '32', '35', '40', '42')
#cluster.sels = c('6', '24', '20', '7')
#cluster.sels = c('13', '17', '18', '11', '25', '10')
# cluster.sels = c('1', '9', '12', '19')
# #ids.sel = c('MSxaaa')
# #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
#cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
#cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]

cluster.sels = c('28', '52', '31')
cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                 # & is.na(seurat.obj$manual.annot.ids)
                                  ]

DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
        cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
        pt.size = 2, label.size = 5) + NoLegend()

#
# if(length(ids.sel)>0){
#   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
#                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
#                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
# }else{
#   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
# }

# select cells with ids
##########################################
index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

ids.current = names(table(seurat.obj$manual.annot.ids[index.pharynx.cells]))
ids.current = ids.current[grep('MSxaap|MSaaapp|MSxapa', ids.current, invert = TRUE)]
#ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
#ids.sels = ids.current
ids.sels = setdiff(ids.current,
                   c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp',
                     'MSxaaa', 'MSpaaapp', 'MSpaaapa', 'MSaappa', 'MSpaaaap.sure',
                      "MSaaaaa",  "MSaaaap", "MSpaaaa", "MSpaaap"))

#ids.sels = c("MSaaaaa", "MSaaaap", "MSpaaaa", "MSpaaap", "MSxaaa.like.or.others")
ids.left = setdiff(ids.current, ids.sels)
# print(ids.left)
# nchar(ids.left)

cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

##########################################
# subset seurat object with selected cells
##########################################
cat(length(cells.sels), ' cells selected to annotate \n')

sub.obj = subset(seurat.obj, cells = cells.sels)

#sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

##########################################
# find new set of variable genes and redo pca and umap
##########################################
Explore.umap.parameters.for.BWMcells = FALSE
if(Explore.umap.parameters.for.BWMcells){

  source.my.script('scRNA_functions.R')
  require(tictoc)
  tic()
  test.umap.params.for.BWM.cells(sub.obj,
                                 pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                 group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                 nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                 n.neighbors.sampling = c(10, 30, 50),
                                 min.dist.sampling = c(0.01, 0.1)
  )
  toc()

}

nfeatures = 2000;
sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
#cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(sub.obj, ndims = 50)

nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
n.neighbors = 10;
min.dist = 0.01; spread = 1
sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                   spread = spread, n.neighbors = n.neighbors,
                   min.dist = min.dist, verbose = TRUE)
if(!Refine.annotated.ids){
  DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
    NoLegend()
}else{
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
    NoLegend()
}

p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
  NoLegend()

#p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
  NoLegend() + ggtitle('scamp prediction')

p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
  NoLegend() + ggtitle('seurat prediction')

p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
  NoLegend()

p0 + p2
p0 + p1
p1 + p2

p2 + p3

##########################################
# redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
##########################################
FindClusters_subclusters = function(sub.obj, resolution = 1.0)
{
  sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
  return(sub.obj$seurat_clusters)
}
sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
  ggtitle('manual.ids') + NoLegend()
p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
             label.size = 8, na.value = "gray", combine = TRUE)
p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

(p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration, '.pdf'),
                         width = 18, height = 16)

jj.missing = which(is.na(sub.obj$manual.annot.ids))
table(sub.obj$seurat_clusters_split[jj.missing])

VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
        group.by = 'seurat_clusters_split')
#dev.off()
##########################################
# check the counts of predicted ids for newly split clusters
##########################################

idents.sel = c('1', '3', '4')
Idents(sub.obj) = sub.obj$seurat_clusters_split
counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
if(exists('idents.sel')){
  if(length(setdiff(colnames(counts), idents.sel)) > 0){
    counts = counts[, !is.na(match(colnames(counts), idents.sel))]
    counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
    counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
    counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
    counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
    counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
  }
}

# features of early stage
features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                  'ceh-36', # MSxaap
                  'tbx-2', 'ceh-27') # MSxaaa
# features of lineage MSxapa
features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                  'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                  'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                  'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                  'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                  'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                  'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                  'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

)
# features of MSxaap lineage
features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                  'sptf-1', 'ceh-22',  #MSaaapaa
                  'lin-12', 'pax-1', # MSxaapap
                  'ngn-1', 'ces-1', # MSpaapaa
                  'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                  'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                  'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                  'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
)
# features of MSxaaa
features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                  'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                  #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                  'ceh-13', # MSaaaap not MSpaaap
                  'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                  'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                  'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                  'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                  'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                  'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                  'igcm-4', 'F54E2.2', 'K04G2.12',
                  'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                  'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                  'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                  )

features.sels = c('pal-1'
                  )
FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

# update the annotation with marker genes
cluster.assingment = list(
  c('2', 'MSxapappp.likely'),
  c('6', 'MSxpapaa')

)

cat(length(cluster.assingment), 'clusters assigned \n')
cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
if(nb.unassigned > 1){
  cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
}

# check info in JM data for specific lineage
markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
#markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
#load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

ids.sel = c('MSppaapp'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)


##########################################
# update the manual annotation if good marker genes or mapped ids were found
##########################################
# add.missing.annotation.for.pharynx.left.over = FALSE
# if(add.missing.annotation.for.pharynx.left.over){
#   jj = which(is.na(sub.obj$manual.annot.ids))
#
#   for(j in jj){
#     #j = 1;
#     cat('cell : ',  j, '\n')
#     index.cluster = sub.obj$seurat_clusters_split[j]
#     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
#     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
#   }
#
#   mm = match(colnames(sub.obj), colnames(seurat.obj))
#   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
#
# }

# sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
# sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
# sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
# sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

for(n in 1:length(cluster.assingment)){
  cluster.index = cluster.assingment[[n]][1];
  id2assign =  cluster.assingment[[n]][2];
  cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
  cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
  sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
  seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
  # if(is.na(id2assign)) {
  #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
  # }else{
  #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
  # }
}

#seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
        pt.size = 2) + NoLegend()

DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
        na.value = "gray") +
  ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
  scale_colour_hue(drop = FALSE) + NoLegend()


saveRDS(seurat.obj, file = RDS2save)


##########################################
  # Main aim:
  # revise again some clusters: cluster.sels = c('28', '52', '31', '50', '17', '43')
  #
  # Notes:
  #
  # cluster 43, 17 were annotated as before, nothing changes;
  # cluster 50 were annotated as likely one
  ##########################################
  GR.iteration = 3 # RG (revison global)
  Refine.annotated.ids = FALSE;

  resDir = paste0("results/", version.analysis, '/annoted_pharynx')
  if(!dir.exists(resDir)){dir.create(resDir)}
  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}


  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                 '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
  RDS2save =  paste0(RdataDir,
                     'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                     '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

  seurat.obj = readRDS(file = RDSsaved)

  #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
  #    ' pharynx cells left to annotate \n')

  cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
  cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

  #pdf(pdfname, width=18, height = 10)
  #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  ##########################################
  # select subset of cells
  ##########################################

  # select cells with cluster index
  ##########################################
  jj = which(is.na(seurat.obj$manual.annot.ids))
  print(table(seurat.obj$seurat_clusters[jj]))

  #cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('6', '24', '20', '7')
  #cluster.sels = c('13', '17', '18', '11', '25', '10')
  # cluster.sels = c('1', '9', '12', '19')
  # #ids.sel = c('MSxaaa')
  # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
  #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]

  cluster.sels = c('28', '52', '31', '50', '17', '43')
  cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                   # & is.na(seurat.obj$manual.annot.ids)
                                    ]
  table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
          cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
          pt.size = 2, label.size = 5) + NoLegend()

  #
  # if(length(ids.sel)>0){
  #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
  #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
  #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
  # }else{
  #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
  # }

  # select cells with ids
  ##########################################
  #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
  #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

  clean.id.names = FALSE # clean id names
  if(clean.id.names){
    jj = which(seurat.obj$manual.annot.ids == "MSpaaaap.sure")
    seurat.obj$manual.annot.ids[jj] = "MSpaaaap"
    jj = which(seurat.obj$manual.annot.ids == "MSaaaappp.MSxapaapp")
    seurat.obj$manual.annot.ids[jj] = "MSaaaappp/MSxapaapp"

    jj = which(seurat.obj$manual.annot.ids == "MSpaaappp.MSxapappa")
    seurat.obj$manual.annot.ids[jj] = "MSpaaappp/MSxapappa"

    jj = which(seurat.obj$manual.annot.ids == 'unknown_MSxpppaa_MSxppppa_later')
    seurat.obj$manual.annot.ids[jj] = "MSxpppaa_MSxppppa_later_unknown"
  }

  ids.current = names(table(seurat.obj$manual.annot.ids))
  ids.current = ids.current[grep('MSxaap|MSaaapp|MSxapa', ids.current, invert = TRUE)]
  #ids.sels = c('MSxapa', 'MSxapap', 'MSxapaa', 'MSxapaap/MSxapaaa')
  #ids.sels = ids.current
  ids.sels = setdiff(ids.current,
                     c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp',
                       'MSxaaa', 'MSpaaapp', 'MSpaaapa', 'MSaappa', 'MSpaaaap.sure',
                        "MSaaaaa",  "MSaaaap", "MSpaaaa", "MSpaaap"))

  #ids.sels = c("MSaaaaa", "MSaaaap", "MSpaaaa", "MSpaaap", "MSxaaa.like.or.others")
  ids.left = setdiff(ids.current, ids.sels)
  # print(ids.left)
  # nchar(ids.left)

  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

  ##########################################
  # subset seurat object with selected cells
  ##########################################
  cat(length(cells.sels), ' cells selected to annotate \n')

  sub.obj = subset(seurat.obj, cells = cells.sels)

  #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  Explore.umap.parameters.for.BWMcells = FALSE
  if(Explore.umap.parameters.for.BWMcells){

    source.my.script('scRNA_functions.R')
    require(tictoc)
    tic()
    test.umap.params.for.BWM.cells(sub.obj,
                                   pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                   group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                   nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                   n.neighbors.sampling = c(10, 30, 50),
                                   min.dist.sampling = c(0.01, 0.1)
    )
    toc()

  }

  nfeatures = 2000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 10;
  min.dist = 0.01; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist, verbose = TRUE)
  if(!Refine.annotated.ids){
    DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()
  }else{
    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()
  }

  p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend()

  #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
  p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend() + ggtitle('scamp prediction')

  p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
    NoLegend() + ggtitle('seurat prediction')

  p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
    NoLegend()

  p0 + p2
  p0 + p1
  p1 + p2

  p2 + p3

  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 1.0)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
    ggtitle('manual.ids') + NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
               label.size = 8, na.value = "gray", combine = TRUE)
  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

  (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration, '.pdf'),
                           width = 18, height = 16)

  #jj.missing = which(is.na(sub.obj$manual.annot.ids))
  #table(sub.obj$seurat_clusters_split[jj.missing])

  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'seurat_clusters_split')
  #dev.off()
  ##########################################
  # check the counts of predicted ids for newly split clusters
  ##########################################

  idents.sel = c('7', '4', '5', '2')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
  counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
  if(exists('idents.sel')){
    if(length(setdiff(colnames(counts), idents.sel)) > 0){
      counts = counts[, !is.na(match(colnames(counts), idents.sel))]
      counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
      counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
      counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
      counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
      counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
    }
  }

  # features of early stage
  features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                    'ceh-36', # MSxaap
                    'tbx-2', 'ceh-27') # MSxaaa
  # features of lineage MSxapa
  features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                    'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                    'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                    'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                    'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                    'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                    'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                    'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

  )
  # features of MSxaap lineage
  features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                    'sptf-1', 'ceh-22',  #MSaaapaa
                    'lin-12', 'pax-1', # MSxaapap
                    'ngn-1', 'ces-1', # MSpaapaa
                    'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                    'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                    'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                    'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
  )
  # features of MSxaaa
  features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                    'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                    #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                    'ceh-13', # MSaaaap not MSpaaap
                    'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                    'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                    'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                    'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                    'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                    'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                    'igcm-4', 'F54E2.2', 'K04G2.12',
                    'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                    'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                    'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                    )

  features.sels = c('pal-1', 'unc-39'
                    )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  # update the annotation with marker genes
  cluster.assingment = list(
    c('2', 'MSxpapaa'),
    c('4', 'MSxpapaa'),
    c('5', 'MSxpapaa'),
    c('7', 'MSxpaap.MSppaapp.likely')

  )

  cat(length(cluster.assingment), 'clusters assigned \n')
  cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
  nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
  if(nb.unassigned > 1){
    cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
  }

  # check info in JM data for specific lineage
  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
  #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

  ids.sel = c('MSppaapp'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  # add.missing.annotation.for.pharynx.left.over = FALSE
  # if(add.missing.annotation.for.pharynx.left.over){
  #   jj = which(is.na(sub.obj$manual.annot.ids))
  #
  #   for(j in jj){
  #     #j = 1;
  #     cat('cell : ',  j, '\n')
  #     index.cluster = sub.obj$seurat_clusters_split[j]
  #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
  #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
  #   }
  #
  #   mm = match(colnames(sub.obj), colnames(seurat.obj))
  #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
  #
  # }

  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    # if(is.na(id2assign)) {
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
    # }else{
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
    # }
  }

  #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
          pt.size = 2) + NoLegend()

  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE) + NoLegend()


  saveRDS(seurat.obj, file = RDS2save)


  ##########################################
    # Main aim:
    # clean a bit the id names
    # last time revise the clusters 31, 52, 28
    #
    # Notes:
    # most of them were reput NA except MSxapappp.likely, MSxpapaa, MSxpaaaa and so on
    #
    ##########################################
    GR.iteration = 4 # RG (revison global)
    Refine.annotated.ids = TRUE

    resDir = paste0("results/", version.analysis, '/annoted_pharynx')
    if(!dir.exists(resDir)){dir.create(resDir)}
    if(Refine.annotated.ids){by.group = 'manual.annot.ids';
    }else{by.group = 'seurat_clusters'}


    RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                   '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
    RDS2save =  paste0(RdataDir,
                       'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                       '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

    seurat.obj = readRDS(file = RDSsaved)

    #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
    DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
            na.value = "gray") +
      ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
      scale_colour_hue(drop = FALSE) +
      NoLegend()

    DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
            na.value = "gray") +
      ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
      scale_colour_hue(drop = FALSE) +
      NoLegend()

    #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
    #    ' pharynx cells left to annotate \n')

    cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
    cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

    #pdf(pdfname, width=18, height = 10)
    #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    ##########################################
    # select subset of cells
    ##########################################

    # select cells with cluster index
    ##########################################
    jj = which(is.na(seurat.obj$manual.annot.ids))
    print(table(seurat.obj$seurat_clusters[jj]))

    #cluster.sels = c('29', '32', '35', '40', '42')
    #cluster.sels = c('6', '24', '20', '7')
    #cluster.sels = c('13', '17', '18', '11', '25', '10')
    # cluster.sels = c('1', '9', '12', '19')
    # #ids.sel = c('MSxaaa')
    # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
    #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
    #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]

    cluster.sels = c('28', '52', '31')
    cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                     # & is.na(seurat.obj$manual.annot.ids)
                                      ]
    table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
    DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
            cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
            pt.size = 2, label.size = 5) + NoLegend()

    #
    # if(length(ids.sel)>0){
    #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
    #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
    #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
    # }else{
    #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
    # }

    # select cells with ids
    ##########################################
    #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
    #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

    clean.id.names = FALSE # clean id names
    if(clean.id.names){
      jj = which(seurat.obj$manual.annot.ids == "MSpaaaap.sure")
      seurat.obj$manual.annot.ids[jj] = "MSpaaaap"
      jj = which(seurat.obj$manual.annot.ids == "MSaaaappp.MSxapaapp")
      seurat.obj$manual.annot.ids[jj] = "MSaaaappp/MSxapaapp"

      jj = which(seurat.obj$manual.annot.ids == "MSpaaappp.MSxapappa")
      seurat.obj$manual.annot.ids[jj] = "MSpaaappp/MSxapappa"

      jj = which(seurat.obj$manual.annot.ids == 'unknown_MSxpppaa_MSxppppa_later')
      seurat.obj$manual.annot.ids[jj] = "MSxpppaa_MSxppppa_later_unknown"

      jj = which(seurat.obj$manual.annot.ids == 'unknown.MSxpppaa')
      seurat.obj$manual.annot.ids[jj] = NA

      jj = which(seurat.obj$manual.annot.ids == 'MSxapapp.pharynx')
      seurat.obj$manual.annot.ids[jj] = NA

      jj = which(seurat.obj$manual.annot.ids == 'likely.nonBWM_origCluster_17')
      seurat.obj$manual.annot.ids[jj] = NA

      jj = which(seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31')
      seurat.obj$manual.annot.ids[jj] = NA

    }

    ids.current = names(table(seurat.obj$manual.annot.ids))
    #ids.sels = ids.current
    ids.sels = setdiff(ids.current,
                       c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp',
                         'MSxaaa', 'MSpaaapp', 'MSpaaapa', 'MSaappa', 'MSpaaaap.sure',
                          "MSaaaaa",  "MSaaaap", "MSpaaaa", "MSpaaap"))


    ids.sels = c('MSxppap/MSxpaaa/MSxpaaap', 'MSxppap', 'MSxpaaa', 'MSxpaaap',
                 "MSxppapp/MSxpappp"
                 )

    ids.left = setdiff(ids.current, ids.sels)
    # print(ids.left)
    # nchar(ids.left)

    cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

    ##########################################
    # subset seurat object with selected cells
    ##########################################
    cat(length(cells.sels), ' cells selected to annotate \n')

    sub.obj = subset(seurat.obj, cells = cells.sels)

    #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
    sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

    DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

    ##########################################
    # find new set of variable genes and redo pca and umap
    ##########################################
    Explore.umap.parameters.for.BWMcells = FALSE
    if(Explore.umap.parameters.for.BWMcells){

      source.my.script('scRNA_functions.R')
      require(tictoc)
      tic()
      test.umap.params.for.BWM.cells(sub.obj,
                                     pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                     group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                     nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                     n.neighbors.sampling = c(10, 30, 50),
                                     min.dist.sampling = c(0.01, 0.1)
      )
      toc()

    }

    nfeatures = 2000;
    sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
    #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
    sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
    sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
    ElbowPlot(sub.obj, ndims = 50)

    nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
    n.neighbors = 20;
    min.dist = 0.01; spread = 1
    sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                       spread = spread, n.neighbors = n.neighbors,
                       min.dist = min.dist, verbose = TRUE)
    if(Refine.annotated.ids){
      DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
        NoLegend()

    }else{
      DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
        NoLegend()

      p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
        NoLegend()

      #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
      p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
        NoLegend() + ggtitle('scamp prediction')

      p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
        NoLegend() + ggtitle('seurat prediction')

      p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
        NoLegend()

      p0 + p2
      p0 + p1
      p1 + p2
      p2 + p3

    }



    ##########################################
    # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
    ##########################################
    FindClusters_subclusters = function(sub.obj, resolution = 1.0)
    {
      sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
      return(sub.obj$seurat_clusters)
    }
    sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
    sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
    DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

    p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
      ggtitle('manual.ids') + NoLegend()
    p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
                 label.size = 8, na.value = "gray", combine = TRUE)
    p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

    (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration, '.pdf'),
                             width = 18, height = 16)

    #jj.missing = which(is.na(sub.obj$manual.annot.ids))
    #table(sub.obj$seurat_clusters_split[jj.missing])

    VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
            group.by = 'seurat_clusters_split')
    #dev.off()
    ##########################################
    # check the counts of predicted ids for newly split clusters
    ##########################################

    idents.sel = c('9', '6', '3', '10')
    Idents(sub.obj) = sub.obj$seurat_clusters_split
    counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
    counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
    counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
    if(exists('idents.sel')){
      if(length(setdiff(colnames(counts), idents.sel)) > 0){
        counts = counts[, !is.na(match(colnames(counts), idents.sel))]
        counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
        counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
        counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
        counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
        counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
      }
    }

    # features of early stage
    features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                      'ceh-36', # MSxaap
                      'tbx-2', 'ceh-27') # MSxaaa
    # features of lineage MSxapa
    features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                      'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                      'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                      'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                      'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                      'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                      'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                      'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

    )
    # features of MSxaap lineage
    features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                      'sptf-1', 'ceh-22',  #MSaaapaa
                      'lin-12', 'pax-1', # MSxaapap
                      'ngn-1', 'ces-1', # MSpaapaa
                      'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                      'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                      'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                      'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
    )
    # features of MSxaaa
    features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                      'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                      #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                      'ceh-13', # MSaaaap not MSpaaap
                      'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                      'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                      'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                      'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                      'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                      'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                      'igcm-4', 'F54E2.2', 'K04G2.12',
                      'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                      'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                      'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                      )

    features.sels = c('maph-1.3', 'shc-2', 'K09G1.1', 'tbx-2', 'twk-31', 'zig-6', 'mec-2', 'stn-2'
                      )
    FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

    # update the annotation with marker genes
    cluster.assingment = list(
      c('0', NA),
      c('2', 'MSxapappp.likely'),
      c('1', NA),
      c('3', NA),
      c('4', NA)
    )

    cat(length(cluster.assingment), 'clusters assigned \n')
    cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
    nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
    if(nb.unassigned > 1){
      cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
    }

    # check info in JM data for specific lineage
    markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
    #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
    #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
    source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

    ids.sel = c('MSppaapp'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
    #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

    ##########################################
    # update the manual annotation if good marker genes or mapped ids were found
    ##########################################
    # add.missing.annotation.for.pharynx.left.over = FALSE
    # if(add.missing.annotation.for.pharynx.left.over){
    #   jj = which(is.na(sub.obj$manual.annot.ids))
    #
    #   for(j in jj){
    #     #j = 1;
    #     cat('cell : ',  j, '\n')
    #     index.cluster = sub.obj$seurat_clusters_split[j]
    #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
    #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
    #   }
    #
    #   mm = match(colnames(sub.obj), colnames(seurat.obj))
    #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
    #
    # }

    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
    # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

    for(n in 1:length(cluster.assingment)){
      cluster.index = cluster.assingment[[n]][1];
      id2assign =  cluster.assingment[[n]][2];
      cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
      cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
      sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
      seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
      # if(is.na(id2assign)) {
      #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
      # }else{
      #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
      # }
    }

    #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
            pt.size = 2) + NoLegend()


    DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
            na.value = "gray") +
      ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
      scale_colour_hue(drop = FALSE) + NoLegend()


    saveRDS(seurat.obj, file = RDS2save)




 ##########################################
 # Main aim:
 # there are two BWM ids MSxppapp/MSxpappp and MSxppap/MSxpaaa/MSxpaaap
 # especiall the latter needs to resolve
 #
 # Notes:
 # at the end I did not change the annotation for MSxppap/MSxpaaa/MSxpaaap, because it is not clear and other cell ids seem to be
 # quite well annotated; making changes with newly split clusters would not improve too much but add cluster-based bias.
 ##########################################
 GR.iteration = 5 # RG (revison global)
 Refine.annotated.ids = TRUE

 resDir = paste0("results/", version.analysis, '/annoted_pharynx')
 if(!dir.exists(resDir)){dir.create(resDir)}
 if(Refine.annotated.ids){by.group = 'manual.annot.ids';
 }else{by.group = 'seurat_clusters'}


 RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
 RDS2save =  paste0(RdataDir,
                    'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                    '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

 seurat.obj = readRDS(file = RDSsaved)

 #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
   scale_colour_hue(drop = FALSE) +
   NoLegend()

 DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
   scale_colour_hue(drop = FALSE) +
   NoLegend()

 #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
 #    ' pharynx cells left to annotate \n')

 cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
 cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

 #pdf(pdfname, width=18, height = 10)
 #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
 ##########################################
 # select subset of cells
 ##########################################

 # select cells with cluster index
 ##########################################
 jj = which(is.na(seurat.obj$manual.annot.ids))
 print(table(seurat.obj$seurat_clusters[jj]))

 #cluster.sels = c('29', '32', '35', '40', '42')
 #cluster.sels = c('6', '24', '20', '7')
 #cluster.sels = c('13', '17', '18', '11', '25', '10')
 # cluster.sels = c('1', '9', '12', '19')
 # #ids.sel = c('MSxaaa')
 # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
 #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
 #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
 cluster.sels = c('28', '52', '31')
 cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                  # & is.na(seurat.obj$manual.annot.ids)
                                   ]
 table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
         cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
         pt.size = 2, label.size = 5) + NoLegend()

 #
 # if(length(ids.sel)>0){
 #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
 #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
 #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
 # }else{
 #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
 # }

 # select cells with ids
 ##########################################
 #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
 #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

 clean.id.names = FALSE # clean id names
 if(clean.id.names){
   # jj = which(seurat.obj$manual.annot.ids == "MSpaaaap.sure")
   # seurat.obj$manual.annot.ids[jj] = "MSpaaaap"
   # jj = which(seurat.obj$manual.annot.ids == "MSaaaappp.MSxapaapp")
   # seurat.obj$manual.annot.ids[jj] = "MSaaaappp/MSxapaapp"
   #
   # jj = which(seurat.obj$manual.annot.ids == "MSpaaappp.MSxapappa")
   # seurat.obj$manual.annot.ids[jj] = "MSpaaappp/MSxapappa"
   #
   # jj = which(seurat.obj$manual.annot.ids == 'unknown_MSxpppaa_MSxppppa_later')
   # seurat.obj$manual.annot.ids[jj] = "MSxpppaa_MSxppppa_later_unknown"
   #
   # jj = which(seurat.obj$manual.annot.ids == 'unknown.MSxpppaa')
   # seurat.obj$manual.annot.ids[jj] = NA
   #
   # jj = which(seurat.obj$manual.annot.ids == 'MSxapapp.pharynx')
   # seurat.obj$manual.annot.ids[jj] = NA
   #
   # jj = which(seurat.obj$manual.annot.ids == 'likely.nonBWM_origCluster_17')
   # seurat.obj$manual.annot.ids[jj] = NA
   #
   # jj = which(seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31')
   # seurat.obj$manual.annot.ids[jj] = NA

   # jj = which(seurat.obj$manual.annot.ids == 'MSxpppaa_MSxppppa_later_unknown')
   # seurat.obj$manual.annot.ids[jj] = 'MSxpppaa_MSxppppa_later_like'
   #
   # jj = which(seurat.obj$manual.annot.ids == 'MSxpaap.MSppaapp.likely')
   # seurat.obj$manual.annot.ids[jj] = 'MSxpaap.MSppaapp.like'
   #
   # jj = which(seurat.obj$manual.annot.ids == 'MSxapappp.likely')
   # seurat.obj$manual.annot.ids[jj] = 'MSxapappp.like'

   #jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")
   #seurat.obj$manual.annot.ids[jj] = 'mixture_transitionToTerminal'

 }

 jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")

 ids.current = names(table(seurat.obj$manual.annot.ids))
 #ids.sels = ids.current
 ids.sels = setdiff(ids.current,
                    c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp',
                      'MSxaaa', 'MSpaaapp', 'MSpaaapa', 'MSaappa', 'MSpaaaap.sure',
                       "MSaaaaa",  "MSaaaap", "MSpaaaa", "MSpaaap"))


 ids.sels = c('MSxppap/MSxpaaa/MSxpaaap',
              'MSxpppp', 'MSxpppa', 'MSxppap', 'MSxpapp'
              #'MSxpapa', 'MSxpaaa'
              #'MSxpaaap'
              #"MSxppapp/MSxpappp"
              )

 ids.left = setdiff(ids.current, ids.sels)
 # print(ids.left)
 # nchar(ids.left)

 cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

 ##########################################
 # subset seurat object with selected cells
 ##########################################
 cat(length(cells.sels), ' cells selected to annotate \n')

 sub.obj = subset(seurat.obj, cells = cells.sels)

 #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
 sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

 DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

 ##########################################
 # find new set of variable genes and redo pca and umap
 ##########################################
 Explore.umap.parameters.for.BWMcells = FALSE
 if(Explore.umap.parameters.for.BWMcells){

   source.my.script('scRNA_functions.R')
   require(tictoc)
   tic()
   test.umap.params.for.BWM.cells(sub.obj,
                                  pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                  group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                  nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                  n.neighbors.sampling = c(10, 30, 50),
                                  min.dist.sampling = c(0.01, 0.1)
   )
   toc()

 }

 nfeatures = 1000;
 sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
 #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
 sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
 sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
 ElbowPlot(sub.obj, ndims = 50)

 nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
 n.neighbors = 20;
 min.dist = 0.01; spread = 1
 sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                    spread = spread, n.neighbors = n.neighbors,
                    min.dist = min.dist, verbose = TRUE)
 DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
   NoLegend()


 if(Refine.annotated.ids){

   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()

 }else{
   DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()

   p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()

   #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
   p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend() + ggtitle('scamp prediction')

   p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
     NoLegend() + ggtitle('seurat prediction')

   p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()

   p0 + p2
   p0 + p1
   p1 + p2
   p2 + p3

 }

 ##########################################
 # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
 ##########################################
 FindClusters_subclusters = function(sub.obj, resolution = 1.0)
 {
   sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
   return(sub.obj$seurat_clusters)
 }
 sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:5, compute.SNN = TRUE)
 sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
 DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

 p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
   ggtitle('manual.ids') + NoLegend()
 p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
              label.size = 8, na.value = "gray", combine = TRUE)
 p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

 (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration, '.pdf'),
                          width = 18, height = 16)

 #jj.missing = which(is.na(sub.obj$manual.annot.ids))
 #table(sub.obj$seurat_clusters_split[jj.missing])

 VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
         group.by = 'seurat_clusters_split')
 #dev.off()
 ##########################################
 # check the counts of predicted ids for newly split clusters
 ##########################################

 idents.sel = c('0', '6', '3', '9')
 Idents(sub.obj) = sub.obj$seurat_clusters_split
 counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
 counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
 counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
 if(exists('idents.sel')){
   if(length(setdiff(colnames(counts), idents.sel)) > 0){
     counts = counts[, !is.na(match(colnames(counts), idents.sel))]
     counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
     counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
     counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
     counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
     counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
   }
 }

 # features of early stage
 features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                   'ceh-36', # MSxaap
                   'tbx-2', 'ceh-27') # MSxaaa
 # features of lineage MSxapa
 features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                   'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                   'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                   'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                   'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                   'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                   'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                   'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

 )
 # features of MSxaap lineage
 features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                   'sptf-1', 'ceh-22',  #MSaaapaa
                   'lin-12', 'pax-1', # MSxaapap
                   'ngn-1', 'ces-1', # MSpaapaa
                   'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                   'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                   'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                   'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
 )
 # features of MSxaaa
 features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                   'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                   #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                   'ceh-13', # MSaaaap not MSpaaap
                   'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                   'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                   'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                   'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                   'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                   'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                   'igcm-4', 'F54E2.2', 'K04G2.12',
                   'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                   'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                   'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                   )

 features.sels = c('zig-8', 'pat-9', 'F19C6.4', 'tbx-8', 'her-1', 'ins-2'
                   )
 FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

 # update the annotation with marker genes
 cluster.assingment = list(
   c('4', 'MSxpaaa'),

   #c('4', ) # not MSxpaaa; not MSxpppa;

 )

 cat(length(cluster.assingment), 'clusters assigned \n')
 cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
 nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
 if(nb.unassigned > 1){
   cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
 }

 # check info in JM data for specific lineage
 markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
 #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
 #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
 source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

 ids.sel = c('MSppaapp'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
 #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

 ##########################################
 # update the manual annotation if good marker genes or mapped ids were found
 ##########################################
 # add.missing.annotation.for.pharynx.left.over = FALSE
 # if(add.missing.annotation.for.pharynx.left.over){
 #   jj = which(is.na(sub.obj$manual.annot.ids))
 #
 #   for(j in jj){
 #     #j = 1;
 #     cat('cell : ',  j, '\n')
 #     index.cluster = sub.obj$seurat_clusters_split[j]
 #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
 #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
 #   }
 #
 #   mm = match(colnames(sub.obj), colnames(seurat.obj))
 #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
 #
 # }

 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

 for(n in 1:length(cluster.assingment)){
   cluster.index = cluster.assingment[[n]][1];
   id2assign =  cluster.assingment[[n]][2];
   cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
   cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
   sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
   seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
   # if(is.na(id2assign)) {
   #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
   # }else{
   #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
   # }
 }

 #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

 DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
         pt.size = 2) + NoLegend()


 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
   scale_colour_hue(drop = FALSE) + NoLegend()


 saveRDS(seurat.obj, file = RDS2save)




 ##########################################
 # Main aim:
 # try to solve MSxppap/MSxpaaa/MSxpaaap and mixture ids
 #
 #
 # Notes:
 # all mixture annotation were solved for the first time and requires comfirmation
 #
 ##########################################
 GR.iteration = 6 # RG (revison global)
 Refine.annotated.ids = TRUE

 resDir = paste0("results/", version.analysis, '/annoted_pharynx')
 if(!dir.exists(resDir)){dir.create(resDir)}
 if(Refine.annotated.ids){by.group = 'manual.annot.ids';
 }else{by.group = 'seurat_clusters'}


 RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
 RDS2save =  paste0(RdataDir,
                    'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                    '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

 seurat.obj = readRDS(file = RDSsaved)

 #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
   scale_colour_hue(drop = FALSE) +
   NoLegend()

 DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
   scale_colour_hue(drop = FALSE) +
   NoLegend()

 #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
 #    ' pharynx cells left to annotate \n')

 cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
 cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

 #pdf(pdfname, width=18, height = 10)
 #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
 ##########################################
 # select subset of cells
 ##########################################

 # select cells with cluster index
 ##########################################
 jj = which(is.na(seurat.obj$manual.annot.ids))
 print(table(seurat.obj$seurat_clusters[jj]))

 #cluster.sels = c('29', '32', '35', '40', '42')
 #cluster.sels = c('6', '24', '20', '7')
 #cluster.sels = c('13', '17', '18', '11', '25', '10')
 # cluster.sels = c('1', '9', '12', '19')
 # #ids.sel = c('MSxaaa')
 # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
 #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
 #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
 cluster.sels = c('28', '52', '31')
 cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                  # & is.na(seurat.obj$manual.annot.ids)
                                   ]
 table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
         cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
         pt.size = 2, label.size = 5) + NoLegend()

 #
 # if(length(ids.sel)>0){
 #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
 #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
 #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
 # }else{
 #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
 # }

 # select cells with ids
 ##########################################
 #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
 #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

 clean.id.names = FALSE # clean id names
 if(clean.id.names){
   # jj = which(seurat.obj$manual.annot.ids == "MSpaaaap.sure")
   # seurat.obj$manual.annot.ids[jj] = "MSpaaaap"
   # jj = which(seurat.obj$manual.annot.ids == "MSaaaappp.MSxapaapp")
   # seurat.obj$manual.annot.ids[jj] = "MSaaaappp/MSxapaapp"
   #
   # jj = which(seurat.obj$manual.annot.ids == "MSpaaappp.MSxapappa")
   # seurat.obj$manual.annot.ids[jj] = "MSpaaappp/MSxapappa"
   #
   # jj = which(seurat.obj$manual.annot.ids == 'unknown_MSxpppaa_MSxppppa_later')
   # seurat.obj$manual.annot.ids[jj] = "MSxpppaa_MSxppppa_later_unknown"
   #
   # jj = which(seurat.obj$manual.annot.ids == 'unknown.MSxpppaa')
   # seurat.obj$manual.annot.ids[jj] = NA
   #
   # jj = which(seurat.obj$manual.annot.ids == 'MSxapapp.pharynx')
   # seurat.obj$manual.annot.ids[jj] = NA
   #
   # jj = which(seurat.obj$manual.annot.ids == 'likely.nonBWM_origCluster_17')
   # seurat.obj$manual.annot.ids[jj] = NA
   #
   # jj = which(seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31')
   # seurat.obj$manual.annot.ids[jj] = NA

   # jj = which(seurat.obj$manual.annot.ids == 'MSxpppaa_MSxppppa_later_unknown')
   # seurat.obj$manual.annot.ids[jj] = 'MSxpppaa_MSxppppa_later_like'
   #
   # jj = which(seurat.obj$manual.annot.ids == 'MSxpaap.MSppaapp.likely')
   # seurat.obj$manual.annot.ids[jj] = 'MSxpaap.MSppaapp.like'
   #
   # jj = which(seurat.obj$manual.annot.ids == 'MSxapappp.likely')
   # seurat.obj$manual.annot.ids[jj] = 'MSxapappp.like'

   #jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")
   #seurat.obj$manual.annot.ids[jj] = 'mixture_transitionToTerminal'

 }


 ids.current = names(table(seurat.obj$manual.annot.ids))
 #ids.sels = ids.current
 ids.sels = setdiff(ids.current,
                    c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp',
                      'MSxaaa', 'MSpaaapp', 'MSpaaapa', 'MSaappa', 'MSpaaaap.sure',
                       "MSaaaaa",  "MSaaaap", "MSpaaaa", "MSpaaap"))


 ids.sels = c("mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap", "mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa",
              "mixture_transitionToTerminal"
              )

 ids.left = setdiff(ids.current, ids.sels)
 # print(ids.left)
 # nchar(ids.left)

 cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

 ##########################################
 # subset seurat object with selected cells
 ##########################################
 cat(length(cells.sels), ' cells selected to annotate \n')

 sub.obj = subset(seurat.obj, cells = cells.sels)

 #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
 sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

 DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

 ##########################################
 # find new set of variable genes and redo pca and umap
 ##########################################
 Explore.umap.parameters.for.BWMcells = FALSE
 if(Explore.umap.parameters.for.BWMcells){

   source.my.script('scRNA_functions.R')
   require(tictoc)
   tic()
   test.umap.params.for.BWM.cells(sub.obj,
                                  pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                  group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                  nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                  n.neighbors.sampling = c(10, 30, 50),
                                  min.dist.sampling = c(0.01, 0.1)
   )
   toc()

 }

 nfeatures = 500;
 sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
 #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
 sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
 sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
 ElbowPlot(sub.obj, ndims = 50)

 nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
 n.neighbors = 10;
 min.dist = 0.01; spread = 1
 sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                    spread = spread, n.neighbors = n.neighbors,
                    min.dist = min.dist, verbose = TRUE)
 DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
   NoLegend()


 if(Refine.annotated.ids){

   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()

 }else{
   DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()

   p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()

   #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
   p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend() + ggtitle('scamp prediction')

   p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
     NoLegend() + ggtitle('seurat prediction')

   p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()

   p0 + p2
   p0 + p1
   p1 + p2
   p2 + p3

 }

 ##########################################
 # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
 ##########################################
 FindClusters_subclusters = function(sub.obj, resolution = 1.0)
 {
   sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
   return(sub.obj$seurat_clusters)
 }
 sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
 sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.5)
 DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

 p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
   ggtitle('manual.ids') + NoLegend()
 p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
              label.size = 8, na.value = "gray", combine = TRUE)
 p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

 (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration, '.pdf'),
                          width = 18, height = 16)

 #jj.missing = which(is.na(sub.obj$manual.annot.ids))
 #table(sub.obj$seurat_clusters_split[jj.missing])

 VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
         group.by = 'seurat_clusters_split')
 #dev.off()
 ##########################################
 # check the counts of predicted ids for newly split clusters
 ##########################################
 idents.sel = c('7', '1', '9')
 Idents(sub.obj) = sub.obj$seurat_clusters_split
 counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
 counts.seurat = table(sub.obj$predicted.ids.seurat, sub.obj$seurat_clusters_split)
 counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
 if(exists('idents.sel')){
   if(length(setdiff(colnames(counts), idents.sel)) > 0){
     counts = counts[, !is.na(match(colnames(counts), idents.sel))]
     counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
     counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
     counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
     counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
     counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
   }
 }

 # features of early stage
 features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                   'ceh-36', # MSxaap
                   'tbx-2', 'ceh-27') # MSxaaa
 # features of lineage MSxapa
 features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                   'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                   'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                   'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                   'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                   'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                   'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                   'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

 )
 # features of MSxaap lineage
 features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                   'sptf-1', 'ceh-22',  #MSaaapaa
                   'lin-12', 'pax-1', # MSxaapap
                   'ngn-1', 'ces-1', # MSpaapaa
                   'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                   'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                   'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                   'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
 )
 # features of MSxaaa
 features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                   'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                   #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                   'ceh-13', # MSaaaap not MSpaaap
                   'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                   'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                   'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                   'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                   'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                   'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                   'igcm-4', 'F54E2.2', 'K04G2.12',
                   'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                   'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                   'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                   )

 features.sels = c('unc-62','gana-1','clec-264', 'maph-1.3', 'zig-6', 'maph-1.2',
                   'mec-2', 'twk-31', 'stn-2', 'ZC449.5',  'shc-2', 'ham-1', 'ceh-34',
                   'K09G1.1',
                   'ceh-13',  'tbx-2', 'D1086.12', 'hot-1', 'sul-1', 'zig-7', 'hmg-1.1'
                   #'gsnl-2', 'abts-1'
                   )
 FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

 # update the annotation with marker genes
 cluster.assingment = list(
   c('6', 'MSxppapp_new'),
   c('4', 'MSxppapp_new'),
   c('0', 'MSxpaaap_new'),
   c('5', 'MSxpaaap_new'),
   c('2', 'MSxpppaa_new'),
   c('8', 'MSxpppaa_new'),
   c('3', 'MSxpppaa_new'),
   c('1', 'MSxapppax_new'),
   c('7', 'MSxpappp_new'),
   c('9', 'MSxpapap_new')
 )

 cat(length(cluster.assingment), 'clusters assigned \n')
 cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
 nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
 if(nb.unassigned >= 1){
   cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
 }

 # check info in JM data for specific lineage
 markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
 #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
 #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
 source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

 ids.sel = c('MSxpaaap'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
 #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

 ##########################################
 # update the manual annotation if good marker genes or mapped ids were found
 ##########################################
 # add.missing.annotation.for.pharynx.left.over = FALSE
 # if(add.missing.annotation.for.pharynx.left.over){
 #   jj = which(is.na(sub.obj$manual.annot.ids))
 #
 #   for(j in jj){
 #     #j = 1;
 #     cat('cell : ',  j, '\n')
 #     index.cluster = sub.obj$seurat_clusters_split[j]
 #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
 #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
 #   }
 #
 #   mm = match(colnames(sub.obj), colnames(seurat.obj))
 #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
 #
 # }

 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

 for(n in 1:length(cluster.assingment)){
   cluster.index = cluster.assingment[[n]][1];
   id2assign =  cluster.assingment[[n]][2];
   cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
   cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
   sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
   seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
   # if(is.na(id2assign)) {
   #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
   # }else{
   #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
   # }
 }

 #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

 DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
         pt.size = 2) + NoLegend()


 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
   scale_colour_hue(drop = FALSE) + NoLegend()


 saveRDS(seurat.obj, file = RDS2save)

 ##########################################
 # Main aim:
 # here we solved the mixed terminal cells ids.sels = c("MSxppapp/MSxpappp")
 #
 #
 # Notes:
 # all mixture annotation were solved for the first time and requires comfirmation
 #
 ##########################################
 GR.iteration = 7 # RG (revison global)
 Refine.annotated.ids = TRUE

 resDir = paste0("results/", version.analysis, '/annoted_pharynx')
 if(!dir.exists(resDir)){dir.create(resDir)}
 if(Refine.annotated.ids){by.group = 'manual.annot.ids';
 }else{by.group = 'seurat_clusters'}


 RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
 RDS2save =  paste0(RdataDir,
                    'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                    '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

 seurat.obj = readRDS(file = RDSsaved)

 #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
   scale_colour_hue(drop = FALSE) +
   NoLegend()

 DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
   scale_colour_hue(drop = FALSE) +
   NoLegend()

 #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
 #    ' pharynx cells left to annotate \n')

 cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
 cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

 #pdf(pdfname, width=18, height = 10)
 #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
 ##########################################
 # select subset of cells
 ##########################################

 # select cells with cluster index
 ##########################################
 jj = which(is.na(seurat.obj$manual.annot.ids))
 print(table(seurat.obj$seurat_clusters[jj]))

 #cluster.sels = c('29', '32', '35', '40', '42')
 #cluster.sels = c('6', '24', '20', '7')
 #cluster.sels = c('13', '17', '18', '11', '25', '10')
 # cluster.sels = c('1', '9', '12', '19')
 # #ids.sel = c('MSxaaa')
 # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
 #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
 #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
 cluster.sels = c('28', '52', '31')
 cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                  # & is.na(seurat.obj$manual.annot.ids)
                                   ]
 table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
         cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
         pt.size = 2, label.size = 5) + NoLegend()

 #
 # if(length(ids.sel)>0){
 #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
 #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
 #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
 # }else{
 #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
 # }

 # select cells with ids
 ##########################################
 #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
 #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

 clean.id.names = FALSE # clean id names
 if(clean.id.names){
   # jj = which(seurat.obj$manual.annot.ids == "MSpaaaap.sure")
   # seurat.obj$manual.annot.ids[jj] = "MSpaaaap"
   # jj = which(seurat.obj$manual.annot.ids == "MSaaaappp.MSxapaapp")
   # seurat.obj$manual.annot.ids[jj] = "MSaaaappp/MSxapaapp"
   #
   # jj = which(seurat.obj$manual.annot.ids == "MSpaaappp.MSxapappa")
   # seurat.obj$manual.annot.ids[jj] = "MSpaaappp/MSxapappa"
   #
   # jj = which(seurat.obj$manual.annot.ids == 'unknown_MSxpppaa_MSxppppa_later')
   # seurat.obj$manual.annot.ids[jj] = "MSxpppaa_MSxppppa_later_unknown"
   #
   # jj = which(seurat.obj$manual.annot.ids == 'unknown.MSxpppaa')
   # seurat.obj$manual.annot.ids[jj] = NA
   #
   # jj = which(seurat.obj$manual.annot.ids == 'MSxapapp.pharynx')
   # seurat.obj$manual.annot.ids[jj] = NA
   #
   # jj = which(seurat.obj$manual.annot.ids == 'likely.nonBWM_origCluster_17')
   # seurat.obj$manual.annot.ids[jj] = NA
   #
   # jj = which(seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31')
   # seurat.obj$manual.annot.ids[jj] = NA

   # jj = which(seurat.obj$manual.annot.ids == 'MSxpppaa_MSxppppa_later_unknown')
   # seurat.obj$manual.annot.ids[jj] = 'MSxpppaa_MSxppppa_later_like'
   #
   # jj = which(seurat.obj$manual.annot.ids == 'MSxpaap.MSppaapp.likely')
   # seurat.obj$manual.annot.ids[jj] = 'MSxpaap.MSppaapp.like'
   #
   # jj = which(seurat.obj$manual.annot.ids == 'MSxapappp.likely')
   # seurat.obj$manual.annot.ids[jj] = 'MSxapappp.like'

   #jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")
   #seurat.obj$manual.annot.ids[jj] = 'mixture_transitionToTerminal'

 }


 ids.current = names(table(seurat.obj$manual.annot.ids))
 #ids.sels = ids.current
 ids.sels = setdiff(ids.current,
                    c('MSxaa', 'MSxap', 'MSpaapaa', 'MSaaapaa', 'MSaaapapp',
                      'MSxaaa', 'MSpaaapp', 'MSpaaapa', 'MSaappa', 'MSpaaaap.sure',
                       "MSaaaaa",  "MSaaaap", "MSpaaaa", "MSpaaap"))

 #ids.sels = c("mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap", "mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa",
 #             "mixture_transitionToTerminal")
 ids.sels = c("MSxppapp/MSxpappp")

 ids.left = setdiff(ids.current, ids.sels)
 # print(ids.left)
 # nchar(ids.left)

 cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

 ##########################################
 # subset seurat object with selected cells
 ##########################################
 cat(length(cells.sels), ' cells selected to annotate \n')

 sub.obj = subset(seurat.obj, cells = cells.sels)

 #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
 sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

 DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

 ##########################################
 # find new set of variable genes and redo pca and umap
 ##########################################
 Explore.umap.parameters.for.BWMcells = FALSE
 if(Explore.umap.parameters.for.BWMcells){

   source.my.script('scRNA_functions.R')
   require(tictoc)
   tic()
   test.umap.params.for.BWM.cells(sub.obj,
                                  pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                  group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                  nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                  n.neighbors.sampling = c(10, 30, 50),
                                  min.dist.sampling = c(0.01, 0.1)
   )
   toc()

 }

 nfeatures = 5000;
 sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
 #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
 sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
 sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
 ElbowPlot(sub.obj, ndims = 50)

 nb.pcs = 5 # nb of pcs depends on the considered clusters or ids
 n.neighbors = 10;
 min.dist = 0.01; spread = 1
 sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                    spread = spread, n.neighbors = n.neighbors,
                    min.dist = min.dist, verbose = TRUE)
 DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
   NoLegend()


 if(Refine.annotated.ids){

   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()

 }else{
   DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()

   p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()

   #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
   p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend() + ggtitle('scamp prediction')

   p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
     NoLegend() + ggtitle('seurat prediction')

   p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()

   p0 + p2
   p0 + p1
   p1 + p2
   p2 + p3

 }

 ##########################################
 # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
 ##########################################
 FindClusters_subclusters = function(sub.obj, resolution = 1.0)
 {
   sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
   return(sub.obj$seurat_clusters)
 }
 sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:5, compute.SNN = TRUE)
 sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
 DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

 p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
   ggtitle('manual.ids') + NoLegend()
 p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
              label.size = 8, na.value = "gray", combine = TRUE)
 p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

 (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration, '.pdf'),
                          width = 18, height = 16)

 #jj.missing = which(is.na(sub.obj$manual.annot.ids))
 #table(sub.obj$seurat_clusters_split[jj.missing])

 VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
         group.by = 'seurat_clusters_split')
 #dev.off()
 ##########################################
 # check the counts of predicted ids for newly split clusters
 ##########################################
 idents.sel = c('0', '1', '2', '3', '4')
 Idents(sub.obj) = sub.obj$seurat_clusters_split
 counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
 counts.seurat = table(sub.obj$predicted.ids.seurat.keep, sub.obj$seurat_clusters_split)
 counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
 if(exists('idents.sel')){
   if(length(setdiff(colnames(counts), idents.sel)) > 0){
     counts = counts[, !is.na(match(colnames(counts), idents.sel))]
     counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
     counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
     counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
     counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
     counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
   }
 }

 # features of early stage
 features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                   'ceh-36', # MSxaap
                   'tbx-2', 'ceh-27') # MSxaaa
 # features of lineage MSxapa
 features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                   'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                   'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                   'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                   'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                   'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                   'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                   'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

 )
 # features of MSxaap lineage
 features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                   'sptf-1', 'ceh-22',  #MSaaapaa
                   'lin-12', 'pax-1', # MSxaapap
                   'ngn-1', 'ces-1', # MSpaapaa
                   'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                   'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                   'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                   'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
 )
 # features of MSxaaa
 features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                   'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                   #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                   'ceh-13', # MSaaaap not MSpaaap
                   'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                   'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                   'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                   'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                   'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                   'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                   'igcm-4', 'F54E2.2', 'K04G2.12',
                   'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                   'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                   'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                   )

 features.sels = c('unc-62','gana-1','clec-264', 'maph-1.3', 'zig-6', 'maph-1.2',
                   'mec-2', 'twk-31', 'stn-2', 'ZC449.5',  'shc-2', 'ham-1', 'ceh-34',
                   'K09G1.1', 'T25G12.11', 'shc-2',
                   'ceh-13',  'tbx-2', 'D1086.12', 'sul-1'

                   #'gsnl-2', 'abts-1'
                   )
 FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

 # update the annotation with marker genes
 cluster.assingment = list(
   c('3', 'MSxpappp_new2'),
   c('1', 'MSxpappp_new2'),
   c('0', 'MSxppapp_new2'),
   c('2', 'MSxppapp_new2'),
   c('4', 'MSxppapp_new2')
 )

 cat(length(cluster.assingment), 'clusters assigned \n')
 cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
 nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
 if(nb.unassigned >= 1){
   cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
 }

 # check info in JM data for specific lineage
 markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
 #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
 #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
 source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

 ids.sel = c('MSxpaaap'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
 #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

 ##########################################
 # update the manual annotation if good marker genes or mapped ids were found
 ##########################################
 # add.missing.annotation.for.pharynx.left.over = FALSE
 # if(add.missing.annotation.for.pharynx.left.over){
 #   jj = which(is.na(sub.obj$manual.annot.ids))
 #
 #   for(j in jj){
 #     #j = 1;
 #     cat('cell : ',  j, '\n')
 #     index.cluster = sub.obj$seurat_clusters_split[j]
 #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
 #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
 #   }
 #
 #   mm = match(colnames(sub.obj), colnames(seurat.obj))
 #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
 #
 # }

 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

 for(n in 1:length(cluster.assingment)){
   cluster.index = cluster.assingment[[n]][1];
   id2assign =  cluster.assingment[[n]][2];
   cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
   cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
   sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
   seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
   # if(is.na(id2assign)) {
   #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
   # }else{
   #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
   # }
 }

 #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

 DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
         pt.size = 2) + NoLegend()


 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
   scale_colour_hue(drop = FALSE) + NoLegend()


 saveRDS(seurat.obj, file = RDS2save)





 ##########################################
   # Main aim:
   # in this iteration, BWM terminal cells were revised to confirm newly resolved ones
   # here MSxpaaap and MSxpapap were focused
   # Notes:
   # at the end, the annotation of terminal cells using terminal cells is still somewhat confusing, because the presence/absence of
   # some marker genes were inconsistant, either because of the cluster cells were mixture of cell ids or there were inaccuracies in
   # the marker genes
   ##########################################
   GR.iteration = 8 # RG (revison global)
   Refine.annotated.ids = TRUE

   resDir = paste0("results/", version.analysis, '/annoted_pharynx')
   if(!dir.exists(resDir)){dir.create(resDir)}
   if(Refine.annotated.ids){by.group = 'manual.annot.ids';
   }else{by.group = 'seurat_clusters'}


   RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                  '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
   RDS2save =  paste0(RdataDir,
                      'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                      '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

   seurat.obj = readRDS(file = RDSsaved)

   #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
   DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
     scale_colour_hue(drop = FALSE) +
     NoLegend()

   DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
     scale_colour_hue(drop = FALSE) +
     NoLegend()

   #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
   #    ' pharynx cells left to annotate \n')

   cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
   cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

   #pdf(pdfname, width=18, height = 10)
   #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
   ##########################################
   # select subset of cells
   ##########################################

   # select cells with cluster index
   ##########################################
   jj = which(is.na(seurat.obj$manual.annot.ids))
   print(table(seurat.obj$seurat_clusters[jj]))

   #cluster.sels = c('29', '32', '35', '40', '42')
   #cluster.sels = c('6', '24', '20', '7')
   #cluster.sels = c('13', '17', '18', '11', '25', '10')
   # cluster.sels = c('1', '9', '12', '19')
   # #ids.sel = c('MSxaaa')
   # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
   #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
   #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
   cluster.sels = c('28', '52', '31')
   cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                    # & is.na(seurat.obj$manual.annot.ids)
                                     ]
   table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
   DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
           cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
           pt.size = 2, label.size = 5) + NoLegend()

   #
   # if(length(ids.sel)>0){
   #   cells.sels = unique(colnames(seurat.obj)[(!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels)) |
   #                                            !is.na(match(seurat.obj$manual.annot.ids, ids.sel))) &
   #                                              is.na(match(seurat.obj$manual.annot.ids, ids.excl))])
   # }else{
   #   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))])
   # }

   # select cells with ids
   ##########################################
   #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
   #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')

   clean.id.names = FALSE # clean id names
   if(clean.id.names){
     # jj = which(seurat.obj$manual.annot.ids == "MSpaaaap.sure")
     # seurat.obj$manual.annot.ids[jj] = "MSpaaaap"
     # jj = which(seurat.obj$manual.annot.ids == "MSaaaappp.MSxapaapp")
     # seurat.obj$manual.annot.ids[jj] = "MSaaaappp/MSxapaapp"
     #
     # jj = which(seurat.obj$manual.annot.ids == "MSpaaappp.MSxapappa")
     # seurat.obj$manual.annot.ids[jj] = "MSpaaappp/MSxapappa"
     #
     # jj = which(seurat.obj$manual.annot.ids == 'unknown_MSxpppaa_MSxppppa_later')
     # seurat.obj$manual.annot.ids[jj] = "MSxpppaa_MSxppppa_later_unknown"
     #
     # jj = which(seurat.obj$manual.annot.ids == 'unknown.MSxpppaa')
     # seurat.obj$manual.annot.ids[jj] = NA
     #
     # jj = which(seurat.obj$manual.annot.ids == 'MSxapapp.pharynx')
     # seurat.obj$manual.annot.ids[jj] = NA
     #
     # jj = which(seurat.obj$manual.annot.ids == 'likely.nonBWM_origCluster_17')
     # seurat.obj$manual.annot.ids[jj] = NA
     #
     # jj = which(seurat.obj$manual.annot.ids == 'likely_nonBWM_origCluster_31')
     # seurat.obj$manual.annot.ids[jj] = NA

     # jj = which(seurat.obj$manual.annot.ids == 'MSxpppaa_MSxppppa_later_unknown')
     # seurat.obj$manual.annot.ids[jj] = 'MSxpppaa_MSxppppa_later_like'
     #
     # jj = which(seurat.obj$manual.annot.ids == 'MSxpaap.MSppaapp.likely')
     # seurat.obj$manual.annot.ids[jj] = 'MSxpaap.MSppaapp.like'
     #
     # jj = which(seurat.obj$manual.annot.ids == 'MSxapappp.likely')
     # seurat.obj$manual.annot.ids[jj] = 'MSxapappp.like'

     #jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")
     #seurat.obj$manual.annot.ids[jj] = 'mixture_transitionToTerminal'

   }


   ids.current = names(table(seurat.obj$manual.annot.ids))
   #ids.sels = ids.current

   ids.current[grep('_new', ids.current)]

   #ids.sels = c("mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap", "mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa",
   #             "mixture_transitionToTerminal")
   ids.sels = c(
               #'MSxppapp', 'MSxpappp', "MSxpappp_new", "MSxpappp_new2", "MSxppapp_new",  "MSxppapp_new2",
               'MSxpapap', 'MSxpaaap', "MSxpaaap_new","MSxpapap_new")

   ids.left = setdiff(ids.current, ids.sels)
   # print(ids.left)
   # nchar(ids.left)

   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

   ##########################################
   # subset seurat object with selected cells
   ##########################################
   cat(length(cells.sels), ' cells selected to annotate \n')

   sub.obj = subset(seurat.obj, cells = cells.sels)

   #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
   sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

   DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

   ##########################################
   # find new set of variable genes and redo pca and umap
   ##########################################
   Explore.umap.parameters.for.BWMcells = FALSE
   if(Explore.umap.parameters.for.BWMcells){

     source.my.script('scRNA_functions.R')
     require(tictoc)
     tic()
     test.umap.params.for.BWM.cells(sub.obj,
                                    pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                    group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                    nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                    n.neighbors.sampling = c(10, 30, 50),
                                    min.dist.sampling = c(0.01, 0.1)
     )
     toc()

   }

   nfeatures = 5000;
   sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
   #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
   sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
   sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
   ElbowPlot(sub.obj, ndims = 50)

   nb.pcs = 5 # nb of pcs depends on the considered clusters or ids
   n.neighbors = 10;
   min.dist = 0.01; spread = 1
   sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                      spread = spread, n.neighbors = n.neighbors,
                      min.dist = min.dist, verbose = TRUE)
   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()


   if(Refine.annotated.ids){

     DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
       NoLegend()

   }else{
     DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
       NoLegend()

     p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
       NoLegend()

     #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
     p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
       NoLegend() + ggtitle('scamp prediction')

     p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
       NoLegend() + ggtitle('seurat prediction')

     p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
       NoLegend()

     p0 + p2
     p0 + p1
     p1 + p2
     p2 + p3

   }

   ##########################################
   # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
   ##########################################
   FindClusters_subclusters = function(sub.obj, resolution = 1.0)
   {
     sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
     return(sub.obj$seurat_clusters)
   }
   sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:5, compute.SNN = TRUE)
   sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.3)
   DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

   p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
     ggtitle('manual.ids') + NoLegend()
   p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
                label.size = 8, na.value = "gray", combine = TRUE)
   p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

   (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration,
                                   '_MSxpaaap_MSxpapap.pdf'),
                            width = 18, height = 16)

   #jj.missing = which(is.na(sub.obj$manual.annot.ids))
   #table(sub.obj$seurat_clusters_split[jj.missing])

   VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
           group.by = 'seurat_clusters_split')
   #dev.off()
   ##########################################
   # check the counts of predicted ids for newly split clusters
   ##########################################
   idents.sel = c('0', '1', '2', '3', '4')
   Idents(sub.obj) = sub.obj$seurat_clusters_split
   counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
   counts.seurat = table(sub.obj$predicted.ids.seurat.keep, sub.obj$seurat_clusters_split)
   counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
   if(exists('idents.sel')){
     if(length(setdiff(colnames(counts), idents.sel)) > 0){
       counts = counts[, !is.na(match(colnames(counts), idents.sel))]
       counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
       counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
       counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
       counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
       counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
     }
   }

   # features of early stage
   features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                     'ceh-36', # MSxaap
                     'tbx-2', 'ceh-27') # MSxaaa
   # features of lineage MSxapa
   features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                     'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                     'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                     'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                     'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                     'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                     'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                     'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

   )
   # features of MSxaap lineage
   features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                     'sptf-1', 'ceh-22',  #MSaaapaa
                     'lin-12', 'pax-1', # MSxaapap
                     'ngn-1', 'ces-1', # MSpaapaa
                     'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                     'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                     'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                     'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
   )
   # features of MSxaaa
   features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                     'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                     #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                     'ceh-13', # MSaaaap not MSpaaap
                     'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                     'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                     'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                     'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                     'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                     'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                     'igcm-4', 'F54E2.2', 'K04G2.12',
                     'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                     'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                     'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                     )

   features.sels = c('unc-62','gana-1','clec-264',
                     'maph-1.3', 'zig-6', 'maph-1.2',
                     'mec-2', 'twk-31', 'stn-2', 'ZC449.5',  'shc-2', 'ham-1', 'ceh-34',
                     'K09G1.1', 'T25G12.11',
                     'ceh-13',  'tbx-2', 'D1086.12', 'sul-1',
                     'Y37E3.30', 'hnd-1', 'Y9C2UA.1'

                     #'gsnl-2', 'abts-1'
                     )
   FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

   # update the annotation with marker genes
   cluster.assingment = list(
     c('3', 'MSxpapap'),
     c('2', 'MSxpapap'),
     c('0', 'MSxpaaap'),
     c('4', 'MSxpaaap'),
     c('1', 'MSxpappp')
   )

   cat(length(cluster.assingment), 'clusters assigned \n')
   cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
   nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
   if(nb.unassigned >= 1){
     cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
   }

   # check info in JM data for specific lineage
   markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
   #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
   #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
   source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

   ids.sel = c('MSxpaaap'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
   #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

   ##########################################
   # update the manual annotation if good marker genes or mapped ids were found
   ##########################################
   # add.missing.annotation.for.pharynx.left.over = FALSE
   # if(add.missing.annotation.for.pharynx.left.over){
   #   jj = which(is.na(sub.obj$manual.annot.ids))
   #
   #   for(j in jj){
   #     #j = 1;
   #     cat('cell : ',  j, '\n')
   #     index.cluster = sub.obj$seurat_clusters_split[j]
   #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
   #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
   #   }
   #
   #   mm = match(colnames(sub.obj), colnames(seurat.obj))
   #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
   #
   # }

   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

   for(n in 1:length(cluster.assingment)){
     cluster.index = cluster.assingment[[n]][1];
     id2assign =  cluster.assingment[[n]][2];
     cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
     cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
     sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
     seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
     # if(is.na(id2assign)) {
     #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
     # }else{
     #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
     # }
   }

   #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
           pt.size = 2) + NoLegend()


   DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
     scale_colour_hue(drop = FALSE) + NoLegend()


   saveRDS(seurat.obj, file = RDS2save)



   ##########################################
  # Main aim:
  # MSxppapp and MSxpappp
  #
  # Notes:
  #
  #
  ##########################################
  GR.iteration = 9 # RG (revison global)
  Refine.annotated.ids = TRUE

  resDir = paste0("results/", version.analysis, '/annoted_pharynx')
  if(!dir.exists(resDir)){dir.create(resDir)}
  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}


  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                 '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
  RDS2save =  paste0(RdataDir,
                     'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                     '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

  seurat.obj = readRDS(file = RDSsaved)

  #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
  #    ' pharynx cells left to annotate \n')

  #cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
  #cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

  #pdf(pdfname, width=18, height = 10)
  #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  ##########################################
  # select subset of cells
  ##########################################

  # select cells with cluster index
  ##########################################
  jj = which(is.na(seurat.obj$manual.annot.ids))
  print(table(seurat.obj$seurat_clusters[jj]))

  #cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('6', '24', '20', '7')
  #cluster.sels = c('13', '17', '18', '11', '25', '10')
  # cluster.sels = c('1', '9', '12', '19')
  # #ids.sel = c('MSxaaa')
  # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
  #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
  cluster.sels = c('28', '52', '31')
  cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                   # & is.na(seurat.obj$manual.annot.ids)
                                    ]
  table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
          cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
          pt.size = 2, label.size = 5) + NoLegend()


  # select cells with ids
  ##########################################
  #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
  #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')
  clean.id.names = FALSE # clean id names
  if(clean.id.names){
    #jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")
    #seurat.obj$manual.annot.ids[jj] = 'mixture_transitionToTerminal'

  }


  ids.current = names(table(seurat.obj$manual.annot.ids))
  #ids.sels = ids.current

  ids.current[grep('_new', ids.current)]

  #ids.sels = c("mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap", "mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa",
  #             "mixture_transitionToTerminal")
  ids.sels = c(
              'MSxppapp', 'MSxpappp', "MSxpappp_new", "MSxpappp_new2", "MSxppapp_new",  "MSxppapp_new2"
              #'MSxpapap', 'MSxpaaap', "MSxpaaap_new","MSxpapap_new"
              )

  ids.left = setdiff(ids.current, ids.sels)
  # print(ids.left)
  # nchar(ids.left)

  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

  ##########################################
  # subset seurat object with selected cells
  ##########################################
  cat(length(cells.sels), ' cells selected to annotate \n')

  sub.obj = subset(seurat.obj, cells = cells.sels)

  #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  Explore.umap.parameters.for.BWMcells = FALSE
  if(Explore.umap.parameters.for.BWMcells){

    source.my.script('scRNA_functions.R')
    require(tictoc)
    tic()
    test.umap.params.for.BWM.cells(sub.obj,
                                   pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                   group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                   nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                   n.neighbors.sampling = c(10, 30, 50),
                                   min.dist.sampling = c(0.01, 0.1)
    )
    toc()

  }

  nfeatures = 5000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 5 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 10;
  min.dist = 0.01; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist, verbose = TRUE)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
    NoLegend()


  if(Refine.annotated.ids){

    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()

  }else{
    DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()

    p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend()

    #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
    p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend() + ggtitle('scamp prediction')

    p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
      NoLegend() + ggtitle('seurat prediction')

    p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend()

    p0 + p2
    p0 + p1
    p1 + p2
    p2 + p3

  }

  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 1.0)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:5, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
    ggtitle('manual.ids') + NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
               label.size = 8, na.value = "gray", combine = TRUE)
  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

  (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration,
                                  '_MSxppapp_MSxpappp.pdf'),
                           width = 18, height = 16)

  #jj.missing = which(is.na(sub.obj$manual.annot.ids))
  #table(sub.obj$seurat_clusters_split[jj.missing])

  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'seurat_clusters_split')
  #dev.off()
  ##########################################
  # check the counts of predicted ids for newly split clusters
  ##########################################
  #idents.sel = c('0', '1', '2', '3', '4')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
  counts.seurat = table(sub.obj$predicted.ids.seurat.keep, sub.obj$seurat_clusters_split)
  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
  if(exists('idents.sel')){
    if(length(setdiff(colnames(counts), idents.sel)) > 0){
      counts = counts[, !is.na(match(colnames(counts), idents.sel))]
      counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
      counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
      counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
      counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
      counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
    }
  }

  # features of early stage
  features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                    'ceh-36', # MSxaap
                    'tbx-2', 'ceh-27') # MSxaaa
  # features of lineage MSxapa
  features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                    'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                    'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                    'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                    'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                    'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                    'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                    'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

  )
  # features of MSxaap lineage
  features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                    'sptf-1', 'ceh-22',  #MSaaapaa
                    'lin-12', 'pax-1', # MSxaapap
                    'ngn-1', 'ces-1', # MSpaapaa
                    'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                    'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                    'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                    'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
  )
  # features of MSxaaa
  features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                    'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                    #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                    'ceh-13', # MSaaaap not MSpaaap
                    'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                    'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                    'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                    'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                    'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                    'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                    'igcm-4', 'F54E2.2', 'K04G2.12',
                    'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                    'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                    'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                    )

  features.sels = c('unc-62','gana-1','clec-264',
                    'maph-1.3', 'zig-6', 'maph-1.2',
                    'mec-2', 'twk-31', 'stn-2', 'ZC449.5',  'shc-2', 'ham-1', 'ceh-34',
                    'K09G1.1', 'T25G12.11',
                    'ceh-13',  'tbx-2', 'D1086.12', 'sul-1',
                    'Y37E3.30', 'hnd-1', 'Y9C2UA.1'

                    #'gsnl-2', 'abts-1'
                    )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  # update the annotation with marker genes
  cluster.assingment = list(
    c('3', 'MSxppapp'),
    c('8', 'MSxppapp'),
    c('5', 'MSxppapp'),
    c('6', 'MSxppapp'),
    c('1', 'MSxppapp'),

    c('7', 'MSxpappp'),
    c('0', 'MSxpappp'),
    c('2', 'MSxpappp'),
    c('4', 'MSxpappp'),
    c('9', 'MSxpappp')
  )

  cat(length(cluster.assingment), 'clusters assigned \n')
  cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
  nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
  if(nb.unassigned >= 1){
    cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
  }

  # check info in JM data for specific lineage
  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
  #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

  ids.sel = c('MSxpaaap'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  # add.missing.annotation.for.pharynx.left.over = FALSE
  # if(add.missing.annotation.for.pharynx.left.over){
  #   jj = which(is.na(sub.obj$manual.annot.ids))
  #
  #   for(j in jj){
  #     #j = 1;
  #     cat('cell : ',  j, '\n')
  #     index.cluster = sub.obj$seurat_clusters_split[j]
  #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
  #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
  #   }
  #
  #   mm = match(colnames(sub.obj), colnames(seurat.obj))
  #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
  #
  # }

  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    # if(is.na(id2assign)) {
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
    # }else{
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
    # }
  }

  #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
          pt.size = 2) + NoLegend()


  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE) + NoLegend()


  saveRDS(seurat.obj, file = RDS2save)


  GR.iteration = 10 # RG (revison global)
  Refine.annotated.ids = TRUE

  resDir = paste0("results/", version.analysis, '/annoted_pharynx')
  if(!dir.exists(resDir)){dir.create(resDir)}
  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}


  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                 '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
  RDS2save =  paste0(RdataDir,
                     'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                     '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

  seurat.obj = readRDS(file = RDSsaved)

  #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
  #    ' pharynx cells left to annotate \n')

  #cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
  #cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

  #pdf(pdfname, width=18, height = 10)
  #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  ##########################################
  # select subset of cells
  ##########################################

  # select cells with cluster index
  ##########################################
  jj = which(is.na(seurat.obj$manual.annot.ids))
  print(table(seurat.obj$seurat_clusters[jj]))

  #cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('6', '24', '20', '7')
  #cluster.sels = c('13', '17', '18', '11', '25', '10')
  # cluster.sels = c('1', '9', '12', '19')
  # #ids.sel = c('MSxaaa')
  # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
  #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
  cluster.sels = c('28', '52', '31')
  cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                   # & is.na(seurat.obj$manual.annot.ids)
                                    ]
  table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
          cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
          pt.size = 2, label.size = 5) + NoLegend()


  # select cells with ids
  ##########################################
  #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
  #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')
  clean.id.names = FALSE # clean id names
  if(clean.id.names){
    #jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")
    #seurat.obj$manual.annot.ids[jj] = 'mixture_transitionToTerminal'

  }


  ids.current = names(table(seurat.obj$manual.annot.ids))

  ids.current[grep('_new', ids.current)]

  #ids.sels = c("mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap", "mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa",
  #             "mixture_transitionToTerminal")
  ids.sels = c(
              "MSxpppaa_new", 'MSxpppaa','MSxpappa'
              #'MSxppapp', 'MSxpappp', "MSxpappp_new", "MSxpappp_new2", "MSxppapp_new",  "MSxppapp_new2"
              #'MSxpapap', 'MSxpaaap', "MSxpaaap_new","MSxpapap_new"

              )

  ids.left = setdiff(ids.current, ids.sels)
  # print(ids.left)
  # nchar(ids.left)

  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

  ##########################################
  # subset seurat object with selected cells
  ##########################################
  cat(length(cells.sels), ' cells selected to annotate \n')

  sub.obj = subset(seurat.obj, cells = cells.sels)

  #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  Explore.umap.parameters.for.BWMcells = FALSE
  if(Explore.umap.parameters.for.BWMcells){

    source.my.script('scRNA_functions.R')
    require(tictoc)
    tic()
    test.umap.params.for.BWM.cells(sub.obj,
                                   pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                   group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                   nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                   n.neighbors.sampling = c(10, 30, 50),
                                   min.dist.sampling = c(0.01, 0.1)
    )
    toc()

  }

  nfeatures = 5000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 10;
  min.dist = 0.01; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist, verbose = TRUE)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
    NoLegend()


  if(Refine.annotated.ids){

    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()

  }else{
    DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()

    p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend()

    #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
    p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend() + ggtitle('scamp prediction')

    p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
      NoLegend() + ggtitle('seurat prediction')

    p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend()

    p0 + p2
    p0 + p1
    p1 + p2
    p2 + p3

  }

  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 1.0)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:5, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
    ggtitle('manual.ids') + NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
               label.size = 8, na.value = "gray", combine = TRUE)
  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

  (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration,
                                  '_MSxpppaa_MSxpappa.pdf'),
                           width = 18, height = 16)


  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'seurat_clusters_split')
  #dev.off()
  ##########################################
  # check the counts of predicted ids for newly split clusters
  ##########################################
  #idents.sel = c('0', '1', '2', '3', '4')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
  counts.seurat = table(sub.obj$predicted.ids.seurat.keep, sub.obj$seurat_clusters_split)
  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
  if(exists('idents.sel')){
    if(length(setdiff(colnames(counts), idents.sel)) > 0){
      counts = counts[, !is.na(match(colnames(counts), idents.sel))]
      counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
      counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
      counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
      counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
      counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
    }
  }

  # features of early stage
  features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                    'ceh-36', # MSxaap
                    'tbx-2', 'ceh-27') # MSxaaa
  # features of lineage MSxapa
  features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                    'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                    'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                    'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                    'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                    'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                    'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                    'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

  )
  # features of MSxaap lineage
  features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                    'sptf-1', 'ceh-22',  #MSaaapaa
                    'lin-12', 'pax-1', # MSxaapap
                    'ngn-1', 'ces-1', # MSpaapaa
                    'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                    'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                    'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                    'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
  )
  # features of MSxaaa
  features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                    'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                    #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                    'ceh-13', # MSaaaap not MSpaaap
                    'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                    'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                    'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                    'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                    'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                    'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                    'igcm-4', 'F54E2.2', 'K04G2.12',
                    'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                    'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                    'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                    )

  features.sels = c('unc-62','gana-1','clec-264',
                    #'maph-1.3', 'zig-6', 'maph-1.2',
                    #'mec-2', 'twk-31',  'ZC449.5',  'shc-2', 'ham-1', 'ceh-34',
                    #'K09G1.1', 'T25G12.11',
                    'stn-2', 'igcm-4', 'E01G4.5', 'fbxc-24', 'fbxb-88', 'D1086.12', 'Y9C2UA.1', # MSxpappa
                    #'Y37E3.30', 'hnd-1', 'Y9C2UA.1'
                    'sul-1', 'tbx-2', 'hot-1', 'hmg-1.1', 'bgal-1', 'zig-7',

                    'gsnl-2', 'abts-1'
                    )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  # update the annotation with marker genes
  cluster.assingment = list(
    c('3', 'MSxpappa'),
    c('1', 'MSxpappa'),
    c('0', 'MSxpppaa'),
    c('2', 'MSxpppaa'),
    c('4', 'MSxpppaa')
  )

  cat(length(cluster.assingment), 'clusters assigned \n')
  cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
  nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
  if(nb.unassigned >= 1){
    cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
  }

  # check info in JM data for specific lineage
  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
  #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

  ids.sel = c('MSxpaaap'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  # add.missing.annotation.for.pharynx.left.over = FALSE
  # if(add.missing.annotation.for.pharynx.left.over){
  #   jj = which(is.na(sub.obj$manual.annot.ids))
  #
  #   for(j in jj){
  #     #j = 1;
  #     cat('cell : ',  j, '\n')
  #     index.cluster = sub.obj$seurat_clusters_split[j]
  #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
  #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
  #   }
  #
  #   mm = match(colnames(sub.obj), colnames(seurat.obj))
  #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
  #
  # }

  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    # if(is.na(id2assign)) {
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
    # }else{
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
    # }
  }

  #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
          pt.size = 2) + NoLegend()


  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE) + NoLegend()


  saveRDS(seurat.obj, file = RDS2save)



  ##########################################
  # Main aim:
  # terminal cells but start with something easily to resolve
  # MSxpppp and MSxpppa
  #
  # Notes:
  # here we annotate MSxppppx and MSxpppax
  #
  ##########################################
  GR.iteration = 11 # RG (revison global)
  Refine.annotated.ids = TRUE

  resDir = paste0("results/", version.analysis, '/annoted_pharynx')
  if(!dir.exists(resDir)){dir.create(resDir)}
  if(Refine.annotated.ids){by.group = 'manual.annot.ids';
  }else{by.group = 'seurat_clusters'}


  RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                 '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
  RDS2save =  paste0(RdataDir,
                     'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                     '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

  seurat.obj = readRDS(file = RDSsaved)

  #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
    scale_colour_hue(drop = FALSE) +
    NoLegend()

  #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
  #    ' pharynx cells left to annotate \n')

  #cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
  #cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

  #pdf(pdfname, width=18, height = 10)
  #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  ##########################################
  # select subset of cells
  ##########################################

  # select cells with cluster index
  ##########################################
  jj = which(is.na(seurat.obj$manual.annot.ids))
  print(table(seurat.obj$seurat_clusters[jj]))

  #cluster.sels = c('29', '32', '35', '40', '42')
  #cluster.sels = c('6', '24', '20', '7')
  #cluster.sels = c('13', '17', '18', '11', '25', '10')
  # cluster.sels = c('1', '9', '12', '19')
  # #ids.sel = c('MSxaaa')
  # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
  #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
  #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
  cluster.sels = c('28', '52', '31')
  cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                   # & is.na(seurat.obj$manual.annot.ids)
                                    ]
  table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
          cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
          pt.size = 2, label.size = 5) + NoLegend()


  # select cells with ids
  ##########################################
  #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
  #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')
  clean.id.names = FALSE # clean id names
  if(clean.id.names){
    #jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")
    #seurat.obj$manual.annot.ids[jj] = 'mixture_transitionToTerminal'

  }


  ids.current = names(table(seurat.obj$manual.annot.ids))

  ids.current[grep('_new', ids.current)]

  #ids.sels = c("mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap", "mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa",
  #             "mixture_transitionToTerminal")
  ids.sels = c(
              #"MSxpppaa_new", 'MSxpppaa','MSxpappa'
              'MSxpppp', 'MSxpppa'
              #'MSxppapp', 'MSxpappp', "MSxpappp_new", "MSxpappp_new2", "MSxppapp_new",  "MSxppapp_new2"
              #'MSxpapap', 'MSxpaaap', "MSxpaaap_new","MSxpapap_new"

              )

  ids.left = setdiff(ids.current, ids.sels)
  # print(ids.left)
  # nchar(ids.left)

  cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

  ##########################################
  # subset seurat object with selected cells
  ##########################################
  cat(length(cells.sels), ' cells selected to annotate \n')

  sub.obj = subset(seurat.obj, cells = cells.sels)

  #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
  sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

  DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

  ##########################################
  # find new set of variable genes and redo pca and umap
  ##########################################
  Explore.umap.parameters.for.BWMcells = FALSE
  if(Explore.umap.parameters.for.BWMcells){

    source.my.script('scRNA_functions.R')
    require(tictoc)
    tic()
    test.umap.params.for.BWM.cells(sub.obj,
                                   pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                   group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                   nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                   n.neighbors.sampling = c(10, 30, 50),
                                   min.dist.sampling = c(0.01, 0.1)
    )
    toc()

  }

  nfeatures = 3000;
  sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
  #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
  sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
  sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
  ElbowPlot(sub.obj, ndims = 50)

  nb.pcs = 10 # nb of pcs depends on the considered clusters or ids
  n.neighbors = 30;
  min.dist = 0.01; spread = 1
  sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                     spread = spread, n.neighbors = n.neighbors,
                     min.dist = min.dist, verbose = TRUE)
  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
    NoLegend()

  DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
    NoLegend()

  if(Refine.annotated.ids){

    DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()

  }else{
    DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
      NoLegend()

    p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend()

    #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
    p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend() + ggtitle('scamp prediction')

    p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
      NoLegend() + ggtitle('seurat prediction')

    p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
      NoLegend()

    p0 + p2
    p0 + p1
    p1 + p2
    p2 + p3

  }

  ##########################################
  # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
  ##########################################
  FindClusters_subclusters = function(sub.obj, resolution = 1.0)
  {
    sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
    return(sub.obj$seurat_clusters)
  }
  sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:5, compute.SNN = TRUE)
  sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 1.0)
  DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

  p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
    ggtitle('manual.ids') + NoLegend()
  p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
               label.size = 8, na.value = "gray", combine = TRUE)
  p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

  (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration,
                                  '_MSxpppp_MSxpppa.pdf'),
                           width = 18, height = 16)


  VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
          group.by = 'seurat_clusters_split')
  #dev.off()
  ##########################################
  # check the counts of predicted ids for newly split clusters
  ##########################################
  idents.sel = c('0', '1', '3', '4', '5', '2', '6')
  Idents(sub.obj) = sub.obj$seurat_clusters_split
  counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
  counts.seurat = table(sub.obj$predicted.ids.seurat.keep, sub.obj$seurat_clusters_split)
  counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
  if(exists('idents.sel')){
    if(length(setdiff(colnames(counts), idents.sel)) > 0){
      counts = counts[, !is.na(match(colnames(counts), idents.sel))]
      counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
      counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
      counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
      counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
      counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
    }
  }

  # features of early stage
  features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                    'ceh-36', # MSxaap
                    'tbx-2', 'ceh-27') # MSxaaa
  # features of lineage MSxapa
  features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                    'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                    'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                    'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                    'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                    'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                    'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                    'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

  )
  # features of MSxaap lineage
  features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                    'sptf-1', 'ceh-22',  #MSaaapaa
                    'lin-12', 'pax-1', # MSxaapap
                    'ngn-1', 'ces-1', # MSpaapaa
                    'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                    'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                    'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                    'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
  )
  # features of MSxaaa
  features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                    'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                    #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                    'ceh-13', # MSaaaap not MSpaaap
                    'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                    'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                    'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                    'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                    'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                    'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                    'igcm-4', 'F54E2.2', 'K04G2.12',
                    'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                    'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                    'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                    )

  features.sels = c('unc-62','gana-1','clec-264', 'lin-39',
                    #'maph-1.3', 'zig-6', 'maph-1.2',
                    #'mec-2', 'twk-31',  'ZC449.5',  'shc-2', 'ham-1', 'ceh-34',
                    #'K09G1.1', 'T25G12.11',
                    'stn-2', 'igcm-4', 'E01G4.5', 'fbxc-24', 'fbxb-88', 'D1086.12', 'Y9C2UA.1', # MSxpappa
                    #'Y37E3.30', 'hnd-1', 'Y9C2UA.1'
                    'sul-1', 'tbx-2', 'hot-1', 'hmg-1.1', 'bgal-1', 'zig-7',
                    'camt-1', 'ctg-1', 'irx-1',
                     'abts-1'
                    )
  FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

  # update the annotation with marker genes
  cluster.assingment = list(
    c('0', 'MSxpppa'),
    c('5', 'MSxpppax'),
    c('2', 'MSxpppax'),
    c('4', 'MSxpppp'),
    c('3', 'MSxpppp'),
    c('1', 'MSxppppx'),
    c('6', 'MSxppppx')

  )

  cat(length(cluster.assingment), 'clusters assigned \n')
  cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
  nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
  if(nb.unassigned >= 1){
    cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
  }

  # check info in JM data for specific lineage
  markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
  #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
  #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
  source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

  ids.sel = c('MSxpppax'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
  #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

  ##########################################
  # update the manual annotation if good marker genes or mapped ids were found
  ##########################################
  # add.missing.annotation.for.pharynx.left.over = FALSE
  # if(add.missing.annotation.for.pharynx.left.over){
  #   jj = which(is.na(sub.obj$manual.annot.ids))
  #
  #   for(j in jj){
  #     #j = 1;
  #     cat('cell : ',  j, '\n')
  #     index.cluster = sub.obj$seurat_clusters_split[j]
  #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
  #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
  #   }
  #
  #   mm = match(colnames(sub.obj), colnames(seurat.obj))
  #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
  #
  # }

  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
  # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

  for(n in 1:length(cluster.assingment)){
    cluster.index = cluster.assingment[[n]][1];
    id2assign =  cluster.assingment[[n]][2];
    cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
    cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
    sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
    seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
    # if(is.na(id2assign)) {
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
    # }else{
    #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
    # }
  }

  #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

  DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
          pt.size = 2) + NoLegend()


  DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
          na.value = "gray") +
    ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
    scale_colour_hue(drop = FALSE) + NoLegend()


  saveRDS(seurat.obj, file = RDS2save)



  ##########################################
 # Main aim:
 # terminal cells but start with something easily to resolve
 # MSxpppp and MSxpppa
 #
 # Notes:
 # for the convergence trajectory, the annotation is more satisfying.
 #
 ##########################################
 GR.iteration = 12 # RG (revison global)
 Refine.annotated.ids = TRUE

 resDir = paste0("results/", version.analysis, '/annoted_pharynx')
 if(!dir.exists(resDir)){dir.create(resDir)}
 if(Refine.annotated.ids){by.group = 'manual.annot.ids';
 }else{by.group = 'seurat_clusters'}


 RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
 RDS2save =  paste0(RdataDir,
                    'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                    '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

 seurat.obj = readRDS(file = RDSsaved)

 #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
   scale_colour_hue(drop = FALSE) +
   NoLegend()

 DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
   scale_colour_hue(drop = FALSE) +
   NoLegend()

 #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
 #    ' pharynx cells left to annotate \n')

 #cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
 #cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

 #pdf(pdfname, width=18, height = 10)
 #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
 ##########################################
 # select subset of cells
 ##########################################

 # select cells with cluster index
 ##########################################
 jj = which(is.na(seurat.obj$manual.annot.ids))
 print(table(seurat.obj$seurat_clusters[jj]))

 #cluster.sels = c('29', '32', '35', '40', '42')
 #cluster.sels = c('6', '24', '20', '7')
 #cluster.sels = c('13', '17', '18', '11', '25', '10')
 # cluster.sels = c('1', '9', '12', '19')
 # #ids.sel = c('MSxaaa')
 # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
 #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
 #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
 cluster.sels = c('28', '52', '31')
 cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                  # & is.na(seurat.obj$manual.annot.ids)
                                   ]
 table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
         cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
         pt.size = 2, label.size = 5) + NoLegend()


 # select cells with ids
 ##########################################
 #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
 #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')
 clean.id.names = FALSE # clean id names
 if(clean.id.names){
   #jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")
   #seurat.obj$manual.annot.ids[jj] = 'mixture_transitionToTerminal'

 }


 ids.current = names(table(seurat.obj$manual.annot.ids))

 ids.current[grep('_new', ids.current)]

 #ids.sels = c("mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap", "mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa",
 #             "mixture_transitionToTerminal")
 ids.sels = c(
             #"MSxpppaa_new", 'MSxpppaa','MSxpappa'
             #'MSxpppp', 'MSxpppa',
             "MSxapppax_new", 'MSxapppp', 'MSxappppx', 'MSxapppa', 'MSxpppax'
             #'MSxppapp', 'MSxpappp', "MSxpappp_new", "MSxpappp_new2", "MSxppapp_new",  "MSxppapp_new2"
             #'MSxpapap', 'MSxpaaap', "MSxpaaap_new","MSxpapap_new"

             )

 ids.left = setdiff(ids.current, ids.sels)
 # print(ids.left)
 # nchar(ids.left)

 cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

 ##########################################
 # subset seurat object with selected cells
 ##########################################
 cat(length(cells.sels), ' cells selected to annotate \n')

 sub.obj = subset(seurat.obj, cells = cells.sels)

 #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
 sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

 DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

 ##########################################
 # find new set of variable genes and redo pca and umap
 ##########################################
 Explore.umap.parameters.for.BWMcells = FALSE
 if(Explore.umap.parameters.for.BWMcells){

   source.my.script('scRNA_functions.R')
   require(tictoc)
   tic()
   test.umap.params.for.BWM.cells(sub.obj,
                                  pdfname = 'UMAP.param.TEST_whole_pharynx.pdf',
                                  group.by = 'predicted.ids.seurat', with_legend = FALSE,
                                  nfeatures.sampling = c(1000, 3000, 5000), nb.pcs.sampling = c(10, 20, 30),
                                  n.neighbors.sampling = c(10, 30, 50),
                                  min.dist.sampling = c(0.01, 0.1)
   )
   toc()

 }

 nfeatures = 5000;
 sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
 #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
 sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
 sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
 ElbowPlot(sub.obj, ndims = 50)

 nb.pcs = 5 # nb of pcs depends on the considered clusters or ids
 n.neighbors = 20;
 min.dist = 0.01; spread = 1
 sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                    spread = spread, n.neighbors = n.neighbors,
                    min.dist = min.dist, verbose = TRUE)
 DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
   NoLegend()

 #DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
 #  NoLegend()

 if(Refine.annotated.ids){

   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()

 }else{
   DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()

   p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()

   #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
   p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend() + ggtitle('scamp prediction')

   p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
     NoLegend() + ggtitle('seurat prediction')

   p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
     NoLegend()

   p0 + p2
   p0 + p1
   p1 + p2
   p2 + p3

 }

 ##########################################
 # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
 ##########################################
 FindClusters_subclusters = function(sub.obj, resolution = 1.0)
 {
   sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
   return(sub.obj$seurat_clusters)
 }
 sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:5, compute.SNN = TRUE)
 sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
 DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

 p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
   ggtitle('manual.ids') + NoLegend()
 p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
              label.size = 8, na.value = "gray", combine = TRUE)
 p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

 (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration,
                                 '_MSxapppxx.pdf'),
                          width = 18, height = 16)


 VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
         group.by = 'seurat_clusters_split')
 #dev.off()
 ##########################################
 # check the counts of predicted ids for newly split clusters
 ##########################################
 idents.sel = c('0', '1', '3', '4', '5', '2', '6')
 Idents(sub.obj) = sub.obj$seurat_clusters_split
 counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
 counts.seurat = table(sub.obj$predicted.ids.seurat.keep, sub.obj$seurat_clusters_split)
 counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
 if(exists('idents.sel')){
   if(length(setdiff(colnames(counts), idents.sel)) > 0){
     counts = counts[, !is.na(match(colnames(counts), idents.sel))]
     counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
     counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
     counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
     counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
     counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
   }
 }

 # features of early stage
 features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                   'ceh-36', # MSxaap
                   'tbx-2', 'ceh-27') # MSxaaa
 # features of lineage MSxapa
 features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                   'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                   'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                   'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                   'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                   'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                   'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                   'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

 )
 # features of MSxaap lineage
 features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                   'sptf-1', 'ceh-22',  #MSaaapaa
                   'lin-12', 'pax-1', # MSxaapap
                   'ngn-1', 'ces-1', # MSpaapaa
                   'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                   'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                   'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                   'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
 )
 # features of MSxaaa
 features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                   'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                   #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                   'ceh-13', # MSaaaap not MSpaaap
                   'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                   'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                   'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                   'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                   'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                   'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                   'igcm-4', 'F54E2.2', 'K04G2.12',
                   'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                   'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                   'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                   )

 features.sels = c('unc-62','gana-1','clec-264', 'lin-39',
                   #'maph-1.3', 'zig-6', 'maph-1.2',
                   #'mec-2', 'twk-31',  'ZC449.5',  'shc-2', 'ham-1', 'ceh-34',
                   #'K09G1.1', 'T25G12.11',
                   'stn-2', 'igcm-4', 'E01G4.5', 'fbxc-24', 'fbxb-88', 'D1086.12', 'Y9C2UA.1', # MSxpappa
                   #'Y37E3.30', 'hnd-1', 'Y9C2UA.1'
                   'sul-1', 'tbx-2', 'hot-1', 'hmg-1.1', 'bgal-1', 'zig-7',
                   #'camt-1', 'ctg-1', 'irx-1',
                    'abts-1'
                   )
 FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

 # update the annotation with marker genes
 cluster.assingment = list(
   c('4', 'MSxapppp'),
   c('1', 'MSxapppp'),
   c('6', 'MSxappppx'),
   c('3', 'MSxappppx'),

   c('2', 'MSxapppa'),
   c('5', 'MSxapppa'),
   c('0', 'MSxapppax')

 )

 cat(length(cluster.assingment), 'clusters assigned \n')
 cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
 nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
 if(nb.unassigned >= 1){
   cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
 }

 # check info in JM data for specific lineage
 markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
 #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
 #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
 source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

 ids.sel = c('MSxpppax'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
 #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

 ##########################################
 # update the manual annotation if good marker genes or mapped ids were found
 ##########################################
 # add.missing.annotation.for.pharynx.left.over = FALSE
 # if(add.missing.annotation.for.pharynx.left.over){
 #   jj = which(is.na(sub.obj$manual.annot.ids))
 #
 #   for(j in jj){
 #     #j = 1;
 #     cat('cell : ',  j, '\n')
 #     index.cluster = sub.obj$seurat_clusters_split[j]
 #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
 #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
 #   }
 #
 #   mm = match(colnames(sub.obj), colnames(seurat.obj))
 #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
 #
 # }

 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
 # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

 for(n in 1:length(cluster.assingment)){
   cluster.index = cluster.assingment[[n]][1];
   id2assign =  cluster.assingment[[n]][2];
   cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
   cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
   sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
   seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
   # if(is.na(id2assign)) {
   #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
   # }else{
   #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
   # }
 }

 #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

 DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
         pt.size = 2) + NoLegend()


 DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
         na.value = "gray") +
   ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
   scale_colour_hue(drop = FALSE) + NoLegend()


 saveRDS(seurat.obj, file = RDS2save)



 ##########################################
   # Main aim:
   # here try to revise MSxppppp and MSxpppap
   #
   #
   # Notes:
   # it is not very clear how to resolve those terminal cells.
   # But I will stop here
   ##########################################
   GR.iteration = 13 # RG (revison global)
   Refine.annotated.ids = TRUE

   resDir = paste0("results/", version.analysis, '/annoted_pharynx')
   if(!dir.exists(resDir)){dir.create(resDir)}
   if(Refine.annotated.ids){by.group = 'manual.annot.ids';
   }else{by.group = 'seurat_clusters'}


   RDSsaved = paste0(RdataDir, 'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                                  '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration -1, '.rds')
   RDS2save =  paste0(RdataDir,
                      'processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot',
                      '_cleanedBWM_and_Pharynx_iteration_GR', GR.iteration, '.rds')

   seurat.obj = readRDS(file = RDSsaved)

   #pdfname = paste0(resDir, "/Manual_cluster_annotation_BDW_test_MSxp_lineage_iteration_", nb.iteration, ".pdf")
   DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
     scale_colour_hue(drop = FALSE) +
     NoLegend()

   DimPlot(seurat.obj, group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10_manual.annoted.IDs_BWM_Pharynx")) +
     scale_colour_hue(drop = FALSE) +
     NoLegend()

   #cat(length(which(is.na(seurat.obj$manual.annot.ids) & !is.na(seurat.obj$seurat_clusters_pharynx))),
   #    ' pharynx cells left to annotate \n')

   #cat(length(which(is.na(seurat.obj$manual.annot.ids))), ' total cells not annotated \n')
   #cat('~', length(which(is.na(seurat.obj$manual.annot.ids))) - 272, 'non-early-stage cells missing annotation \n')

   #pdf(pdfname, width=18, height = 10)
   #par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
   ##########################################
   # select subset of cells
   ##########################################

   # select cells with cluster index
   ##########################################
   jj = which(is.na(seurat.obj$manual.annot.ids))
   print(table(seurat.obj$seurat_clusters[jj]))

   #cluster.sels = c('29', '32', '35', '40', '42')
   #cluster.sels = c('6', '24', '20', '7')
   #cluster.sels = c('13', '17', '18', '11', '25', '10')
   # cluster.sels = c('1', '9', '12', '19')
   # #ids.sel = c('MSxaaa')
   # #ids.excl = c('MSxapp', 'MSxppa', 'MSpappa')
   #cluster.sels = c('3', '4', '8', '22', '0', '21', '23')
   #cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters_pharynx, cluster.sels))]
   cluster.sels = c('28', '52', '31')
   cells.sels = colnames(seurat.obj)[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))
                                    # & is.na(seurat.obj$manual.annot.ids)
                                     ]
   table(seurat.obj$manual.annot.ids[!is.na(match(seurat.obj$seurat_clusters, cluster.sels))], useNA = 'ifany')
   DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE,
           cells.highlight = cells.sels, cols.highlight = 'red', sizes.highlight = 1,
           pt.size = 2, label.size = 5) + NoLegend()


   # select cells with ids
   ##########################################
   #index.pharynx.cells = which(!is.na(seurat.obj$seurat_clusters_pharynx))
   #table(seurat.obj$manual.annot.ids[index.pharynx.cells], useNA = 'ifany')
   clean.id.names = FALSE # clean id names
   if(clean.id.names){
     #jj = which(seurat.obj$manual.annot.ids == "MSxppap/MSxpaaa/MSxpaaap")
     #seurat.obj$manual.annot.ids[jj] = 'mixture_transitionToTerminal'

   }


   ids.current = names(table(seurat.obj$manual.annot.ids))
   #ids.current[grep('_new', ids.current)]
   #ids.sels = c("mixture_MSxpaaap.MSxppapp.MSxpappp.MSxpapap", "mixture_MSxppppp.MSxppppa.MSxpppap.MSxpppaa.MSxpappa",
   #             "mixture_transitionToTerminal")
   ids.sels = c(
               #"MSxpppaa_new", 'MSxpppaa','MSxpappa'
               #'MSxpppp', 'MSxpppa',
               #"MSxapppax_new", 'MSxapppp', 'MSxappppx', 'MSxapppa', 'MSxpppax'
               #'MSxppapp', 'MSxpappp', "MSxpappp_new", "MSxpappp_new2", "MSxppapp_new",  "MSxppapp_new2"
               #'MSxpapap', 'MSxpaaap', "MSxpaaap_new","MSxpapap_new"
                'MSxppppp', 'MSxpppap'
                #'MSxppppa',  'MSxpppaa', 'MSxpappa'
               )

   ids.left = setdiff(ids.current, ids.sels)
   # print(ids.left)
   # nchar(ids.left)

   cells.sels = unique(colnames(seurat.obj)[!is.na(match(seurat.obj$manual.annot.ids, ids.sels))])

   ##########################################
   # subset seurat object with selected cells
   ##########################################
   cat(length(cells.sels), ' cells selected to annotate \n')

   sub.obj = subset(seurat.obj, cells = cells.sels)

   #sub.obj$seurat_clusters = as.integer(as.character(sub.obj$seurat_clusters_pharynx))
   sub.obj$timingEst = as.numeric(as.character(sub.obj$timingEst))

   DimPlot(sub.obj, reduction = 'umap', label = TRUE, group.by = by.group) + NoLegend()

   ##########################################
   # find new set of variable genes and redo pca and umap
   ##########################################
   Explore.umap.parameters.for.BWMcells = FALSE
   if(Explore.umap.parameters.for.BWMcells){

     source.my.script('scRNA_functions.R')
     require(tictoc)
     tic()
     test.umap.params.for.BWM.cells(sub.obj,
                                    pdfname = 'UMAP.param.TEST_BWM_MSxpppxx.pdf',
                                    group.by = 'manual.annot.ids', with_legend = FALSE,
                                    nfeatures.sampling = c(1000, 3000, 5000, 8000), nb.pcs.sampling = c(5, 10, 30),
                                    n.neighbors.sampling = c(10, 30),
                                    min.dist.sampling = c(0.01)
     )
     toc()

   }

   nfeatures = 5000;
   sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures)
   #cat('nb of variableFeatures excluding timer genes : ', length(VariableFeatures(sub.obj)), '\n')
   sub.obj = ScaleData(sub.obj, features = rownames(sub.obj))
   sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, weight.by.var = FALSE)
   ElbowPlot(sub.obj, ndims = 50)

   nb.pcs = 5 # nb of pcs depends on the considered clusters or ids
   n.neighbors = 10;
   min.dist = 0.01; spread = 1
   sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = c(1:nb.pcs),
                      spread = spread, n.neighbors = n.neighbors,
                      min.dist = min.dist, verbose = TRUE)
   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
     NoLegend()

   #DimPlot(sub.obj, group.by = 'seurat_clusters_split', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
   #  NoLegend()

   if(Refine.annotated.ids){

     DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
       NoLegend()

   }else{
     DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 6, pt.size = 2.0, repel = TRUE) +
       NoLegend()

     p0 = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
       NoLegend()

     #p1 = DimPlot(sub.obj, group.by = 'pred.ids.seurat.keep.bwm.all', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)
     p1 = DimPlot(sub.obj, group.by = 'predicted.ids.scmap', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
       NoLegend() + ggtitle('scamp prediction')

     p2 = DimPlot(sub.obj, group.by = 'predicted.ids.seurat', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE)  +
       NoLegend() + ggtitle('seurat prediction')

     p3 = DimPlot(sub.obj, group.by = 'seurat_clusters', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE) +
       NoLegend()

     p0 + p2
     p0 + p1
     p1 + p2
     p2 + p3

   }

   ##########################################
   # redo the clustering using seurat FindCluster (SLM alogrithm) after testing k-mean from RaceID
   ##########################################
   FindClusters_subclusters = function(sub.obj, resolution = 1.0)
   {
     sub.obj <- FindClusters(sub.obj, resolution = resolution, algorithm = 3)
     return(sub.obj$seurat_clusters)
   }
   sub.obj <- FindNeighbors(object = sub.obj, reduction = "pca", k.param = 10, dims = 1:10, compute.SNN = TRUE)
   sub.obj$seurat_clusters_split = FindClusters_subclusters(sub.obj, resolution = 0.5)
   DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5)

   p1  = DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 8, repel = TRUE,  pt.size = 3) +
     ggtitle('manual.ids') + NoLegend()
   p2 = DimPlot(sub.obj, group.by = "seurat_clusters_split", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 3,
                label.size = 8, na.value = "gray", combine = TRUE)
   p4 = VlnPlot(sub.obj, features = c('timingEst'), ncol = 1, group.by = 'seurat_clusters_split') + NoLegend()

   (p2 + p4) / p1  + ggsave(paste0(resDir, '/splitcluster_timingEstimation_manual.IDs_iteration_GR', GR.iteration,
                                   '_MSxppppp_MSxpppap.pdf'),
                            width = 18, height = 16)


   VlnPlot(sub.obj, features = c("FSC_log2", "BSC_log2"), ncol = 2,
           group.by = 'seurat_clusters_split')
   #dev.off()
   ##########################################
   # check the counts of predicted ids for newly split clusters
   ##########################################
   idents.sel = c('0', '1', '3', '4', '5', '2', '6')
   Idents(sub.obj) = sub.obj$seurat_clusters_split
   counts = table(sub.obj$predicted.ids.scmap, sub.obj$seurat_clusters_split)
   counts.seurat = table(sub.obj$predicted.ids.seurat.keep, sub.obj$seurat_clusters_split)
   counts.annot = table(sub.obj$manual.annot.ids, sub.obj$seurat_clusters_split)
   if(exists('idents.sel')){
     if(length(setdiff(colnames(counts), idents.sel)) > 0){
       counts = counts[, !is.na(match(colnames(counts), idents.sel))]
       counts = counts[apply(as.matrix(counts), 1, sum) >0, ]
       counts.seurat = counts.seurat[, !is.na(match(colnames(counts.seurat), idents.sel))]
       counts.seurat = counts.seurat[apply(as.matrix(counts.seurat), 1, sum)>0, ]
       counts.annot = counts.annot[, !is.na(match(colnames(counts.annot), idents.sel))]
       counts.annot = counts.annot[apply(as.matrix(counts.annot), 1, sum) >0, ]
     }
   }

   # features of early stage
   features.sels = c('pha-4', 'hnd-1', 'cft-1', 'alr-1', 'irx-1', # MSxaa
                     'ceh-36', # MSxaap
                     'tbx-2', 'ceh-27') # MSxaaa
   # features of lineage MSxapa
   features.sels = c('dmd-4', 'cnd-1', 'swt-3', 'ceh-16', # MSxapa
                     'cnd-1', 'asp-4', 'dmd-4', 'cgef-1', # MSxapap
                     'dmd-4', 'K10D3.6', 'ZK829.9', 'srd-32', # MSxapaa
                     'srd-32', 'gst-20', 'sdz-21', 'ceh-13', 'hlh-6', 'F54E2.2', #MSxapaap
                     'ceh-22', 'ZK829.9', 'srd-32', 'ceh-13', 'clec-258', 'C07C7.1', # MSxapaaa
                     'ceh-22', 'spp-7', 'tnc-2', 'F54E2.2', 'K04G2.12', # MSaaaappp/MSxapaapp
                     'nfki-1', 'unc-62', 'ser-2', 'tnc-2', 'irx-1', 'fem-1', # MSpaaappp/MSxapappa
                     'hlh-6', 'ces-1', 'Y51H7C.10', 'Y62F5A.9', # MSxapaapa

   )
   # features of MSxaap lineage
   features.sels = c('lin-12', 'ceh-36', 'tbx-2', 'tbx-7', 'ceh-34', 'ceh-32', 'clec-258', # MSxaapa
                     'sptf-1', 'ceh-22',  #MSaaapaa
                     'lin-12', 'pax-1', # MSxaapap
                     'ngn-1', 'ces-1', # MSpaapaa
                     'fos-1', 'ceh-36', 'ceh-34', 'irx-1',  # MSaaapp
                     'jun-1', 'ttx-1', 'irx-1', # MSaaappp
                     'pax-1', 'lin-12', 'ttx-1', 'agr-1', 'irx-1', 'cwn-2', # MSxaapap and MSxaapapa
                     'ref-1', 'aff-1', 'K04G2.12', 'unc-129' # MSaaapapp
   )
   # features of MSxaaa
   features.sels = c('tbx-2', 'ceh-27', # MSxaaa
                     'ceh-32', # MSaaaaa and MSpaaaa not MSxaaa and the other two daughters
                     #'let-381', 'fem-1', 'hphd-1', 'F26B1.1', # only MSaaaa
                     'ceh-13', # MSaaaap not MSpaaap
                     'ceh-27', 'spi-1', 'fem-1', 'unc-62', 'tbx-2', #MSpaaapp
                     'lim-7', 'dod-6', 'ces-1', # MSpaaapa
                     'let-381', 'F26B1.1', 'unc-30', 'ceh-32', 'fem-1', 'ceh-34', 'igcm-4', 'hphd-1', # MSaaaaa and daughers
                     'ceh-27', 'hlh-3', 'spi-1', 'K04G2.12', 'F54E2.2', 'igcm-4', 'fem-1', # MSaaaapa and MSaaaaap
                     'mnp-1', 'ceh-13', 'fkh-2', 'shc-2', # MSaaaaap
                     'dod-6', 'ces-1', 'ceh-53', # MSaaaaapa (probably missing)
                     'igcm-4', 'F54E2.2', 'K04G2.12',
                     'ceh-32', 'ceh-34', 'ceh-27', 'fem-1', 'spi-1', 'C45G7.4', # MSpaaaaa
                     'dod-6', 'ceh-53', 'hlh-3', 'ceh-34', 'fkh-2', 'hlh-6', # MSaaaaapa
                     'let-381', 'F26B1.1', 'unc-30' # MSaaaaaax/MSxpaaaax
                     )

   features.sels = c('unc-62','gana-1','clec-264', 'lin-39',
                     #'maph-1.3', 'zig-6', 'maph-1.2',
                     #'mec-2', 'twk-31',  'ZC449.5',  'shc-2', 'ham-1', 'ceh-34',
                     #'K09G1.1', 'T25G12.11',
                     'stn-2', 'igcm-4', 'E01G4.5', 'fbxc-24', 'fbxb-88', 'D1086.12', 'Y9C2UA.1', # MSxpappa
                     #'Y37E3.30', 'hnd-1', 'Y9C2UA.1'
                     'sul-1', 'tbx-2', 'hot-1', 'hmg-1.1', 'bgal-1', 'zig-7',
                     #'camt-1', 'ctg-1', 'irx-1',
                      'abts-1',
                     'F37H8.5', 'B0379.1', 'hmg-1.1', 'ZK512.1'
                     )
   FeaturePlot(sub.obj, reduction = 'umap', features = features.sels)

   # update the annotation with marker genes
   cluster.assingment = list(
     c('1', 'MSxppppp'),
     c('2', 'MSxppppp'),
     c('0', 'MSxpppap')

   )

   cat(length(cluster.assingment), 'clusters assigned \n')
   cat(length(levels(sub.obj$seurat_clusters_split)), ' split clusters \n')
   nb.unassigned = length(levels(sub.obj$seurat_clusters_split)) - length(cluster.assingment)
   if(nb.unassigned >= 1){
     cat(' Error  : ', nb.unassigned, 'clusters unassigned \n')
   }

   # check info in JM data for specific lineage
   markers.JM = read.xlsx('data/Supplementary_Tables_190611.xlsx', sheet=  4, startRow = 8, colNames = TRUE)
   #markers.JM = readRDS(file = paste0(RdataDir, 'BWM_markerGenes_JM.rds'))
   #load(file = paste0(RdataDir, 'Seurat.object_JM_BWM_data_markers.Rdata'))
   source.my.script('scRNA_cluster_annotation_utilityFunctions.R')

   ids.sel = c('MSxpppax'); find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers.JM)
   #find.markerGenes.used.in.JM.scRNAseq(ids = ids.sel, markers = markers)

   ##########################################
   # update the manual annotation if good marker genes or mapped ids were found
   ##########################################
   # add.missing.annotation.for.pharynx.left.over = FALSE
   # if(add.missing.annotation.for.pharynx.left.over){
   #   jj = which(is.na(sub.obj$manual.annot.ids))
   #
   #   for(j in jj){
   #     #j = 1;
   #     cat('cell : ',  j, '\n')
   #     index.cluster = sub.obj$seurat_clusters_split[j]
   #     ids.clusters = table(sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == index.cluster)])
   #     sub.obj$manual.annot.ids[j] = names(ids.clusters)[which.max(ids.clusters)]
   #   }
   #
   #   mm = match(colnames(sub.obj), colnames(seurat.obj))
   #   seurat.obj$manual.annot.ids[mm] = sub.obj$manual.annot.ids
   #
   # }

   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'mixture_MSxppapp_MSxpappp')] = 'MSxppapp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxpappp')] = 'MSxppapp/MSxpappp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids == 'likely_MSxppapp')] = 'MSxppapp/MSxpappp'
   # sub.obj$manual.annot.ids[which(sub.obj$manual.annot.ids.5 != 'MSxppapp')] = 'MSxppapp/MSxpappp'

   for(n in 1:length(cluster.assingment)){
     cluster.index = cluster.assingment[[n]][1];
     id2assign =  cluster.assingment[[n]][2];
     cat('cluster ', cluster.index, 'assinged to ', id2assign, '\n')
     cells = colnames(sub.obj)[which(sub.obj$seurat_clusters_split == cluster.index)]
     sub.obj$manual.annot.ids[which(sub.obj$seurat_clusters_split == cluster.index)] = id2assign
     seurat.obj$manual.annot.ids[match(cells, colnames(seurat.obj))] = id2assign
     # if(is.na(id2assign)) {
     #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = NA
     # }else{
     #   seurat.obj$BWM.cells[match(cells, colnames(seurat.obj))] = 'BWM'
     # }
   }

   #seurat.obj$manual.annot.ids[match(colnames(sub.obj), colnames(seurat.obj))] = sub.obj$manual.annot.ids

   DimPlot(sub.obj, group.by = 'manual.annot.ids', reduction = 'umap', label = TRUE, label.size = 5, repel = TRUE,
           pt.size = 2) + NoLegend()


   DimPlot(seurat.obj, group.by = "manual.annot.ids", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 5,
           na.value = "gray") +
     ggtitle(paste0("Seurat_clustering_SLM_resolution3_3000variableFeatures_20pca_k10")) +
     scale_colour_hue(drop = FALSE) + NoLegend()


   saveRDS(seurat.obj, file = RDS2save)
   
