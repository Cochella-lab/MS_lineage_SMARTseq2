##########################################################################
##########################################################################
# Project: Aleks' scRNA-seq MS project
# Script purpose: functions for clustering  
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Nov 22 15:38:57 2019
##########################################################################
##########################################################################
test.umap.params = function(seurat.obj, pdfname = paste0(resDir, '/umap_params_test.pdf'), 
                            nb.pcs.sampling = seq(20, 50, by = 10),
                            n.neighbors.sampling = seq(10, 50, by = 10),
                            min.dist.sampling = c(0.01, 0.05, seq(0.1, 0.5, by = 0.1))
)
{
  #pdfname = paste0(resDir, "/cluster_annotations/test_umap_params_for_seurat.obj.pdf")
  pdf(pdfname, width=12, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #nb.pcs = ncol(seurat.obj[['pca']]);
  #n.neighbors = 50; min.dist = 0.01;
  
  for(nb.pcs in unique(c(nb.pcs.sampling, ncol(seurat.obj[['pca']]))))
  {
    if(nb.pcs <= ncol(seurat.obj[['pca']])){
      for(n.neighbors in n.neighbors.sampling)
      {
        for(min.dist in min.dist.sampling)
        {
          cat('nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, ', min.dist - ', min.dist, '\n')
          seurat.obj <- RunUMAP(object = seurat.obj, reduction = 'pca', 
                                dims = 1:nb.pcs, 
                                n.neighbors = n.neighbors, min.dist = min.dist)
          p1 =  DimPlot(seurat.obj, reduction = 'umap', label = TRUE, pt.size = 1, label.size = 6) + 
            NoLegend() + 
            ggtitle(paste0('nb.pcs - ', nb.pcs, '; n.neighbors - ', n.neighbors, ', min.dist - ', min.dist))
          plot(p1)
        }
      }
    }
  }
  
  dev.off()
  
}

########################################################
########################################################
# Section : label transferring from Murray dataset
# 1) seurat
# 2) scmap, svm, RF
# 3) utility functions
########################################################
########################################################
seurat.transfer.labels.from.Murray.scRNA.to.scRNA = function(seurat.obj)
{
  # seurat.obj = ms;
  # process and import Parker et al. data
  process.import.Murray.scRNA()
  ee <- FindVariableFeatures(
    object = ee,
    nfeatures = 3000
  )
  
  ee <- ScaleData(object = ee)
  Idents(ee) = ee$lineage
  
  # Here, we process the gene activity matrix 
  # in order to find anchors between cells in the scATAC-seq dataset 
  # and the scRNA-seq dataset.
  Idents(seurat.obj) = seurat.obj$SCT_snn_res.12
  
  DimPlot(seurat.obj, reduction = "umap", label = TRUE, pt.size = 2,  label.size = 5, repel = FALSE) + NoLegend()
  
  DefaultAssay(seurat.obj) <- 'RNA'
  #nb.variableFeatures = 5000
  #seurat.obj <- FindVariableFeatures(seurat.obj, nfeatures = nb.variableFeatures)
  #seurat.obj <- NormalizeData(seurat.obj)
  #seurat.obj <- ScaleData(seurat.obj)
  
  #seurat.obj <- RunPCA(seurat.obj, npcs = 50, verbose = FALSE, reduction.name = 'pca.ga')
  
  transfer.anchors <- FindTransferAnchors(
    reference = ee,
    query = seurat.obj,
    features = unique(VariableFeatures(seurat.obj)),
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
  
  cat('nb of cells in query and in reference as anchors : ', 
      length(unique(transfer.anchors@anchors[, 1])), '--',  length(unique(transfer.anchors@anchors[, 2])), '\n')
  
  
  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    #refdata = Idents(tintori),
    refdata = Idents(ee),
    #refdata = as.vector(Idents(seurat.obj)),
    #weight.reduction = seurat.obj[['pca']],
    weight.reduction = seurat.obj[['pca']],
    dims = 1:30,
    k.weight = 50
  )
  
  seurat.obj <- AddMetaData(object = seurat.obj, metadata = predicted.labels)
  #nb.pcs = 50; n.neighbors = 30; min.dist = 0.4;
  #seurat.obj <- RunUMAP(object = seurat.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
  #                      n.neighbors = n.neighbors, min.dist = min.dist)
  
  DimPlot(seurat.obj, group.by = "predicted.id", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5) + 
    ggtitle("transferred labels") +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  table(seurat.obj$prediction.score.max > 0.5)
  
  hist(seurat.obj$prediction.score.max)
  abline(v = 0.5, col = "red")
  
  seurat.obj.filtered <- subset(seurat.obj, subset = prediction.score.max > 0.5)
  
  # to make the colors match
  seurat.obj.filtered$predicted.id <- factor(seurat.obj.filtered$predicted.id, levels = levels(ee))  
  DimPlot(seurat.obj.filtered, group.by = "predicted.id", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5) + 
  ggtitle("transferred labels with probability threshod = 0.5") +
  scale_colour_hue(drop = FALSE) + 
  NoLegend()
  
  DimPlot(ms_correlation_idents, reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5) + 
    ggtitle("aleks correlation assignment") +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  ##########################################
  # make summary of cluster-to-predicted label mapping
  ##########################################
  res = data.frame(clusters = Idents(seurat.obj), predicted.labels, stringsAsFactors = FALSE)
  
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

reference.based.cluster.annotation = function(seurat.obj, method = 'scmap', 
                                              nb.features.scmap = 500, threshold.scmap = 0.7,
                                              threshold.svm = 0.5, threshold.rf = 0.5)
{
  # seurat.obj = ms;
  # method = 'scmap';  nb.features.scmap = 500; threshold.scmap = 0.7; nb.features.svm = 1500; threshold.svm = 0.5;threshold.rf = 0.5;
  
  ##########################################
  # first step: project aleks cells to the reference using scmap
  # transfer labels with stringent threshold
  ##########################################
  ## process aleks data for scmap
  library(SingleCellExperiment)
  library(scmap)
  sce = Seurat::as.SingleCellExperiment(seurat.obj)
  sce <- sce[!duplicated(rownames(sce)), ]
  rowData(sce)$feature_symbol <- rownames(sce)
  counts(sce) = as.matrix(counts(sce)) # sce object converted from seurat object was using spare matrix
  logcounts(sce) = as.matrix(logcounts(sce))
  
  ## import and process Murray datfa for scmap as reference
  ee = process.import.Murray.scRNA()
  ee = Seurat::as.SingleCellExperiment(ee)
  counts(ee) = as.matrix(counts(ee))
  logcounts(ee) = as.matrix(logcounts(ee))
  rowData(ee)$feature_symbol <- rownames(ee)
  ee$cell_type1 = ee$lineage
  
  ## feature selection for scmap
  ee <- selectFeatures(ee, suppress_plot = FALSE, n_features = nb.features.scmap)
  table(rowData(ee)$scmap_features)
  #as.character(unique(ee$cell_type1))
  
  if(method == 'scmap'){
    ee_ref = indexCluster(ee)
    
    #head(metadata(ee_ref)$scmap_cluster_index)
    #heatmap(as.matrix(metadata(ee_ref)$scmap_cluster_index))
    
    scmapCluster_results <- scmapCluster(
      projection = sce, 
      index_list = list(
        murray = metadata(ee_ref)$scmap_cluster_index
      ),
      threshold = threshold.scmap
    )
    
    length(scmapCluster_results$scmap_cluster_labs)
    length(scmapCluster_results$combined_labs)
    ident.murray = unique(ee$lineage)
    ident.projection = unique(scmapCluster_results$scmap_cluster_labs)
    ident.missed = ident.murray[which(is.na(match(ident.murray, ident.projection)))]
    print(ident.missed)
    #head(scmapCluster_results$scmap_cluster_labs)
    #head(scmapCluster_results$scmap_cluster_siml)
    
    hist(scmapCluster_results$scmap_cluster_siml, breaks = 100)
    abline(v = threshold.scmap, col = 'red')
    head(scmapCluster_results$combined_labs)
    
    predicted.id = scmapCluster_results$scmap_cluster_labs
    predicted.id[which(predicted.id == 'unassigned')] = NA
    
    cat('nb of assigned cells :',  length(predicted.id[!is.na(predicted.id)]), '\n')
    cat('percent of assigned cells: ', length(predicted.id[!is.na(predicted.id)])/length(predicted.id), '\n')
    
    seurat.obj$predicted.id.scmap = predicted.id
    p1 = DimPlot(seurat.obj, group.by = "predicted.id.scmap", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5,
            na.value = "gray") + 
      ggtitle(paste0("projection into Murray data with scmap (nfeature = ", nb.features.scmap,", threshold = ", 
                     threshold.scmap, ")")) +
      scale_colour_hue(drop = FALSE) + 
      NoLegend()
    
  }
  
  ##########################################
  # step 2: predict unassigned cells using assgined cells with rf and svm
  ##########################################
  library(Seurat)
  seurat.obj = Seurat::FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = nb.features.svm)
  
  sce2 = Seurat::as.SingleCellExperiment(seurat.obj)
  counts(sce2) = as.matrix(counts(sce2))
  logcounts(sce2) = as.matrix(logcounts(sce2))
  rowData(sce2)$feature_symbol <- rownames(sce2)
  
  ## split the Murray data into train and test
  features.sels = rowData(sce2)$vst.variable
  
  index.train = which(!is.na(sce2$predicted.id.scmap))
  index.test = which(is.na(sce2$predicted.id.scmap))
  
  train = sce2[features.sels, index.train]
  test = sce2[features.sels, index.test]
  y <- as.factor(train$predicted.id.scmap)
  
  train = logcounts(train)
  test = logcounts(test)
  train = t(train)
  train <- as.data.frame(train)
  test = t(test)
  test <- as.data.frame(test)
  rownames(train) <- NULL
  rownames(test) <- NULL
  
  ## run rf prediction
  rf.res = run.classifier.rf(train, y, test, ntree = 200, param.tuning = FALSE)
  rf.res = data.frame(rf.res, index.test, stringsAsFactors = FALSE)
  cat('RF nb of cells with probability > 0.5 :', length(which(rf.res$prob>0.5)), '\n')
  cat('RF nb of cells with probability > 0.7 :', length(which(rf.res$prob>0.7)), '\n')
  
  ## run svm prediction
  svm.res = run.classifier.svm(train, y, test, cost = 0.1, param.tuning = FALSE)
  svm.res = data.frame(svm.res, index.test, stringsAsFactors = FALSE)
  
  cat('SVM nb of cells with probability > 0.5 :', length(which(svm.res$prob>0.5)), '\n')
  cat('SVM nb of cells with probability > 0.7 :', length(which(svm.res$prob>0.7)), '\n')
  
  rf.filter = rf.res[which(rf.res$prob > threshold.rf), ]
  seurat.obj$predicted.id.rf = NA;
  seurat.obj$predicted.id.rf[rf.filter$index] = as.character(rf.filter$label)
  
  svm.filter = svm.res[which(svm.res$prob > threshold.svm), ]
  seurat.obj$predicted.id.svm = NA;
  seurat.obj$predicted.id.svm[svm.filter$index] = as.character(svm.filter$label)
  
  #seurat.obj$predicted.id.scmap.svm = seurat.obj$predicted.id.scmap
  #seurat.obj$predicted.id.scmap.svm[svm.filter$index] = as.character(svm.filter$label)
  
  #predicted.id = seurat.obj$predicted.id.scmap.svm
  #cat('after svm classification nb of assigned cells :',  length(predicted.id[!is.na(predicted.id)]), '\n')
  #cat('after svm classification percent of assigned cells: ', length(predicted.id[!is.na(predicted.id)])/length(predicted.id), '\n')
  
  #Idents(seurat.obj) = seurat.obj$predicted.id.scmap.svm
  p2 = DimPlot(seurat.obj, group.by = "predicted.id.rf", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("projection into Murray data with rf (nfeature = ", nb.features.svm,", threshold = ", 
                   threshold.rf, ")")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
  p3 = DimPlot(seurat.obj, group.by = "predicted.id.svm", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5,
          na.value = "gray") + 
    ggtitle(paste0("projection into Murray data with svm (nfeature = ", nb.features.scmap,", threshold = ", 
                   threshold.svm, ")")) +
    scale_colour_hue(drop = FALSE) + 
    NoLegend()
  
   CombinePlots(list(p2, p3))
   
   seurat.obj = annotate.clusters.using.predicted.id(seurat.obj)
   
}

##########################################
# utility functions for the function reference.based.cluster.annotation.scmap
##########################################
annotate.clusters.using.predicted.id = function(seurat.obj)
{
  ms <- FindNeighbors(object = ms, reduction = "mnn", k.param = 20, dims = 1:20)
  ms <- FindClusters(ms, resolution = 12, algorithm = 3)
  
  cluster.ids = unique(seurat.obj$seurat_clusters)
  cluster.ids = cluster.ids[order(cluster.ids)]
  
  pdfname = paste0(resDir, "/annotate_clusters_using_predicted.ids_v1.pdf")
  pdf(pdfname, width=10, height = 8)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(n in 1:length(cluster.ids))
  {
    # n = 2
    cat('cluster ', cluster.ids[n], '\n')
    DimPlot(seurat.obj, 
            cells.highlight = colnames(seurat.obj)[sels],
            group.by = "seurat_clusters", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 4,
            sizes.highlight = 2,
            na.value = "gray") + 
      ggtitle(paste0("cluster ", cluster.ids[n])) +
      NoLegend()
    
    sels = which(seurat.obj$seurat_clusters == cluster.ids[n])
    p1 = DimPlot(seurat.obj, 
            cells.highlight = colnames(seurat.obj)[sels],
            group.by = "predicted.id.scmap", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 1, label.size = 4,
            sizes.highlight = 2,
            na.value = "gray") + 
      ggtitle(paste0("cluster ", cluster.ids[n])) +
      NoLegend()
    
    plot(p1)
    
    pred.ids = seurat.obj$predicted.id.scmap[sels]
    pred.ids[which(is.na(pred.ids))] = 'unassigned'
    par(mar=c(5,8,4,2)) # increase y-axis margin.
    barplot(table(pred.ids)/length(pred.ids), las = 2, horiz = TRUE, xlim =c(0, 1), main = paste0("cluster ", cluster.ids[n]))
    
  }
  
  dev.off()
  
  
}

test.classifier.for.Murray.data = function(method = c('scmap', 'rf', 'gbm', 'svm'))
{
  ##########################################
  # here we test classical predition methods using Johm Murray's data 
  ##########################################
  library(SingleCellExperiment)
  library(irr)
  
  ee = process.import.Murray.scRNA()
  ee = Seurat::as.SingleCellExperiment(ee)
  counts(ee) = as.matrix(counts(ee))
  logcounts(ee) = as.matrix(logcounts(ee))
  rowData(ee)$feature_symbol <- rownames(ee)
  ee$cell_type1 = ee$lineage
  
  ## split the Murray data into train and test
  set.seed(101)
  index.train = sample(1:ncol(ee), 8000)
  train = ee[, index.train]
  test = ee[, -index.train]
  
  # test scmap
  library(scmap)
  nb.features = 500
  threshold = 0.7
  
  train <- selectFeatures(train, suppress_plot = FALSE, n_features = nb.features)
  
  table(rowData(train)$scmap_features)
  as.character(unique(train$cell_type1))
  
  train = indexCluster(train)
  head(metadata(train)$scmap_cluster_index)
  
  scmapCluster_results <- scmapCluster(
    projection = test, 
    index_list = list(
      murray = metadata(train)$scmap_cluster_index
    ),
    threshold = threshold
  )
  
  head(scmapCluster_results$scmap_cluster_labs)
  length(scmapCluster_results$scmap_cluster_labs)
  
  head(scmapCluster_results$scmap_cluster_siml)
  
  hist(scmapCluster_results$scmap_cluster_siml, breaks = 40)
  abline(v = threshold, col = 'red')
  head(scmapCluster_results$combined_labs)
  predicted.id = scmapCluster_results$scmap_cluster_labs
  predicted.id[which(predicted.id == 'unassigned')] = NA
  cat(length(predicted.id[!is.na(predicted.id)]), length(predicted.id[!is.na(predicted.id)])/length(predicted.id), '\n')
  
  kappa0 = kappa2(data.frame(predicted.id, test$lineage), weight = c("unweighted"), sort.levels = FALSE)
  kappa2(data.frame(predicted.id, test$lineage)[!is.na(predicted.id), ], weight = c("unweighted"), sort.levels = FALSE)
  
  ## prepare train and test matrix for rf and gbm
  train.all.feature = train; test.all.feature = test
  table(rowData(train.all.feature)$scmap_features)
  
  sels = rowData(train.all.feature)$scmap_features
  train = train.all.feature[sels, ]
  test = test.all.feature[sels, ]
  y <- as.factor(train$lineage)
  y.test = as.factor(test$lineage)
  train = logcounts(train)
  test = logcounts(test)
  train = t(train)
  train <- as.data.frame(train)
  test = t(test)
  test <- as.data.frame(test)
  rownames(train) <- NULL
  rownames(test) <- NULL
  
  # run random forest
  #res.rf = run.classifier.rf()
  # run svm
  
}

run.classifier.rf = function(train, y, test, ntree = 200, param.tuning = FALSE)
{
  #https://cran.r-project.org/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html
  # http://uc-r.github.io/random_forests # the one I am referring to 
  #library(tree)
  library(randomForest)
  library(ranger)
  library(stats)
  library(mlbench)
  library(caret)
  library(tictoc)
  
  if(!param.tuning){
    #train = scale(train)
    #study = scale(study)
    #ntree = 100
    #mtry = 46
    # x = train
    tic()
    train_rf <- randomForest::randomForest(x = train, y = y, ntree = ntree, keep.forest = TRUE, 
                                           importance = FALSE)
    toc()
    #plot(train_rf)
    
    # tic()
    # m1 = ranger::ranger(x = train, y = y, num.trees = ntree, mtry = mtry, write.forest = TRUE, classification = TRUE, 
    #                     local.importance = TRUE)
    # toc()
    
  }else{
    
    ## tune mtry parameter with tuneRF from randomForest package
    m2 <- tuneRF(
      x          = train,
      y          = y,
      ntreeTry   = ntree,
      mtryStart  = floor(sqrt(ncol(x))),
      stepFactor = 1.5,
      improve    = 0.01,
      trace      = FALSE      # to not show real-time progress 
    )
    print(m2)
    
    # hyperparameter grid search
    hyper_grid <- expand.grid(
      mtry       = seq(20, 60, by = 5),
      node_size  = seq(1, 10, by = 2),
      #sample_size = c(.55, .632, .70, .80),
      sample_size = c(.632),
      prediction_err   = 0
    )
    
    # total number of combinations
    nrow(hyper_grid)
    ## [1] 96
    
    for(i in 1:nrow(hyper_grid)) {
      # train model
      model <- ranger(
        x = train, 
        y = y,
        num.trees       = ntree,
        mtry            = hyper_grid$mtry[i],
        min.node.size   = hyper_grid$node_size[i],
        sample.fraction = hyper_grid$sample_size[i],
        seed            = 123
      )
      
      # add OOB error to grid
      hyper_grid$prediction_err[i] <- model$prediction.error
    }
    
    hyper_grid %>% 
      dplyr::arrange(prediction_err) %>%
      head(10)
    
    index.optimal = which.min(hyper_grid$prediction_err)
    tic()
    train.optimal = ranger::ranger(x = train, y = y, 
                                   num.trees = ntree, 
                                   mtry = hyper_grid$mtry[index.optimal], 
                                   min.node.size = hyper_grid$node_size[index.optimal],
                                   sample.fraction = hyper_grid$sample_size[index.optimal], 
                                   write.forest = TRUE, 
                                   probability = TRUE,
                                   classification = TRUE, 
                                   local.importance = TRUE)
    toc()
    
    saveRDS(train.optimal, file = paste0(RdataDir, 'MurrayData_classifier_test_RF_optimalParam.rds'))
    rf.fit = readRDS(file = paste0(RdataDir, 'MurrayData_classifier_test_RF_optimalParam.rds'))
    #pred_test <-stats::predict(rf.fit, test)
    
    err.test = mean(pred_test==factor(y.test, levels = levels(pred_test)))
  }
  
  Prediction <- stats::predict(train_rf, test, type = "prob")
  #Prediction <- stats::predict(m1, test, type = "prob")
  #Prediction <- predict(train_rf, test, type = "response")$predictions
  rf.res = data.frame(label = apply(Prediction, 1, function(x) colnames(Prediction)[which.max(x)]),
                      prob = apply(Prediction, 1, function(x) x[which.max(x)]))
  #hist(rf.res$prob)
  #kappa2(data.frame(rf.res$label, y.test), weight = c("unweighted"), sort.levels = FALSE)
  return(rf.res)
  
}

run.classifier.svm = function(train, y, test, cost = 1, param.tuning = FALSE)
{
  library(e1071)
  library(stats)
  library(tidyverse)    # data manipulation and visualization
  library(kernlab)      # SVM methodology
  library(RColorBrewer)
  library(tictoc)
  
  if(param.tuning){
    ##########################################
    # tune function in e1071 is extremely slow and not used here
    ##########################################
    # tune.out <- tune(svm, train.x = train, train.y = y, kernel = "linear",
    #                  ranges = list(cost = c(0.1, 1, 5)))
    # 
    # # extract the best model
    # svmfit <- tune.out$best.model
    
    ##########################################
    # here parallel method is used
    # https://www.r-bloggers.com/improve-svm-tuning-through-parallelism/
    ##########################################
    pkgs <- c('foreach', 'doParallel')
    lapply(pkgs, require, character.only = T)
    registerDoParallel(cores = 4)
    
    #set.seed(2016)
    #df2$fold <- caret::createFolds(1:nrow(df2), k = 4, list = FALSE)
    
    #gamma <- c(1, 2)
    parms <- expand.grid(
      cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100) 
    )
    
    ### LOOP THROUGH PARAMETER VALUES ###
    result <- foreach(i = 1:nrow(parms), .combine = rbind) %do% {
      c <- parms$cost[i]
      
      svmfit = e1071::svm(train, y, kernel = 'linear', cost = c, type = "C-classification", 
                          scale = FALSE, probability = TRUE)
      pred <- predict(svmfit, test, probability = TRUE)
      Prediction = attr(pred, "probabilities")
      svm.res = data.frame(label = apply(Prediction, 1, function(x) colnames(Prediction)[which.max(x)]),
                           prob = apply(Prediction, 1, function(x) x[which.max(x)]))
      
      pred_train <-predict(svmfit,train)
      err.train = mean(pred_train==y)
      
      pred_test <-predict(svmfit,test)
      err.test = mean(pred_test==factor(y.test, levels = levels(pred_test)))
      ckappa = kappa2(data.frame(pred_test, y.test))$value
      #print(kappa2(data.frame(svm.res$label, y.test), weight = c("unweighted"), sort.levels = FALSE), '\n')
      
      ### CALCULATE SVM PERFORMANCE ###
      #roc <- pROC::roc(as.factor(out$y), out$prob) 
      data.frame(cost = parms[i, ], err.train, err.test, ckappa)
    }
    
  }else{
    
    tic()
    svmfit <- e1071::svm(train, y, kernel = 'linear',  cost = cost, scale = FALSE, probability = TRUE)
    toc()
    
  }
  
  pred <- predict(svmfit, test, probability = TRUE)
  pred = attr(pred, "probabilities")
  svm.res = data.frame(label = apply(pred, 1, function(x) colnames(pred)[which.max(x)]),
                       prob = apply(pred, 1, function(x) x[which.max(x)]), stringsAsFactors = FALSE)
  
  return(svm.res)
  
}

run.classifier.gbm = function(train, test)
{
  # http://uc-r.github.io/gbm_regression
  library(rsample)
  library(gbm)          # basic implementation
  library(xgboost)      # a faster implementation of gbm
  library(caret)        # an aggregator package for performing many machine learning models
  library(ggplot2)      # model visualization
  #library(lime) 
  library(tictoc)
  
  data = data.frame(y, train)
  
  ## test gbm package
  set.seed(123)
  tic()
  gbm.fit <- gbm(
    formula = y ~ .,
    distribution = "multinomial",
    data = data,
    n.trees = 1000,
    interaction.depth = 1,
    shrinkage = 0.001,
    cv.folds = 4,
    n.cores = NULL, # will use all cores by default
    verbose = FALSE
  )  
  toc()
  # print results
  print(gbm.fit)
  
  ## test xgboost
  set.seed(123)
  params <- list(booster = "gbtree", objective = "multi:softprob", num_class = length(unique(y)), eval_metric = "mlogloss")
  xgb.data = xgb.DMatrix(as.matrix(data[, -1]), label = data[, 1], missing=median)
  xgb.fit1 <- xgb.cv(
    params = params,
    data = xgb.data,
    nrounds = 100,
    eta=0.2,
    nfold = 4,
    verbose = 0               # silent,
  )
  
  
}

process.import.Murray.scRNA = function()
{
  library(VisCello.celegans)
  eset = readRDS(file = paste0('data/Parker_et_al_dataSet_afterFiltering_89701cell.rds'))
  pmeda = data.frame(pData(eset))
  
  ## select the cells for MS lineages
  kk = grep('^MS', pmeda$lineage)
  kk1 = which(pmeda$lineage == '28_cell_or_earlier'| pmeda$lineage == 'ABaxx'| pmeda$lineage == 'Cx'|
                pmeda$lineage == 'Dx'|pmeda$lineage == 'Dxa'|pmeda$lineage == 'Exx')
  
  kk = unique(c(kk, kk1))
  cat('nb of cell in reference -- ', length(kk), '\n')
  cat('nb of cell states in reference -- ', length(unique(pmeda$lineage[kk])), '\n')
  
  ee = CreateSeuratObject(counts = eset@assayData$exprs[,kk], assay = 'RNA', meta.data = pmeda[kk, ])
  ee@assays$RNA@data = eset@assayData$norm_exprs[,kk]
  
  return(ee)
  
}

reference.based.cell.projection.rf.svm = function()
{
  # prepare train and test tables
  #ee <- selectFeatures(ee, suppress_plot = FALSE, n_features = 500)
  #table(rowData(ee)$scmap_features)
  features.sel = rownames(ee)[rowData(ee)$scmap_features]
  features.sel = features.sel[!is.na(match(features.sel, rownames(sce)))]
  ee.sel = ee[match(features.sel, rownames(ee)),]
  sce.sel = sce[match(features.sel, rownames(sce))]
  
  train = logcounts(ee.sel)
  study = logcounts(sce.sel)
  train = t(train)
  train <- as.data.frame(train)
  y <- as.factor(ee.sel$lineage)
  study = t(study)
  study <- as.data.frame(study)
  rownames(train) <- NULL
  rownames(study) <- NULL
  
  if(method == 'rf'){
    rf.res = run.classifier.rf(train, study)
    
    pred.prob.threshod = 0.5
    
    
    
    predicted.id = rf.res$label
    predicted.id[which(rf.res$prob < pred.prob.threshod)] = NA
    seurat.obj$predicted.id = predicted.id
    
    DimPlot(seurat.obj, group.by = "predicted.id", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5,
            na.value = "grey10") + 
      ggtitle(paste0("projection with RF (threshold = ", pred.prob.threshod, ")")) +
      scale_colour_hue(drop = FALSE) + 
      NoLegend()
  }
  
  if(method == 'svm'){
    svm.res = run.classifier.svm(train, study)
    predicted.id = svm.res$label
    predicted.id[which(svm.res$prob < 0.5)] = NA
    seurat.obj$predicted.id = predicted.id
    
    DimPlot(seurat.obj, group.by = "predicted.id", reduction = 'umap', label = TRUE, repel = TRUE, pt.size = 2, label.size = 5,
            na.value = "grey10") + 
      ggtitle("projection with SVM ") +
      scale_colour_hue(drop = FALSE) + 
      NoLegend()
  }
  
}

test.scmap.similarity = function(ee_ref, sce)
{
  ee.sel = as.matrix(metadata(ee_ref)$scmap_cluster_index)
  features.sel = rownames(ee.sel)
  features.sel = features.sel[!is.na(match(features.sel, rownames(sce)))]
  
  ee.sel = ee.sel[match(features.sel, rownames(ee.sel)),]
  sce.sel = sce[match(features.sel, rownames(sce))]
  
  train = ee.sel
  study = logcounts(sce.sel)
  
  library(lsa)
  for(n in 1:nrow(study))
  {
    # n = 1
    xx1 = cor(study[,n], train, method = 'pearson')
    kk1 = which.max(xx1)
    xx1 = xx1[which.max(xx1)]
    xx2 = cor(study[,n], train, method = 'spearman')
    kk2 = which.max(xx2)
    xx2 = xx2[which.max(xx2)]
    xx3= lsa::cosine(study[,n], train)
    kk3 = which.max(xx3)
    xx3 = xx3[which.max(xx3)]
    
    cat(colnames(study)[n], ': ', xx1, xx2, xx3, colnames(train)[c(kk1, kk2, kk3)],  '\n',
        colnames(sce)[n], ': scamp -- ', scmapCluster_results$scmap_cluster_siml[n], '\n')
    n = n + 1
    #readline() 
    #break()
  }
}


##########################################
# Aleks' cell-state assignment using correlation 
##########################################
cell.state.assignment.correlation.based.aleks = function()
{
  dir.processed.data = "/Users/jiwang/workspace/imp/scRNAseq_MS_lineage_dev/results_aleks/results/all_batches_202005/Rdata/"
  load(paste0(dir.processed.data, "all_batches_QCed_cells_genes_filtered_timingEst_Normed_bc_Seurat.Rdata"))
  
  library(dplyr)
  #library(dplyr)
  #setwd("/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/results_aleks/results/all_batches_202005/Rdata/")
  #load("./all_batches_QCed_cells_genes_filtered_timingEst_Normed_bc_Seurat.Rdata")
  
  ms <- FindNeighbors(object = ms, reduction = "mnn", k.param = 20, dims = 1:20)
  ms <- FindClusters(ms, resolution = 12, algorithm = 3)
  
  DimPlot(object = ms, cells = colnames(ms), group.by = 'ident', label = TRUE, pt.size = 4)
  all_ms.markers <- FindAllMarkers(ms, min.pct = 0.25, logfc.threshold = 0.25)
  all_ms.top10 <- all_ms.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  
  lin_sc_expr_190602 <- readRDS("data/lin_sc_expr_190602.rds")
  MS_names <- lin_sc_expr_190602[,grepl(pattern = "MS.", x = colnames(lin_sc_expr_190602))]
  #lin_sc_expr_190602[,MS_names]
  #load(file = paste0(dir.processed.data, 'ms_correlation_idents.rds'))
  
  combined_MS_names <- data.frame(matrix(NA, dim(MS_names)[1], 1))
  used_cells <- c()
  
  for (cell in colnames(MS_names)){
    if (cell %in% used_cells){next}
    indx <- (MS_names[1:100,] == MS_names[1:100,cell])[1,]
    indx_name <- colnames(MS_names)[indx]
    tmp_df <- data.frame(MS_names[,indx_name[1]])
    colnames(tmp_df) <- paste(indx_name, collapse = "/")
    combined_MS_names <- cbind(combined_MS_names,tmp_df)
    used_cells <- c(used_cells, indx_name)
  }
  
  JM_data <- combined_MS_names[,-1]
  my_MS <- as.matrix(ms[["SCT"]]@data)
  
  indexes <- match(all_ms.top10$gene, rownames(my_MS))
  my_MS <- my_MS[indexes,]
  indexes <- match(rownames(my_MS), rownames(JM_data))
  mathed_gene_names <- rownames(JM_data)[indexes[!is.na(indexes)]]
  JM_data <- JM_data[mathed_gene_names,]
  my_MS <- my_MS[mathed_gene_names,]
  
  cor_vec <- c()
  
  combined_df <- data.frame(row.names = 1:(dim(JM_data)[2]))
  combined_df_corv <- data.frame(row.names = 1:(dim(JM_data)[2]))
  
  library(progress)
  pb <- progress_bar$new(
    format = " progress [:bar] :percent eta: :eta",
    total = length(colnames(my_MS)))
  
  for (a in colnames(my_MS)){
    for (i in 1:dim(JM_data)[2]){
      cor_vec[i] <- cor(my_MS[,a], JM_data[,i])
    }
    combined_df[,a] <- colnames(JM_data)[order(cor_vec, decreasing = T)]
    combined_df_corv[,a] <- cor_vec[order(cor_vec, decreasing = T)]
    pb$tick()
  }
  
  
  View(combined_df_corv[1:10,1:40])
  cell.names <- colnames(combined_df)
  
  Idents(ms, cells = cell.names) <- combined_df[1,]
  
  pb <- progress_bar$new(
    format = " progress [:bar] :percent eta: :eta",
    total = length(levels(Idents(ms))))
  
  pdfname = paste0(paste0(resDir, "/idents_correlation_aleks.pdf"))
  pdf(pdfname, width=16, height = 16)
  
  for (idntt in levels(Idents(ms))){
    
    #pdfname = paste0(paste0("~/Documents/plots/", which(levels(Idents(test_ms)) == idntt), ".pdf"))
    #pdf(pdfname, width=8, height = 8)
    
    cells_to_show <- list(c(WhichCells(ms, idents = idntt)))
    names(cells_to_show) <- idntt
    print(DimPlot(ms, reduction = "umap", label = F, cells.highlight = cells_to_show))
    #dev.off()
    pb$tick()
    
  }
  
  dev.off()
  
  ms_correlation_idents <- ms
  
  saveRDS(ms_correlation_idents, file = paste0(RdataDir,  'ms_correlation_idents.rds'))
  saveRDS(combined_df, file = paste0(RdataDir,  'cell.states_assignment_to_reference.JM_correlation_idents.rds'))
  saveRDS(combined_df_corv, file = paste0(RdataDir,  'cell.states.correlation_assignment_to_reference.JM_correlation_idents.rds'))
  
}

########################################################
# Section : Clustering section by integrating various informations: 
# gene expression, fac info, estimated timing and velocity 
########################################################
########################################################
cluster.by.integration.gene.expression.fac.timingEst.velocity = function()
{
  library(scater)
  library(scran)
  library(scRNA.seq.funcs)
  library(matrixStats)
  #library(M3Drop)
  library(RColorBrewer)
  library(SingleCellExperiment)
  
  load(file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_bcMNN.Rdata')) 
  
  Seurat.clustering = TRUE
  Test.scran.clustering = FALSE
  fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
  
  ## double check the umap and tsnet visulization and 
  set.seed(100)
  sce <- runUMAP(sce, use_dimred="MNN", n_dimred = 15)
  p1 = plotUMAP(sce, ncomponents = 2, colour_by="batches", size_by = "FSC_log2", point_size= 0.01) + ggtitle("Corrected-umap") 
  plot(p1)
  
  sce = runTSNE(sce, use_dimred="MNN", perplexity = 20, n_dimred = 15)
  p2 = plotTSNE(sce, ncomponents = 2, colour_by="batches", size_by = "FSC_log2", point_size= 0.01) + ggtitle("Corrected-tSNE") 
  plot(p2)
  
  p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  multiplot(p1, p2, cols = 2)
  
  
  p1 = plotTSNE(sce, colour_by="pha-4", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  p2 = plotTSNE(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
    fontsize + ggtitle("MNN corrected")
  multiplot(p1, p2, cols = 2)
  
  ##########################################
  # (test clustering method) but Seurat method is used at the end 
  # https://master.bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/work-1-reads.html/
  # DE analysis is tested with scran but will be considered to be replaced by MAST or Seurat
  # DE analysis (or marker gene finding) following the cluster analysis
  ##########################################
  if(Seurat.clustering)
  { 
    require(Seurat)
    library(cowplot)
    
    check.individualExample.geneMarker = FALSE
    ntops = 3 # nb of top gene markers
    resolutions = c(0.4, 0.6, 0.8, seq(1.0, 4.0, by = 0.5))
    
    for(rr in resolutions){
      
      rr = 1.2
      cat("--- resolution is :", rr, "---\n")
      
      pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_batchCorrected_clustering.Seurat_geneMarkers.scran_resolution_", 
                       rr, ".pdf")
      pdf(pdfname, width=22, height = 18)
      par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
      
      ##########################################
      # test graph-based Louvain algorithm 
      ##########################################
      pbmc = as.Seurat(sce)
      #pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
      #pbmc <- ScaleData(pbmc, features = rownames(pbmc))
      # pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
      pbmc = FindNeighbors(object = pbmc, reduction = "MNN", k.param = 20, dims = 1:20)
      pbmc = FindClusters(pbmc, resolution = rr, algorithm = 3)
      
      # note that you can set `label = TRUE` or use the LabelClusters function to help label
      # individual clusters
      #pbmc <- Seurat::RunUMAP(pbmc, dims = 1:10, reduction = "MNN", reduction.key = "umap", n.neighbors = 15, repulsion.strength = 2)
      #DimPlot(pbmc, reduction = "umap")
      #FeaturePlot(pbmc, features = c("pha-4", "hnd-1"), reduction = "umap")
      # FeaturePlot(pbmc, features = c("pha-4", "hnd-1"), reduction = "UMAP")
      sce$cluster_seurat <- factor(pbmc@active.ident)
      sce$cluster <- factor(pbmc@active.ident)
      
      plotUMAP(sce, colour_by="cluster", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      
      plotTSNE(sce, colour_by="cluster", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      
      p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      multiplot(p1, p2, cols = 2)
      
      p1 = plotTSNE(sce, colour_by="pha-4", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      p2 = plotTSNE(sce, colour_by="hnd-1", size_by = "FSC_log2") + 
        fontsize + ggtitle("Seurat clustering")
      multiplot(p1, p2, cols = 2)
      
      my.clusters = as.numeric(as.character(sce$cluster_seurat))
      cat(table(my.clusters), "\n")
      
      ## run the find markers and then collect markers for each clusters (scran) 
      # https://bioconductor.org/packages/3.10/workflows/vignettes/simpleSingleCell/inst/doc/de.html#blocking-on-uninteresting-factors-of-variation
      #design <- model.matrix( ~ sce$batches)
      #design <- design[,-1,drop=FALSE]
      markers <- findMarkers(sce, my.clusters, block=sce$batches, direction="up")
      top.markers = c()
      
      for(n in unique(my.clusters)){
        #n = 0
        marker.set <- markers[[as.character(n)]]
        #marker.set <- markers[["1"]]
        #head(marker.set, 5)
        top.markers <- c(top.markers, rownames(marker.set)[marker.set$Top <= ntops])  
      }
      
      top.markers = unique(top.markers)
      cat(length(top.markers), " gene markers found by scran \n")
      
      plotHeatmap(sce, features=top.markers,
                  columns=order(sce$cluster_seurat), 
                  colour_columns_by=c("cluster",  "FSC_log2", "batches"),
                  cluster_cols=FALSE, show_colnames = FALSE,
                  center=TRUE, symmetric=TRUE, zlim=c(-5, 5))
      
      if(check.individualExample.geneMarker){
        for(n in 1:length(top.markers)) {
          xx = plotTSNE(sce, colour_by = top.markers[n]) 
          plot(xx)
        }
      }
      
      dev.off()
      
    }
    
  }
  
  ##########################################
  # here select subset of whole dataset to redo clustering
  # and marker gene finding
  # this is the beginning of iterative clustering, 
  # the gene marker finding will be a criterion to stop iteration
  ##########################################
  Select.Early.timePoints = FALSE
  if(Select.Early.timePoints){
    
    pbmc = as.Seurat(sce)
    #Seurat::DimPlot(pbmc, dims = c(1, 2), reduction = "MNN")
    pbmc = FindNeighbors(object = pbmc, reduction = "MNN", k.param = 10, dims = 1:10)
    pbmc = FindClusters(pbmc, resolution = 2, algorithm = 3)
    sce$cluster_seurat <- factor(pbmc@active.ident)
    sce$cluster <- factor(pbmc@active.ident)
    
    xx = table(sce$cluster,sce$Batch)
    #colnames(xx) = sce$Batch
    cluster4early = rownames(xx)[which(xx[, 1]>=5|xx[,2]>=5)]
    
    mm = match(sce$cluster, factor(cluster4early))
    
    sels = which(!is.na(mm))
    
    cat("estimated cell in early stage -- ", length(cluster4early), 
        "clusters with", length(sels),  "cells\n")
    
    sce.sel = sce[, sels ]
    # here rerun the clustering for clusters or stages
    pcs <- reducedDim(sce.sel, "MNN")
    my.dist <- dist(pcs[, c(1:20)])
    my.tree <- hclust(my.dist, method="ward.D2")
    
    #hist(my.tree)
    library(dynamicTreeCut)
    my.clusters <- unname(cutreeDynamic(my.tree,
                                        distM=as.matrix(my.dist), 
                                        minClusterSize=5, verbose=0))
    
    table(my.clusters, sce.sel$Batch)
    
    sce.sel$cluster <- factor(my.clusters)
    
    set.seed(100)
    sce.sel <- runTSNE(sce.sel,  use_dimred = "MNN", perplexity = 20, n_dimred = 20)
    
    set.seed(100)
    sce.sel = runUMAP(sce.sel, use_dimred="MNN", perplexity = 20, n_dimred = 20)
    
    plotTSNE(sce.sel, colour_by="cluster", size_by = "total_features_by_counts") + 
      fontsize + ggtitle("scran -- hcluster")
    plotUMAP(sce.sel, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
      fontsize + ggtitle("scran -- hcluster")
    
    design <- model.matrix( ~ sce.sel$Batch)
    design <- design[,-1,drop=FALSE]
    
    markers <- findMarkers(sce.sel, sce.sel$cluster, design=design, direction = 'any')
    
    ntops = 5;
    top.markers = c()
    
    for(n in unique(my.clusters)){
      #n = 0
      marker.set <- markers[[as.character(n)]]
      #marker.set <- markers[["1"]]
      #head(marker.set, 5)
      top.markers <- c(top.markers, rownames(marker.set)[marker.set$Top <= ntops])  
    }
    top.markers = unique(top.markers)
    
    pdfname = paste0(resDir, "/scRNAseq_QCed_filtered_normalized_batchCorrected_clustering_markerGenes_earlyTimepoint_tSNEexample.pdf")
    pdf(pdfname, width=12, height = 10)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    plotTSNE(sce.sel, colour_by="cluster", size_by = "total_features_by_counts") + 
      fontsize + ggtitle("scran -- hcluster (tSNE)")
    plotUMAP(sce.sel, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
      fontsize + ggtitle("scran -- hcluster (UMAP)")
    
    plotHeatmap(sce.sel, features=top.markers,
                columns=order(sce.sel$cluster), 
                colour_columns_by=c("cluster"),
                cluster_cols=FALSE, show_colnames = FALSE,
                center=TRUE, symmetric=TRUE, zlim = c(-5, 5))
    
    for(n in 1:length(top.markers)) {
      xx = plotTSNE(sce.sel, colour_by = top.markers[n]) 
      plot(xx)
    }
    
    dev.off()
    
  }
}


########################################################
########################################################
# Section :
# test code and not used anymore 
########################################################
########################################################
test.clustering.method.scran = function(sce){
  if(Test.scran.clustering){
    ##########################################
    # test grapha method from scran
    ##########################################
    snn.gr <- buildSNNGraph(sce, use.dimred="MNN")
    clusters <- igraph::cluster_walktrap(snn.gr)
    table(clusters$membership, sce$Batch)
    
    sce$cluster <- factor(clusters$membership)
    plotTSNE(sce, colour_by="cluster") + ggtitle("scran -- graph based clustering")
    plotUMAP(sce, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
      fontsize + ggtitle("scran -- graph based clustering")
    
    ##########################################
    # test hierachy clustering from scran
    ##########################################
    pcs <- reducedDim(sce, "MNN")
    my.dist <- dist(pcs)
    my.tree <- hclust(my.dist, method="ward.D2")
    
    #hist(my.tree)
    library(dynamicTreeCut)
    my.clusters <- unname(cutreeDynamic(my.tree, cutHeight = 3, 
                                        distM=as.matrix(my.dist), 
                                        minClusterSize=5, verbose=0))
    
    table(my.clusters, sce$Batch)
    
    sce$cluster <- factor(my.clusters)
    
    #plotTSNE(sce, colour_by="cluster", size_by = "total_features_by_counts") + fontsize + ggtitle("scran -- hcluster")
    plotUMAP(sce, colour_by="cluster", size_by = "total_features_by_counts", shape_by = "Batch") + 
      fontsize + ggtitle("scran -- hcluster")
    
    #plotDiffusionMap(sce, colour_by="cluster", size_by = "total_features_by_counts") + fontsize
    library(cluster)
    clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
    sil <- silhouette(my.clusters, dist = my.dist)
    sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
    sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
    plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
         border=sil.cols, col=sil.cols, do.col.sort=FALSE) 
    
  }
}

##########################################
# test seurat for clustering 
# original codes in https://satijalab.org/seurat/pbmc3k_tutorial.html
##########################################
Test.Seurat.workflow = FALSE
if(Test.Seurat.workflow){
  library(Seurat)
  library(dplyr)
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata'))
  
  #pbmc = as.seurat(from = sce)
  #pbmc = Convert(from = sce, to = 'seurat', raw.data.slot = "counts",data.slot = "logcounts")
  pbmc <- CreateSeuratObject(raw.data = counts(sce), min.cells = 0, min.genes = 200, 
                             project = "seurat_test")
  pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", 
                        scale.factor = 10000)
  
  pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                            x.low.cutoff = 0.5, y.cutoff = 0.5)
  
  length(x = pbmc@var.genes)
  
  pbmc <- ScaleData(object = pbmc)
  
  pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, 
                 genes.print = 5)
  
  PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
  VizPCA(object = pbmc, pcs.use = 1:2)
  
  PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
  
  pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
  
  PCHeatmap(object = pbmc, pc.use = 1:4, cells.use = 80, do.balanced = TRUE, label.columns = FALSE)
  
  pbmc <- JackStraw(object = pbmc, num.replicate = 100, display.progress = FALSE)
  JackStrawPlot(object = pbmc, PCs = 1:12)
  
  PCElbowPlot(object = pbmc)
  
  pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, 
                       resolution = 2, print.output = 0, save.SNN = TRUE)
  
  PrintFindClustersParams(object = pbmc)
  #PrintCalcParams(pbmc, calculation = "FindClusters")
  
  pbmc <- RunTSNE(object = pbmc, dims.use = 1:10, do.fast = TRUE, perplexity=10, eta=2000)
  TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 3.0)
  
  ## test phate plot
  pbmc_phate <- RunPHATE(object = pbmc, gamma=1, npca = 20, k=5, seed.use = 10, plot.optimal.t = TRUE, t=14)
  # Plot results
  DimPlot(object = pbmc_phate, reduction.use = 'phate', pt.size = 2.0)
  #DimPlot(object = pbmc_phate, reduction.use = 'pca')
}

##################################################
# test Hamberg's single-cell RNA seq analysis, clustering part
# original code https://hemberg-lab.github.io/scRNA.seq.course/biological-analysis.html
##################################################
### test SC3
TEST.Hamberg.workflow.clustering = FALSE
if(TEST.Hamberg.workflow.clustering)
{
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE.Rdata')) 
  library(SC3)
  #library()
  expr_matrix =  exp(logcounts(sce))
  Brennecke_HVG <- BrenneckeGetVariableGenes(
    expr_mat = expr_matrix,
    spikes = NA,
    fdr = 0.2,
    minBiolDisp = 0.2
  )
  #HVG_genes <- Brennecke_HVG
  
  ## another method to identify by Kolodziejczyk AA, Kim JK, Tsang JCH et al. (2015)
  assay(sce, "normcounts") <- exp(logcounts(sce))
  means <- rowMeans(normcounts(sce))
  cv2 <- apply(normcounts(sce), 1, var)/means^2
  dm.stat <- DM(means, cv2)
  #head(dm.stat)
  
  DM_HVG = names(dm.stat)[which(dm.stat>0.3)]
  #DM_HVG = DM_HVG[which()]
  #dev.off()
  
  sce.HVG.Brenneck = sce[rownames(sce)%in%Brennecke_HVG, ]
  sce.HVG.DM = sce[rownames(sce)%in%DM_HVG, ]
  
  save(sce, sce.HVG.Brenneck, sce.HVG.DM, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  #HVG_genes <- Brennecke_HVG$Gene
  #load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_HVGsels.Rdata'))
  
  sce.sels = sce.HVG.Brenneck
  #sce.sels = sce.HVG.DM
  sce = sce.sels
  # define feature names in feature_symbol column
  rowData(sce)$feature_symbol <- rownames(sce)
  # remove features with duplicated names
  sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
  
  plotPCA(
    sce,
    run_args = list(exprs_values = "logcounts"),
    colour_by = "total_counts",
    size_by = "total_features_by_counts"
  )
  
  ### test SC3 clustering method
  sce = sc3_estimate_k(sce)
  metadata(sce)$sc3$k_estimation
  
  sce <- sc3(sce.sels, gene_filter = FALSE, ks = 2:30, biology = TRUE, n_cores = 6)
  #rowData(sce)$feature_symbol
  #sce <- sc3(sce, ks = 2, gene_filter = TRUE, biology = TRUE)
  #deng <- plotTSNE(deng, rand_seed = 1, return_SCE = TRUE)
  
  col_data <- colData(sce)
  head(col_data[ , grep("sc3_", colnames(col_data))])
  plotPCA(
    sce, 
    colour_by = "sc3_10_clusters", 
    size_by = "sc3_10_log2_outlier_score"
  )
  
  plotPCA(
    sce, 
    colour_by = "sc3_3_clusters",
    shape_by = "celltypes",
    size_by = "total_features_by_counts"
  )
  
  row_data <- rowData(sce)
  head(row_data[ , grep("sc3_", colnames(row_data))])
  
  plotFeatureData(
    sce, 
    aes(
      x = sc3_3_markers_clusts, 
      y = sc3_3_markers_auroc, 
      colour = sc3_3_markers_padj
    )
  )
  
  set.seed(1)
  plotTSNE(sce, colour_by="sc3_6_clusters", size_by = "total_features_by_counts",
           run_args = list(perplexity=20)) 
  
  #sc3_plot_consensus(sce, k = 3)
  sc3_plot_consensus(
    sce, k = 3, 
    show_pdata = c(
      "celltypes", 
      "sc3_3_clusters", 
      "log10_total_features_by_counts",
      "log10_total_counts",
      
      "sc3_3_log2_outlier_score"
    )
  )
  
  sc3_plot_cluster_stability(sce, k = 6)
  
  sc3_plot_de_genes(sce, k = 3, p.val = 0.2)
  
  sc3_plot_markers(sce, k = 3)
  
}

##########################################
# # configuration of umap
# https://github.com/cran/umap/blob/master/R/umap.R
##########################################
config.umap = function(n_neighbors=15,
                       n_components=2,
                       metric="euclidean",
                       n_epochs=200,
                       input="data",
                       init="spectral",
                       min_dist=0.5,
                       set_op_mix_ratio=1,
                       local_connectivity=1,
                       bandwidth=1.0,
                       alpha=1,
                       gamma=1.0,
                       negative_sample_rate=5,
                       a=NA,
                       b=NA,
                       spread=1,
                       random_state=NA,
                       transform_state=NA,
                       knn_repeats=1,
                       verbose=TRUE,
                       umap_learn_args = NA)
{
  umap.defaults = list(
    n_neighbors = n_neighbors,
    n_components = n_components,
    metric = metric,
    n_epochs = n_epochs,
    input = input,
    init = init,
    min_dist = min_dist,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    bandwidth = bandwidth,
    alpha = alpha,
    gamma = gamma,
    negative_sample_rate = negative_sample_rate,
    a = a,
    b = b,
    spread = spread,
    random_state = random_state,
    transform_state = transform_state,
    knn_repeats = knn_repeats,
    verbose = verbose,
    umap_learn_args = umap_learn_args
  )
  
  class(umap.defaults) = "umap.config"
  
  return(umap.defaults)
  
}




