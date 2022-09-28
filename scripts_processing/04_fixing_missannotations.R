#Figuring out the problem with MSxa cluster
#Correcting some small bugs in the dataset


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



setwd(source.path)

library(Seurat)
library(ggplot2)
library(dplyr)
library(slingshot)



######## fixing NHR-67 CONTAMINATION in MSXA
Idents(ms) <- ms$manual.annot.ids

MSxa_MSxapp <- subset(ms, idents =c('MSxa','MSxapp','MSxap'))

table(MSxa_MSxapp$flowcell.lane)

DoHeatmap(MSxa_MSxapp, features = 'nhr-67', group.by = 'ident')
View(MSxa_MSxapp)

assay_data <- GetAssayData(MSxa_MSxapp[,MSxa_MSxapp@active.ident == 'MSxa'], slot = "scale.data", assay = 'RNA')
MSa_contaminants <- colnames(assay_data)[assay_data[row.names(assay_data) == 'nhr-67',] > 0]


MSxa_MSxapp[,colnames(MSxa_MSxapp) %in% MSa_contaminants]$index.well

subset(ms, idents=c('MSxapp'))$index.well

MSxa_MSxapp_dataset <- data.frame()
isContaminant.temp <- rownames(MSxa_MSxapp_dataset) %in% MSa_contaminants
isContaminant.temp <- isContaminant
isContaminant[isContaminant.temp] <- 'Contaminant'
isContaminant[!isContaminant.temp] <- 'Not'
MSxa_MSxapp_dataset <- cbind(MSxa_MSxapp$samples, as.numeric(MSxa_MSxapp$FSC_log2), as.numeric(MSxa_MSxapp$BSC_log2), MSxa_MSxapp$manual.annot.ids, isContaminant)
colnames(MSxa_MSxapp_dataset) <- c('name','FSC_log2','BSC_log2','ID','contaminant')
MSxa_MSxapp_dataset <- as.data.frame(MSxa_MSxapp_dataset)

MSxa_MSxapp_dataset$FSC_log2 <- as.numeric(MSxa_MSxapp_dataset$FSC_log2)
MSxa_MSxapp_dataset$BSC_log2 <- as.numeric(MSxa_MSxapp_dataset$BSC_log2)

library(tidyr)
#pivot_longer(MSxa_MSxapp_dataset, cols = c('FSC_log2', 'BSC_log2'), names_to = 'facs', values_to = 'facs_value')

ggplot(as.data.frame(MSxa_MSxapp_dataset), aes(x=FSC_log2, y=BSC_log2, color=ID, shape=contaminant)) + 
  geom_point()+
  scale_shape_manual(values=c(16, 5))+
  scale_color_manual(values=c('#007c92','#f3ff82','#ac14d5'))+
  scale_size_manual(values=c(0.5))



######## DEALING WITH another contamination

Idents(ms) <- ms$manual.annot.ids
MSxa_MSxapp <- subset(ms, idents =c('MSx','MSxa','MSxaa'))
table(MSxa_MSxapp$flowcell.lane)
DoHeatmap(MSxa_MSxapp, features = 'ceh-51', group.by = 'ident')
View(MSxa_MSxapp)
assay_data <- GetAssayData(MSxa_MSxapp[,MSxa_MSxapp@active.ident == 'MSx'], slot = "scale.data", assay = 'RNA')
MSa_contaminants <- colnames(assay_data)[assay_data[row.names(assay_data) == 'ceh-51',] < 1.6]
ms_filtered <- ms[,!colnames(ms) %in% MSa_contaminants]
ms <- ms_filtered
saveRDS(ms, file = '~/R_projects/BWMs/processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_cleanedBWM_and_Pharynx_iteration_GR15_fixed.rds')


# FIXING MSpaaapp
MSpaaapp <- subset(ms, idents ='MSpaaapp')
levels(MSpaaapp) <- 'MSpaaapp'
DoHeatmap(MSpaaapp, features = c('myo-3', 'spp-15', 'unc-15', "emb-9", 'lev-11'))

MSpaaapp <- FindNeighbors(object = MSpaaapp, reduction = "MNN", k.param = 20, dims = 1:20)
MSpaaapp <- FindClusters(MSpaaapp, resolution = 0.7, algorithm = 3)
DimPlot(MSpaaapp, group.by = 'ident', label = T, reduction = 'UMAP_MNN', pt.size = 5)

HVGs <- FindAllMarkers(MSpaaapp, min.pct = 0.35)
top.markers <- c()
ntops <-  30
for(n in unique(HVGs$cluster)){
  #n = 0
  marker.set <- subset(HVGs, cluster == n)
  #marker.set <- markers[["1"]]
  #head(marker.set, 5)
  #top.markers <- c(top.markers, rownames(marker.set)[1:ntops])  
  top.markers <- c(top.markers, marker.set$gene[1:ntops])  
}
top.markers <- unique(top.markers)
DoHeatmap(MSpaaapp, features = c(top.markers,'myo-3', 'spp-15', 'unc-15', "emb-9", 'lev-11'))



Idents(ms, cells = colnames(MSpaaapp[,MSpaaapp@active.ident == '0'])) <- 'unidentified_MSxapp-derived_BWMs'
Idents(ms, cells = colnames(MSpaaapp[,MSpaaapp@active.ident == '1'])) <- 'MSaaapapp_pm8_potentially'

saveRDS(ms, file = '~/R_projects/BWMs/processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_cleanedBWM_and_Pharynx_iteration_GR15_fixed.rds')


# FIXING annotaions of the MSxpappp and MSxpappa based on expression of the clec-264 expression

former_MSxpappa <- Cells(subset(x = ms, idents = 'MSxpappa'))
former_MSxpappp <- Cells(subset(x = ms, idents = 'MSxpappp'))

Idents(ms, cells = former_MSxpappa) <- 'MSxpappp'
Idents(ms, cells = former_MSxpappp) <- 'MSxpappa'

ms_sbst <- subset(ms, idents = c('MSxpappp', 'MSxpappa', "MSxppppp" ))

DoHeatmap(ms_sbst, features = c('clec-264', 'zig-6', 'T08H10.1', 'lin-39'))
saveRDS(ms, file = '~/R_projects/BWMs/processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_cleanedBWM_and_Pharynx_iteration_GR15_fixed.rds')





