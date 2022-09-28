##########################################################################
##########################################################################
# Project: Figures for the paper

##########################################################################
##########################################################################

install.packages("devtools")
devtools::install_github("cole-trapnell-lab/monocle-release@develop")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')


#Function for soursing functions
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
library(viridisLite)
library(viridis)
library(patchwork)
library(princurve)
library(ggpubr)
library(RColorBrewer)
library(slingshot)
library(ggsignif)
library(ggpubr)
library(tidyr)
library(sccore)
library(pheatmap)
library(psych)


ms <- readRDS('~/R_projects/BWMs/processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_cleanedBWM_and_Pharynx_iteration_GR15_fixed.rds')
ms_backup <- ms



Idents(ms) <- ms$manual.annot.ids
all_ids <- c(levels(Idents(ms)))[-74]
ms_new <- subset(ms, idents = all_ids)
ms_new@active.ident <- factor(x = ms_new$manual.annot.ids, levels =  all_ids)
ms_new <- subset(ms, idents =c("MSxpp","MSxapa","MSxap", "MSxp",  "MSxpa", "MSxaa",
                               "MSxaaa","MSxa","MSaaaaaa","MSx","MSxppa","MSxapap",
                               "MSxpaa","MSxppaa","MSaaaaap","MSaaapp", "MSxppp","MSxpapp",
                               "MSxpaaa", "MSxapp","MSaaaap", "MSxppppx","MSxpapa", "MSpaapaa",                    
                               "MSxppap", "MSxpppa", "MSxaap","MSaaapaa","MSxpppp", "MSpappa",
                               "MSxapaa", "MSxappp", "MSxaapapx","MSxpap","MSaappa", "MSaaaaaaa",                 
                               "MSxaapa", "MSpaaaa", "MSxapapa.MSxapapaa","MSpaaap", "MSxapppa","MSxpaaaa",
                               "MSaaaapp","MSxapaap","MSxpapap","MSxapaaa","MSapaapp","MSxpapaa",
                               "MSaaaaapa","MSpaaaaa","MSpaaapa","MSaaaappp/MSxapaapp","MSxppapp","MSxpppap",
                               "MSxppppa","MSxapppp","MSxppaaa","MSxaapap","MSxpappp","MSpaaaap",
                               "MSpaaapp","MSxapapp","MSxpppaa","MSxapppax","MSxppppp","MSxpaaap",
                               "MSaaaapa","MSxppaap","MSxappppx","MSxapappp.like","MSxapaapa","MSxpappa",
                               "MSaaappp","MSpaaappp/MSxapappa","MSxpaap.MSppaapp.like","MSxpppaa_MSxppppa_later_like"))

ms_new$manual.annot.ids



################
## Figure 1 
## Panel C - Timing Estimation with the corr. data
################

levels(ms$timingEst) <- c('<30', levels(ms$timingEst))
tmp <- which(as.numeric(as.character(ms$timingEst)) <= 30)
ms$timingEst[tmp] <- '<30'
levels(ms$timingEst) <- c(levels(ms$timingEst), '>400')
ms$timingEst[as.numeric(as.character(ms$timingEst)) >= 400] <- '>400'
DimPlot(ms, reduction = 'umap', group.by = 'timingEst', label = F, cols = viridis(30))
#export as PNG 700 x 550

################
## Figure 1 
## Panel D - Timing Estimation with FACS
################

rcorel <- matrix(c(ms$FSC_log2, ms$BSC_log2),ncol = 2)
plot(rcorel) # Just to see the FACS data on the plot

fit <- principal_curve(rcorel, maxit = 20)
plot(fit) # just the representation of the fitted curve
lines(fit)
points(fit)
whiskers(rcorel, fit$s)

# just to see how timeEst corellates with lambda
rcorellation <- data.frame(x= as.numeric(as.character(ms$timingEst)), y =fit$lambda)
ggscatter(rcorellation, x = "x", y = "y", size = 0.1,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "estimated time", ylab = "cell size")

ms$facs_groups <- fit$lambda
fsc_factor <- cut(ms$facs_groups, breaks = seq(max(ms$facs_groups), min(ms$facs_groups), length.out = 31))
levels(fsc_factor) <- rev(levels(fsc_factor))
ms$facs_groups <- fsc_factor
levels(ms$facs_groups) <- paste0('Group ', 1:30)
DimPlot(ms, reduction = 'umap', group.by = 'facs_groups', label = F, cols = viridis(30))
#export as PNG 700 x 550

################
## Figure 2
## Panel B,C 
################

#nb of cells is necessary to just create an empty plot. Export 550 x 500 px - so that the size of the font is bigger
FeaturePlot(ms, features = c('pha-4', 'hnd-1', 'unc-120', 'nb.cells'), cols = viridis(14),reduction = 'umap', keep.scale = "all", combine = T, slot = 'data')

#remove cells with weird name
ms_subset <- subset(ms, idents = 'weirdos.possible.doublets', invert = TRUE)
DimPlot(ms_subset, reduction = 'umap', group.by = 'ident', label = T, cols = turbo(100), repel = TRUE) + NoLegend()

################
## Figure 2
## Panel D,E,F 
################

#the list of all TFs in C elegans
all_TFs <- c("cnd-1","cep-1","egl-46","ham-2","eor-1","C17E4.6","F28C6.1","F28C6.2","K06A1.1","Y62E10A.17",
             "hif-1","mxl-1","tag-295","Y53H1A.5","unc-120","mef-2","C09G9.7","Y51H4A.17","tbx-38","tbx-31","tbx-18","tbx-2","tbx-11","mls-1","tbx-8","tbx-9","tbx-41","mab-9","tbx-34","tbx-37","tbx-30","Y59E9AR.5","tbx-33","tbx-39","tbx-40","tbx-35","tbx-7",
             "tbx-32","tbx-36","C30G4.7","T18D3.7","Y48C3A.12","fkh-6","fkh-10","fkh-3","fkh-4","C34B4.2","unc-130","fkh-5","let-381","pha-4","fkh-8","fkh-9","lin-31","fkh-2","T27A8.2","pes-1","daf-16","Y53C10A.3",
             "dpl-1","efl-1","efl-2","F49E12.6","grl-25","fkh-7","K04C1.3","ekl-2","C18G1.2","elt-2","elt-4","elt-6","egl-18","end-3","elt-3","med-2","med-1","C16A3.7","ZK1067.2","nhr-150","nhr-63","nhr-50","nhr-152","nhr-94","nhr-155","nhr-257","nhr-72",
             "nhr-30","nhr-258","nhr-74","nhr-73","nhr-43","nhr-162","nhr-107","nhr-163","nhr-42","nhr-161","nhr-139","nhr-140","nhr-164","nhr-165","nhr-81","nhr-75","nhr-261","nhr-167","nhr-169","nhr-170","nhr-171","nhr-172","nhr-137","nhr-89","nhr-174","nhr-264","nhr-176",
             "F16B12.8","nhr-177","nhr-179","nhr-117","nhr-178","nhr-45","nhr-27","nhr-180","nhr-15","nhr-108","nhr-78","nhr-54","nhr-183","nhr-182","nhr-82","nhr-265","nhr-20","nhr-39","nhr-126","nhr-96","nhr-134","nhr-18","nhr-103","nhr-128","nhr-56","nhr-133",
             "nhr-184","nhr-185","nhr-186","nhr-187","nhr-188","nhr-189","nhr-83","nhr-7","nhr-266","F57A10.6","nhr-193","nhr-195","nhr-143","nhr-267","nhr-51","nhr-199","nhr-53","nhr-196","nhr-268","nhr-71","nhr-99","nhr-98","nhr-123","nhr-125","nhr-12","nhr-269","dpr-1","nhr-132","nhr-210","nhr-58","nhr-131","nhr-106",
             "nhr-211","nhr-212","nhr-55","nhr-130","nhr-271","nhr-57","nhr-213","nhr-102","nhr-84","nhr-272","nhr-214","nhr-215","nhr-216","nhr-217","nhr-109","nhr-218","nhr-77","odr-7","nhr-44","nhr-220","srw-49","T26E4.16","nhr-223","nhr-79","nhr-59","nhr-115","nhr-227","nhr-228","nhr-229","nhr-230","nhr-231","nhr-232","tag-122","nhr-234","nhr-95","nhr-92","nhr-122","nhr-110","nhr-275","nhr-13","nhr-241",
             "nhr-112","nhr-276","nhr-243","nhr-277","nhr-245","nhr-113","nhr-246","nhr-248","nhr-9","nhr-250","nhr-90","nhr-251","nhr-278","nhr-252","nhr-253","nhr-254","nhr-255","nhr-256","nhr-219","nhr-135",
             "C01G8.9","cfi-1","rbr-2","C27B7.4","isw-1","tag-347","C44F1.2","ahr-1","F23F1.1","F40F9.7","dro-1","T26A5.8","W10D9.4","Y53F4B.3","ceh-26","mbf-1","T01C1.3","B0261.1",
             "F10E7.11","R151.8","T04D1.4","T07F8.4","ekl-4","psa-1","gei-8","D1081.8","dnj-11","F54F2.9","gei-11","C01B12.2","F53H4.5","F22D6.2","rabx-5","F25H8.6","K09A11.1","Y105C5A.15","Y39B6A.12","Y47D3B.9","ZK673.4",
             "srab-2","ada-2","F53H10.2","B0035.1","C08B11.3","die-1","F55H12.6","F26H11.2","H20J04.2","tag-350","F15E6.1","athp-1","Y116A8C.22","hmg-12","atf-7","C27D6.4","C34D1.5","C48E7.11","F17A9.3",
             "F23F12.9","F29G9.4","atf-6","mgl-2","F57B10.1","R07H5.10","xbp-1","atf-5","skn-1","T27F2.4","sknr-1","zip-3","W08E12.1","Y51H4A.4","Y75B8A.29","ZC376.7","ces-2","atf-2","lin-28",
             "cey-1","cey-2","cey-3","cey-4","sox-1","hmg-3","sox-3","F47G4.6","sox-2","hmg-4","egl-13 ","gei-3","W02D9.3","hmg-1.1","Y71H2AM.17","hmg-5","hmg-1.2",
             "unc-3","rnt-1","daf-19","C28H8.9","C26C6.1","aha-1","T24H10.7","C02F12.5","C07E3.5","C26E1.3","ceh-7","F34D6.2","F42G2.6","F45C12.2","F45C12.3","R03E1.4","T13C5.4","ZC376.4","R04A9.5","Y66A7A.5","ZK1193.5","ceh-48","ceh-49","ceh-38",
             "ceh-41","ceh-21","ceh-39","ceh-44","F54A5.1","C07E3.6","lin-39","egl-5","mab-5","ceh-16","pal-1","ceh-19","ceh-12","vab-7","R06F6.6","ceh-13","php-3","nob-1","Y80D3A.3","lim-7","ttx-3","mec-3",
             "lim-6","lin-11","lim-4","ceh-5","ceh-2","ceh-43","ceh-31","ceh-30","mls-2","C50H2.6","ceh-1","ceh-22","tab-1","ceh-27","ceh-24","ceh-28","pha-2","cog-1","vab-15","ceh-9","ceh-23",
             "unc-86","F28H6.2","ceh-6","ceh-18","Y38E10A.6","F45C12.15","C18B12.3","ceh-37","C40D2.4","ceh-17","unc-4","unc-42","R08B4.2","ceh-10","ttx-1","ceh-45","ceh-8","vab-3","pax-3","Y53C12C.1","unc-30",
             "ceh-34","ceh-33","ceh-35","ceh-32","irx-1","C49C3.5","ceh-40","ceh-20","unc-62","ceh-60","C12D12.5","spr-1","C24A1.2","C33A11.4","lin-1","ets-5","C50A2.4","C52B9.2","F19F10.1",
             "F19F10.5","ets-4","ast-1","hsf-1","zag-1","ZC123.3","dmd-4","mab-23","C34D1.1","C34D1.2","dmd-5","F13G11.1","K08B12.2","dmd-3","Y67D8A.3","T22H9.4","mab-3","end-1","egr-1","Y74C9A.4","nhr-10","nhr-23",
             "nhr-17","nhr-76","nhr-35","nhr-67","nhr-28","nhr-153","nhr-154","nhr-136","nhr-124","nhr-47","nhr-120","nhr-31","nhr-100","nhr-138","nhr-2","nhr-260","nhr-64","nhr-46","nhr-6","nhr-166","nhr-129","nhr-168",
             "nhr-173","nhr-19","nhr-121","nhr-262","nhr-116","nhr-175","daf-12","nhr-25","nhr-118","nhr-21","nhr-141","nhr-4","nhr-8","nhr-181","sex-1","nhr-37","nhr-142","nhr-111","nhr-191","unc-55",
             "fax-1","nhr-60","nhr-192","nhr-34","nhr-3","nhr-80","nhr-68","nhr-101","nhr-97","nhr-38","nhr-22","nhr-52","nhr-197","nhr-198","nhr-88","nhr-32","nhr-49","nhr-119","nhr-204",
             "nhr-205","nhr-206","nhr-207","nhr-208","nhr-209","nhr-1","nhr-104","nhr-270","nhr-14","nhr-40","nhr-279","nhr-26","nhr-66","nhr-273","nhr-16","nhr-127","nhr-69","nhr-225","nhr-226","nhr-61",
             "nhr-85","nhr-91","nhr-65","nhr-145","nhr-86","nhr-274","nhr-146","nhr-87","nhr-114","nhr-70","nhr-239","nhr-62","nhr-242","nhr-5","nhr-11","nhr-247","nhr-33","nhr-48","B0019.2","C05D10.1","C34F6.9","F16B12.6","T04G9.1","T10D4.6","T11A5.1","T23G5.6","Y39E4B.2","gak-1","C38D4.3",
             "Y18D10A.1","lin-15b","F09G2.9","F13C5.2","hmg-11","pqn-75","hlh-1","hlh-8","cky-1","hlh-27","hlh-25","hlh-14","hlh-12","hlh-15","hnd-1","hlh-16","hlh-17","F38C2.8","mxl-2","mxl-3","hlh-13","hlh-19","hlh-11","hlh-2","mdl-1","T01D3.2",
             "hlh-4","lin-32","hlh-6","mml-1","hlh-3","W02C12.3","Y105C5B.29","Y39A3CR.6","sbp-1","lin-22","ngn-1","hlh-10","hlh-26","hlh-28","hlh-29","ref-1","lpd-2","D1005.3","F17C11.1","zip-2","zip-4","zip-1","lfi-1",
             "F23B12.7","grh-1","F59H6.6","ZC204.2","ceh-14","C09G12.1","ceh-36","pop-1","lag-1","T05C1.4","C01G12.1","C06A5.4","C46F11.3","F36D1.1","F57G4.6","T07C12.11","Y106G6H.4","ZC416.1","H20J04.3","sma-4 ","sma-2","F21A10.2","pqn-47",
             "egl-38","F21D12.5","F48B9.5","pax-2","pax-1","R13.2","F26H9.2","Y51H4A.19","egl-44","wrm-1","sdc-2","bar-1","hmp-2","Y65B4BR.5","ham-1","F57A8.1","T28D9.9","Y73B3A.5","dac-1","2L52.1",
             "C01F6.9","C02F5.12","C06E2.1","pqn-21","C38D4.7","C46E10.8","C52E12.1","ztf-3","D1046.2","D2030.7","F11A10.2","lin-26","lir-1","lir-2","F21D5.9","F21G4.5","F23A7.6","F23C8.4","F27D4.6","F37B4.10","lir-3",
             "F40G9.14","F44D12.10","F45C12.1","F45H11.1","F49E8.2","F52B5.7","F53B7.2","F55B11.4","K05F1.5","K10D2.3","sea-2","K11D12.12","K11H3.4","K12H6.12","M04G12.4","M163.2","R02D3.7","R05D3.3","R07E5.5",
             "T06G6.5","zim-2","zim-3","him-8","T07G12.13","T08G5.7","T09F3.1","ztf-4","eea-1","T22C8.3","T22C8.4","T23F11.4","tlp-1","T24C4.7","W04B5.2","W04D2.4","Y37F4.6","ztf-5",
             "Y47G6A.7","Y48A6C.1","Y48A6C.3","Y48C3A.4","Y48E1B.7","Y51H1A.6","Y52B11A.9","Y53F4B.5","Y53H1A.2","Y54G2A.20","Y56A3A.18","Y57A10A.31","Y82E9BR.17","Y87G2A.3","ZK177.3","ZK185.1","ZK652.6","ZK686.4","F33H1.4","lin-36","Y48G8AL.10","spr-4",
             "lin-13","B0310.2","C06E1.8","C08G9.2","C55C2.1","F21A9.2","F36F12.8","F44E2.7","F47H4.1","sdc-1","F52E4.8","K09H9.7","K10B3.5","R10E4.11","R144.3","dnj-17","T09A5.12","T20F7.1","T20H4.2",
             "W03F9.2","ztf-6","Y17G7B.22","Y22D7AL.16","Y48G8AL.9","Y48G9A.11","Y56A3A.28","Y66D12A.12","Y6G8.3","Y71A12B.8","Y82E9BR.1","Y95B8A.7","ZK892.7","odd-1","C04F5.9","C10A4.8","C16A3.4","C27C12.2","odd-2","ref-2","unc-98","ztf-2","F13H6.1",
             "F26A10.2","pqm-1","F47E1.3","F53F8.1","mua-1","F56F11.3","F57C9.4","K11D2.4","T22C8.5","T27B1.2","W02D7.6","Y37E11B.1","Y39B6A.46","Y40B1A.4","Y59A8B.13","Y67H2A.10","Y73F8A.33","Y79H2A.3","ZC328.2","ZC395.8","ZK546.5","ZK686.5","F54C4.3",
             "C46E10.9","C52E12.6","che-1","F10B5.3","odd-3","ces-1","ztf-7","ztf-1","F56D1.1","H16D19.3","K02D7.2","egl-43","T07D10.3","ZK337.2","ZK867.1","C28G1.4","F37D6.2","pag-3","lsy-2","lsl-1","mep-1","R06C7.9",
             "lin-29","Y38H8A.5","tra-1","Y5F2A.4","B0250.4","C01B7.1","C09F5.3","C27A12.2","tag-146","F26F4.8","F53B3.1","F58G1.2","M03D4.4","T05G11.1","Y111B2A.10",
             "Y54E10BR.8","ehn-3","spr-3","F12E12.5","sem-4","F35H8.3","sma-9","R08E3.4","hbl-1","F39B2.1","Y55F3AM.14","B0336.3","C33H5.17","tag-331","T11G6.8","W05B10.2","W05H7.4","Y60A9.2","uaf-2","mex-6","oma-1","C34D10.2","C35D6.4","F27D4.4","moe-3","F38B7.1","F38C2.5",
             "F38C2.7","pos-1","K02H8.1","T02E1.3","T26A8.4","mex-5","mex-1","Y116A8C.17","Y116A8C.19","Y116A8C.20","pie-1","Y53G8AR.9","Y57G11C.25","Y60A9.3","oma-2","Y55F3AM.6","AC3.10","C17D12.1","C43H6.7","D2021.2","F09B12.2","F33D11.12","F59C6.2","H32C10.3","K02G10.1","M18.8","R13F6.5","T22E7.2",
             "Y39E4B.7","Y47H9C.2","ZK757.4","C26E6.2","peb-1","Y11D7A.12","Y11D7A.13","elt-1","egl-27","R13.1","gei-17","nhr-148","nhr-149","nhr-147","nhr-105","nhr-156","nhr-157","nhr-158",
             "nhr-159","srt-58","nhr-36","nhr-190","nhr-201","nhr-202","nhr-203","nhr-222","nhr-221","nhr-41","nhr-237","nhr-238","B0336.7","cdc-14","F49E10.5","crh-1","daf-3","tag-68","smad-1","sma-3","nfi-1","lin-48")


##############
# Checking HVGs between MSxa and MSxp
## Panel D
################
subset_ms_temp <- subset(ms_new, idents =c("MSxa","MSxp"))
levels(subset_ms_temp) <- c("MSxa","MSxp")
HVGs <- FindAllMarkers(subset_ms_temp, min.pct = 0.5,  assay = "RNA")
top.markers <- c()
ntops <-  200
for(n in unique(HVGs$cluster)){
  #n = 0
  marker.set <- subset(HVGs, cluster == n)
  #marker.set <- markers[["1"]]
  #head(marker.set, 5)
  #top.markers <- c(top.markers, rownames(marker.set)[1:ntops])  
  top.markers <- c(top.markers, marker.set$gene[1:ntops])  
}
top.markers <- unique(top.markers)
top.markers_TFs <- top.markers[top.markers %in% all_TFs]
#subset_ms_temp <- AverageExpression(subset_ms_temp, return.seurat = TRUE, verbose = TRUE)

#DoHeatmap(subset_ms_temp, features = all_TFs, label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()

top.markers_TFs <- top.markers[top.markers %in% all_TFs]
DoHeatmap(subset_ms_temp, features = top.markers[top.markers %in% all_TFs], label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()

#selecting HVGs
top.markers_TFs_ordered <- c("C34F6.9", "K09H9.7", "nhr-57", "F27D4.4",
                             "B0310.2",  "fkh-4", "ceh-36", "ceh-49", "egl-18", "ceh-36",  "fkh-3",
                             'ngn-1', "tbx-11", 'hnd-1','pha-4', "ceh-51",  'nhr-67')

DoHeatmap(subset_ms_temp, features = c(top.markers_TFs_ordered, 'ceh-51'), label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()

##########
## selecting clusters and predicting pseudotime
#MSx-MSxa 
sling_ms_early <- subset(ms_new, idents = c('MSx', "MSxa", "MSxaa"))
levels(sling_ms_early) <- c('MSx', "MSxa", "MSxaa")

sling_ms_early <- NormalizeData(object = sling_ms_early, verbose = T)
sling_ms_early <- FindVariableFeatures(sling_ms_early, nfeatures = 2000)
sling_ms_early <- ScaleData(sling_ms_early, features = rownames(sling_ms_early))

sling_ms_early <- RunPCA(object = sling_ms_early, features = VariableFeatures(sling_ms_early), verbose = FALSE)

nb.pcs = 20; n.neighbors = 10; min.dist = 0.3;
sling_ms_early <- RunUMAP(object = sling_ms_early, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)

DimPlot(sling_ms_early, reduction = "umap", group.by = 'ident', label = T)
DimPlot(sling_ms_early, reduction = "umap", group.by = 'timingEst', label = T) 

colnames(subset(ms_new, idents = c('MSpaaapp')))

DimPlot(sling_ms_early, reduction = "pca", group.by = 'ident', label = T)

sling_ms_early <- as.SingleCellExperiment(sling_ms_early)


## Visualize all lineages with clusters
# Not included in the paper, but it's nice to see how exactly the pseudotime axis is alligned
sling <- slingshot(sling_ms_early, clusterLabels = sling_ms_early@colData@listData[["manual.annot.ids"]], reducedDim = 'UMAP', start.clus = 'MSx', end.clus = c('MSxaa'))

pallete_length <- length(levels(sling$ident))
colmathcer <- data.frame(levels(sling$ident), brewer.pal(pallete_length, name = "Dark2"))
col_vector <- brewer.pal(pallete_length, name = "Dark2")[match(sling$ident, colmathcer[,1])]

plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
lines(SlingshotDataSet(sling), lwd=2, type = 'lineages')
lines(SlingshotDataSet(sling), lwd=2)
legend('bottomright',legend = levels(sling$ident), col = brewer.pal(pallete_length, name = "Dark2"), fill= brewer.pal(pallete_length, name = "Dark2"))

## Visualize all lineages with timingEst
# Not included in the paper, but it's nice to see how exactly the pseudotime axis is alligned
sling$timingEst <- droplevels(sling$timingEst)

pallete_length <- length(levels(sling$timingEst))
colmathcer <- data.frame(levels(sling$timingEst), plasma(pallete_length))
col_vector <- plasma(pallete_length)[match(sling$timingEst, colmathcer[,1])]
plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
legend('bottomleft',legend = levels(sling$timingEst), fill= plasma(pallete_length))
lines(SlingshotDataSet(sling), lwd=2)

#Visualise clusters
sling_ms_early$slingPseudotime_1 <- slingPseudotime(sling, na=FALSE)[,1]
sling_ms_early <- as.Seurat(sling_ms_early)
sling_ms_early <- subset(ms_new, idents = c('MSx', "MSxa", "MSxaa"))
levels(sling_ms_early) <- c('MSx', "MSxa", "MSxaa")

DoHeatmap(sling_ms_early, features = c( top.markers_TFs_ordered,"ceh-51"), cells = colnames(sling)[order(sling$slingPseudotime_1, decreasing = T)])

##########
## selecting clusters and predicting pseudotime
#MSx-MSxp
sling_ms_early <- subset(ms_new, idents = c('MSx', "MSxp", "MSxpp"))
levels(sling_ms_early) <- c('MSx', "MSxp", "MSxpp")
sling_ms_early <- NormalizeData(object = sling_ms_early, verbose = T)
sling_ms_early <- FindVariableFeatures(sling_ms_early, nfeatures = 2000)
sling_ms_early <- ScaleData(sling_ms_early, features = rownames(sling_ms_early))
sling_ms_early <- RunPCA(object = sling_ms_early, features = VariableFeatures(sling_ms_early), verbose = FALSE)
nb.pcs = 20; n.neighbors = 10; min.dist = 0.3;
sling_ms_early <- RunUMAP(object = sling_ms_early, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(sling_ms_early, reduction = "umap", group.by = 'ident', label = T)
DimPlot(sling_ms_early, reduction = "umap", group.by = 'timingEst', label = T)
DimPlot(sling_ms_early, reduction = "pca", group.by = 'ident', label = T)
sling_ms_early <- as.SingleCellExperiment(sling_ms_early)
## Visualize all lineages with clusters
# Not included in the paper, but it's nice to see how exactly the pseudotime axis is alligned
sling <- slingshot(sling_ms_early, clusterLabels = sling_ms_early@colData@listData[["manual.annot.ids"]], reducedDim = 'UMAP', start.clus = 'MSx', end.clus = c('MSxpp'))

pallete_length <- length(levels(sling$ident))
colmathcer <- data.frame(levels(sling$ident), brewer.pal(pallete_length, name = "Dark2"))
col_vector <- brewer.pal(pallete_length, name = "Dark2")[match(sling$ident, colmathcer[,1])]
plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
lines(SlingshotDataSet(sling), lwd=2, type = 'lineages')
lines(SlingshotDataSet(sling), lwd=2)
legend('bottomright',legend = levels(sling$ident), col = brewer.pal(pallete_length, name = "Dark2"), fill= brewer.pal(pallete_length, name = "Dark2"))

## Visualize all lineages with timingEst
# Not included in the paper, but it's nice to see how exactly the pseudotime axis is alligned
sling$timingEst <- droplevels(sling$timingEst)
pallete_length <- length(levels(sling$timingEst))
colmathcer <- data.frame(levels(sling$timingEst), plasma(pallete_length))
col_vector <- plasma(pallete_length)[match(sling$timingEst, colmathcer[,1])]
plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
legend('bottomleft',legend = levels(sling$timingEst), fill= plasma(pallete_length))
lines(SlingshotDataSet(sling), lwd=2)

#Visualise clusters
sling_ms_early$slingPseudotime_1 <- slingPseudotime(sling, na=FALSE)[,1]
sling_ms_early <- as.Seurat(sling_ms_early)
sling_ms_early <- subset(ms_new, idents = c( 'MSx',"MSxp", "MSxpp"))
levels(sling_ms_early) <- c('MSx',"MSxp",  "MSxpp")
DoHeatmap(sling_ms_early, features =  c(top.markers_TFs_ordered,"ceh-51"), cells = colnames(sling)[order(sling$slingPseudotime_1, decreasing = F)])



##############
# Checking HVGs between MSxa and MSxp
## Panel E
################

subset_ms_temp <- subset(ms_new, idents =c("MSxaa","MSxap"))
levels(subset_ms_temp) <- c("MSxaa","MSxap")


HVGs <- FindAllMarkers(subset_ms_temp, min.pct = 0.5,  assay = "RNA")
top.markers <- c()
ntops <-  500
for(n in unique(HVGs$cluster)){
  #n = 0
  marker.set <- subset(HVGs, cluster == n)
  #marker.set <- markers[["1"]]
  #head(marker.set, 5)
  #top.markers <- c(top.markers, rownames(marker.set)[1:ntops])  
  top.markers <- c(top.markers, marker.set$gene[1:ntops])  
}
top.markers <- unique(top.markers)
top.markers_TFs <- top.markers[top.markers %in% all_TFs]
#subset_ms_temp <- AverageExpression(subset_ms_temp, return.seurat = TRUE, verbose = TRUE)

#DoHeatmap(subset_ms_temp, features = all_TFs, label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()

top.markers_TFs <- top.markers[top.markers %in% all_TFs]
DoHeatmap(subset_ms_temp, features = top.markers[top.markers %in% all_TFs], label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()


top.markers_TFs_ordered <- c('Y53G8AR.9', 'ceh-39', 'egl-43', 'nhr-78','spr-3',  'odd-2',
                             'hnd-1','pha-4', "ceh-51",  'nhr-67')
# to see the HVGs
DoHeatmap(subset_ms_temp, features = c(top.markers_TFs_ordered, 'ceh-51'), label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()

##########
## selecting clusters and predicting pseudotime

#MSax-MSxaa
sling_ms_early <- subset(ms_new, idents = c("MSxaa","MSxa","MSxaaa"))
levels(sling_ms_early) <- c("MSxaa","MSxa","MSxaaa")

sling_ms_early <- NormalizeData(object = sling_ms_early, verbose = T)
sling_ms_early <- FindVariableFeatures(sling_ms_early, nfeatures = 2000)
sling_ms_early <- ScaleData(sling_ms_early, features = rownames(sling_ms_early))

sling_ms_early <- RunPCA(object = sling_ms_early, features = VariableFeatures(sling_ms_early), verbose = FALSE)

nb.pcs = 20; n.neighbors = 10; min.dist = 0.3;
sling_ms_early <- RunUMAP(object = sling_ms_early, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)

DimPlot(sling_ms_early, reduction = "umap", group.by = 'ident', label = T)
DimPlot(sling_ms_early, reduction = "umap", group.by = 'timingEst', label = T)

DimPlot(sling_ms_early, reduction = "pca", group.by = 'ident', label = T)

sling_ms_early <- as.SingleCellExperiment(sling_ms_early)
## Visualize all lineages with clusters

sling <- slingshot(sling_ms_early, clusterLabels = sling_ms_early@colData@listData[["manual.annot.ids"]], reducedDim = 'UMAP', start.clus = 'MSxa', end.clus = c('MSxaaa'))

library(RColorBrewer)
library(viridis)
pallete_length <- length(levels(sling$ident))
colmathcer <- data.frame(levels(sling$ident), brewer.pal(pallete_length, name = "Dark2"))
col_vector <- brewer.pal(pallete_length, name = "Dark2")[match(sling$ident, colmathcer[,1])]

plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
lines(SlingshotDataSet(sling), lwd=2, type = 'lineages')
lines(SlingshotDataSet(sling), lwd=2)
legend('bottomright',legend = levels(sling$ident), col = brewer.pal(pallete_length, name = "Dark2"), fill= brewer.pal(pallete_length, name = "Dark2"))

## Visualize all lineages with timingEst
sling$timingEst <- droplevels(sling$timingEst)
pallete_length <- length(levels(sling$timingEst))
colmathcer <- data.frame(levels(sling$timingEst), plasma(pallete_length))
col_vector <- plasma(pallete_length)[match(sling$timingEst, colmathcer[,1])]
plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
legend('bottomleft',legend = levels(sling$timingEst), fill= plasma(pallete_length))
lines(SlingshotDataSet(sling), lwd=2)


#Visualise clusters
sling_ms_early$slingPseudotime_1 <- slingPseudotime(sling, na=FALSE)[,1]
sling_ms_early <- as.Seurat(sling_ms_early)
sling_ms_early <- subset(ms_new, idents = c("MSxa","MSxaa", 'MSxaaa'))
levels(sling_ms_early) <- c("MSxa","MSxaa",'MSxaaa')
DoHeatmap(sling_ms_early, features = c( top.markers_TFs_ordered), cells = colnames(sling)[order(sling$slingPseudotime_1, decreasing = F)])




##########
## selecting clusters and predicting pseudotime
#MSxa-MSxap
sling_ms_early <- subset(ms_new, idents = c("MSxap","MSxa","MSxapp"))
levels(sling_ms_early) <- c("MSxap","MSxa","MSxapp")
sling_ms_early <- NormalizeData(object = sling_ms_early, verbose = T)
sling_ms_early <- FindVariableFeatures(sling_ms_early, nfeatures = 2000)
sling_ms_early <- ScaleData(sling_ms_early, features = rownames(sling_ms_early))
sling_ms_early <- RunPCA(object = sling_ms_early, features = VariableFeatures(sling_ms_early), verbose = FALSE)

nb.pcs = 20; n.neighbors = 10; min.dist = 0.3;
sling_ms_early <- RunUMAP(object = sling_ms_early, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(sling_ms_early, reduction = "umap", group.by = 'ident', label = T)
DimPlot(sling_ms_early, reduction = "umap", group.by = 'timingEst', label = T)
DimPlot(sling_ms_early, reduction = "pca", group.by = 'ident', label = T)
sling_ms_early <- as.SingleCellExperiment(sling_ms_early)
## Visualize all lineages with clusters
sling <- slingshot(sling_ms_early, clusterLabels = sling_ms_early@colData@listData[["manual.annot.ids"]], reducedDim = 'UMAP', start.clus = 'MSxa', end.clus = c('MSxapp'))
library(RColorBrewer)
library(viridis)
pallete_length <- length(levels(sling$ident))
colmathcer <- data.frame(levels(sling$ident), brewer.pal(pallete_length, name = "Dark2"))
col_vector <- brewer.pal(pallete_length, name = "Dark2")[match(sling$ident, colmathcer[,1])]
plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
lines(SlingshotDataSet(sling), lwd=2, type = 'lineages')
lines(SlingshotDataSet(sling), lwd=2)
legend('bottomright',legend = levels(sling$ident), col = brewer.pal(pallete_length, name = "Dark2"), fill= brewer.pal(pallete_length, name = "Dark2"))

## Visualize all lineages with timingEst
sling$timingEst <- droplevels(sling$timingEst)
pallete_length <- length(levels(sling$timingEst))
colmathcer <- data.frame(levels(sling$timingEst), plasma(pallete_length))
col_vector <- plasma(pallete_length)[match(sling$timingEst, colmathcer[,1])]
plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
legend('bottomleft',legend = levels(sling$timingEst), fill= plasma(pallete_length))
lines(SlingshotDataSet(sling), lwd=2)


#Visualise clusters
sling_ms_early$slingPseudotime_1 <- slingPseudotime(sling, na=FALSE)[,1]
sling_ms_early <- as.Seurat(sling_ms_early)
sling_ms_early <- subset(ms_new, idents = c("MSxa","MSxap","MSxapp"))
levels(sling_ms_early) <- c("MSxa","MSxap","MSxapp")
DoHeatmap(sling_ms_early, features = c( top.markers_TFs_ordered), cells = colnames(sling)[order(sling$slingPseudotime_1, decreasing = F)])


##############
# Checking HVGs between MSxa and MSxp
## Panel F
################
subset_ms_temp <- subset(ms_new, idents =c("MSxapa","MSxapp"))
levels(subset_ms_temp) <- c("MSxapa","MSxapp")

HVGs <- FindAllMarkers(subset_ms_temp, min.pct = 0.5,  assay = "RNA")
top.markers <- c()
ntops <-  500
for(n in unique(HVGs$cluster)){
  #n = 0
  marker.set <- subset(HVGs, cluster == n)
  #marker.set <- markers[["1"]]
  #head(marker.set, 5)
  #top.markers <- c(top.markers, rownames(marker.set)[1:ntops])  
  top.markers <- c(top.markers, marker.set$gene[1:ntops])  
}
top.markers <- unique(top.markers)
top.markers_TFs <- top.markers[top.markers %in% all_TFs]
#subset_ms_temp <- AverageExpression(subset_ms_temp, return.seurat = TRUE, verbose = TRUE)
#DoHeatmap(subset_ms_temp, features = all_TFs, label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()

top.markers_TFs <- top.markers[top.markers %in% all_TFs]
DoHeatmap(subset_ms_temp, features = top.markers[top.markers %in% all_TFs], label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()


top.markers_TFs_ordered <- c('dmd-4', 'irx-1', 'fkh-7', 'ceh-34', 'zip-2', 'tbx-7', 'ham-1', 'egl-44', 'spb-1', 'tab-1', "unc-62",
                             'hnd-1','pha-4', "ceh-51",  'nhr-67')

DoHeatmap(subset_ms_temp, features = c(top.markers_TFs_ordered, 'ceh-51'), label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()

############
### selecting individual branches for plotting
##########

"MSxap" "MSxapa"
sling_ms_early <- subset(ms_new, idents = c("MSxap","MSxapa","MSxapaa"))
levels(sling_ms_early) <- c("MSxap","MSxapa","MSxapaa")

sling_ms_early <- NormalizeData(object = sling_ms_early, verbose = T)
sling_ms_early <- FindVariableFeatures(sling_ms_early, nfeatures = 2000)
sling_ms_early <- ScaleData(sling_ms_early, features = rownames(sling_ms_early))
sling_ms_early <- RunPCA(object = sling_ms_early, features = VariableFeatures(sling_ms_early), verbose = FALSE)

nb.pcs = 20; n.neighbors = 10; min.dist = 0.3;
sling_ms_early <- RunUMAP(object = sling_ms_early, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(sling_ms_early, reduction = "umap", group.by = 'ident', label = T)
DimPlot(sling_ms_early, reduction = "umap", group.by = 'timingEst', label = T)
DimPlot(sling_ms_early, reduction = "pca", group.by = 'ident', label = T)
sling_ms_early <- as.SingleCellExperiment(sling_ms_early)
## Visualize all lineages with clusters
sling <- slingshot(sling_ms_early, clusterLabels = sling_ms_early@colData@listData[["manual.annot.ids"]], reducedDim = 'UMAP', start.clus = 'MSxap', end.clus = c('MSxapaa'))
library(RColorBrewer)
library(viridis)
pallete_length <- length(levels(sling$ident))
colmathcer <- data.frame(levels(sling$ident), brewer.pal(pallete_length, name = "Dark2"))
col_vector <- brewer.pal(pallete_length, name = "Dark2")[match(sling$ident, colmathcer[,1])]
plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
lines(SlingshotDataSet(sling), lwd=2, type = 'lineages')
lines(SlingshotDataSet(sling), lwd=2)
legend('bottomright',legend = levels(sling$ident), col = brewer.pal(pallete_length, name = "Dark2"), fill= brewer.pal(pallete_length, name = "Dark2"))

## Visualize all lineages with timingEst
sling$timingEst <- droplevels(sling$timingEst)
pallete_length <- length(levels(sling$timingEst))
colmathcer <- data.frame(levels(sling$timingEst), plasma(pallete_length))
col_vector <- plasma(pallete_length)[match(sling$timingEst, colmathcer[,1])]
plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
legend('bottomleft',legend = levels(sling$timingEst), fill= plasma(pallete_length))
lines(SlingshotDataSet(sling), lwd=2)
#Visualise clusters
sling_ms_early$slingPseudotime_1 <- slingPseudotime(sling, na=FALSE)[,1]
sling_ms_early <- as.Seurat(sling_ms_early)
sling_ms_early <- subset(ms_new, idents = c("MSxap","MSxapa","MSxapaa"))
levels(sling_ms_early) <- c("MSxap","MSxapa","MSxapaa")
DoHeatmap(sling_ms_early, features = c( top.markers_TFs_ordered), cells = colnames(sling)[order(sling$slingPseudotime_1, decreasing = F)])


#MSxap-MSxapp
sling_ms_early <- subset(ms_new, idents = c("MSxap","MSxapp","MSxappp"))
levels(sling_ms_early) <- c("MSxap","MSxapp","MSxappp")
sling_ms_early <- NormalizeData(object = sling_ms_early, verbose = T)
sling_ms_early <- FindVariableFeatures(sling_ms_early, nfeatures = 2000)
sling_ms_early <- ScaleData(sling_ms_early, features = rownames(sling_ms_early))
sling_ms_early <- RunPCA(object = sling_ms_early, features = VariableFeatures(sling_ms_early), verbose = FALSE)
nb.pcs = 20; n.neighbors = 10; min.dist = 0.3;
sling_ms_early <- RunUMAP(object = sling_ms_early, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
DimPlot(sling_ms_early, reduction = "umap", group.by = 'ident', label = T)
DimPlot(sling_ms_early, reduction = "umap", group.by = 'timingEst', label = T)
DimPlot(sling_ms_early, reduction = "pca", group.by = 'ident', label = T)
sling_ms_early <- as.SingleCellExperiment(sling_ms_early)
## Visualize all lineages with clusters
sling <- slingshot(sling_ms_early, clusterLabels = sling_ms_early@colData@listData[["manual.annot.ids"]], reducedDim = 'UMAP', start.clus = 'MSxap', end.clus = c('MSxappp'))
library(RColorBrewer)
library(viridis)
pallete_length <- length(levels(sling$ident))
colmathcer <- data.frame(levels(sling$ident), brewer.pal(pallete_length, name = "Dark2"))
col_vector <- brewer.pal(pallete_length, name = "Dark2")[match(sling$ident, colmathcer[,1])]
plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
lines(SlingshotDataSet(sling), lwd=2, type = 'lineages')
lines(SlingshotDataSet(sling), lwd=2)
legend('bottomright',legend = levels(sling$ident), col = brewer.pal(pallete_length, name = "Dark2"), fill= brewer.pal(pallete_length, name = "Dark2"))

## Visualize all lineages with timingEst
sling$timingEst <- droplevels(sling$timingEst)
pallete_length <- length(levels(sling$timingEst))
colmathcer <- data.frame(levels(sling$timingEst), plasma(pallete_length))
col_vector <- plasma(pallete_length)[match(sling$timingEst, colmathcer[,1])]
plot(reducedDims(sling)$UMAP, col = col_vector, pch=16, asp = 1)
legend('bottomleft',legend = levels(sling$timingEst), fill= plasma(pallete_length))
lines(SlingshotDataSet(sling), lwd=2)


#Visualise clusters
sling_ms_early$slingPseudotime_1 <- slingPseudotime(sling, na=FALSE)[,1]
sling_ms_early <- as.Seurat(sling_ms_early)
sling_ms_early <- subset(ms_new, idents = c("MSxap","MSxapp","MSxappp"))
levels(sling_ms_early) <- c("MSxap","MSxapp","MSxappp")
DoHeatmap(sling_ms_early, features = c( top.markers_TFs_ordered), cells = colnames(sling)[order(sling$slingPseudotime_1, decreasing = F)])




###################
## FIGURE 3a
###########


ms_subset_hnd1 <- subset(ms, idents =c("MSxa","MSxp",'MSxpa', "MSxpp", "MSxap", "MSxaa" , "MSxapa","MSxapp" ))
levels(ms_subset_hnd1) <- c("MSxa","MSxp",'MSxpa', "MSxpp", "MSxap", "MSxaa" , "MSxapa","MSxapp" )
DoHeatmap(ms_subset_hnd1, features = c('pha-4',  'hnd-1', 'hlh-1', 'pat-9', 'nhr-67'), label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()
DoHeatmap(ms_subset_hnd1, features = c('pha-4',  'hnd-1'), label = T, draw.lines = T, lines.width = 6, angle = 90, group.by = "ident") + NoLegend()
Idents(ms_subset_hnd1)
ms_subset_hnd1$manual.annot.ids
VlnPlot(ms_subset_hnd1, "hnd-1", group.by = 'manual.annot.ids') +
  geom_boxplot()

for_bxplpt <- data.frame(FetchData(object = ms_subset_hnd1, vars = c('manual.annot.ids', 'pha-4', 'hnd-1')))
for_bxplpt <- pivot_longer(for_bxplpt, cols = 2:3, names_to = 'gene', values_to = 'expression_level' )
p <- ggboxplot(for_bxplpt,x ="manual.annot.ids", y="expression_level", color= 'gene',  order = c("MSxa","MSxp",'MSxpa', "MSxpp", "MSxap", "MSxaa" , "MSxapa","MSxapp" ))
p


############################
### FIGURE 4 ####

## Plot notch components for MSxxx stage
ms_subset_hnd1 <- subset(ms, idents =c('MSxaa', "MSxap", "MSxpa", "MSxpp"))
levels(ms_subset_hnd1) <- c('MSxaa', "MSxap", "MSxpa", "MSxpp")
DoHeatmap(ms_subset_hnd1, features = c('pha-4',  'hnd-1', 'lag-2', 'glp-1', 'lin-12', 'lag-1', 'ref-1'), label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()



############################
### FIGURE 5 ####
# panel a

ms_new_termianl_BWMs <- subset(ms_new, idents = c('MSxappppx','MSxapppax',
                                                  'MSxppppp', 'MSxppppa', 
                                                  'MSxpppap', 'MSxpppaa',
                                                  'MSxppapp','MSxpapap',
                                                  'MSxpappp','MSxpappa',
                                                  'MSxpaaap'))

ms_new_termianl_BWMs$ident <- ms_new_termianl_BWMs@active.ident
levels(ms_new_termianl_BWMs) <- c('MSxappppx','MSxapppax',
                                  'MSxppppp', 'MSxppppa', 
                                  'MSxpppap', 'MSxpppaa',
                                  'MSxppapp','MSxpapap',
                                  'MSxpappp','MSxpappa',
                                  'MSxpaaap')

selected_genes <- c("maph-1.3","maph-1.2","ensh-1","maph-1.2","egl-13","zig-6","unc-52",  
                    "pde-1","far-1","lbp-3","dhs-29","W03G11.3","cah-5","bgal-1","R04B5.6", 
                    "gana-1","pck-1","mrps-7","unc-49","M153.4","F53C11.9","srp-1","fbxa-220", 
                    "srh-163","cwn-1","cwn-2","unc-129","egl-15","sfrp-1","clec-264","lin-39", 
                    "eya-1","E04F6.10","F47B10.8","nra-2","sre-23","str-148","lron-12","T28D6.5", 
                    "K09G1.1","C15C7.6","frpr-8","kvs-5","ttr-16","D1086.12","strl-1","B0379.1", 
                    "ZC374.2","F52H2.7","lntl-1","F35B12.10")
DoHeatmap(ms_new_termianl_BWMs, features = selected_genes, draw.lines = T, lines.width = 2,label = T, angle = 90, group.by = "ident")


# panel b 
# REQUIRES INSTALLING EQTL DATASET
# Used data :
# https://github.com/eyalbenda/worm_sceQTL/blob/master/study_data/backgroundCorrectedMonocle3_V5.rds
# Load RDS first as eqtls

eqtls <- readRDS('~/Downloads/backgroundCorrectedMonocle3_V5.rds')

eqtlsbwms <- eqtls[,eqtls@colData@listData[["broad_tissue"]]=="Body Wall Muscle"]
eqtlsbwms_selected <- choose_cells(eqtlsbwms) # remove some scattered cells
plot_cells(eqtlsbwms_selected, color_cells_by = 'cluster')


marker_test_res <- top_markers(eqtlsbwms_selected, group_cells_by="cluster", cores=8)

selected_markers <- unique(marker_test_res$gene_short_name)


top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_markers_mine <- c('lin-39', 'clec-264', 'lbp-2', 'mig-1', 'mig-6', 'cwn-1', 'cwn-2', 
                               'B0379.1', 'tlp-1', 'maph-1.2', 'T08H10.1', "sul-2", 'gana-1', 'zig-6' ,
                               "col-111", 'col-118', "F41D9.1", "Y47D3B.6",'zmp-1','C03B1.1','ceh-13','ham-1')

top_specific_markers_mine <- c('lin-39', 'clec-264', 'lbp-2', 'mig-1', 'mig-6', 'cwn-1', 'cwn-2', 
                               'B0379.1', 'tlp-1',  "sul-2", 'gana-1', 'zig-6' ,
                               "col-111", 'col-118', 'zmp-1','ceh-13','ham-1')

plot_cells(eqtlsbwms_selected, genes=top_specific_markers_mine, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)


plot_cells(eqtlsbwms_selected, genes=selected_markers, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

sbst_clls <- c("MSxapppax","MSxappppx","MSxpaaap","MSxpapap",  #,"MSxppppx","MSxpaap.MSppaapp.like", "MSxapppa","MSxapppp",
               "MSxpappa","MSxpappp", "MSxppapp","MSxpppaa","MSxpppap", "MSxppppa", "MSxppppp")
ms_sbst <- subset(ms, idents = sbst_clls)
ms_sbst@active.ident
levels(ms_sbst) <- sbst_clls 
HVGs <- FindAllMarkers(ms_sbst, min.pct = 0.35)

ntops <- 10 
top.markers <- c()

for(n in unique(HVGs$cluster)){
  #n = 0
  marker.set <- subset(HVGs, cluster == n)
  #marker.set <- markers[["1"]]
  #head(marker.set, 5)
  #top.markers <- c(top.markers, rownames(marker.set)[1:ntops])  
  top.markers <- c(top.markers, marker.set$gene[1:ntops])  
}

top.markers = unique(top.markers)

selected_markers_from_terminal_cells <- top.markers 
exclude_genes <- c('nid-1', 'maph-1.3', 'eya-1', 'F59D12.1', 'srbc-5', 'B0554.4', "F01G12.1", 'C15C7.6', 'gst-7', 'M153.4', 'frpr-8', 'C56G2.4', 'maph-1.2', 'F09F9.4', 'pnc-1',
                   'T25G12.11', 'sto-1', 'F36H5.4', 'W10D9.1', 'C08F1.10', 'Y41D4A.3', 'T28D6.5', 'pxn-2', 'pat-10', 'lev-11', 'sre-23', 'egl-13', 'F35B12.10',
                   'str-148', 'egl-15', 'arrd-15', 'ctbp-1', 'unc-24', ' W04A8.4', 'fkb-4', "F15B10.3")

selected_markers_from_terminal_cells <- selected_markers_from_terminal_cells[!selected_markers_from_terminal_cells%in%exclude_genes]

#This is used for the supplementary figure
plot_cells(eqtlsbwms_selected, genes= unique(c(top_specific_markers_mine, selected_markers_from_terminal_cells)), 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)

#This is used for the supplementary figure
DoHeatmap(ms_sbst, features = unique(c(top_specific_markers_mine, selected_markers_from_terminal_cells)), label = T, draw.lines = T, lines.width = 4, angle = 90, group.by = "ident") + NoLegend() #+ scale_fill_gradientn(colours = c("#000004",viridis(3)))


#This is used for the main figure figure
small_genes_panel <- c( 'unc-54','lbp-3', 'zig-6', 'clec-264')
plot_cells(eqtlsbwms_selected, genes= small_genes_panel, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)




#############################################################
## Supplemental figures
#############################################################

#######
## Supplemental figure 1
#######

## Requires data set from John Murray's paper.
# https://github.com/qinzhu/VisCello.celegans

eset <- readRDS("###########/VisCello.celegans/inst/app/data/eset.rds")
cells_of_interest <- c(grep(pattern = 'MSx.{1,100}$',unique(eset@phenoData@data[["lineage"]]), value = T))
idents <-  unique(c("MSx", "MSxa", "MSxp",
                    "MSxaa", "MSxap", "MSxpa", "MSxpp", 
                    'MSxaaa', 'MSxaap', 'MSxapa', 'MSxapp', 
                    'MSxpaa', 'MSxpap', 'MSxppa', 'MSxppp',
                    'MSpaaaa', 'MSaaaap',"MSpaaap", 'MSxaapa', 'MSaaapp',
                    'MSpaaaa','MSpaaap','MSaaaap','MSxaapa','MSaaapp','MSxapaa','MSxapap','MSaappa', 'MSpappa', 'MSxappp',
                    'MSxpaaa', 'MSxpapa', 'MSxpapp', 'MSxppaa', 'MSxppap', 'MSxpppa', 'MSxpppp',
                    'MSpaaaaa', 'MSpaaaap','MSpaaapa','MSpaaapp', 'MSpaapaa', 'MSxaapap', 'MSaaaaaa', 'MSaaaaap', 'MSaaaapa', 'MSaaaapp', 'MSaaapaa',
                    'MSxapaaa', 'MSxapaap', 'MSxapapp', 'MSxapppa', 'MSxapppp',
                    'MSxpaaaa', 'MSxpaaap', 'MSxpapaa', 'MSxpapap', 'MSxpappa', 'MSxpappp', 'MSxppaaa', 'MSxppaap','MSxppapp',
                    'MSxpppaa', 'MSxpppap', 'MSxppppa', 'MSxppppp',
                    "MSpappax","MSxaapapx", "MSaaaaaaa", "MSapaapp", "MSaaaaapa", "MSxapppax", "MSxappppx","MSxapaapa","MSaaappp"))
setdiff(idents, cells_of_interest)
setdiff(cells_of_interest,idents)

cell_list <- intersect(idents, cells_of_interest)
cell_order <- unique(c("MSx", "MSxa", "MSxp",
                       "MSxaa", "MSxap", "MSxpa", "MSxpp", 
                       'MSxaaa', 'MSxaap', 'MSxapa', 'MSxapp', 
                       'MSxpaa', 'MSxpap', 'MSxppa', 'MSxppp',
                       'MSpaaaa', 'MSaaaap',"MSpaaap", 'MSxaapa', 'MSaaapp',
                       'MSpaaaa','MSpaaap','MSaaaap','MSxaapa','MSaaapp','MSxapaa','MSxapap','MSaappa', 'MSpappa', 'MSxappp',
                       'MSxpaaa', 'MSxpapa', 'MSxpapp', 'MSxppaa', 'MSxppap', 'MSxpppa', 'MSxpppp',
                       'MSpaaaaa', 'MSpaaaap','MSpaaapa','MSpaaapp', 'MSpaapaa', 'MSxaapap', 'MSaaaaaa', 'MSaaaaap', 'MSaaaapa', 'MSaaaapp', 'MSaaapaa',
                       'MSxapaaa', 'MSxapaap', 'MSxapapp', 'MSxapppa', 'MSxapppp',
                       'MSxpaaaa', 'MSxpaaap', 'MSxpapaa', 'MSxpapap', 'MSxpappa', 'MSxpappp', 'MSxppaaa', 'MSxppaap','MSxppapp',
                       'MSxpppaa', 'MSxpppap', 'MSxppppa', 'MSxppppp',
                       "MSpappax","MSxaapapx", "MSaaaaaaa", "MSapaapp", "MSaaaaapa", "MSxapppax", "MSxappppx","MSxapaapa","MSaaappp"))

cell_order%in%cell_list
cell_list%in%cell_order

cell_count_df <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c('dataset_owner','cell_name', 'cell_id', 'n_features'))))

for (cell_name in cell_list){
  print(cell_name)
  temp_df <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c('dataset_owner','cell_name', 'cell_id', 'n_features'))))
  cell_ids <- which(eset@phenoData@data[["lineage"]] %in% cell_name)
  for (ID in cell_ids){
    temp_n_feature <- sum(eset@assayData[["norm_exprs"]][,ID] != 0)
    temp_df <- rbind.data.frame( list('JM', cell_name, ID, temp_n_feature), temp_df, make.row.names = F)
  }
  colnames(temp_df) <- c('dataset_owner','cell_name', 'cell_id', 'n_features')
  colnames(cell_count_df) <- c('dataset_owner','cell_name', 'cell_id', 'n_features')
  cell_count_df <- rbind.data.frame( temp_df,cell_count_df, make.row.names = F)
}

jm_data <- cell_count_df


#cell_count_df <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c('dataset_owner','cell_name', 'cell_id', 'n_features'))))

head(ms@meta.data[["nFeature_SCT"]][4037])
cell_count_df <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c('dataset_owner','cell_name', 'cell_id', 'n_features'))))

for (cell_name in cell_list){
  print(cell_name)
  temp_df <- data.frame(matrix(ncol=4,nrow=0, dimnames=list(NULL, c('dataset_owner','cell_name', 'cell_id', 'n_features'))))
  cell_ids <- which(ms@active.ident %in% cell_name)
  for (ID in cell_ids){
    temp_n_feature <- ms@meta.data[["nFeature_SCT"]][ID]
    temp_df <- rbind.data.frame( list('LC', cell_name, ID, temp_n_feature), temp_df, make.row.names = F)
  }
  colnames(temp_df) <- c('dataset_owner','cell_name', 'cell_id', 'n_features')
  colnames(cell_count_df) <- c('dataset_owner','cell_name', 'cell_id', 'n_features')
  cell_count_df <- rbind.data.frame( temp_df,cell_count_df, make.row.names = F)
}


lc_data <- cell_count_df

library(ggplot2)


ggplot(rbind(jm_data, lc_data), aes(x = factor(cell_name, level = cell_order[cell_order%in%cell_list]), y = n_features)) +
  theme_linedraw()+
  geom_boxplot(aes(color = dataset_owner), position = position_dodge(0.8)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                          colour = "lightgrey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_size(guide="none")+
  scale_color_discrete(name="Data sets:",
                       breaks=c("JM", "LC"),
                       labels=c("Packer et. al. 2019", "Current study"))+
  labs(title = "Number of genes detected per cell identity", x = "Cell identity name", y = "Number of genes detected")



ggplot(rbind(jm_data, lc_data), aes(x = dataset_owner, y = n_features)) +
  geom_boxplot(aes(fill = dataset_owner), position = position_dodge(0.8)) + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(size = 0.5, linetype = 'solid',
                                          colour = "lightgrey"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightgrey"))+
  scale_fill_discrete(name="Data sets:",
                      breaks=c("JM", "LC"),
                      labels=c("Packer et. al. 2019", "Current study"))+
  scale_y_continuous(breaks=seq(0,6000,500))





######################################
# Supplementary figure 2
######################################

BWMS <- c("MSxa","MSxp", 
          "MSxap", "MSxpa", "MSxpp", 
          "MSxapp","MSxpaa", "MSxpap", "MSxppa", "MSxppp",
          "MSpappa", "MSxappp", "MSxpaaa",  "MSxpapa", "MSxpapp", "MSxppap","MSxpppa", "MSxpppp","MSxppppx",
          "MSxapppa", "MSxapppp", "MSxpapap","MSxpaaap",  "MSxpappa", "MSxpappp", "MSxppapp", "MSxpppaa", "MSxpppap", 
          "MSxppppa", "MSxppppp","MSxapppax", "MSxappppx")
ms_BWMs <- subset(ms, idents = BWMS)

levels(ms_BWMs) <- BWMS

DoHeatmap(ms_BWMs, features = c('pha-4', 'hnd-1', 'nhr-67', 'T08H10.1', 'zig-6', 'clec-264', 
                                'F41D9.2', 'lbp-2', 'tlp-1', 'sul-2', 'col-111'), label = T, draw.lines = T, lines.width = 2,angle = 90, raster = F) + NoLegend()



##########################################################################################################################################################
### upplementary figure 3 FBXs genes
##########################################################################################################################################################

fbxs_list <- rownames(ms_new)
ms_fbx <- subset(ms, idents =unique(c("MSx", "MSxa", "MSxp",
                                      "MSxaa", "MSxap", "MSxpa", "MSxpp", 
                                      'MSxaaa', 'MSxaap', 'MSxapa', 'MSxapp', 
                                      'MSxpaa', 'MSxpap', 'MSxppa', 'MSxppp',
                                      'MSpaaaa', 'MSaaaap',"MSpaaap", 'MSxaapa', 'MSaaapp',
                                      'MSpaaaa','MSpaaap','MSaaaap','MSxaapa','MSaaapp','MSxapaa','MSxapap','MSaappa', 'MSpappa', 'MSxappp',
                                      'MSxpaaa', 'MSxpapa', 'MSxpapp', 'MSxppaa', 'MSxppap', 'MSxpppa', 'MSxpppp',
                                      'MSpaaaaa', 'MSpaaaap','MSpaaapa','MSpaaapp', 'MSpaapaa', 'MSxaapap', 'MSaaaaaa', 'MSaaaaap', 'MSaaaapa', 'MSaaaapp', 'MSaaapaa',
                                      'MSxapaaa', 'MSxapaap', 'MSxapapp', 'MSxapppa', 'MSxapppp',
                                      'MSxpaaaa', 'MSxpaaap', 'MSxpapaa', 'MSxpapap', 'MSxpappa', 'MSxpappp', 'MSxppaaa', 'MSxppaap','MSxppapp',
                                      'MSxpppaa', 'MSxpppap', 'MSxppppa', 'MSxppppp',
                                      "MSxaapapx", "MSaaaaaaa", "MSapaapp", "MSaaaaapa", "MSxapppax", "MSxappppx","MSxapaapa","MSaaappp")))
levels(ms_fbx) <- unique(c("MSx", "MSxa", "MSxp",
                           "MSxaa", "MSxap", "MSxpa", "MSxpp", 
                           'MSxaaa', 'MSxaap', 'MSxapa', 'MSxapp', 
                           'MSxpaa', 'MSxpap', 'MSxppa', 'MSxppp',
                           'MSpaaaa', 'MSaaaap',"MSpaaap", 'MSxaapa', 'MSaaapp',
                           'MSpaaaa','MSpaaap','MSaaaap','MSxaapa','MSaaapp','MSxapaa','MSxapap','MSaappa', 'MSpappa', 'MSxappp',
                           'MSxpaaa', 'MSxpapa', 'MSxpapp', 'MSxppaa', 'MSxppap', 'MSxpppa', 'MSxpppp',
                           'MSpaaaaa', 'MSpaaaap','MSpaaapa','MSpaaapp', 'MSpaapaa', 'MSxaapap', 'MSaaaaaa', 'MSaaaaap', 'MSaaaapa', 'MSaaaapp', 'MSaaapaa',
                           'MSxapaaa', 'MSxapaap', 'MSxapapp', 'MSxapppa', 'MSxapppp',
                           'MSxpaaaa', 'MSxpaaap', 'MSxpapaa', 'MSxpapap', 'MSxpappa', 'MSxpappp', 'MSxppaaa', 'MSxppaap','MSxppapp',
                           'MSxpppaa', 'MSxpppap', 'MSxppppa', 'MSxppppp',
                           "MSxaapapx", "MSaaaaaaa", "MSapaapp", "MSaaaaapa", "MSxapppax", "MSxappppx","MSxapaapa","MSaaappp"))

ms_fbx$ident <- ms_fbx$manual.annot.ids


#setdiff(levels(ms), levels(ms_fbx))

additional_Fboxes <- c("B0294.3", "BATH-28", "C07D10.1", "C14B1.3", 'C17C3.6',"C27F2.7", 'dre-1', 'F3139.3', 'F43C9.1',
                       'F56F10.2', 'fog-2', 'fsn-1', 'ftr-1', 'K02A6.3', 'K06A9.2', 'lin-23', 'mec-15', 'mfb-1', 'R10E8.8', 
                       'sel-10', 'skpt-1', 'T03F6.4', 'T05F1.13', 'T10C6.15', 'T28B11.1', 'T28B4.1', 'Y37H2A.7', "Y72A10A.1", 'Y9C9A.8', 'ZK328.6')
ms_fbx_aver <- AverageExpression(ms_fbx, features = c(additional_Fboxes, fbxs_list[grepl("fbx[a,b,c]*",fbxs_list)]), return.seurat = F, slot = "scale.data" )

ms_fbx_aver@active.ident

ms_fbx_aver <- ms_fbx_aver$RNA

# Clusterrrize the expression
hc <- hclust(dist(ms_fbx_aver))
plot(hc)
hc$labels[hc$order]
fbx_selected_genes <- hc$labels[hc$order] 

fbx_remove <- c('dre-1',"mfb-1", "C07D10.1", "T28B4.1", "K06A9.2", "F56F10.2","fbxc-57","fbxa-46","fbxa-137","fbxb-101","fbxa83","fbxa-67","fbxa-23","fbxa-90","fbxc-54","fbxc-35","fbxa-128","fbxb-87","fbxa-153","fbxa-10","fbxa-218","fbxa-11","fbxa-97","fbxa-93","fbxa-114","fbxc-23","fbxa-64","fbxa-167","fbxl-1","fbxa-146","fbxc-33","fbxa-13","fbxc-44","fbxa-108","fbxa-215","fbxa-149","fbxa-124","fbxa-200","fbxc-53","fbxc-20")

fbx_selected_genes <- fbx_selected_genes[! fbx_selected_genes %in% fbx_remove]


#DoHeatmap(ms_fbx, features = fbxs_list[grepl("fbx[a,b,c]*",fbxs_list)], label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()
#DoHeatmap(ms_fbx, features = hc$labels[hc$order], label = T, draw.lines = T, lines.width = 2, angle = 90, group.by = "ident") + NoLegend()
ms_fbx_aver <- AverageExpression(ms_fbx, features = c(additional_Fboxes, fbxs_list[grepl("fbx[a,b,c]*",fbxs_list)]), return.seurat = T, slot = "scale.data" )

DoHeatmap(ms_fbx_aver, features = fbx_selected_genes, label = T, draw.lines = F, angle = 90, raster = F) + NoLegend()



##########################################################################################################################################################
### Calculating JHU
##########################################################################################################################################################
BWMS <- c("MSxa","MSxp", 
          "MSxap", "MSxpa", "MSxpp", 
          "MSxapp","MSxpaa", "MSxpap", "MSxppa", "MSxppp",
          "MSpappa", "MSxappp", "MSxpaaa",  "MSxpapa", "MSxpapp", "MSxppap","MSxpppa", "MSxpppp",#"MSxppppx",
          "MSxapppa", "MSxapppp", "MSxpapap","MSxpaaap",  "MSxpappa", "MSxpappp", "MSxppapp", "MSxpppaa", "MSxpppap", 
          "MSxppppa", "MSxppppp","MSxapppax", "MSxappppx")

ms_all_BWMs <- subset(ms, idents =BWMS)
levels(ms_all_BWMs) <- BWMS

# Using variable features for reducing the gene space in the averaged gene expression table. In principle it doesn't change the final result much and not really necessary.
use_Var_features <-  F
if (use_Var_features){
  ms_all_BWMs <- FindVariableFeatures(ms_all_BWMs, nfeatures = 5000)
  ave_ms_subset <- data.frame(AverageExpression(ms_all_BWMs, return.seurat = FALSE,  group.by = "manual.annot.ids",  
                                                add.ident = NULL, slot = "data",  verbose = TRUE, features = VariableFeatures(ms_all_BWMs)))
}else{
  ave_ms_subset <- data.frame(AverageExpression(ms_all_BWMs, return.seurat = FALSE,  group.by = "manual.annot.ids",  
                                                add.ident = NULL, slot = "data",  verbose = TRUE))
}

ave_ms_subset <- ave_ms_subset[,paste0('RNA.',BWMS)]
colnames(ave_ms_subset) <- BWMS

JSH_dist_matrix <- jsDist(as.matrix(data.frame(ave_ms_subset)))
JSH_dist_matrix <- as.data.frame(JSH_dist_matrix)
colnames(JSH_dist_matrix) <- colnames(ave_ms_subset)
rownames(JSH_dist_matrix) <- colnames(ave_ms_subset)
JSH_dist_matrix <- as.matrix(JSH_dist_matrix/10^6)
JSH_dist_matrix <- round(x = JSH_dist_matrix, digits = 2)

#JSH_dist_matrix <- cbind(rownames(JSH_dist_matrix),JSH_dist_matrix)

#JSH_dist_matrix <- pivot_longer(as.data.frame(JSH_dist_matrix), cols = 2:35, names_to = 'cell_ident', values_to = 'JSHD')
#colnames(JSH_dist_matrix) <- c('cell_ident1', 'cell_ident2', 'JSHD')
#JSH_dist_matrix$JSHD <- as.numeric(JSH_dist_matrix$JSHD)
#ggplot(JSH_dist_matrix, aes(cell_ident1, cell_ident2, fill= JSHD)) + 
#geom_tile()


#install.packages('pheatmap') # if not installed already


pheatmap(JSH_dist_matrix, 
         display_numbers = T, 
         color = viridis(50), 
         cluster_rows = F, 
         cluster_cols = F, 
         fontsize_number = 8,
         legend = T,
         number_color = 'white')

##############################
# daughter - mother JSHD box plot
##############################
MSxx <- c("MSxa","MSxp") 
MSxxx <- c("MSxap", "MSxpa", "MSxpp") 
MSxxxx <- c("MSxapp","MSxpaa", "MSxpap", "MSxppa", "MSxppp")
MSxxxxx <- c("MSpappa", "MSxappp", "MSxpaaa",  "MSxpapa", "MSxpapp", "MSxppap","MSxpppa", "MSxpppp")
MSxxxxxx <- c("MSxapppa", "MSxapppp", "MSxpapap","MSxpaaap",  "MSxpappa", "MSxpappp", "MSxppapp", "MSxpppaa", "MSxpppap","MSxppppa", "MSxppppp","MSxapppax", "MSxappppx")

MSxxx_values <- JSH_dist_matrix[MSxxx, MSxxx][upper.tri(JSH_dist_matrix[MSxxx, MSxxx])]
MSxxx_to_MSxx_values <- as.vector(JSH_dist_matrix[MSxxx, MSxx])
list_1 <- list(values = c(MSxxx_values, MSxxx_to_MSxx_values), 
               cell_division = c("MSxxx"),
               group = c(rep("MSxxx", length(MSxxx_values)), rep("MSxxx_to_MSxx", length(MSxxx_to_MSxx_values)) ) )
as.data.frame(list_1)

MSxxxx_values <- JSH_dist_matrix[MSxxxx, MSxxxx][upper.tri(JSH_dist_matrix[MSxxxx, MSxxxx])]
MSxxxx_to_MSxxx_values <- as.vector(JSH_dist_matrix[MSxxxx, MSxxx])
list_2 <- list(values = c(MSxxxx_values, MSxxxx_to_MSxxx_values), 
               cell_division = c("MSxxxx"),
               group = c(rep("MSxxxx", length(MSxxxx_values)), rep("MSxxxx_to_MSxxx", length(MSxxxx_to_MSxxx_values)) ) )
as.data.frame(list_2)

MSxxxxx_values <- JSH_dist_matrix[MSxxxxx, MSxxxxx][upper.tri(JSH_dist_matrix[MSxxxxx, MSxxxxx])]
MSxxxxx_to_MSxxxx_values <- as.vector(JSH_dist_matrix[MSxxxxx, MSxxxx])
list_3 <- list(values = c(MSxxxxx_values, MSxxxxx_to_MSxxxx_values), 
               cell_division = c("MSxxxxx"),
               group = c(rep("MSxxxxx", length(MSxxxxx_values)), rep("MSxxxxx_to_MSxxxx", length(MSxxxxx_to_MSxxxx_values)) ) )
as.data.frame(list_3)


MSxxxxxx_values <- JSH_dist_matrix[MSxxxxxx, MSxxxxxx][upper.tri(JSH_dist_matrix[MSxxxxxx, MSxxxxxx])]
MSxxxxxx_to_MSxxxxx_values <- as.vector(JSH_dist_matrix[MSxxxxxx, MSxxxxx])
list_4 <- list(values = c(MSxxxxxx_values, MSxxxxxx_to_MSxxxxx_values), 
               cell_division = c("MSxxxxxx"),
               group = c(rep("MSxxxxxx", length(MSxxxxxx_values)), rep("MSxxxxxx_to_MSxxxxx", length(MSxxxxxx_to_MSxxxxx_values)) ) )
as.data.frame(list_4)

data_to_plot <- rbind(as.data.frame(list_1),
                      as.data.frame(list_2),
                      as.data.frame(list_3),
                      as.data.frame(list_4))

ggplot(data_to_plot, aes(x=cell_division, y=values, fill=group)) + 
  geom_boxplot()
