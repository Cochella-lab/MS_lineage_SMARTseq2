##########################################################################
##########################################################################
# Project: Aleks' MS lineage project
# Script purpose: processing John Murray's scRNA-seq 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Mar 19 11:10:31 2020
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

version.DATA = 'scRNA_Murray_2019'
version.analysis =  paste0(version.DATA, '_20200319')
dataDir = paste0("../data/gene_counts/")
resDir = paste0("results/", version.analysis)
tabDir = paste0("results/", version.analysis, "/tables/")

RdataDir = paste0("results/", version.analysis, "/Rdata/")
if(!dir.exists("results/")){dir.create("results/")}
if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

# source.my.script('scATAC_functions.R')
########################################################
########################################################
# Section : processing the Murray's scRNA-seq raw data
# 
########################################################
########################################################
process.scRNAseq.for.early.embryo.packer.et.al = function()
{
  Install.VisCello.celegans = FALSE
  if(Install.VisCello.celegans){
    ##########################################
    # details see https://github.com/qinzhu/VisCello.celegans
    ##########################################
    devtools::install_local("../VisCello.celegans", force=T)
    packageurl <- "https://cran.r-project.org/src/contrib/Archive/tidytree/tidytree_0.2.6.tar.gz"
    install.packages(packageurl, repos=NULL, type="source")
    library(VisCello.celegans)
    cello()
  }
  Check.Cello.DataSet = FALSE
  if(Check.Cello.DataSet){
    ##########################################
    # more details in 
    # https://github.com/qinzhu/VisCello.celegans/blob/master/inst/app/global.R
    ##########################################
    cello.data.path = "../VisCello.celegans//inst/app/data/"
    
    eset = readRDS(paste0(cello.data.path, 'eset.rds'))
    saveRDS(eset, file =  paste0(RdataDir, 'cello_Parker_et_al_allData.rds'))
    #eset <- readRDS("data/eset.rds")
    
    ## clusters with coordinates in reduced dimensions (PCA, UMAP)
    clist <- readRDS(paste0(cello.data.path, "clist.rds"))
    
    elist <- readRDS(paste0(cello.data.path, "elist.rds"))
    ct_tbl <-  readRDS("data/s6_tbl.rds")
    lin_tbl <-  readRDS("data/s7_tbl.rds")
    tree_tbl <- as_tibble(readRDS("data/lineage_tree_tbl.rds"))
    lin_sc_expr <- readRDS("data/lin_sc_expr_190602.rds")
    lin.expanded.list <- readRDS("data/lin_expanded_list_0602.rds")
    avail_nodes <- readRDS("data/avail_nodes.rds")
    cell_type_markers <- read.xlsx("data/Supplementary_Tables_190611.xlsx",sheet=1, startRow=4)
    lineage_markers <-  read.xlsx("data/Supplementary_Tables_190611.xlsx",sheet=4, startRow=7)
    
  }else{
    eset = readRDS(file = paste0(RdataDir, 'cello_Parker_et_al_allData.rds'))
    pmeda = data.frame(pData(eset))
    
    sels = which(pmeda$embryo.time.bin == '< 100')
    
    pp = pmeda[sels, ]
     
  }
}


process.scRNAseq.for.early.embryo.Tintori.et.al = function(start.from.raw.counts = FALSE)
{
  if(!start.from.raw.counts)
  {
    cello.data.path = "../Celegans.Tintori.Cello/"
    
    eset = readRDS(paste0(cello.data.path, 'eset.rds'))
    saveRDS(eset, file =  paste0(RdataDir, 'VisCello_Tintori_et_al.rds'))
    
    pmeda = data.frame(pData(eset))
    
    library(Seurat)
    sels = which(pmeda$Usable.Quality. == 'Yes')
    
    eset = eset[, sels]
    
    ee = CreateSeuratObject(counts = eset@assayData$exprs, assay = 'RNA', meta.data = pmeda)
    ee@assays$RNA@data = eset@assayData$norm_exprs
    
    FeatureScatter(ee, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    
    # ee <- subset(ee, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
    #al = subset(ee, cells = Usable.Quality. == 'Yes')
    
    al <- FindVariableFeatures(object = ee, nfeatures = 2000)
    al <- ScaleData(object = al)
    al <- RunPCA(object = al, npcs = 50, verbose = FALSE)
    
    Idents(al) = al$lineage
    al <- RunUMAP(object = al, dims = 1:30, n.neighbors = 10, min.dist = 0.3)
    
    DimPlot(al, label = TRUE, 
            pt.size = 2, label.size = 3, repel = TRUE)
    
    #ee = al[, which(al$Usable.Quality. == 'Yes')]
    saveRDS(al, file =  paste0(RdataDir, 'Tintori_et_al_highQualtiyCells.rds'))
  }else{
    ##########################################
    # process the count data
    ##########################################
    dataDir = paste0("../Celegans.Tintori.Cello/nf_output")
    library(data.table)
    # first concatenate the data
    xlist <-list.files(path=dataDir, pattern = "*merged_gene_counts.txt", full.names = TRUE)
    aa = NULL
    
    for(n in 1:length(xlist)){
      # n = 1
      cat(n, '\t')
      cat(xlist[n], '\n')
      
      if(n==1){
        aa = data.frame(fread(xlist[n], header=TRUE, sep="\t", stringsAsFactors=FALSE), stringsAsFactors = FALSE)
      }else{
        test = data.frame(fread(xlist[n], header=TRUE, sep="\t", stringsAsFactors=FALSE), stringsAsFactors = FALSE)
        mm = match(aa$ENSEMBL_ID, test$ENSEMBL_ID)
        aa = data.frame(aa, test[mm, -grep('ENSEMBL_ID', colnames(test))])
      }
    }
    
    colnames(aa)[1] = 'gene';
    if(length(grep("out_gene.featureCounts.txt", colnames(aa)))>0) {
      aa = aa[, -grep("out_gene.featureCounts.txt", colnames(aa))]
    }
    
    library(SingleCellExperiment)
    library(scater)
    library(scRNA.seq.funcs)
    library(scran)
    options(stringsAsFactors = FALSE)
    
    source.my.script("scRNAseq_functions.R")
    annot = read.csv(file = "data/annotations/BioMart_WBcel235_noFilters.csv",
                     header = TRUE)
    counts = convertGeneNames.forCountTable(aa = aa, annot = annot)
    gg.Mt = find.particular.geneSet("Mt", annot = annot)
    gg.ribo = find.particular.geneSet("ribo", annot = annot)
    
    ##########################################
    # preprocess metadata 
    ##########################################
    design = read.delim(file = '../Celegans.Tintori.Cello/nf_output/PRJNA312176.txt', sep = '\t', header= TRUE)
    
    eset = readRDS(paste0('../Celegans.Tintori.Cello/eset.rds'))
    pmeda = data.frame(pData(eset))
    
    mm = match(design$sample_title, pmeda$sample)
    design = data.frame(design, pmeda, stringsAsFactors = FALSE)
    
    kk = match(colnames(counts), design$run_accession)
    design = design[kk, ]
    colnames(counts) = design$sample
    rownames(design) = design$sample
    design = data.frame(design)
    
    ## using scran and scater for QC and normalization 
    
    
    
    library(dplyr)
    library(Seurat)
    library(patchwork)
    #pbmc <- CreateSeuratObject(counts = counts, project = "Tintori", min.cells = 3, min.features = 200, 
                               meta.data = design)
    #pbmc
    
  }
  
  
}




