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
    devtools::install_local("../VisCello.celegans", force=T)
    packageurl <- "https://cran.r-project.org/src/contrib/Archive/tidytree/tidytree_0.2.6.tar.gz"
    install.packages(packageurl, repos=NULL, type="source")
    library(VisCello.celegans)
    cello()
    
    cello.data.path = "/Volumes/groups/cochella/jiwang/Projects/Aleks/scRNAseq_published_dataSets/VisCello.celegans"
    cello = readRDS(paste0(cello.data.path, '/inst/app/data/eset.rds'))
    saveRDS(cello, file =  paste0(RdataDir, 'cello_Parker_et_al_allData.rds'))
    
    eset <- readRDS("data/eset.rds")
    clist <- readRDS("data/clist.rds")
    elist <- readRDS("data/elist.rds")
    ct_tbl <-  readRDS("data/s6_tbl.rds")
    lin_tbl <-  readRDS("data/s7_tbl.rds")
    tree_tbl <- as_tibble(readRDS("data/lineage_tree_tbl.rds"))
    lin_sc_expr <- readRDS("data/lin_sc_expr_190602.rds")
    lin.expanded.list <- readRDS("data/lin_expanded_list_0602.rds")
    avail_nodes <- readRDS("data/avail_nodes.rds")
    cell_type_markers <- read.xlsx("data/Supplementary_Tables_190611.xlsx",sheet=1, startRow=4)
    lineage_markers <-  read.xlsx("data/Supplementary_Tables_190611.xlsx",sheet=4, startRow=7)
    
  }else{
    est = readRDS(file = paste0(RdataDir, 'cello_Parker_et_al_allData.rds'))
  }
  
}



