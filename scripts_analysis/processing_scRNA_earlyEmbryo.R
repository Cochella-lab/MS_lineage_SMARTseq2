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
    library(SingleCellExperiment)
    library(scater)
    options(stringsAsFactors = FALSE)
    
    #load(file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_design_technicalRepMerged.Rdata'))
    sce <- SingleCellExperiment(assays = list(counts = counts),
                                colData = as.data.frame(design),
                                rowData = data.frame(gene_names = rownames(counts), feature_symbol = rownames(counts)))
    
    keep_feature <- rowSums(counts(sce) > 0) > 0
    sce <- sce[keep_feature, ]
    
    #is.spike <- grepl("^ERCC", rownames(sce))
    is.mito <- rownames(sce) %in% gg.Mt;
    is.ribo <- rownames(sce) %in% gg.ribo;
    summary(is.mito)
    summary(is.ribo)
    
    sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito, Ribo=is.ribo))
    
    
    pdfname = paste0(resDir, "/Tintori.et.al_scRNAseq_QCs_cells_filterting.pdf")
    pdf(pdfname, width=12, height = 6)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    # some general statistics for each request and lane
    plotColData(sce, y="pct_counts_Ribo", x="total_counts") + ggtitle("% of rRNA contamination")
    plotColData(sce, y="pct_counts_Mt", x="total_counts") + ggtitle("% of Mt")
    
    plotColData(sce, y="log10_total_counts", x="total_counts") + ggtitle("total nb of reads mapped to transcripts")
    
    plotColData(sce, y="total_features_by_counts", x="total_counts") + ggtitle("total nb of genes")
    
    plotColData(sce, y="Number.of.ERCC.reads", x="total_counts") + ggtitle("total nb of genes")
    
    plotColData(sce,
                x = "log10_total_counts",
                y = "log10_total_features_by_counts",
                #colour_by = "percent_mapped",
                colour_by = "Usable.Quality.",
                size_by = "pct_counts_Mt"
    ) + scale_x_continuous(limits=c(4, 7)) +
      scale_y_continuous(limits = c(2.5, 4.1)) +
      geom_hline(yintercept=log10(c(500, 1000, 5000)) , linetype="dashed", color = "darkgray", size=0.5) +
      geom_vline(xintercept = c(4:6), linetype="dotted", color = "black", size=0.5)
    
    dev.off()
    
    
    threshod.total.counts.per.cell = 10^5
    threshod.nb.detected.genes.per.cell = 750 # Adjust the number for each batch
    
    filter_by_total_counts <- (sce$total_counts > threshod.total.counts.per.cell)
    table(filter_by_total_counts)
    filter_by_expr_features <- (sce$total_features_by_counts > threshod.nb.detected.genes.per.cell)
    table(filter_by_expr_features)
    filter_by_MT = sce$pct_counts_Mt < 10 # Adjust the number for each batch
    table(filter_by_MT)
    filter_by_Ribo = sce$pct_counts_Ribo < 50
    
    sce$use <- (
      filter_by_expr_features & # sufficient features (genes)
        filter_by_total_counts & # sufficient molecules counted
        filter_by_Ribo & # sufficient endogenous RNA
        filter_by_MT # remove cells with unusual number of reads in MT genes
    )
    table(sce$use)
    
    sce = sce[, sce$use]
    
    
    pdfname = paste0(resDir, "/Tintori.et.al_scRNAseq_QCs_genes_filterting_", version.analysis, ".pdf")
    pdf(pdfname, width=16, height = 10)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    #plotQC(reads, type = "highest-expression", n=20)
    fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
    plotHighestExprs(sce, n=50) + fontsize
    
    #fontsize <- theme(axis.text=element_text(size=16), axis.title=element_text(size=16))
    #plotHighestExprs(sce, n=30) + fontsize
    
    ave.counts <- calcAverage(sce)
    hist(log10(ave.counts), breaks=100, main="", col="grey80",
         xlab=expression(Log[10]~"average count"))
    
    num.cells <- nexprs(sce, byrow=TRUE)
    smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
                  xlab=expression(Log[10]~"average count"))
    dev.off()
    
    genes.to.keep <- num.cells > 5 & ave.counts >= 1     # detected in >= 5 cells, ave.counts >=5 but not too high
    summary(genes.to.keep)
    
    # remove mt and ribo genes
    genes.to.keep = genes.to.keep & ! rownames(sce) %in% gg.Mt & ! rownames(sce) %in% gg.ribo
    summary(genes.to.keep)
    
    sce <- sce[genes.to.keep, ]
    
    ##########################################
    #  normalization
    ##########################################
    plotColData(sce,
                x = "log10_total_counts",
                y = "log10_total_features_by_counts",
                colour_by = "Usable.Quality.",
                #colour_by = "",
                size_by = "pct_counts_Mt"
    ) + scale_x_continuous(limits=c(4, 7)) +
      scale_y_continuous(limits = c(2.5, 4.1)) +
      geom_hline(yintercept=log10(c(500, 1000, 5000)) , linetype="dashed", color = "darkgray", size=0.5) +
      geom_vline(xintercept = c(4:6), linetype="dotted", color = "black", size=0.5)
    
    sce$library.size = apply(counts(sce), 2, sum)
    qclust <- quickCluster(sce)
    sce <- computeSumFactors(sce, clusters = qclust)
    sce <- logNormCounts(sce, log = TRUE, pseudo_count = 1)
    
    par(mfrow = c(1, 2))
    plot(sizeFactors(sce), sce$library.size/1e6, log="xy", ylab="Library size (millions)", xlab="Size factor")
    plot(sizeFactors(sce), 1/sce$Number.of.ERCC.reads*10^6,  log = 'xy', 
         ylab="1/ERCC", xlab="Size factor")
    
    library(Seurat)
    library(ggplot2)
    
    ms = as.Seurat(sce, counts = 'counts', data = 'logcounts', assay = "RNA") # scran normalized data were kept in Seurat
    ms = ms[, which(ms$Usable.Quality. == 'Yes')]
    
    nfeatures = 2000
    ms <- FindVariableFeatures(ms, selection.method = "vst", nfeatures = nfeatures)
    
    #top10 <- head(VariableFeatures(ms), 10) # Identify the 10 most highly variable genes
    #plot1 <- VariableFeaturePlot(ms)
    #plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    #CombinePlots(plots = list(plot1, plot2)) # plot variable features with and without labels
    
    ms = ScaleData(ms, features = rownames(ms))
    
    # new normalization from Seurat
    # tried regress out the pct_counts_Mt but works less well
    #ms <- SCTransform(object = ms, variable.features.n = nfeatures) 
    ms <- RunPCA(object = ms, features = VariableFeatures(ms), verbose = FALSE)
    ElbowPlot(ms)
    
    nb.pcs = 30; n.neighbors = 20; min.dist = 0.2;
    ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
    DimPlot(ms, reduction = "umap", group.by = 'lineage',  label = TRUE, pt.size = 2, label.size = 5, repel = TRUE) + 
    scale_colour_hue(drop = FALSE) + ggtitle('Tintori (size factor)')
    
  }
  
  
}




