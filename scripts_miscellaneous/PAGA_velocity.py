# -*- coding: utf-8 -*-
"""
Spyder Editor
scanpy and spyder editor both were installed under conda env scanpy
This is a script file for testing scanpy, in particular for paga and velocity analysis
the current code was modified based on the example from 
https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html
"""
import numpy as np
#import umap
from matplotlib import rcParams
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

#%matplotlib inline 


sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
results_file = './MS_BWM_scanpy_visulization_paga.h5ad'
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor='white')  # low dpi (dots per inch) yields small inline figures

# =============================================================================
# 
# =============================================================================
"""
import seurat object to scanpy
code was adapted from https://github.com/LuckyMD/Code_snippets/blob/master/Seurat_to_anndata.ipynb
here not used anymore, because the r-seurat makes problems for scanpy 
"""
# import anndata2ri
# from rpy2.robjects import r
# # Activate the anndata2ri conversion between SingleCellExperiment and AnnData
# anndata2ri.activate()

# #Loading the rpy2 extension enables cell magic to be used
# #This runs R code in jupyter notebook cells
# %load_ext rpy2.ipython
# sc.settings.verbosity = 3
# sc.logging.print_header()
# R
# #library(Seurat)
# suppressPackageStartupMessages(library(Seurat))
# # Load PBMC dataset from Seurat tutorial
# savedRDS = '../results_jiwang/results/all_batches_202008/Rdata/processed_cells_scran.normalized_reference.based.annotation.scmap.seurat_ManualClusterAnnot_36.rds'
# bwm.data <- readRDS(file = savedRDS)
# # Initialize the Seurat object with the raw (non-normalized data).
# #pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# read the loom file saved from Seurat and convert it to Anndata for scanpy (this takes some time)
#adata = sc.read_loom('seuratObj_BWM_manual_cellIds_iteration_36.loom', sparse=True, 
#                     cleanup=False, X_name='RNA', obs_names='CellID', var_names='Gene', dtype='float32')

#adata.write('Adata_from_seurat_BWM.h5ad')
adata = sc.read('Adata_from_seurat_BWM.h5ad')

sc.pl.scatter(adata, x='FSC_log2', y='timingEst')

#sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
#sc.tl.pca(adata, svd_solver='arpack')
#sc.pp.neighbors(adata, n_neighbors=30, svd_solver='arpack)
#sc.pl.pca(adata, color='pha-4')
#sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors = 50, n_pcs = 30, use_rep='pca_cell_embeddings')
#sc$tl$leiden(adata, resolution = 1.0)
sc.tl.draw_graph(adata)

#sc.pl.draw_graph(adata)

sc.pl.draw_graph(adata, color='manual_ids', legend_loc='on data') # using clustering from Seurat

#sc.tl.louvain(adata, resolution=6)
#sc.pl.draw_graph(adata, color='louvain', legend_loc='on data')

sc.tl.paga(adata, groups='manual_ids')
#sc.tl.paga(adata)
sc.pl.paga(adata)

sc.pl.paga(adata, threshold=0.8, layout='fa')
#sc.pl.paga(adata, threshold=0.7, layout='rt', frameon=True, fontsize=10, node_size_scale=0.5, max_edge_width=10)
sc.pl.paga(adata,
    threshold=0.9,           
    solid_edges='connectivities',
    #dashed_edges='connectivities', c
    layout='fa',
    node_size_scale=3,
    node_size_power=1,
    max_edge_width=2,
    frameon=False,
    fontsize=20
    )

#plt.ion()
#plt.show()
adata.write(results_file)

adata.obs['manual_ids'].cat.categories

sc.tl.draw_graph(adata, init_pos='paga')

sc.tl.draw_graph(adata, init_pos='paga', layout='fr')

X = adata.obsm['X_draw_graph_fr'].copy()
adata.obsm['X_draw_graph_fr'] = X.copy()
# adata.obsm['X_draw_graph_fr'][adata.obs['clusters'] == '0', 1] -= 400
adata.obsm['X_draw_graph_fr'][adata.obs['manual_ids'] == '0', 0] += 500
# adata.obsm['X_draw_graph_fr'][adata.obs['clusters'] == '1', 1] -= 1000
adata.obsm['X_draw_graph_fr'][adata.obs['manual_ids'] == '1', 0] -= 500

sc.pl.draw_graph(adata, color='manual_ids', edges=False, size=40, legend_loc='on data')

# umap now working 
sc.tl.umap(adata)
sc.pl.umap(adata, color = 'manual_ids')

sc.tl.umap(adata, init_pos=sc.tl._utils.get_init_pos_from_paga(adata))
#sc.tl.umap(adata, init_pos='paga')
sc.pl.umap(adata, color = 'manual_ids')


sc.pl.paga_compare(
  adata,
  basis = "umap",
  threshold = 0.15,
  edge_width_scale = 0.5,
  save = False
)

"""
Run scVI for imputation before gene dynamics 
original code from https://www.scvi-tools.org/en/stable/tutorials/basic_tutorial.html
working directory is :
    /Users/jiwang/workspace/imp/scRNAseq_MS_lineage_dev/results_jiwang/results/all_batches_202008/scVI
"""
import sys
IN_COLAB = "google.colab" in sys.modules

def allow_notebook_for_test():
    print("Testing the basic tutorial notebook")

show_plot = True
test_mode = False
save_path = "./"

if not test_mode:
    save_path = "../../data"

if IN_COLAB:
    !pip install --quiet scvi[notebooks]

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import scvi



##################################
##
## test scVelo
##################################
import scvelo as scv
scv.logging.print_version()

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.set_figure_params('scvelo')  # for beautified visualization

adata =scv.read_loom('merged.loom') ## extremely slow
adata.write('Adata_for_RNAvelocity.h5ad')

# show proportions of spliced/unspliced abundances
scv.utils.show_proportions(adata)
adata

scv.pp.filter_genes(adata, min_shared_counts=10)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=3000)
scv.pp.log1p(adata)

scv.pp.filter_and_normalize(adata, min_shared_counts=30, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata)
scv.tl.velocity_confidence(adata)

#scv.tl.velocity_graph(adata)
scv.tl.umap(adata)

scv.pl.velocity_embedding_stream(adata, basis='umap')

scv.pl.velocity_embedding(adata, basis='umap', arrow_length=1.2, arrow_size=1.2, dpi=150)
scv.pl.velocity_embedding_grid(adata, color='pha-4', layer=['velocity', 'spliced'], arrow_size=1.5)

scv.pl.velocity_graph(adata)

#############
# how to combine the paga and velocity 
# some advice from the author 
# https://github.com/theislab/paga/issues/11
#############
scv.tl.velocity_graph(adata)

sc.tl.louvain(adata, resolution=6)
sc.tl.paga(adata)
sc.tl.paga(adata, use_rna_velocity=True)

scv.pl.paga(adata, threshold=0.01, transitions='transitions_confidence')

scv.pl.paga()
sc.tl.draw_graph(adata, init_pos='paga', layout='fr', maxiter=50)
X = adata.obsm['X_draw_graph_fr'].copy()
adata.obsm['X_draw_graph_fr'] = X.copy()
# adata.obsm['X_draw_graph_fr'][adata.obs['clusters'] == '0', 1] -= 400
adata.obsm['X_draw_graph_fr'][adata.obs['louvain'] == '0', 0] += 500
# adata.obsm['X_draw_graph_fr'][adata.obs['clusters'] == '1', 1] -= 1000
adata.obsm['X_draw_graph_fr'][adata.obs['louvain'] == '1', 0] -= 500
sc.pl.draw_graph(adata, color='louvain', edges=True, size=40, legend_loc='on data')