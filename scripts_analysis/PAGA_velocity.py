# -*- coding: utf-8 -*-
"""
Spyder Editor
This is a temporary script file.
"""
import scanpy as sc
import numpy as np
from matplotlib import rcParams
#import matplotlib.pyplot as pl

import matplotlib.pyplot as plt
#%matplotlib inline 

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, color_map='viridis')  # low dpi (dots per inch) yields small inline figures
sc.logging.print_versions()

adata = sc.read_loom('Seurat_tmp.loom', sparse=True, cleanup=False, X_name='SCT', obs_names='CellID', var_names='Gene', dtype='float32')

sc.pp.neighbors(adata, n_neighbors = 10, n_pcs = 20, use_rep='pca_cell_embeddings')
#sc$tl$leiden(adata, resolution = 1.0)
sc.tl.draw_graph(adata)
sc.pl.draw_graph(adata)

sc.pl.draw_graph(adata, color='ClusterID', legend_loc='on data')

sc.tl.paga(adata, groups='ClusterName')

#plt.figure()
sc.pl.paga(adata, threshold=0.7, layout='fa', frameon=False, fontsize=4, node_size_scale=0.3, max_edge_width=0.2, save = 'paga_test.pdf')
#plt.ion()
#plt.show()

sc.tl.umap(adata, init_pos = "paga")


sc$pl$paga_compare(
  adata,
  basis = "umap",
  threshold = 0.15,
  edge_width_scale = 0.5,
  save = TRUE
)


import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy.api as sc

#print("hello world")
#pl.plot(arange(5))

sc.set_figure_params(dpi=100)  # low dpi (dots per inch) yields small inline figures
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()

results_file = '../results/paga_test/planaria_extended.h5ad'

paga_plot_params = dict(
    legend_fontsize=5,
    solid_edges='confidence_tree',
    dashed_edges='confidence',
    root='neoblast 1',
    layout='rt_circular',
    node_size_scale=0.5,
    node_size_power=0.9,
    max_edge_width=0.7,
    fontsize=3.5)

adata = sc.read('../results/paga_test/R_pca_seurat.txt')
clusters = pd.read_csv('../results/paga_test/R_annotation.txt', header=None)
adata.obs['clusters'] = clusters[0].values
adata.uns['iroot'] = 6  # root cell (the first neoblast in the file)

sc.utils.sanitize_anndata(adata)

adata.obs['clusters'].cat.reorder_categories([
    'early epidermal progenitors', 'activated early epidermal progenitors',
    'epidermal neoblasts', 'epidermis', 'epidermis DVb',
    'epidermis DVb neoblast', 'glia', 'phagocytes', 'goblet cells',
    'psd+ cells', 'gut progenitors', 'late epidermal progenitors 1',
    'late epidermal progenitors 2', 'muscle body', 'muscle pharynx',
    'muscle progenitors', 'neoblast 1', 'neoblast 2', 'neoblast 3',
    'neoblast 4', 'neoblast 5', 'neoblast 6', 'neoblast 7', 'neoblast 8',
    'neoblast 9', 'neoblast 10', 'neoblast 11', 'neoblast 12',
    'neoblast 13', 'ChAT neurons 1', 'ChAT neurons 2', 'GABA neurons',
    'otf+ cells 1', 'otf+ cells 2', 'spp-11+ neurons', 'npp-18+ neurons',
    'cav-1+ neurons', 'neural progenitors', 'pharynx cell type progenitors',
    'pgrn+ parenchymal cells', 'ldlrr-1+ parenchymal cells',
    'psap+ parenchymal cells', 'aqp+ parenchymal cells',
    'parenchymal progenitors', 'pharynx cell type', 'pigment',
    'protonephridia', 'secretory 1', 'secretory 2', 'secretory 3',
    'secretory 4'], inplace=True)


## run tsne extremely SLOW !!!
#sc.tl.tsne(adata)

colors = pd.read_csv('../results/paga_test/colors_dataset.txt', header=None, sep='\t')
# transform to dict where keys are cluster names
colors = {k: c for k, c in colors.values}
adata.uns['clusters_colors'] = [colors[clus] for clus in adata.obs['clusters'].cat.categories]

sc.pl.tsne(adata)
sc.pl.tsne(adata, color='red', legend_loc='on data', legend_fontsize=5, show=False)
#pl.plot()

adata.write(results_file)


