#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  5 17:49:29 2019

@author: jingkui.wang
"""

####################################
# test the analysis for output of STITCH 
# https://nbviewer.jupyter.org/github/theislab/paga/blob/master/zebrafish/zebrafish.ipynb
####################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
from scipy.sparse import csr_matrix
import scanpy as sc

sc.set_figure_params(dpi=200, color_map='viridis')  # low dpi (dots per inch) yields small inline figures
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
results_file = './zebrafish_graph/zebrafish.h5ad'

var = pd.read_csv('./zebrafish_graph/output_Wagner2018/genes.txt', index_col=0, header=None)

cluster_indices = np.genfromtxt('./zebrafish_graph/output_Wagner2018/cell_IDs.txt', dtype=int)
clusters = pd.Categorical(cluster_indices)

clusters.categories = ['?'] + np.genfromtxt('./zebrafish_graph/output_Wagner2018/cell_IDs_names.txt', dtype=str, delimiter=',').tolist()
obs = pd.DataFrame({'clusters': clusters})


cell_cell_edges = np.genfromtxt('./zebrafish_graph/output_Wagner2018/edges.csv', dtype=int, delimiter=';')
#cell_cell_edges -= 1  # our indexing starts with 0 as common in Python, C etc.

from scipy.sparse import coo_matrix
rows = cell_cell_edges[:, 0]
cols = cell_cell_edges[:, 1]
length = len(obs)
ones = np.ones(len(rows), np.uint32)
connectivities = coo_matrix((ones, (rows, cols)), shape=(length, length))
# make sure it's symmsetric
connectivities = connectivities + connectivities.T

uns = {'neighbors': {'connectivities': connectivities.tocsr()}}

adata = sc.AnnData(obs=obs, var=var, uns=uns)
# remove the outliers
adata = adata[adata.obs['clusters'] != '?']

adata

adata.write(results_file)

# Roughly embedding the graph
adata = sc.read(results_file)

sc.tl.draw_graph(adata, maxiter=100)
sc.pl.draw_graph(adata)

#sc.pl.draw_graph(adata, color='ClusterName', edges=True, size=40, legend_loc='on data')


#sc.pl.draw_graph(adata, color=['clusters_coarse', 'time'], legend_loc='on data', legend_fontsize=4)

#sc.tl.umap(adata)



