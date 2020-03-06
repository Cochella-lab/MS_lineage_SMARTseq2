#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 11:46:34 2019

@author: jingkui.wang
"""

import scvelo as scv
import glob
from matplotlib import pyplot as plt
import sys
from pathlib import Path
import pandas as pd

def scv_velocity_on_embedding(f_embedding, PATH, f_cell_cycle, fout_figure):
    print('Params:')
    print([f_embedding, PATH, f_cell_cycle, fout_figure])
    print('Loading loom...')
    #path = Path(PATH)
    #print(path.parent+'loom')
    f_loom = glob.glob(PATH+'/../loom/'+'*.loom')[0]
    print(f_loom)
    adata_raw = scv.read(f_loom, sparse=True, cache=True)

    print('Loading cell cycle genes...')
    cell_cycle_genes = pd.read_csv(f_cell_cycle).gene.ravel()
    
    print('Loading embedding...')
    emb_csv = scv.read(f_embedding) #'d2_ccrm.umap.csv'
    ###emb_csv.obs.index = ['possorted_genome_bam_A6HDS:'+x.rstrip('-1')+'x' for x in emb_csv.obs.index]
    emb_csv.obs.index = [x.rstrip('-1') for x in emb_csv.obs.index]
    adata_raw.obs.index = [x.rstrip('x').split(':')[-1] for x in adata_raw.obs.index]

    print('Sample names in embedding:')
    print(emb_csv.obs.index[1:10])
    print('Sample names in loom:')
    print(adata_raw.obs.index[1:10])
    print('before filtering:')
    print(adata_raw.shape)
    
    shared_cells = list(set(emb_csv.obs.index).intersection( set(adata_raw.obs.index) ))
    print('Shared cells:')
    print(len(shared_cells))

    #emb_csv.obs.index
    adata = adata_raw[shared_cells].copy()
    adata.obsm['X_umap'] = emb_csv[shared_cells].X.copy()
    print('after cell filtering:')
    print(adata.shape)

    adata = adata[:, [x not in cell_cycle_genes for x in adata.var.index]]
    print('after removal of cell cycle genes:')
    print(adata.shape)

    #sc.pp.recipe_seurat(adata)

    #scv.utils.show_proportions(adata)
    #scv.utils.cleanup(adata, clean='all')

    #adata

    print('Velocity...')
    #scv.pp.filter_and_normalize(adata, min_counts=20, min_counts_u=10, n_top_genes=3000)
    scv.pp.moments(adata, n_pcs=3, n_neighbors=30)

    scv.tl.velocity(adata)
    scv.tl.velocity_graph(adata)

    print('Clustering...')
    scv.tl.louvain(adata)
    #sc.pl.umap(adata)

    print('Plot...')
    print(fout_figure)
    #scv.pl.velocity_embedding_grid(adata, legend_loc='on data', basis='umap', dpi=200)
    #scv.pl.velocity_embedding_stream(adata, legend_loc='on data', basis='umap', dpi=200)
    scv.pl.velocity_embedding_stream(adata, legend_loc='on data', basis='umap', dpi=400)#, save=fout_figure+'.stream.png')
    plt.savefig(fout_figure+'.stream.png')


    #!scv.pl.velocity_embedding_grid(adata, color=['Xbp1'], legend_loc='on data', basis='umap', dpi=200)
    scv.pl.velocity_embedding_grid(adata, legend_loc='on data', basis='umap', dpi=400)#, save=fout_figure)
    plt.savefig(fout_figure)


#adata_d2 = scv.read("/home/artem/BCells/timecourse/HTMFJAFXY_all/aligned/97051/loom/possorted_genome_bam_A6HDS.loom", sparse=True, cache=True)
#adata_d3 = scv.read("./10x73640/outs/loom/possorted_genome_bam_U2IDE.loom", sparse=True, cache=True)

if __name__=='__main__':
    print(sys.argv)
    if len(sys.argv) == 5:
        _, f_embedding, path, f_cell_cycle, fout_figure = sys.argv
        scv_velocity_on_embedding(f_embedding, path, f_cell_cycle, fout_figure)