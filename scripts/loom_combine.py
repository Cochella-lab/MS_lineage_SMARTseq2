#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 13:45:39 2020
the version of loompy should be the consistent with velocyto; otherwise error returned
@author: jingkui.wang
"""
import loompy
import glob 


files = glob.glob("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S117008_R9533/LOOMS/*.loom")
loompy.combine(files, 
               "/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S117008_R9533/LOOMS/S117008_R9533_merged.loom")

files = glob.glob("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S117007_R9533/LOOMS/*.loom")
loompy.combine(files, 
               "/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S117007_R9533/LOOMS/S117007_R9533_merged.loom")

files = glob.glob("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S117009_R9533/LOOMS/*.loom")
loompy.combine(files, 
               "/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S117009_R9533/LOOMS/S117009_R9533_merged.loom")


# folder S124890_R9968
files = glob.glob("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S124890_R9968/LOOMS/*.loom")
loompy.combine(files, 
               "/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S124890_R9968/LOOMS/S124890_R9968_merged.loom")

# folder S124889_R9968
files = glob.glob("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S124889_R9968/LOOMS/*.loom")
loompy.combine(files, 
               "/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S124889_R9968/LOOMS/S124889_R9968_merged.loom")

# folder S124891_R9968
files = glob.glob("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S124891_R9968/LOOMS/*.loom")
loompy.combine(files, 
               "/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/S124891_R9968/LOOMS/S124891_R9968_merged.loom")


# here combine merged loom files from all folders
files = glob.glob("/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/*/LOOMS/*merged.loom")
loompy.combine(files, '/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/raw_ngs_data/all_merged.loom')

loompy.combine(files, '/Users/jiwang/workspace/imp/scRNAseq_MS_lineage_dev/data/velocyto_all_merged.loom')

#dsmerged = loompy.connect('/Users/jiwang/workspace/imp/scRNAseq_MS_lineage_dev/data/velocyto_all_merged.loom')

#ds0.shape
#dsmerged.shape
