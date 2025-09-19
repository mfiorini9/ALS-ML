## run this in Narval

salloc -A def-tdurcan --time=0-8 -c 1 --mem=200g
salloc -A def-sfarhan --time=0-2 -c 1 --mem=10g

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Fibro 

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

nano generalizable_HVG_and_NMF_Fibro_ALS_jumbo.sh

#!/bin/bash  
#SBATCH --account=def-sfarhan
#SBATCH --time=02-00:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=200g          # memory per cor
#SBATCH --job-name=generalizable_HVG_and_NMF_Fibro_ALS_jumbo
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/generalizable_HVG_and_NMF_Fibro_ALS_jumbo.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano generalizable_HVG_and_NMF_Fibro_ALS_jumbo.py


################################################################## clean version ------ USE THIS ONE

import os
import sys

# If you truly need a custom site-packages path, keep this. Otherwise remove it.
PYTHONENV_PATH = "/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA/lib/python3.8/site-packages"
if PYTHONENV_PATH not in sys.path:
    sys.path.insert(0, PYTHONENV_PATH)

# Core numerics / data
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.sparse import csr_matrix, vstack
from scipy import optimize
from scipy.spatial import procrustes
from scipy.linalg import orthogonal_procrustes
import scipy.linalg
import scipy.sparse as sp

# Plotting
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# Scanpy / AnnData ecosystem
import anndata as ad
import scanpy as sc
import scanpy.external as sce
import umap  # optional; remove if unused

# Machine learning (scikit-learn)
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA, NMF
from sklearn.model_selection import train_test_split, GroupKFold
from sklearn.metrics import balanced_accuracy_score, classification_report
from sklearn.preprocessing import OneHotEncoder, MaxAbsScaler
from sklearn.neighbors import KNeighborsClassifier

# Torch
import torch
from torch import nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

# Topic models / factorization (optional; keep only if used)
from cnmf import cNMF, Preprocess
from scETM import scETM, UnsupervisedTrainer, evaluate, prepare_for_transfer

################################################################################################ 

###################################
# Parameters
###################################
par_prep = "CombatSeq"
remove = ["SFTLD", "C9FTLD"]
remove_li = ['C9-FTD']
par_keep_cell_type =  "Fibro"


###################################
# Cell Specific parameters -- Fibro
###################################
par_ann_data_Pineda_BA4 = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_BA4_{par_keep_cell_type}_int.h5ad"
par_ann_data_Pineda_BA9 = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_BA9_{par_keep_cell_type}_int.h5ad"

par_ann_data_Li_BA4 = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_combat_BA4_{par_keep_cell_type}_int.h5ad"
par_ann_data_Li_BA9 = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_combat_BA9_{par_keep_cell_type}_int.h5ad"

par_ann_data_Limone_BA4 = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_combat_BA4_{par_keep_cell_type}_int.h5ad"

###################################
# Load information
###################################
print(par_keep_cell_type)  

#Load LIME-selected genes
lime_file_BA4 = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_All ALS_BA4_{par_keep_cell_type}_0_abs_case_control_narval_2.csv'
lime_file_BA9 = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_All ALS_BA9_{par_keep_cell_type}_0_abs_case_control_narval_2.csv'

features_BA4 = pd.read_csv(lime_file_BA4)["gene"].tolist()
features_BA9 = pd.read_csv(lime_file_BA9)["gene"].tolist()

features_HVG = list(set(features_BA4 + features_BA9))

len(features_HVG)

# Load model report
report_file = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_{par_prep}_SALS_BA4_{par_keep_cell_type}_narval_2.csv"
report_df = pd.read_csv(report_file)
batch_size = int(report_df.at[0, 'batch_size'])
learning_rate = report_df.at[0, 'learning_rate']

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
############################################################################################### Create Merged Pineda object to compute 10000 HVGs

###################################
# Pineda BA4
###################################
adata_pineda_BA4 = sc.read_h5ad(par_ann_data_Pineda_BA4)
adata_pineda_BA4 = adata_pineda_BA4[~adata_pineda_BA4.obs['Group'].isin(remove)]
adata_pineda_BA4 = adata_pineda_BA4[adata_pineda_BA4.obs['CellType'] == par_keep_cell_type]
adata_pineda_BA4 = adata_pineda_BA4[adata_pineda_BA4.obs['Region'] == 'BA4']
adata_pineda_BA4.obs_names = [f"Cell_{i:d}" for i in range(adata_pineda_BA4.n_obs)]
sc.pp.filter_cells(adata_pineda_BA4, min_genes=200) ## <------ moved here from line 272

set(adata_pineda_BA4.obs['Group'])

# Map disease status
mapping = {'C9ALS': 1, 'SALS': 1, 'PN': 0}
adata_pineda_BA4.obs['Group'] = adata_pineda_BA4.obs['Group'].map(mapping)

## Add dataset column
adata_pineda_BA4.obs['Dataset'] = "Pineda"

## Add comprehensive donor column
adata_pineda_BA4.obs['Donor_comp'] = adata_pineda_BA4.obs["Donor"].astype(str) + "_" + adata_pineda_BA4.obs["Region"].astype(str)

set(adata_pineda_BA4.obs['Group'])
set(adata_pineda_BA4.obs['CellType'])
set(adata_pineda_BA4.obs['Region'])
set(adata_pineda_BA4.obs['Dataset'])
set(adata_pineda_BA4.obs['Donor_comp'])


###################################
# Pineda BA9
###################################
adata_pineda_BA9 = sc.read_h5ad(par_ann_data_Pineda_BA9)
adata_pineda_BA9 = adata_pineda_BA9[~adata_pineda_BA9.obs['Group'].isin(remove)]
adata_pineda_BA9 = adata_pineda_BA9[adata_pineda_BA9.obs['CellType'] == par_keep_cell_type]
adata_pineda_BA9 = adata_pineda_BA9[adata_pineda_BA9.obs['Region'] == 'BA9']
adata_pineda_BA9.obs_names = [f"Cell_{i:d}" for i in range(adata_pineda_BA9.n_obs)]
sc.pp.filter_cells(adata_pineda_BA9, min_genes=200) ## <------ moved here from line 272

set(adata_pineda_BA9.obs['Group'])

# Map disease status
mapping = {'C9ALS': 1, 'SALS': 1, 'PN': 0}
adata_pineda_BA9.obs['Group'] = adata_pineda_BA9.obs['Group'].map(mapping)

## Add dataset column
adata_pineda_BA9.obs['Dataset'] = "Pineda"

## Add comprehensive donor column
adata_pineda_BA9.obs['Donor_comp'] = adata_pineda_BA9.obs["Donor"].astype(str) + "_" + adata_pineda_BA9.obs["Region"].astype(str)

set(adata_pineda_BA9.obs['Group'])
set(adata_pineda_BA9.obs['CellType'])
set(adata_pineda_BA9.obs['Region'])
set(adata_pineda_BA9.obs['Dataset'])
set(adata_pineda_BA9.obs['Donor_comp'])

###################################
# Overlapping genes across datasets
###################################
common_genes = set(adata_pineda_BA4.var_names)
len(common_genes)

# Intersect with each of the others
common_genes &= set(adata_pineda_BA9.var_names)
len(common_genes)


keep_genes = common_genes

###################################
# Create a unified object
###################################
keep_genes_1 = [g for g in adata_pineda_BA4.var_names if g in keep_genes]
adata_pineda_BA4 = adata_pineda_BA4[:, keep_genes_1]
len(adata_pineda_BA4.var_names)

keep_genes_1 = [g for g in adata_pineda_BA9.var_names if g in keep_genes]
adata_pineda_BA9 = adata_pineda_BA9[:, keep_genes_1]
len(adata_pineda_BA9.var_names)

# Combine X matrices
X_combined = vstack([adata_pineda_BA4.X, adata_pineda_BA9.X])
X_combined.shape

# Combine obs dataframes 'CellType', 'Sample_ID', 'Group', 'Region', 'Donor',
adata_pineda_BA4.obs = adata_pineda_BA4.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_pineda_BA4.obs.columns.tolist())

adata_pineda_BA9.obs = adata_pineda_BA9.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_pineda_BA9.obs.columns.tolist())

obs_combined = pd.concat([adata_pineda_BA4.obs, adata_pineda_BA9.obs])
obs_combined.shape[0] == X_combined.shape[0]

var_combined = adata_pineda_BA4.var.copy()

# Create new AnnData
adata_combined = ad.AnnData(X=X_combined, obs=obs_combined, var=var_combined)
dupes = adata_combined.obs_names[adata_combined.obs_names.duplicated()]
print(dupes)
adata_combined.obs_names_make_unique()

###################################
# Preprocess combined data
################################### 
sc.pp.filter_genes(adata_combined, min_cells=3)
sc.pp.normalize_total(adata_combined, target_sum=1e4)
sc.pp.log1p(adata_combined)
sc.pp.highly_variable_genes(adata_combined, n_top_genes=10000, flavor='seurat')

keep_genes_NMF = adata_combined.var_names[adata_combined.var['highly_variable']].tolist()


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
############################################################################################### Start from scratch to process everything NMF

###################################
# Pineda BA4
###################################
adata_pineda_BA4 = sc.read_h5ad(par_ann_data_Pineda_BA4)
adata_pineda_BA4 = adata_pineda_BA4[~adata_pineda_BA4.obs['Group'].isin(remove)]
adata_pineda_BA4 = adata_pineda_BA4[adata_pineda_BA4.obs['CellType'] == par_keep_cell_type]
adata_pineda_BA4 = adata_pineda_BA4[adata_pineda_BA4.obs['Region'] == 'BA4']
adata_pineda_BA4.obs_names = [f"Cell_{i:d}" for i in range(adata_pineda_BA4.n_obs)]
sc.pp.filter_cells(adata_pineda_BA4, min_genes=200) ## <------ moved here from line 272

set(adata_pineda_BA4.obs['Group'])

# Map disease status
mapping = {'C9ALS': 1, 'SALS': 1, 'PN': 0}
adata_pineda_BA4.obs['Group'] = adata_pineda_BA4.obs['Group'].map(mapping)

## Add dataset column
adata_pineda_BA4.obs['Dataset'] = "Pineda"

## Add comprehensive donor column
adata_pineda_BA4.obs['Donor_comp'] = adata_pineda_BA4.obs["Donor"].astype(str) + "_" + adata_pineda_BA4.obs["Region"].astype(str)

set(adata_pineda_BA4.obs['Group'])
set(adata_pineda_BA4.obs['CellType'])
set(adata_pineda_BA4.obs['Region'])
set(adata_pineda_BA4.obs['Dataset'])
set(adata_pineda_BA4.obs['Donor_comp'])


###################################
# Pineda BA9
###################################
adata_pineda_BA9 = sc.read_h5ad(par_ann_data_Pineda_BA9)
adata_pineda_BA9 = adata_pineda_BA9[~adata_pineda_BA9.obs['Group'].isin(remove)]
adata_pineda_BA9 = adata_pineda_BA9[adata_pineda_BA9.obs['CellType'] == par_keep_cell_type]
adata_pineda_BA9 = adata_pineda_BA9[adata_pineda_BA9.obs['Region'] == 'BA9']
adata_pineda_BA9.obs_names = [f"Cell_{i:d}" for i in range(adata_pineda_BA9.n_obs)]
sc.pp.filter_cells(adata_pineda_BA9, min_genes=200) ## <------ moved here from line 272

set(adata_pineda_BA9.obs['Group'])

# Map disease status
mapping = {'C9ALS': 1, 'SALS': 1, 'PN': 0}
adata_pineda_BA9.obs['Group'] = adata_pineda_BA9.obs['Group'].map(mapping)

## Add dataset column
adata_pineda_BA9.obs['Dataset'] = "Pineda"

## Add comprehensive donor column
adata_pineda_BA9.obs['Donor_comp'] = adata_pineda_BA9.obs["Donor"].astype(str) + "_" + adata_pineda_BA9.obs["Region"].astype(str)

set(adata_pineda_BA9.obs['Group'])
set(adata_pineda_BA9.obs['CellType'])
set(adata_pineda_BA9.obs['Region'])
set(adata_pineda_BA9.obs['Dataset'])
set(adata_pineda_BA9.obs['Donor_comp'])


###################################
# Li BA4
###################################
adata_li_BA4 = sc.read_h5ad(par_ann_data_Li_BA4)
adata_li_BA4 = adata_li_BA4[~adata_li_BA4.obs['Group'].isin(remove_li)]
adata_li_BA4 = adata_li_BA4[adata_li_BA4.obs['CellType'] == par_keep_cell_type]
adata_li_BA4 = adata_li_BA4[adata_li_BA4.obs['Region'] == 'motor cortex']
adata_li_BA4.obs_names = [f"Cell_{i:d}" for i in range(adata_li_BA4.n_obs)]
sc.pp.filter_cells(adata_li_BA4, min_genes=200) ## <------ moved here from line 272

set(adata_li_BA4.obs['Group'])

# Map disease status
mapping = {'C9-ALS': 1, 'Control': 0}
adata_li_BA4.obs['Group'] = adata_li_BA4.obs['Group'].map(mapping)

## Add dataset column
adata_li_BA4.obs['Dataset'] = "Li"

## Add proper region column
adata_li_BA4.obs['Region'] = "BA4"

## Add comprehensive donor column
adata_li_BA4.obs['Donor_comp'] = adata_li_BA4.obs["Donor"].astype(str) + "_" + adata_li_BA4.obs["Region"].astype(str)

set(adata_li_BA4.obs['Group'])
set(adata_li_BA4.obs['CellType'])
set(adata_li_BA4.obs['Region'])
set(adata_li_BA4.obs['Dataset'])
set(adata_li_BA4.obs['Donor_comp'])

###################################
# Li BA9
###################################
adata_li_BA9 = sc.read_h5ad(par_ann_data_Li_BA9)
adata_li_BA9 = adata_li_BA9[~adata_li_BA9.obs['Group'].isin(remove_li)]
adata_li_BA9 = adata_li_BA9[adata_li_BA9.obs['CellType'] == par_keep_cell_type]
adata_li_BA9 = adata_li_BA9[adata_li_BA9.obs['Region'] == 'medial frontal cortex']
adata_li_BA9.obs_names = [f"Cell_{i:d}" for i in range(adata_li_BA9.n_obs)]
sc.pp.filter_cells(adata_li_BA9, min_genes=200) ## <------ moved here from line 272

set(adata_li_BA9.obs['Group'])

# Map disease status
mapping = {'C9-ALS': 1, 'Control': 0}
adata_li_BA9.obs['Group'] = adata_li_BA9.obs['Group'].map(mapping)

## Add dataset column
adata_li_BA9.obs['Dataset'] = "Li"

## Add proper region column
adata_li_BA9.obs['Region'] = "BA9"

## Add comprehensive donor column
adata_li_BA9.obs['Donor_comp'] = adata_li_BA9.obs["Donor"].astype(str) + "_" + adata_li_BA9.obs["Region"].astype(str)

set(adata_li_BA9.obs['Group'])
set(adata_li_BA9.obs['CellType'])
set(adata_li_BA9.obs['Region'])
set(adata_li_BA9.obs['Dataset'])
set(adata_li_BA9.obs['Donor_comp'])

###################################
# Limone BA4
###################################
adata_limone_BA4 = sc.read_h5ad(par_ann_data_Limone_BA4)
adata_limone_BA4 = adata_limone_BA4[adata_limone_BA4.obs['CellType'] == par_keep_cell_type]
adata_limone_BA4 = adata_limone_BA4[adata_limone_BA4.obs['Region'] == 'BA4']
adata_limone_BA4.obs_names = [f"Cell_{i:d}" for i in range(adata_limone_BA4.n_obs)]
sc.pp.filter_cells(adata_limone_BA4, min_genes=200) ## <------ moved here from line 272

set(adata_limone_BA4.obs['Group'])

# Map disease status
mapping = {'ALS': 1, 'Control': 0}
adata_limone_BA4.obs['Group'] = adata_limone_BA4.obs['Group'].map(mapping)

## Add dataset column
adata_limone_BA4.obs['Dataset'] = "Limone"

## Add proper region column
adata_limone_BA4.obs['Region'] = "BA4"

## Add comprehensive donor column
adata_limone_BA4.obs['Donor_comp'] = adata_limone_BA4.obs["Donor"].astype(str) + "_" + adata_limone_BA4.obs["Region"].astype(str)

set(adata_limone_BA4.obs['Group'])
set(adata_limone_BA4.obs['CellType'])
set(adata_limone_BA4.obs['Region'])
set(adata_limone_BA4.obs['Dataset'])
set(adata_limone_BA4.obs['Donor_comp'])

###################################
# Overlapping genes across datasets
###################################
common_genes = set(adata_pineda_BA4.var_names)
len(common_genes)

# Intersect with each of the others
common_genes &= set(adata_pineda_BA9.var_names)
len(common_genes)
common_genes &= set(adata_li_BA4.var_names)
len(common_genes)
common_genes &= set(adata_li_BA9.var_names)
len(common_genes)
common_genes &= set(adata_limone_BA4.var_names)
len(common_genes)

features = keep_genes_NMF
keep_genes_NMF = common_genes.intersection(features)
len(keep_genes_NMF) ## USE THIS ONE
len(features)

###################################
# Create a unified object
###################################
keep_genes_1 = [g for g in adata_pineda_BA4.var_names if g in keep_genes_NMF]
adata_pineda_BA4 = adata_pineda_BA4[:, keep_genes_1]
len(adata_pineda_BA4.var_names)

keep_genes_1 = [g for g in adata_pineda_BA9.var_names if g in keep_genes_NMF]
adata_pineda_BA9 = adata_pineda_BA9[:, keep_genes_1]
len(adata_pineda_BA9.var_names)

keep_genes_1 = [g for g in adata_li_BA4.var_names if g in keep_genes_NMF]
adata_li_BA4 = adata_li_BA4[:, keep_genes_1]
len(adata_li_BA4.var_names)

keep_genes_1 = [g for g in adata_li_BA9.var_names if g in keep_genes_NMF]
adata_li_BA9 = adata_li_BA9[:, keep_genes_1]
len(adata_li_BA9.var_names)

keep_genes_1 = [g for g in adata_limone_BA4.var_names if g in keep_genes_NMF]
adata_limone_BA4 = adata_limone_BA4[:, keep_genes_1]
len(adata_limone_BA4.var_names)

# Combine X matrices
X_combined = vstack([adata_pineda_BA4.X, adata_pineda_BA9.X, adata_li_BA4.X, adata_li_BA9.X, adata_limone_BA4.X ])
X_combined.shape

# Combine obs dataframes 'CellType', 'Sample_ID', 'Group', 'Region', 'Donor',
adata_pineda_BA4.obs = adata_pineda_BA4.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_pineda_BA4.obs.columns.tolist())

adata_pineda_BA9.obs = adata_pineda_BA9.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_pineda_BA9.obs.columns.tolist())

adata_li_BA4.obs = adata_li_BA4.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_li_BA4.obs.columns.tolist())

adata_li_BA9.obs = adata_li_BA9.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_li_BA9.obs.columns.tolist())

adata_limone_BA4.obs = adata_limone_BA4.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_limone_BA4.obs.columns.tolist())

obs_combined = pd.concat([adata_pineda_BA4.obs, adata_pineda_BA9.obs, adata_li_BA4.obs, adata_li_BA9.obs, adata_limone_BA4.obs])
obs_combined.shape[0] == X_combined.shape[0]

var_combined = adata_pineda_BA4.var.copy()

# Create new AnnData
adata_combined_NMF = ad.AnnData(X=X_combined, obs=obs_combined, var=var_combined)
dupes = adata_combined_NMF.obs_names[adata_combined_NMF.obs_names.duplicated()]
print(dupes)
adata_combined_NMF.obs_names_make_unique()

###################################
# Preprocess combined data
################################### 
sc.pp.filter_genes(adata_combined_NMF, min_cells=3)
sc.pp.normalize_total(adata_combined_NMF, target_sum=1e4)
sc.pp.log1p(adata_combined_NMF)
adata_combined_NMF.raw = adata_combined_NMF

###################################
# Aply combat normalization
###################################
print(adata_combined_NMF.obs["Dataset"].unique())
adata_combined_NMF.X.max() #8.400470713570888
sc.pp.combat(adata_combined_NMF, key='Dataset', covariates=['Group'])
adata_combined_NMF.X.max() #28.953134240526 --> Changes after combat batch correction
adata_combined_NMF.X.min()

###################################
# Global shift to avoid negative values
###################################
X = adata_combined_NMF.X
if sp.issparse(X):
    X = X.toarray()

min_val = X.min()
if min_val < 0:
    X = X - min_val  

adata_combined_NMF.layers["combat_shifted"] = X  # keep a copy
adata_combined_NMF.X = X

adata_combined_NMF.X.max() 
adata_combined_NMF.X.min()

###################################
# Print UMAP
################################### 
sc.pp.pca(adata_combined_NMF, n_comps=50)
sc.pp.neighbors(adata_combined_NMF)       
sc.tl.umap(adata_combined_NMF)  

sc.pl.umap(adata_combined_NMF, color="Group", show=False)
plt.savefig("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", bbox_inches="tight")

###################################
# Subset back out Pineda from the combine object
################################### 
adata_combined_pineda_NMF = adata_combined_NMF[adata_combined_NMF.obs['Dataset'] == "Pineda"]
adata_combined_pineda_NMF.obs



###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
############################################################################################### Start from scratch to process everything HVG

###################################
# Pineda BA4
###################################
adata_pineda_BA4 = sc.read_h5ad(par_ann_data_Pineda_BA4)
adata_pineda_BA4 = adata_pineda_BA4[~adata_pineda_BA4.obs['Group'].isin(remove)]
adata_pineda_BA4 = adata_pineda_BA4[adata_pineda_BA4.obs['CellType'] == par_keep_cell_type]
adata_pineda_BA4 = adata_pineda_BA4[adata_pineda_BA4.obs['Region'] == 'BA4']
adata_pineda_BA4.obs_names = [f"Cell_{i:d}" for i in range(adata_pineda_BA4.n_obs)]
sc.pp.filter_cells(adata_pineda_BA4, min_genes=200) ## <------ moved here from line 272

set(adata_pineda_BA4.obs['Group'])

# Map disease status
mapping = {'C9ALS': 1, 'SALS': 1, 'PN': 0}
adata_pineda_BA4.obs['Group'] = adata_pineda_BA4.obs['Group'].map(mapping)

## Add dataset column
adata_pineda_BA4.obs['Dataset'] = "Pineda"

## Add comprehensive donor column
adata_pineda_BA4.obs['Donor_comp'] = adata_pineda_BA4.obs["Donor"].astype(str) + "_" + adata_pineda_BA4.obs["Region"].astype(str)

set(adata_pineda_BA4.obs['Group'])
set(adata_pineda_BA4.obs['CellType'])
set(adata_pineda_BA4.obs['Region'])
set(adata_pineda_BA4.obs['Dataset'])
set(adata_pineda_BA4.obs['Donor_comp'])


###################################
# Pineda BA9
###################################
adata_pineda_BA9 = sc.read_h5ad(par_ann_data_Pineda_BA9)
adata_pineda_BA9 = adata_pineda_BA9[~adata_pineda_BA9.obs['Group'].isin(remove)]
adata_pineda_BA9 = adata_pineda_BA9[adata_pineda_BA9.obs['CellType'] == par_keep_cell_type]
adata_pineda_BA9 = adata_pineda_BA9[adata_pineda_BA9.obs['Region'] == 'BA9']
adata_pineda_BA9.obs_names = [f"Cell_{i:d}" for i in range(adata_pineda_BA9.n_obs)]
sc.pp.filter_cells(adata_pineda_BA9, min_genes=200) ## <------ moved here from line 272

set(adata_pineda_BA9.obs['Group'])

# Map disease status
mapping = {'C9ALS': 1, 'SALS': 1, 'PN': 0}
adata_pineda_BA9.obs['Group'] = adata_pineda_BA9.obs['Group'].map(mapping)

## Add dataset column
adata_pineda_BA9.obs['Dataset'] = "Pineda"

## Add comprehensive donor column
adata_pineda_BA9.obs['Donor_comp'] = adata_pineda_BA9.obs["Donor"].astype(str) + "_" + adata_pineda_BA9.obs["Region"].astype(str)

set(adata_pineda_BA9.obs['Group'])
set(adata_pineda_BA9.obs['CellType'])
set(adata_pineda_BA9.obs['Region'])
set(adata_pineda_BA9.obs['Dataset'])
set(adata_pineda_BA9.obs['Donor_comp'])


###################################
# Li BA4
###################################
adata_li_BA4 = sc.read_h5ad(par_ann_data_Li_BA4)
adata_li_BA4 = adata_li_BA4[~adata_li_BA4.obs['Group'].isin(remove_li)]
adata_li_BA4 = adata_li_BA4[adata_li_BA4.obs['CellType'] == par_keep_cell_type]
adata_li_BA4 = adata_li_BA4[adata_li_BA4.obs['Region'] == 'motor cortex']
adata_li_BA4.obs_names = [f"Cell_{i:d}" for i in range(adata_li_BA4.n_obs)]
sc.pp.filter_cells(adata_li_BA4, min_genes=200) ## <------ moved here from line 272

set(adata_li_BA4.obs['Group'])

# Map disease status
mapping = {'C9-ALS': 1, 'Control': 0}
adata_li_BA4.obs['Group'] = adata_li_BA4.obs['Group'].map(mapping)

## Add dataset column
adata_li_BA4.obs['Dataset'] = "Li"

## Add proper region column
adata_li_BA4.obs['Region'] = "BA4"

## Add comprehensive donor column
adata_li_BA4.obs['Donor_comp'] = adata_li_BA4.obs["Donor"].astype(str) + "_" + adata_li_BA4.obs["Region"].astype(str)

set(adata_li_BA4.obs['Group'])
set(adata_li_BA4.obs['CellType'])
set(adata_li_BA4.obs['Region'])
set(adata_li_BA4.obs['Dataset'])
set(adata_li_BA4.obs['Donor_comp'])

###################################
# Li BA9
###################################
adata_li_BA9 = sc.read_h5ad(par_ann_data_Li_BA9)
adata_li_BA9 = adata_li_BA9[~adata_li_BA9.obs['Group'].isin(remove_li)]
adata_li_BA9 = adata_li_BA9[adata_li_BA9.obs['CellType'] == par_keep_cell_type]
adata_li_BA9 = adata_li_BA9[adata_li_BA9.obs['Region'] == 'medial frontal cortex']
adata_li_BA9.obs_names = [f"Cell_{i:d}" for i in range(adata_li_BA9.n_obs)]
sc.pp.filter_cells(adata_li_BA9, min_genes=200) ## <------ moved here from line 272

set(adata_li_BA9.obs['Group'])

# Map disease status
mapping = {'C9-ALS': 1, 'Control': 0}
adata_li_BA9.obs['Group'] = adata_li_BA9.obs['Group'].map(mapping)

## Add dataset column
adata_li_BA9.obs['Dataset'] = "Li"

## Add proper region column
adata_li_BA9.obs['Region'] = "BA9"

## Add comprehensive donor column
adata_li_BA9.obs['Donor_comp'] = adata_li_BA9.obs["Donor"].astype(str) + "_" + adata_li_BA9.obs["Region"].astype(str)

set(adata_li_BA9.obs['Group'])
set(adata_li_BA9.obs['CellType'])
set(adata_li_BA9.obs['Region'])
set(adata_li_BA9.obs['Dataset'])
set(adata_li_BA9.obs['Donor_comp'])

###################################
# Limone BA4
###################################
adata_limone_BA4 = sc.read_h5ad(par_ann_data_Limone_BA4)
adata_limone_BA4 = adata_limone_BA4[adata_limone_BA4.obs['CellType'] == par_keep_cell_type]
adata_limone_BA4 = adata_limone_BA4[adata_limone_BA4.obs['Region'] == 'BA4']
adata_limone_BA4.obs_names = [f"Cell_{i:d}" for i in range(adata_limone_BA4.n_obs)]
sc.pp.filter_cells(adata_limone_BA4, min_genes=200) ## <------ moved here from line 272

set(adata_limone_BA4.obs['Group'])

# Map disease status
mapping = {'ALS': 1, 'Control': 0}
adata_limone_BA4.obs['Group'] = adata_limone_BA4.obs['Group'].map(mapping)

## Add dataset column
adata_limone_BA4.obs['Dataset'] = "Limone"

## Add proper region column
adata_limone_BA4.obs['Region'] = "BA4"

## Add comprehensive donor column
adata_limone_BA4.obs['Donor_comp'] = adata_limone_BA4.obs["Donor"].astype(str) + "_" + adata_limone_BA4.obs["Region"].astype(str)

set(adata_limone_BA4.obs['Group'])
set(adata_limone_BA4.obs['CellType'])
set(adata_limone_BA4.obs['Region'])
set(adata_limone_BA4.obs['Dataset'])
set(adata_limone_BA4.obs['Donor_comp'])

###################################
# Overlapping genes across datasets
###################################
common_genes = set(adata_pineda_BA4.var_names)
len(common_genes)

# Intersect with each of the others
common_genes &= set(adata_pineda_BA9.var_names)
len(common_genes)
common_genes &= set(adata_li_BA4.var_names)
len(common_genes)
common_genes &= set(adata_li_BA9.var_names)
len(common_genes)
common_genes &= set(adata_limone_BA4.var_names)
len(common_genes)

features = features_HVG
keep_genes_HVG = common_genes.intersection(features)
len(keep_genes_HVG) ## USE THIS ONE
len(features)

###################################
# Create a unified object
###################################
keep_genes_1 = [g for g in adata_pineda_BA4.var_names if g in keep_genes_HVG]
adata_pineda_BA4 = adata_pineda_BA4[:, keep_genes_1]
len(adata_pineda_BA4.var_names)

keep_genes_1 = [g for g in adata_pineda_BA9.var_names if g in keep_genes_HVG]
adata_pineda_BA9 = adata_pineda_BA9[:, keep_genes_1]
len(adata_pineda_BA9.var_names)

keep_genes_1 = [g for g in adata_li_BA4.var_names if g in keep_genes_HVG]
adata_li_BA4 = adata_li_BA4[:, keep_genes_1]
len(adata_li_BA4.var_names)

keep_genes_1 = [g for g in adata_li_BA9.var_names if g in keep_genes_HVG]
adata_li_BA9 = adata_li_BA9[:, keep_genes_1]
len(adata_li_BA9.var_names)

keep_genes_1 = [g for g in adata_limone_BA4.var_names if g in keep_genes_HVG]
adata_limone_BA4 = adata_limone_BA4[:, keep_genes_1]
len(adata_limone_BA4.var_names)

# Combine X matrices
X_combined = vstack([adata_pineda_BA4.X, adata_pineda_BA9.X, adata_li_BA4.X, adata_li_BA9.X, adata_limone_BA4.X ])
X_combined.shape

# Combine obs dataframes 'CellType', 'Sample_ID', 'Group', 'Region', 'Donor',
adata_pineda_BA4.obs = adata_pineda_BA4.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_pineda_BA4.obs.columns.tolist())

adata_pineda_BA9.obs = adata_pineda_BA9.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_pineda_BA9.obs.columns.tolist())

adata_li_BA4.obs = adata_li_BA4.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_li_BA4.obs.columns.tolist())

adata_li_BA9.obs = adata_li_BA9.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_li_BA9.obs.columns.tolist())

adata_limone_BA4.obs = adata_limone_BA4.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
print(adata_limone_BA4.obs.columns.tolist())

obs_combined = pd.concat([adata_pineda_BA4.obs, adata_pineda_BA9.obs, adata_li_BA4.obs, adata_li_BA9.obs, adata_limone_BA4.obs])
obs_combined.shape[0] == X_combined.shape[0]

var_combined = adata_pineda_BA4.var.copy()

# Create new AnnData
adata_combined_HVG = ad.AnnData(X=X_combined, obs=obs_combined, var=var_combined)
dupes = adata_combined_HVG.obs_names[adata_combined_HVG.obs_names.duplicated()]
print(dupes)
adata_combined_HVG.obs_names_make_unique()
adata_combined_HVG.X.shape

###################################
# Preprocess combined data
################################### 
sc.pp.filter_genes(adata_combined_HVG, min_cells=3)
sc.pp.normalize_total(adata_combined_HVG, target_sum=1e4)
sc.pp.log1p(adata_combined_HVG)
adata_combined_HVG.raw = adata_combined_HVG

###################################
# Aply combat normalization
###################################
print(adata_combined_HVG.obs["Dataset"].unique())
adata_combined_HVG.X.max() #8.400470713570888
sc.pp.combat(adata_combined_HVG, key='Dataset', covariates=['Group'])
adata_combined_HVG.X.max() #28.953134240526 --> Changes after combat batch correction
adata_combined_HVG.X.min()

###################################
# Global shift to avoid negative values
###################################
X = adata_combined_HVG.X
if sp.issparse(X):
    X = X.toarray()

min_val = X.min()
if min_val < 0:
    X = X - min_val  

adata_combined_HVG.layers["combat_shifted"] = X  # keep a copy
adata_combined_HVG.X = X

adata_combined_HVG.X.max() 
adata_combined_HVG.X.min()

###################################
# Print UMAP
################################### 
sc.pp.pca(adata_combined_HVG, n_comps=50)
sc.pp.neighbors(adata_combined_HVG)       
sc.tl.umap(adata_combined_HVG)  

sc.pl.umap(adata_combined_HVG, color="Group", show=False)
plt.savefig("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", bbox_inches="tight")

###################################
# Subset back out Pineda from the combine object
################################### 
adata_combined_pineda_HVG = adata_combined_HVG[adata_combined_HVG.obs['Dataset'] == "Pineda"]
adata_combined_pineda_HVG.obs

###################################
# Define the models
################################### 
# Main model class
class MainModel(nn.Module):
    def __init__(self, input_size, num_classes):
        super().__init__()
        self.shared = nn.Sequential(
            nn.Linear(input_size, 100), nn.ReLU(),
            nn.Linear(100, 50), nn.ReLU(),
            nn.Linear(50, 25), nn.ReLU()
        )
        self.classifier = nn.Linear(25, num_classes)
    
    def forward(self, x):
        shared = self.shared(x)
        return self.classifier(shared), shared

# Main model class
class DomainClassifier(nn.Module):
    def __init__(self, input_size, num_domains):
        super().__init__()
        self.model = nn.Sequential(
            nn.Linear(input_size, 25), nn.ReLU(),
            nn.Linear(25, num_domains)
        )
    
    def forward(self, x):
        return self.model(x)

# Train function for Pineda LOSO with domain donor
def train_main_LOSO(main_model, domain_model, dataloader, epochs, device):
    main_model.to(device)
    domain_model.to(device)
    optimizer_main = optim.Adam(main_model.parameters(), lr=learning_rate)
    optimizer_domain = optim.Adam(domain_model.parameters(), lr=0.001)
    loss_class = nn.CrossEntropyLoss()
    loss_domain = nn.CrossEntropyLoss()
    
    for epoch in range(epochs):
        main_model.train()
        domain_model.train()
        total_loss_c = 0.0
        total_loss_d = 0.0
        num_batches = 0
        correct = 0
        total = 0
        
        for X, y, d in dataloader:
            X, y, d = X.to(device), y.to(device), d.to(device)
            y_out, shared = main_model(X)
            loss_c = loss_class(y_out, y)
            d_out = domain_model(shared.detach())
            loss_d = loss_domain(d_out, d)
            loss = loss_c - loss_d
    
            optimizer_main.zero_grad()
            loss.backward(retain_graph=True)
            optimizer_main.step()
    
            optimizer_domain.zero_grad()
            loss_d.backward()
            optimizer_domain.step()
     
            total_loss_c += loss_c.item()
            total_loss_d += loss_d.item()
            num_batches += 1
        
            # Accuracy computation
            _, predicted = torch.max(y_out, 1)
            correct += (predicted == y).sum().item()
            total += y.size(0)
        
        avg_loss_c = total_loss_c / num_batches
        avg_loss_d = total_loss_d / num_batches
        accuracy = 100.0 * correct / total
        
        print(f"Epoch {epoch}, Avg Class Loss: {avg_loss_c:.4f}, Avg Domain Loss: {avg_loss_d:.4f}, Accuracy: {accuracy:.2f}%")

# Train function for Entire Pineda with domain donor -- need to return the Main model
def train_main_total(main_model, domain_model, dataloader, epochs, device):
    main_model.to(device)
    domain_model.to(device)
    optimizer_main = optim.Adam(main_model.parameters(), lr=learning_rate)
    optimizer_domain = optim.Adam(domain_model.parameters(), lr=0.001)
    loss_class = nn.CrossEntropyLoss()
    loss_domain = nn.CrossEntropyLoss()
    
    for epoch in range(epochs):
        main_model.train()
        domain_model.train()
        total_loss_c = 0.0
        total_loss_d = 0.0
        num_batches = 0
        correct = 0
        total = 0
        
        for X, y, d in dataloader:
            X, y, d = X.to(device), y.to(device), d.to(device)
            y_out, shared = main_model(X)
            loss_c = loss_class(y_out, y)
            d_out = domain_model(shared.detach())
            loss_d = loss_domain(d_out, d)
            loss = loss_c - loss_d
    
            optimizer_main.zero_grad()
            loss.backward(retain_graph=True)
            optimizer_main.step()
    
            optimizer_domain.zero_grad()
            loss_d.backward()
            optimizer_domain.step()
     
            total_loss_c += loss_c.item()
            total_loss_d += loss_d.item()
            num_batches += 1
        
            # Accuracy computation
            _, predicted = torch.max(y_out, 1)
            correct += (predicted == y).sum().item()
            total += y.size(0)
        
        avg_loss_c = total_loss_c / num_batches
        avg_loss_d = total_loss_d / num_batches
        accuracy = 100.0 * correct / total
        
        print(f"Epoch {epoch}, Avg Class Loss: {avg_loss_c:.4f}, Avg Domain Loss: {avg_loss_d:.4f}, Accuracy: {accuracy:.2f}%")
       
    return main_model

# Train function for transfer with domain dataset
def train_transfer(main_model, domain_model, dataloader, epochs, device):
    main_model.to(device)
    domain_model.to(device)
    optimizer_main = optim.Adam(main_model.parameters(), lr=learning_rate)
    optimizer_domain = optim.Adam(domain_model.parameters(), lr=0.001)
    loss_class = nn.CrossEntropyLoss()
    loss_domain = nn.CrossEntropyLoss()
    
    for epoch in range(epochs):
        main_model.train()
        domain_model.train()
        total_loss_c = 0.0
        total_loss_d = 0.0
        num_batches = 0
        correct = 0
        total = 0
        
        for X, y, d in dataloader:
            X, y, d = X.to(device), y.to(device), d.to(device)
            y_out, shared = main_model(X)
            loss_c = loss_class(y_out, y)
            d_out = domain_model(shared.detach())
            loss_d = loss_domain(d_out, d)
            loss = loss_c - loss_d
    
            optimizer_main.zero_grad()
            loss.backward(retain_graph=True)
            optimizer_main.step()
    
            optimizer_domain.zero_grad()
            loss_d.backward()
            optimizer_domain.step()
     
            total_loss_c += loss_c.item()
            total_loss_d += loss_d.item()
            num_batches += 1
        
            # Accuracy computation
            _, predicted = torch.max(y_out, 1)
            correct += (predicted == y).sum().item()
            total += y.size(0)
        
        avg_loss_c = total_loss_c / num_batches
        avg_loss_d = total_loss_d / num_batches
        accuracy = 100.0 * correct / total
        
        print(f"Epoch {epoch}, Avg Class Loss: {avg_loss_c:.4f}, Avg Domain Loss: {avg_loss_d:.4f}, Accuracy: {accuracy:.2f}%")
                
###################################
# Perform KNN LOSO with Pineda for model validation
###################################
num_genes = 100 + adata_combined_pineda_HVG.shape[1]
data = []


set(adata_combined_pineda_HVG.obs['Donor_comp']) == set(adata_combined_pineda_NMF.obs['Donor_comp'])

# LOSO loop
sample_IDs = set(adata_combined_pineda_HVG.obs['Donor_comp'])
#kNN_threshold_list = [.99, .95, .9, .85, .8]
kNN_threshold = 0.9
print(f"KNN Threshold: {kNN_threshold}")
training_epoch = 25

## Check
adata_combined_pineda_HVG.shape #(2524, 1314)
adata_combined_pineda_NMF.shape #(2524, 4753)

## Balance the counts
group_col = 'Group'
counts = adata_combined_pineda_NMF.obs[group_col].value_counts()
min_cells = counts.min()  # target number per group
keep_indices = []

for grp in counts.index:
    grp_indices = np.where(adata_combined_pineda_NMF.obs[group_col] == grp)[0]
    # Randomly choose min_cells from this group
    sampled = np.random.choice(grp_indices, size=min_cells, replace=False)
    keep_indices.extend(sampled)

# Subset the AnnData object
adata_balanced_NMF = adata_combined_pineda_NMF[keep_indices].copy()
adata_balanced_HVG = adata_combined_pineda_HVG[keep_indices].copy()

# Check
adata_balanced_NMF.obs[group_col].value_counts()
adata_balanced_HVG.obs[group_col].value_counts()

# We should downsample pineda -- takes wayyy too long
if adata_balanced_NMF.n_obs > 20000:
    random_indices = np.random.choice(adata_balanced_NMF.n_obs, size=20000, replace=False)
    adata_sub_NMF = adata_balanced_NMF[random_indices, :].copy()
    adata_sub_HVG = adata_balanced_HVG[random_indices, :].copy()
else:
    adata_sub_NMF = adata_balanced_NMF
    adata_sub_HVG = adata_balanced_HVG

adata_sub_NMF.X.shape
adata_sub_HVG.X.shape


for donor in sample_IDs:
    print(f"Processing donor: {donor}")
    
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### prep NMF
    
    adata_train_NMF_pre = adata_sub_NMF[adata_sub_NMF.obs['Donor_comp'] != donor]
    adata_test_NMF_pre = adata_sub_NMF[adata_sub_NMF.obs['Donor_comp'] == donor]
    
    if adata_test_NMF_pre.n_obs == 0:
        print(f"Skipping donor {donor} — no cells in test set.")
        continue  # Skip to next donor
    
    num_cells = adata_test_NMF_pre.shape[0]
    
    ############################
    ## NMF model for train
    ############################
    model = NMF(n_components=100, init='random', tol=1e-3, random_state=42, max_iter = 10000, verbose=True) 
    ## Fit the model
    model.fit(adata_train_NMF_pre.X)
    W = model.transform(adata_train_NMF_pre.X)
    H = model.components_
    print(W.shape)
    print(H.shape)
    H_df = pd.DataFrame(H, columns=adata_train_NMF_pre.var_names)
    H_df
    W_df = pd.DataFrame(W)
    W_df
    
    ## Check class
    type(W), type(adata_train_NMF_pre.X)
        
    ## Create dataframe
    new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
    
    ## Create new anndata object for method 3
    adata_train_NMF = sc.AnnData(X = W.copy(),
    obs = adata_train_NMF_pre.obs.copy(),
    var = new_var,
    uns = adata_train_NMF_pre.uns.copy(),
    obsm = adata_train_NMF_pre.obsm.copy(),
    raw = adata_train_NMF_pre.raw.copy(),
    dtype = "float32",
    shape = None,
    obsp = adata_train_NMF_pre.obsp.copy(),
    varp = adata_train_NMF_pre.varp
    )
    adata_train_NMF.__dict__['_raw'].__dict__['_var'] = adata_train_NMF_pre.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
        
    ############################
    ## Apply NMF model to test
    ############################
    ## NMF model
    W = model.transform(adata_test_NMF_pre.X)
    H = model.components_
    W.shape
    H.shape
    
    H_df = pd.DataFrame(H, columns=adata_test_NMF_pre.var_names)
    H_df
    
    W_df = pd.DataFrame(W)
    W_df
    
    ## Check class
    type(W), type(adata_test_NMF_pre.X)
    
    ## Create dataframe
    new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
    
    ## Create new anndata object for method 3
    adata_test_NMF = sc.AnnData(X = W.copy(),
    obs = adata_test_NMF_pre.obs.copy(),
    var = new_var,
    uns = adata_test_NMF_pre.uns.copy(),
    obsm = adata_test_NMF_pre.obsm.copy(),
    raw = adata_test_NMF_pre.raw.copy(),
    dtype = "float32",
    shape = None,
    obsp = adata_test_NMF_pre.obsp.copy(),
    varp = adata_test_NMF_pre.varp
    )
    adata_test_NMF.__dict__['_raw'].__dict__['_var'] = adata_test_NMF_pre.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
         
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### prep HVG
    
    adata_train_HVG_pre = adata_sub_HVG[adata_sub_HVG.obs['Donor_comp'] != donor]
    adata_test_HVG_pre = adata_sub_HVG[adata_sub_HVG.obs['Donor_comp'] == donor]
    
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### Train and evaluate the model
    ## Prep NMF for merge
    X_train_NMF = adata_train_NMF.X
    X_train_NMF_sparse = sparse.csr_matrix(X_train_NMF)
    
    X_test_NMF = adata_test_NMF.X
    X_test_NMF_sparse = sparse.csr_matrix(X_test_NMF)
    
    ## Prep HVG for merge
    X_train_HVG = adata_train_HVG_pre.X
    X_test_HVG = adata_test_HVG_pre.X
     
    ## combine
    X_combined_train = sparse.hstack([X_train_NMF_sparse, X_train_HVG]).tocsr()
    X_combined_test = sparse.hstack([X_test_NMF_sparse, X_test_HVG]).tocsr()
     
    print(X_combined_train.shape)   
    print(X_combined_test.shape)
     
    for _ in range(5):
        ## Lets try to visualize the embeddings        
        X_train = torch.FloatTensor(X_combined_train.toarray() if hasattr(X_combined_train, 'toarray') else X_combined_train)
        y_train = torch.LongTensor(adata_train_HVG_pre.obs['Group'].values)
        domains = pd.factorize(adata_train_HVG_pre.obs['Donor'])[0] ## Modified
        d_train = torch.LongTensor(domains)
         
        dataset = TensorDataset(X_train, y_train, d_train)
        input_size = X_train.shape[1]
        num_classes = len(np.unique(y_train))
        num_domains = len(np.unique(domains))
        
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
        
        model_main = MainModel(input_size, num_classes)
        model_domain = DomainClassifier(25, num_domains)
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
         
        train_main_LOSO(model_main, model_domain, loader, training_epoch, device)
         
        # Evaluation
        X_test = torch.FloatTensor(X_combined_test.toarray() if hasattr(X_combined_test, 'toarray') else X_combined_test)
        y_test = adata_test_HVG_pre.obs['Group']
     
        model_main.eval()
        with torch.no_grad():
            _, train_embeddings = model_main(X_train.to(device))
            train_embeddings = train_embeddings.cpu().numpy()
            y_train_np = y_train.numpy()
     
            knn = KNeighborsClassifier(n_neighbors=5)
            knn.fit(train_embeddings, y_train_np)
     
            y_logits, test_embeddings = model_main(X_test.to(device))
            test_embeddings = test_embeddings.cpu().numpy()
     
            test_proba = knn.predict_proba(test_embeddings)
            test_preds = knn.predict(test_embeddings)
            test_confidence = test_proba.max(axis=1)
     
            confidence_threshold = kNN_threshold
            keep_indices = np.where(test_confidence >= confidence_threshold)[0]
     
            y_logits_all, _ = model_main(X_test.to(device))
            y_probs_all = torch.softmax(y_logits_all, dim=1).cpu().numpy()
     
            if len(keep_indices) > 0:
                y_logits_conf = y_logits_all[keep_indices]
                y_labels = torch.argmax(y_logits_conf, dim=1)
                y_true = y_test.iloc[keep_indices].values
     
                # Step 1: Compute class prior from high-confidence predictions
                high_conf_labels = y_labels.cpu().numpy()
                class_counts = np.bincount(high_conf_labels, minlength=num_classes)
                class_prior = class_counts / class_counts.sum()
     
                # Accuracy on high-confidence cells
                acc_high = balanced_accuracy_score(y_true, y_labels.numpy())
            else:
                print("No test cells passed the confidence threshold.")
                y_labels = torch.argmax(y_logits_all, dim=1)
                y_true = y_test.values
                class_prior = np.ones(num_classes) / num_classes  # uniform prior if no confident cells
                acc_high = np.nan
     
            print("High-confidence cells:")
            print(classification_report(y_true, y_labels.numpy()))
     
            # Step 2: Apply prior to low-confidence cells
            low_indices = np.setdiff1d(np.arange(len(X_test)), keep_indices)
            if len(low_indices) > 0:
                y_probs_low = y_probs_all[low_indices]
     
                if y_probs_low.shape[0] == 0 or class_prior.shape[0] != y_probs_low.shape[1]:
                    print("Skipping low-confidence adjustment due to shape mismatch:")
                    print(f"class_prior.shape = {class_prior.shape}")
                    print(f"y_probs_low.shape = {y_probs_low.shape}")
                    y_pred_adjusted = np.array([], dtype=int)
                    acc_low = np.nan
                else:
                    y_probs_adjusted = y_probs_low * class_prior[np.newaxis, :]
                    y_probs_adjusted /= y_probs_adjusted.sum(axis=1, keepdims=True)
     
                    y_pred_adjusted = y_probs_adjusted.argmax(axis=1)
                    y_true_low = y_test.iloc[low_indices].values
     
                    print("Low-confidence (prior-influenced) cells:")
                    print(classification_report(y_true_low, y_pred_adjusted))
     
                    # Accuracy on low-confidence cells
                    acc_low = balanced_accuracy_score(y_true_low, y_pred_adjusted)
            else:
                y_pred_adjusted = np.array([], dtype=int)
                acc_low = np.nan
                print("No low-confidence cells to assign.")
     
            # Accuracy on all test cells
            y_all_pred = np.full(len(X_test), -1, dtype=int)
            if len(keep_indices) > 0:
                y_all_pred[keep_indices] = y_labels.numpy()
     
            if len(low_indices) > 0 and y_pred_adjusted.shape[0] == low_indices.shape[0]:
                y_all_pred[low_indices] = y_pred_adjusted
            else:
                print("Skipping assignment due to shape mismatch:")
                print(f"low_indices.shape = {low_indices.shape}")
                print(f"y_pred_adjusted.shape = {y_pred_adjusted.shape}")
     
            valid_indices = y_all_pred != -1
            if valid_indices.any():
                y_all_true = y_test.values[valid_indices]
                y_all_pred = y_all_pred[valid_indices]
                acc_all = balanced_accuracy_score(y_all_true, y_all_pred)
            else:
                acc_all = np.nan
     
            data.append({
                'prep': par_prep, 'donor': donor,
                'group': 'ALS', 'celltype': par_keep_cell_type,
                'n_genes': num_genes, 'n_cells': len(y_true), 'train_epoch': training_epoch,
                'learning_rate': learning_rate, 'batch_size': batch_size,
                'test_accuracy_high': acc_high, 'test_accuracy_low': acc_low,
                'test_accuracy_all': acc_all, 'method': 'LIME + kNN-filtered + prior-adjusted', 'kNN_thresh': kNN_threshold
            })

pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_HVG_and_NMF_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_{par_prep}_ALS_{par_keep_cell_type}.csv"
results_df.to_csv(out_path, index=False)

###################################
# Train model on entire Pineda dataset
###################################
#### New
adata_combined_pineda_NMF.shape #(2524, 4753)
adata_combined_pineda_HVG.shape #(2524, 1314)

## Balance the counts
group_col = 'Group'
counts = adata_combined_pineda_NMF.obs[group_col].value_counts()
min_cells = counts.min()  # target number per group
keep_indices = []

for grp in counts.index:
    grp_indices = np.where(adata_combined_pineda_NMF.obs[group_col] == grp)[0]
    # Randomly choose min_cells from this group
    sampled = np.random.choice(grp_indices, size=min_cells, replace=False)
    keep_indices.extend(sampled)

# Subset the AnnData object
adata_balanced_NMF_pre = adata_combined_pineda_NMF[keep_indices].copy()
adata_balanced_HVG_pre = adata_combined_pineda_HVG[keep_indices].copy()

# Check
adata_balanced_NMF_pre.obs[group_col].value_counts()
adata_balanced_HVG_pre.obs[group_col].value_counts()


######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
###################################################################################################### Prep NMF -- apply to entire Pineda
model = NMF(n_components=100, init='random', tol=1e-3, random_state=42, max_iter = 10000, verbose=True) 

## Fit the model
model.fit(adata_balanced_NMF_pre.X)
W = model.transform(adata_balanced_NMF_pre.X)
H = model.components_
print(W.shape)
print(H.shape)
H_df = pd.DataFrame(H, columns=adata_balanced_NMF_pre.var_names)
H_df
W_df = pd.DataFrame(W)
W_df

## Check class
type(W), type(adata_balanced_NMF_pre.X)
    
## Create dataframe
new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])

## Create new anndata object
adata_balanced_NMF = sc.AnnData(X = W.copy(),
obs = adata_balanced_NMF_pre.obs.copy(),
var = new_var,
uns = adata_balanced_NMF_pre.uns.copy(),
obsm = adata_balanced_NMF_pre.obsm.copy(),
raw = adata_balanced_NMF_pre.raw.copy(),
dtype = "float32",
shape = None,
obsp = adata_balanced_NMF_pre.obsp.copy(),
varp = adata_balanced_NMF_pre.varp
)
adata_balanced_NMF.__dict__['_raw'].__dict__['_var'] = adata_balanced_NMF_pre.__dict__['_raw'].__dict__['_var'].rename(
    columns={'_index': 'features'})
        
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
###################################################################################################### prep HVG

adata_balanced_HVG = adata_balanced_HVG_pre

######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
###################################################################################################### Train model

## Prep NMF for merge
X_train_NMF = adata_balanced_NMF.X
X_train_NMF_sparse = sparse.csr_matrix(X_train_NMF)

## Prep HVG for merge
X_train_HVG = adata_balanced_HVG.X

## combine
X_combined_train = sparse.hstack([X_train_NMF_sparse, X_train_HVG]).tocsr()
  
print(X_combined_train.shape)   

## Apply model to entire Pineda dataset

X_train = torch.FloatTensor(X_combined_train.toarray() if hasattr(X_combined_train, 'toarray') else X_combined_train)
y_train = torch.LongTensor(adata_balanced_NMF.obs['Group'].values)
domains = pd.factorize(adata_balanced_NMF.obs['Donor'])[0]
d_train = torch.LongTensor(domains)

dataset = TensorDataset(X_train, y_train, d_train)
input_size = X_train.shape[1]
num_classes = len(np.unique(y_train))
num_domains = len(np.unique(domains))

loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

model_main = MainModel(input_size, num_classes)
model_domain = DomainClassifier(25, num_domains)
device = 'cuda' if torch.cuda.is_available() else 'cpu'

main_model_total_pineda = train_main_total(model_main, model_domain, loader, 50, device)


###################################
# Apply LOSO transfer learning for Li and Limone -- Unbalanced dataset counts for transfer learning. 
###################################
data2 = []
len(set(adata_combined_NMF.obs['Donor_comp'])) == len(set(adata_combined_HVG.obs['Donor_comp']))

num_genes = 100 + adata_balanced_HVG.X.shape[1]

sample_IDs_combo = set(
    adata_combined_NMF.obs.loc[adata_combined_NMF.obs['Dataset'].isin(['Li', 'Limone']), 'Donor_comp']
)

#kNN_threshold_list = [.99, .95, .9, .85, .8]
kNN_threshold = 0.9
print(f"KNN Threshold: {kNN_threshold}")
training_epoch = 50

## Balance case and control counts in adata_combined for Pineda dataset
group_col = 'Group'

# Mask for Pineda dataset
mask_pineda = adata_combined_NMF.obs['Dataset'] == "Pineda"

# Count cells per group only within Pineda
counts = adata_combined_NMF.obs.loc[mask_pineda, group_col].value_counts()
min_cells = counts.min()

keep_indices = []

# Balance only Pineda cells
for grp in counts.index:
    grp_indices = np.where((adata_combined_NMF.obs[group_col] == grp) & mask_pineda)[0]
    sampled = np.random.choice(grp_indices, size=min_cells, replace=False)
    keep_indices.extend(sampled)

# Add all non-Pineda cells
non_pineda_indices = np.where(~mask_pineda)[0]
keep_indices.extend(non_pineda_indices)

# Subset the AnnData object
adata_combined_balanced_NMF = adata_combined_NMF[keep_indices].copy()
adata_combined_balanced_HVG = adata_combined_HVG[keep_indices].copy()

# Check
print(adata_combined_balanced_NMF.obs.groupby(['Group', 'Dataset']).size())
print(adata_combined_balanced_HVG.obs.groupby(['Group', 'Dataset']).size())

for donor in sample_IDs_combo:
    print(f"Processing donor: {donor}")
    model_temp = model
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### prep NMF
     
    adata_train_NMF_pre = adata_combined_balanced_NMF[adata_combined_balanced_NMF.obs['Donor_comp'] != donor]
    adata_test_NMF_pre = adata_combined_balanced_NMF[adata_combined_balanced_NMF.obs['Donor_comp'] == donor]
    
    if adata_test_NMF_pre.n_obs == 0:
        print(f"Skipping donor {donor} — no cells in test set.")
        continue  # Skip to next donor
    
    num_cells = adata_test_NMF_pre.shape[0]
     
    ############################
    ## transfer NMF model for train
    ############################
    W = model_temp.transform(adata_train_NMF_pre.X)
    H = model_temp.components_
    print(W.shape)
    print(H.shape)
    H_df = pd.DataFrame(H, columns=adata_train_NMF_pre.var_names)
    H_df
    W_df = pd.DataFrame(W)
    W_df
    
    ## Check class
    type(W), type(adata_train_NMF_pre.X)
        
    ## Create dataframe
    new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
    
    ## Create new anndata object for method 3
    adata_train_NMF = sc.AnnData(X = W.copy(),
    obs = adata_train_NMF_pre.obs.copy(),
    var = new_var,
    uns = adata_train_NMF_pre.uns.copy(),
    obsm = adata_train_NMF_pre.obsm.copy(),
    raw = adata_train_NMF_pre.raw.copy(),
    dtype = "float32",
    shape = None,
    obsp = adata_train_NMF_pre.obsp.copy(),
    varp = adata_train_NMF_pre.varp
    )
    adata_train_NMF.__dict__['_raw'].__dict__['_var'] = adata_train_NMF_pre.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
    
    ############################
    ## Apply NMF model to test
    ############################
    ## NMF model
    W = model_temp.transform(adata_test_NMF_pre.X)
    H = model_temp.components_
    W.shape
    H.shape
    
    H_df = pd.DataFrame(H, columns=adata_test_NMF_pre.var_names)
    H_df
    
    W_df = pd.DataFrame(W)
    W_df
    
    ## Check class
    type(W), type(adata_test_NMF_pre.X)
    
    ## Create dataframe
    new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
    
    ## Create new anndata object for method 3
    adata_test_NMF = sc.AnnData(X = W.copy(),
    obs = adata_test_NMF_pre.obs.copy(),
    var = new_var,
    uns = adata_test_NMF_pre.uns.copy(),
    obsm = adata_test_NMF_pre.obsm.copy(),
    raw = adata_test_NMF_pre.raw.copy(),
    dtype = "float32",
    shape = None,
    obsp = adata_test_NMF_pre.obsp.copy(),
    varp = adata_test_NMF_pre.varp
    )
    adata_test_NMF.__dict__['_raw'].__dict__['_var'] = adata_test_NMF_pre.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
     
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### prep HVG
    
    adata_train_HVG_pre = adata_combined_balanced_HVG[adata_combined_balanced_HVG.obs['Donor_comp'] != donor]
    adata_test_HVG_pre = adata_combined_balanced_HVG[adata_combined_balanced_HVG.obs['Donor_comp'] == donor]
    
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### Train and evaluate the model
    ## Prep NMF for merge
    X_train_NMF = adata_train_NMF.X
    X_train_NMF_sparse = sparse.csr_matrix(X_train_NMF)
    
    X_test_NMF = adata_test_NMF.X
    X_test_NMF_sparse = sparse.csr_matrix(X_test_NMF)
    
    ## Prep HVG for merge
    X_train_HVG = adata_train_HVG_pre.X
    X_test_HVG = adata_test_HVG_pre.X
     
    ## combine
    X_combined_train = sparse.hstack([X_train_NMF_sparse, X_train_HVG]).tocsr()
    X_combined_test = sparse.hstack([X_test_NMF_sparse, X_test_HVG]).tocsr()
     
    print(X_combined_train.shape)   
    print(X_combined_test.shape)
    
    for _ in range(5):
        ## Lets try to visualize the embeddings
        ## Lets try to visualize the embeddings        
        X_train = torch.FloatTensor(X_combined_train.toarray() if hasattr(X_combined_train, 'toarray') else X_combined_train)
        y_train = torch.LongTensor(adata_train_HVG_pre.obs['Group'].values)
        domains = pd.factorize(adata_train_HVG_pre.obs['Dataset'])[0] ## Modified
        d_train = torch.LongTensor(domains)
         
        dataset = TensorDataset(X_train, y_train, d_train)
        input_size = X_train.shape[1]
        num_classes = len(np.unique(y_train))
        num_domains = len(np.unique(domains))
        
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
        
        main_model_transfer = MainModel(input_size, num_classes).to(device)
        domain_model_transfer = DomainClassifier(25, num_domains).to(device)
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        
        # Optional: load pretrained weights if you want transfer learning start
        main_model_transfer.load_state_dict(main_model_total_pineda.state_dict())
    
        # Then fine-tune
        train_transfer(main_model_transfer, domain_model_transfer, loader, training_epoch, device) 
       
        # Evaluation
        X_test = torch.FloatTensor(X_combined_test.toarray() if hasattr(X_combined_test, 'toarray') else X_combined_test)
        y_test = adata_test_HVG_pre.obs['Group']
     
        main_model_transfer.eval()
        with torch.no_grad():
            _, train_embeddings = main_model_transfer(X_train.to(device))
            train_embeddings = train_embeddings.cpu().numpy()
            y_train_np = y_train.numpy()
     
            knn = KNeighborsClassifier(n_neighbors=5)
            knn.fit(train_embeddings, y_train_np)
     
            y_logits, test_embeddings = main_model_transfer(X_test.to(device))
            test_embeddings = test_embeddings.cpu().numpy()
     
            test_proba = knn.predict_proba(test_embeddings)
            test_preds = knn.predict(test_embeddings)
            test_confidence = test_proba.max(axis=1)
     
            confidence_threshold = kNN_threshold
            keep_indices = np.where(test_confidence >= confidence_threshold)[0]
     
            y_logits_all, _ = main_model_transfer(X_test.to(device))
            y_probs_all = torch.softmax(y_logits_all, dim=1).cpu().numpy()
     
            if len(keep_indices) > 0:
                y_logits_conf = y_logits_all[keep_indices]
                y_labels = torch.argmax(y_logits_conf, dim=1)
                y_true = y_test.iloc[keep_indices].values
     
                # Step 1: Compute class prior from high-confidence predictions
                high_conf_labels = y_labels.cpu().numpy()
                class_counts = np.bincount(high_conf_labels, minlength=num_classes)
                class_prior = class_counts / class_counts.sum()
     
                # Accuracy on high-confidence cells
                acc_high = balanced_accuracy_score(y_true, y_labels.numpy())
            else:
                print("No test cells passed the confidence threshold.")
                y_labels = torch.argmax(y_logits_all, dim=1)
                y_true = y_test.values
                class_prior = np.ones(num_classes) / num_classes  # uniform prior if no confident cells
                acc_high = np.nan
     
            print("High-confidence cells:")
            print(classification_report(y_true, y_labels.numpy()))
     
            # Step 2: Apply prior to low-confidence cells
            low_indices = np.setdiff1d(np.arange(len(X_test)), keep_indices)
            if len(low_indices) > 0:
                y_probs_low = y_probs_all[low_indices]
     
                if y_probs_low.shape[0] == 0 or class_prior.shape[0] != y_probs_low.shape[1]:
                    print("Skipping low-confidence adjustment due to shape mismatch:")
                    print(f"class_prior.shape = {class_prior.shape}")
                    print(f"y_probs_low.shape = {y_probs_low.shape}")
                    y_pred_adjusted = np.array([], dtype=int)
                    acc_low = np.nan
                else:
                    y_probs_adjusted = y_probs_low * class_prior[np.newaxis, :]
                    y_probs_adjusted /= y_probs_adjusted.sum(axis=1, keepdims=True)
     
                    y_pred_adjusted = y_probs_adjusted.argmax(axis=1)
                    y_true_low = y_test.iloc[low_indices].values
     
                    print("Low-confidence (prior-influenced) cells:")
                    print(classification_report(y_true_low, y_pred_adjusted))
     
                    # Accuracy on low-confidence cells
                    acc_low = balanced_accuracy_score(y_true_low, y_pred_adjusted)
            else:
                y_pred_adjusted = np.array([], dtype=int)
                acc_low = np.nan
                print("No low-confidence cells to assign.")
     
            # Accuracy on all test cells
            y_all_pred = np.full(len(X_test), -1, dtype=int)
            if len(keep_indices) > 0:
                y_all_pred[keep_indices] = y_labels.numpy()
     
            if len(low_indices) > 0 and y_pred_adjusted.shape[0] == low_indices.shape[0]:
                y_all_pred[low_indices] = y_pred_adjusted
            else:
                print("Skipping assignment due to shape mismatch:")
                print(f"low_indices.shape = {low_indices.shape}")
                print(f"y_pred_adjusted.shape = {y_pred_adjusted.shape}")
     
            valid_indices = y_all_pred != -1
            if valid_indices.any():
                y_all_true = y_test.values[valid_indices]
                y_all_pred = y_all_pred[valid_indices]
                acc_all = balanced_accuracy_score(y_all_true, y_all_pred)
            else:
                acc_all = np.nan
     
            data2.append({
                'donor': donor,
                'group': 'ALS', 'celltype': par_keep_cell_type,
                'n_genes': num_genes, 'n_cells': len(y_true), 'train_epoch': training_epoch,
                'learning_rate': learning_rate, 'batch_size': batch_size,
                'test_accuracy_high': acc_high, 'test_accuracy_low': acc_low,
                'test_accuracy_all': acc_all, 'method': 'LIME + kNN-filtered + prior-adjusted', 'kNN_thresh': kNN_threshold
            })

# Save results
pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data2)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_HVG_and_NMF_Li_Limone_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_{par_prep}_ALS_{par_keep_cell_type}_unbalanced.csv"
results_df.to_csv(out_path, index=False)


###################################
# Apply LOSO transfer learning for Li and Limone -- Balanced dataset counts for transfer learning -- 1X maximum. 
###################################
data4 = []
len(set(adata_combined_NMF.obs['Donor_comp'])) == len(set(adata_combined_HVG.obs['Donor_comp']))
num_genes = 100 + adata_balanced_HVG.X.shape[1]

sample_IDs_combo = set(
    adata_combined_NMF.obs.loc[adata_combined_NMF.obs['Dataset'].isin(['Li', 'Limone']), 'Donor_comp']
)

#kNN_threshold_list = [.99, .95, .9, .85, .8]
kNN_threshold = 0.9
print(f"KNN Threshold: {kNN_threshold}")
training_epoch = 50

## Balance case and control counts in adata_combined for Pineda dataset
group_col = 'Group'

# Mask for Pineda dataset
mask_pineda = adata_combined_NMF.obs['Dataset'] == "Pineda"

# Count cells per group only within Pineda
counts = adata_combined_NMF.obs.loc[mask_pineda, group_col].value_counts()
min_cells = counts.min()

keep_indices = []

# Balance only Pineda cells
for grp in counts.index:
    grp_indices = np.where((adata_combined_NMF.obs[group_col] == grp) & mask_pineda)[0]
    sampled = np.random.choice(grp_indices, size=min_cells, replace=False)
    keep_indices.extend(sampled)

# Add all non-Pineda cells
non_pineda_indices = np.where(~mask_pineda)[0]
keep_indices.extend(non_pineda_indices)

# Subset the AnnData object
adata_combined_balanced_NMF = adata_combined_NMF[keep_indices].copy()
adata_combined_balanced_HVG = adata_combined_HVG[keep_indices].copy()

# Check
print(adata_combined_balanced_NMF.obs.groupby(['Group', 'Dataset']).size())
print(adata_combined_balanced_HVG.obs.groupby(['Group', 'Dataset']).size())

## Prepocessed
adata_combined_balanced_NMF.obs['Dataset'].value_counts()

largest_non_pineda = (
    adata_combined_balanced_NMF.obs['Dataset']
    .value_counts()
    .drop('Pineda')
    .max()
)

pineda_mask = adata_combined_balanced_NMF.obs['Dataset'] == 'Pineda'
pineda_indices = np.where(pineda_mask)[0]

np.random.seed(42)  # for reproducibility
downsampled_indices = np.random.choice(pineda_indices, size=1*largest_non_pineda, replace=False)

other_indices = np.where(adata_combined_balanced_NMF.obs['Dataset'] != 'Pineda')[0]

final_indices = np.concatenate([downsampled_indices, other_indices])

adata_combined_balanced_NMF = adata_combined_balanced_NMF[final_indices].copy()
adata_combined_balanced_HVG = adata_combined_balanced_HVG[final_indices].copy()
print(adata_combined_balanced_NMF.obs['Dataset'].value_counts())
print(adata_combined_balanced_HVG.obs['Dataset'].value_counts())

for donor in sample_IDs_combo:
    print(f"Processing donor: {donor}")
    model_temp = model
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### prep NMF
     
    adata_train_NMF_pre = adata_combined_balanced_NMF[adata_combined_balanced_NMF.obs['Donor_comp'] != donor]
    adata_test_NMF_pre = adata_combined_balanced_NMF[adata_combined_balanced_NMF.obs['Donor_comp'] == donor]
    
    if adata_test_NMF_pre.n_obs == 0:
        print(f"Skipping donor {donor} — no cells in test set.")
        continue  # Skip to next donor
    
    num_cells = adata_test_NMF_pre.shape[0]
     
    ############################
    ## transfer NMF model for train
    ############################
    W = model_temp.transform(adata_train_NMF_pre.X)
    H = model_temp.components_
    print(W.shape)
    print(H.shape)
    H_df = pd.DataFrame(H, columns=adata_train_NMF_pre.var_names)
    H_df
    W_df = pd.DataFrame(W)
    W_df
    
    ## Check class
    type(W), type(adata_train_NMF_pre.X)
        
    ## Create dataframe
    new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
    
    ## Create new anndata object for method 3
    adata_train_NMF = sc.AnnData(X = W.copy(),
    obs = adata_train_NMF_pre.obs.copy(),
    var = new_var,
    uns = adata_train_NMF_pre.uns.copy(),
    obsm = adata_train_NMF_pre.obsm.copy(),
    raw = adata_train_NMF_pre.raw.copy(),
    dtype = "float32",
    shape = None,
    obsp = adata_train_NMF_pre.obsp.copy(),
    varp = adata_train_NMF_pre.varp
    )
    adata_train_NMF.__dict__['_raw'].__dict__['_var'] = adata_train_NMF_pre.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
    
    ############################
    ## Apply NMF model to test
    ############################
    ## NMF model
    W = model_temp.transform(adata_test_NMF_pre.X)
    H = model_temp.components_
    W.shape
    H.shape
    
    H_df = pd.DataFrame(H, columns=adata_test_NMF_pre.var_names)
    H_df
    
    W_df = pd.DataFrame(W)
    W_df
    
    ## Check class
    type(W), type(adata_test_NMF_pre.X)
    
    ## Create dataframe
    new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
    
    ## Create new anndata object for method 3
    adata_test_NMF = sc.AnnData(X = W.copy(),
    obs = adata_test_NMF_pre.obs.copy(),
    var = new_var,
    uns = adata_test_NMF_pre.uns.copy(),
    obsm = adata_test_NMF_pre.obsm.copy(),
    raw = adata_test_NMF_pre.raw.copy(),
    dtype = "float32",
    shape = None,
    obsp = adata_test_NMF_pre.obsp.copy(),
    varp = adata_test_NMF_pre.varp
    )
    adata_test_NMF.__dict__['_raw'].__dict__['_var'] = adata_test_NMF_pre.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
     
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### prep HVG
    
    adata_train_HVG_pre = adata_combined_balanced_HVG[adata_combined_balanced_HVG.obs['Donor_comp'] != donor]
    adata_test_HVG_pre = adata_combined_balanced_HVG[adata_combined_balanced_HVG.obs['Donor_comp'] == donor]
    
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### Train and evaluate the model
    ## Prep NMF for merge
    X_train_NMF = adata_train_NMF.X
    X_train_NMF_sparse = sparse.csr_matrix(X_train_NMF)
    
    X_test_NMF = adata_test_NMF.X
    X_test_NMF_sparse = sparse.csr_matrix(X_test_NMF)
    
    ## Prep HVG for merge
    X_train_HVG = adata_train_HVG_pre.X
    X_test_HVG = adata_test_HVG_pre.X
     
    ## combine
    X_combined_train = sparse.hstack([X_train_NMF_sparse, X_train_HVG]).tocsr()
    X_combined_test = sparse.hstack([X_test_NMF_sparse, X_test_HVG]).tocsr()
     
    print(X_combined_train.shape)   
    print(X_combined_test.shape)
    
    for _ in range(5):
        ## Lets try to visualize the embeddings
        ## Lets try to visualize the embeddings        
        X_train = torch.FloatTensor(X_combined_train.toarray() if hasattr(X_combined_train, 'toarray') else X_combined_train)
        y_train = torch.LongTensor(adata_train_HVG_pre.obs['Group'].values)
        domains = pd.factorize(adata_train_HVG_pre.obs['Dataset'])[0] ## Modified
        d_train = torch.LongTensor(domains)
         
        dataset = TensorDataset(X_train, y_train, d_train)
        input_size = X_train.shape[1]
        num_classes = len(np.unique(y_train))
        num_domains = len(np.unique(domains))
        
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
        
        main_model_transfer = MainModel(input_size, num_classes).to(device)
        domain_model_transfer = DomainClassifier(25, num_domains).to(device)
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        
        # Optional: load pretrained weights if you want transfer learning start
        main_model_transfer.load_state_dict(main_model_total_pineda.state_dict())
    
        # Then fine-tune
        train_transfer(main_model_transfer, domain_model_transfer, loader, training_epoch, device) 
       
        # Evaluation
        X_test = torch.FloatTensor(X_combined_test.toarray() if hasattr(X_combined_test, 'toarray') else X_combined_test)
        y_test = adata_test_HVG_pre.obs['Group']
     
        main_model_transfer.eval()
        with torch.no_grad():
            _, train_embeddings = main_model_transfer(X_train.to(device))
            train_embeddings = train_embeddings.cpu().numpy()
            y_train_np = y_train.numpy()
     
            knn = KNeighborsClassifier(n_neighbors=5)
            knn.fit(train_embeddings, y_train_np)
     
            y_logits, test_embeddings = main_model_transfer(X_test.to(device))
            test_embeddings = test_embeddings.cpu().numpy()
     
            test_proba = knn.predict_proba(test_embeddings)
            test_preds = knn.predict(test_embeddings)
            test_confidence = test_proba.max(axis=1)
     
            confidence_threshold = kNN_threshold
            keep_indices = np.where(test_confidence >= confidence_threshold)[0]
     
            y_logits_all, _ = main_model_transfer(X_test.to(device))
            y_probs_all = torch.softmax(y_logits_all, dim=1).cpu().numpy()
     
            if len(keep_indices) > 0:
                y_logits_conf = y_logits_all[keep_indices]
                y_labels = torch.argmax(y_logits_conf, dim=1)
                y_true = y_test.iloc[keep_indices].values
     
                # Step 1: Compute class prior from high-confidence predictions
                high_conf_labels = y_labels.cpu().numpy()
                class_counts = np.bincount(high_conf_labels, minlength=num_classes)
                class_prior = class_counts / class_counts.sum()
     
                # Accuracy on high-confidence cells
                acc_high = balanced_accuracy_score(y_true, y_labels.numpy())
            else:
                print("No test cells passed the confidence threshold.")
                y_labels = torch.argmax(y_logits_all, dim=1)
                y_true = y_test.values
                class_prior = np.ones(num_classes) / num_classes  # uniform prior if no confident cells
                acc_high = np.nan
     
            print("High-confidence cells:")
            print(classification_report(y_true, y_labels.numpy()))
     
            # Step 2: Apply prior to low-confidence cells
            low_indices = np.setdiff1d(np.arange(len(X_test)), keep_indices)
            if len(low_indices) > 0:
                y_probs_low = y_probs_all[low_indices]
     
                if y_probs_low.shape[0] == 0 or class_prior.shape[0] != y_probs_low.shape[1]:
                    print("Skipping low-confidence adjustment due to shape mismatch:")
                    print(f"class_prior.shape = {class_prior.shape}")
                    print(f"y_probs_low.shape = {y_probs_low.shape}")
                    y_pred_adjusted = np.array([], dtype=int)
                    acc_low = np.nan
                else:
                    y_probs_adjusted = y_probs_low * class_prior[np.newaxis, :]
                    y_probs_adjusted /= y_probs_adjusted.sum(axis=1, keepdims=True)
     
                    y_pred_adjusted = y_probs_adjusted.argmax(axis=1)
                    y_true_low = y_test.iloc[low_indices].values
     
                    print("Low-confidence (prior-influenced) cells:")
                    print(classification_report(y_true_low, y_pred_adjusted))
     
                    # Accuracy on low-confidence cells
                    acc_low = balanced_accuracy_score(y_true_low, y_pred_adjusted)
            else:
                y_pred_adjusted = np.array([], dtype=int)
                acc_low = np.nan
                print("No low-confidence cells to assign.")
     
            # Accuracy on all test cells
            y_all_pred = np.full(len(X_test), -1, dtype=int)
            if len(keep_indices) > 0:
                y_all_pred[keep_indices] = y_labels.numpy()
     
            if len(low_indices) > 0 and y_pred_adjusted.shape[0] == low_indices.shape[0]:
                y_all_pred[low_indices] = y_pred_adjusted
            else:
                print("Skipping assignment due to shape mismatch:")
                print(f"low_indices.shape = {low_indices.shape}")
                print(f"y_pred_adjusted.shape = {y_pred_adjusted.shape}")
     
            valid_indices = y_all_pred != -1
            if valid_indices.any():
                y_all_true = y_test.values[valid_indices]
                y_all_pred = y_all_pred[valid_indices]
                acc_all = balanced_accuracy_score(y_all_true, y_all_pred)
            else:
                acc_all = np.nan
     
            data4.append({
                'donor': donor,
                'group': 'ALS', 'celltype': par_keep_cell_type,
                'n_genes': num_genes, 'n_cells': len(y_true), 'train_epoch': training_epoch,
                'learning_rate': learning_rate, 'batch_size': batch_size,
                'test_accuracy_high': acc_high, 'test_accuracy_low': acc_low,
                'test_accuracy_all': acc_all, 'method': 'LIME + kNN-filtered + prior-adjusted', 'kNN_thresh': kNN_threshold
            })

# Save results
pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data4)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_HVG_and_NMF_Li_Limone_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_{par_prep}_ALS_{par_keep_cell_type}_balanced_1X.csv"
results_df.to_csv(out_path, index=False)



###################################
# Apply LOSO transfer learning for Li and Limone -- Balanced dataset counts for transfer learning -- 2X maximum. 
###################################
data3 = []
len(set(adata_combined_NMF.obs['Donor_comp'])) == len(set(adata_combined_HVG.obs['Donor_comp']))
num_genes = 100 + adata_balanced_HVG.X.shape[1]

sample_IDs_combo = set(
    adata_combined_NMF.obs.loc[adata_combined_NMF.obs['Dataset'].isin(['Li', 'Limone']), 'Donor_comp']
)

#kNN_threshold_list = [.99, .95, .9, .85, .8]
kNN_threshold = 0.9
print(f"KNN Threshold: {kNN_threshold}")
training_epoch = 50

## Balance case and control counts in adata_combined for Pineda dataset
group_col = 'Group'

# Mask for Pineda dataset
mask_pineda = adata_combined_NMF.obs['Dataset'] == "Pineda"

# Count cells per group only within Pineda
counts = adata_combined_NMF.obs.loc[mask_pineda, group_col].value_counts()
min_cells = counts.min()

keep_indices = []

# Balance only Pineda cells
for grp in counts.index:
    grp_indices = np.where((adata_combined_NMF.obs[group_col] == grp) & mask_pineda)[0]
    sampled = np.random.choice(grp_indices, size=min_cells, replace=False)
    keep_indices.extend(sampled)

# Add all non-Pineda cells
non_pineda_indices = np.where(~mask_pineda)[0]
keep_indices.extend(non_pineda_indices)

# Subset the AnnData object
adata_combined_balanced_NMF = adata_combined_NMF[keep_indices].copy()
adata_combined_balanced_HVG = adata_combined_HVG[keep_indices].copy()

# Check
print(adata_combined_balanced_NMF.obs.groupby(['Group', 'Dataset']).size())
print(adata_combined_balanced_HVG.obs.groupby(['Group', 'Dataset']).size())

## Prepocessed
adata_combined_balanced_NMF.obs['Dataset'].value_counts()

largest_non_pineda = (
    adata_combined_balanced_NMF.obs['Dataset']
    .value_counts()
    .drop('Pineda')
    .max()
)

pineda_mask = adata_combined_balanced_NMF.obs['Dataset'] == 'Pineda'
pineda_indices = np.where(pineda_mask)[0]

np.random.seed(42)  # for reproducibility
downsampled_indices = np.random.choice(pineda_indices, size=2*largest_non_pineda, replace=False)

other_indices = np.where(adata_combined_balanced_NMF.obs['Dataset'] != 'Pineda')[0]

final_indices = np.concatenate([downsampled_indices, other_indices])

adata_combined_balanced_NMF = adata_combined_balanced_NMF[final_indices].copy()
adata_combined_balanced_HVG = adata_combined_balanced_HVG[final_indices].copy()
print(adata_combined_balanced_NMF.obs['Dataset'].value_counts())
print(adata_combined_balanced_HVG.obs['Dataset'].value_counts())

for donor in sample_IDs_combo:
    print(f"Processing donor: {donor}")
    model_temp = model
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### prep NMF
     
    adata_train_NMF_pre = adata_combined_balanced_NMF[adata_combined_balanced_NMF.obs['Donor_comp'] != donor]
    adata_test_NMF_pre = adata_combined_balanced_NMF[adata_combined_balanced_NMF.obs['Donor_comp'] == donor]
    
    if adata_test_NMF_pre.n_obs == 0:
        print(f"Skipping donor {donor} — no cells in test set.")
        continue  # Skip to next donor
    
    num_cells = adata_test_NMF_pre.shape[0]
     
    ############################
    ## transfer NMF model for train
    ############################
    W = model_temp.transform(adata_train_NMF_pre.X)
    H = model_temp.components_
    print(W.shape)
    print(H.shape)
    H_df = pd.DataFrame(H, columns=adata_train_NMF_pre.var_names)
    H_df
    W_df = pd.DataFrame(W)
    W_df
    
    ## Check class
    type(W), type(adata_train_NMF_pre.X)
        
    ## Create dataframe
    new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
    
    ## Create new anndata object for method 3
    adata_train_NMF = sc.AnnData(X = W.copy(),
    obs = adata_train_NMF_pre.obs.copy(),
    var = new_var,
    uns = adata_train_NMF_pre.uns.copy(),
    obsm = adata_train_NMF_pre.obsm.copy(),
    raw = adata_train_NMF_pre.raw.copy(),
    dtype = "float32",
    shape = None,
    obsp = adata_train_NMF_pre.obsp.copy(),
    varp = adata_train_NMF_pre.varp
    )
    adata_train_NMF.__dict__['_raw'].__dict__['_var'] = adata_train_NMF_pre.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
    
    ############################
    ## Apply NMF model to test
    ############################
    ## NMF model
    W = model_temp.transform(adata_test_NMF_pre.X)
    H = model_temp.components_
    W.shape
    H.shape
    
    H_df = pd.DataFrame(H, columns=adata_test_NMF_pre.var_names)
    H_df
    
    W_df = pd.DataFrame(W)
    W_df
    
    ## Check class
    type(W), type(adata_test_NMF_pre.X)
    
    ## Create dataframe
    new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
    
    ## Create new anndata object for method 3
    adata_test_NMF = sc.AnnData(X = W.copy(),
    obs = adata_test_NMF_pre.obs.copy(),
    var = new_var,
    uns = adata_test_NMF_pre.uns.copy(),
    obsm = adata_test_NMF_pre.obsm.copy(),
    raw = adata_test_NMF_pre.raw.copy(),
    dtype = "float32",
    shape = None,
    obsp = adata_test_NMF_pre.obsp.copy(),
    varp = adata_test_NMF_pre.varp
    )
    adata_test_NMF.__dict__['_raw'].__dict__['_var'] = adata_test_NMF_pre.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
     
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### prep HVG
    
    adata_train_HVG_pre = adata_combined_balanced_HVG[adata_combined_balanced_HVG.obs['Donor_comp'] != donor]
    adata_test_HVG_pre = adata_combined_balanced_HVG[adata_combined_balanced_HVG.obs['Donor_comp'] == donor]
    
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ######################################################################################################
    ###################################################################################################### Train and evaluate the model
    ## Prep NMF for merge
    X_train_NMF = adata_train_NMF.X
    X_train_NMF_sparse = sparse.csr_matrix(X_train_NMF)
    
    X_test_NMF = adata_test_NMF.X
    X_test_NMF_sparse = sparse.csr_matrix(X_test_NMF)
    
    ## Prep HVG for merge
    X_train_HVG = adata_train_HVG_pre.X
    X_test_HVG = adata_test_HVG_pre.X
     
    ## combine
    X_combined_train = sparse.hstack([X_train_NMF_sparse, X_train_HVG]).tocsr()
    X_combined_test = sparse.hstack([X_test_NMF_sparse, X_test_HVG]).tocsr()
     
    print(X_combined_train.shape)   
    print(X_combined_test.shape)
    
    for _ in range(5):
        ## Lets try to visualize the embeddings
        ## Lets try to visualize the embeddings        
        X_train = torch.FloatTensor(X_combined_train.toarray() if hasattr(X_combined_train, 'toarray') else X_combined_train)
        y_train = torch.LongTensor(adata_train_HVG_pre.obs['Group'].values)
        domains = pd.factorize(adata_train_HVG_pre.obs['Dataset'])[0] ## Modified
        d_train = torch.LongTensor(domains)
         
        dataset = TensorDataset(X_train, y_train, d_train)
        input_size = X_train.shape[1]
        num_classes = len(np.unique(y_train))
        num_domains = len(np.unique(domains))
        
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
        
        main_model_transfer = MainModel(input_size, num_classes).to(device)
        domain_model_transfer = DomainClassifier(25, num_domains).to(device)
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        
        # Optional: load pretrained weights if you want transfer learning start
        main_model_transfer.load_state_dict(main_model_total_pineda.state_dict())
    
        # Then fine-tune
        train_transfer(main_model_transfer, domain_model_transfer, loader, training_epoch, device) 
       
        # Evaluation
        X_test = torch.FloatTensor(X_combined_test.toarray() if hasattr(X_combined_test, 'toarray') else X_combined_test)
        y_test = adata_test_HVG_pre.obs['Group']
     
        main_model_transfer.eval()
        with torch.no_grad():
            _, train_embeddings = main_model_transfer(X_train.to(device))
            train_embeddings = train_embeddings.cpu().numpy()
            y_train_np = y_train.numpy()
     
            knn = KNeighborsClassifier(n_neighbors=5)
            knn.fit(train_embeddings, y_train_np)
     
            y_logits, test_embeddings = main_model_transfer(X_test.to(device))
            test_embeddings = test_embeddings.cpu().numpy()
     
            test_proba = knn.predict_proba(test_embeddings)
            test_preds = knn.predict(test_embeddings)
            test_confidence = test_proba.max(axis=1)
     
            confidence_threshold = kNN_threshold
            keep_indices = np.where(test_confidence >= confidence_threshold)[0]
     
            y_logits_all, _ = main_model_transfer(X_test.to(device))
            y_probs_all = torch.softmax(y_logits_all, dim=1).cpu().numpy()
     
            if len(keep_indices) > 0:
                y_logits_conf = y_logits_all[keep_indices]
                y_labels = torch.argmax(y_logits_conf, dim=1)
                y_true = y_test.iloc[keep_indices].values
     
                # Step 1: Compute class prior from high-confidence predictions
                high_conf_labels = y_labels.cpu().numpy()
                class_counts = np.bincount(high_conf_labels, minlength=num_classes)
                class_prior = class_counts / class_counts.sum()
     
                # Accuracy on high-confidence cells
                acc_high = balanced_accuracy_score(y_true, y_labels.numpy())
            else:
                print("No test cells passed the confidence threshold.")
                y_labels = torch.argmax(y_logits_all, dim=1)
                y_true = y_test.values
                class_prior = np.ones(num_classes) / num_classes  # uniform prior if no confident cells
                acc_high = np.nan
     
            print("High-confidence cells:")
            print(classification_report(y_true, y_labels.numpy()))
     
            # Step 2: Apply prior to low-confidence cells
            low_indices = np.setdiff1d(np.arange(len(X_test)), keep_indices)
            if len(low_indices) > 0:
                y_probs_low = y_probs_all[low_indices]
     
                if y_probs_low.shape[0] == 0 or class_prior.shape[0] != y_probs_low.shape[1]:
                    print("Skipping low-confidence adjustment due to shape mismatch:")
                    print(f"class_prior.shape = {class_prior.shape}")
                    print(f"y_probs_low.shape = {y_probs_low.shape}")
                    y_pred_adjusted = np.array([], dtype=int)
                    acc_low = np.nan
                else:
                    y_probs_adjusted = y_probs_low * class_prior[np.newaxis, :]
                    y_probs_adjusted /= y_probs_adjusted.sum(axis=1, keepdims=True)
     
                    y_pred_adjusted = y_probs_adjusted.argmax(axis=1)
                    y_true_low = y_test.iloc[low_indices].values
     
                    print("Low-confidence (prior-influenced) cells:")
                    print(classification_report(y_true_low, y_pred_adjusted))
     
                    # Accuracy on low-confidence cells
                    acc_low = balanced_accuracy_score(y_true_low, y_pred_adjusted)
            else:
                y_pred_adjusted = np.array([], dtype=int)
                acc_low = np.nan
                print("No low-confidence cells to assign.")
     
            # Accuracy on all test cells
            y_all_pred = np.full(len(X_test), -1, dtype=int)
            if len(keep_indices) > 0:
                y_all_pred[keep_indices] = y_labels.numpy()
     
            if len(low_indices) > 0 and y_pred_adjusted.shape[0] == low_indices.shape[0]:
                y_all_pred[low_indices] = y_pred_adjusted
            else:
                print("Skipping assignment due to shape mismatch:")
                print(f"low_indices.shape = {low_indices.shape}")
                print(f"y_pred_adjusted.shape = {y_pred_adjusted.shape}")
     
            valid_indices = y_all_pred != -1
            if valid_indices.any():
                y_all_true = y_test.values[valid_indices]
                y_all_pred = y_all_pred[valid_indices]
                acc_all = balanced_accuracy_score(y_all_true, y_all_pred)
            else:
                acc_all = np.nan
     
            data3.append({
                'donor': donor,
                'group': 'ALS', 'celltype': par_keep_cell_type,
                'n_genes': num_genes, 'n_cells': len(y_true), 'train_epoch': training_epoch,
                'learning_rate': learning_rate, 'batch_size': batch_size,
                'test_accuracy_high': acc_high, 'test_accuracy_low': acc_low,
                'test_accuracy_all': acc_all, 'method': 'LIME + kNN-filtered + prior-adjusted', 'kNN_thresh': kNN_threshold
            })

# Save results
pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data3)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_HVG_and_NMF_Li_Limone_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_{par_prep}_ALS_{par_keep_cell_type}_balanced_2X.csv"
results_df.to_csv(out_path, index=False)


