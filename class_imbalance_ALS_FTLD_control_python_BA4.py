salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

module load StdEnv/2020 
module load python/3.8.10
module load r/4.3.1

export R_HOME=$(R RHOME) #2

PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  


#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("sva")

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

nano class_imbalance_ALS_FTLD_control_python_BA4.sh

#!/bin/bash  
#SBATCH --account=def-sfarhan
#SBATCH --time=00-05:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=100g          # memory per cor
#SBATCH --job-name=class_imbalance_ALS_FTLD_control_python_BA4
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
module load r/4.3.1

export R_HOME=$(R RHOME)

PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/class_imbalance_ALS_FTLD_control_python_BA4.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano class_imbalance_ALS_FTLD_control_python_BA4.py


################################################################## clean version ------ USE THIS ONE

import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.model_selection import train_test_split
from sklearn.metrics import balanced_accuracy_score, classification_report
from torch.utils.data import DataLoader, TensorDataset
from torch_geometric.data import Data as GeoData
from torch_geometric.nn import GCNConv
from sklearn.preprocessing import MaxAbsScaler
from sklearn.neighbors import KNeighborsClassifier
import matplotlib.pyplot as plt
from scipy.sparse import vstack
import torch

import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Function
from rpy2.robjects import pandas2ri, r
import rpy2.robjects.packages as rpackages

#if not rpackages.isinstalled('sva'):
#    utils = rpackages.importr('utils')
#    utils.install_packages('sva')


import rpy2.robjects as ro
from rpy2.robjects import Formula
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

from scipy.spatial.distance import euclidean

sva = rpackages.importr('sva')
base = rpackages.importr('base')

from itertools import combinations
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

#import pandas as pd
#from rpy2.robjects import pandas2ri
#pandas2ri.activate()

################################################################################################ 

###################################
# Parameters
###################################
par_prep = "CombatSeq"

data = []

cell_types = [
    "L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6",
    "PV", "5HT3aR", "Rosehip", "SOM",
    "Oligo", "Astro", "OPC", "Micro", "Mural", "Endo", "Fibro", "L5"
]
   
for par_keep_cell_type in cell_types:
    print(f"Processing {par_keep_cell_type}")
    
    ###################################
    # Cell Specific parameters -- Rosehip
    ###################################
    par_ann_data_Pineda_BA4 = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_BA4_{par_keep_cell_type}_int.h5ad"
    
    par_ann_data_Li_BA4 = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_combat_BA4_{par_keep_cell_type}_int.h5ad"
    
    par_ann_data_Limone_BA4 = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_combat_BA4_{par_keep_cell_type}_int.h5ad"
    
    ###################################
    # Load information
    ###################################
    print(par_keep_cell_type)  
    
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
    ############################################################################################### Create Merged Pineda object to compute HVGs
    
    ###################################
    # Pineda BA4
    ###################################
    adata_pineda_BA4 = sc.read_h5ad(par_ann_data_Pineda_BA4)
    adata_pineda_BA4 = adata_pineda_BA4[adata_pineda_BA4.obs['CellType'] == par_keep_cell_type]
    adata_pineda_BA4 = adata_pineda_BA4[adata_pineda_BA4.obs['Region'] == 'BA4']
    adata_pineda_BA4.obs_names = [f"Cell_{i:d}" for i in range(adata_pineda_BA4.n_obs)]
    sc.pp.filter_cells(adata_pineda_BA4, min_genes=200) ## <------ moved here from line 272
    
    set(adata_pineda_BA4.obs['Group'])
    
    # Map disease status
    mapping = {'C9ALS': 1, 'SALS': 1, 'SFTLD': 2, 'C9FTLD': 2,  'PN': 0}
    adata_pineda_BA4.obs['Group'] = adata_pineda_BA4.obs['Group'].map(mapping)
    
    ## Add dataset column
    adata_pineda_BA4.obs['Dataset'] = "Pineda"
    
    ## Add comprehensive donor column
    adata_pineda_BA4.obs['Donor_comp'] = adata_pineda_BA4.obs["Donor"].astype(str) + "_" + adata_pineda_BA4.obs["Region"].astype(str) + "_" + adata_pineda_BA4.obs["Group"].astype(str)
    
    set(adata_pineda_BA4.obs['Group'])
    set(adata_pineda_BA4.obs['CellType'])
    set(adata_pineda_BA4.obs['Region'])
    set(adata_pineda_BA4.obs['Dataset'])
    set(adata_pineda_BA4.obs['Donor_comp'])
    
    ###################################
    # Preprocess combined data
    ################################### 
    sc.pp.filter_genes(adata_pineda_BA4, min_cells=3)
    sc.pp.normalize_total(adata_pineda_BA4, target_sum=1e4)
    sc.pp.log1p(adata_pineda_BA4)
    sc.pp.highly_variable_genes(adata_pineda_BA4, flavor='seurat')
    
    features = adata_pineda_BA4.var_names[adata_pineda_BA4.var['highly_variable']].tolist()
    len(features)
    
    ###############################################################################################
    ###############################################################################################
    ###############################################################################################
    ###############################################################################################
    ###############################################################################################
    ############################################################################################### Start from scratch
    
    ###################################
    # Pineda BA4
    ###################################
    adata_pineda_BA4 = sc.read_h5ad(par_ann_data_Pineda_BA4)
    adata_pineda_BA4 = adata_pineda_BA4[adata_pineda_BA4.obs['CellType'] == par_keep_cell_type]
    adata_pineda_BA4 = adata_pineda_BA4[adata_pineda_BA4.obs['Region'] == 'BA4']
    adata_pineda_BA4.obs_names = [f"Cell_{i:d}" for i in range(adata_pineda_BA4.n_obs)]
    sc.pp.filter_cells(adata_pineda_BA4, min_genes=200) ## <------ moved here from line 272
    
    set(adata_pineda_BA4.obs['Group'])
    
    # Map disease status
    mapping = {'C9ALS': 1, 'SALS': 1, 'SFTLD': 2, 'C9FTLD': 2,  'PN': 0}
    adata_pineda_BA4.obs['Group'] = adata_pineda_BA4.obs['Group'].map(mapping)
    
    ## Add dataset column
    adata_pineda_BA4.obs['Dataset'] = "Pineda"
    
    ## Add comprehensive donor column
    adata_pineda_BA4.obs['Donor_comp'] = adata_pineda_BA4.obs["Donor"].astype(str) + "_" + adata_pineda_BA4.obs["Region"].astype(str) + "_" + adata_pineda_BA4.obs["Group"].astype(str)
    
    set(adata_pineda_BA4.obs['Group'])
    set(adata_pineda_BA4.obs['CellType'])
    set(adata_pineda_BA4.obs['Region'])
    set(adata_pineda_BA4.obs['Dataset'])
    set(adata_pineda_BA4.obs['Donor_comp'])
    
    ###################################
    # Li BA4
    ###################################
    adata_li_BA4 = sc.read_h5ad(par_ann_data_Li_BA4)
    adata_li_BA4 = adata_li_BA4[adata_li_BA4.obs['CellType'] == par_keep_cell_type]
    adata_li_BA4 = adata_li_BA4[adata_li_BA4.obs['Region'] == 'motor cortex']
    adata_li_BA4.obs_names = [f"Cell_{i:d}" for i in range(adata_li_BA4.n_obs)]
    sc.pp.filter_cells(adata_li_BA4, min_genes=200) ## <------ moved here from line 272
    
    set(adata_li_BA4.obs['Group'])
    
    # Map disease status
    mapping = {'C9-ALS': 1 'C9-FTD': 2,  'Control': 0}
    adata_li_BA4.obs['Group'] = adata_li_BA4.obs['Group'].map(mapping)
    
    ## Add dataset column
    adata_li_BA4.obs['Dataset'] = "Li"
    
    ## Add proper region column
    adata_li_BA4.obs['Region'] = "BA4"
    
    ## Add comprehensive donor column
    adata_li_BA4.obs['Donor_comp'] = adata_li_BA4.obs["Donor"].astype(str) + "_" + adata_li_BA4.obs["Region"].astype(str) + "_" + adata_li_BA4.obs["Group"].astype(str)
    
    set(adata_li_BA4.obs['Group'])
    set(adata_li_BA4.obs['CellType'])
    set(adata_li_BA4.obs['Region'])
    set(adata_li_BA4.obs['Dataset'])
    set(adata_li_BA4.obs['Donor_comp'])
    
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
    adata_limone_BA4.obs['Donor_comp'] = adata_limone_BA4.obs["Donor"].astype(str) + "_" + adata_limone_BA4.obs["Region"].astype(str) + "_" + adata_limone_BA4.obs["Group"].astype(str)
    
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
    common_genes &= set(adata_li_BA4.var_names)
    len(common_genes)
    common_genes &= set(adata_limone_BA4.var_names)
    len(common_genes)
    
    keep_genes = common_genes.intersection(features)
    len(keep_genes) ## USE THIS ONE
    len(features)
    
    ###################################
    # Create a unified object
    ###################################
    keep_genes_1 = [g for g in adata_pineda_BA4.var_names if g in keep_genes]
    adata_pineda_BA4 = adata_pineda_BA4[:, keep_genes_1]
    len(adata_pineda_BA4.var_names)
    
    keep_genes_1 = [g for g in adata_li_BA4.var_names if g in keep_genes]
    adata_li_BA4 = adata_li_BA4[:, keep_genes_1]
    len(adata_li_BA4.var_names)
    
    keep_genes_1 = [g for g in adata_limone_BA4.var_names if g in keep_genes]
    adata_limone_BA4 = adata_limone_BA4[:, keep_genes_1]
    len(adata_limone_BA4.var_names)
    
    # Combine X matrices
    X_combined = vstack([adata_pineda_BA4.X, adata_li_BA4.X, adata_limone_BA4.X ])
    X_combined.shape
    
    # Combine obs dataframes 'CellType', 'Sample_ID', 'Group', 'Region', 'Donor',
    adata_pineda_BA4.obs = adata_pineda_BA4.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
    print(adata_pineda_BA4.obs.columns.tolist())
    
    adata_li_BA4.obs = adata_li_BA4.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
    print(adata_li_BA4.obs.columns.tolist())
    
    adata_limone_BA4.obs = adata_limone_BA4.obs[['CellType', 'Sample_ID', 'Group', 'Region', 'Donor', 'Dataset', 'Donor_comp', 'nCount_RNA', 'nFeature_RNA']]
    print(adata_limone_BA4.obs.columns.tolist())
    
    obs_combined = pd.concat([adata_pineda_BA4.obs, adata_li_BA4.obs, adata_limone_BA4.obs])
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
    sc.pp.highly_variable_genes(adata_combined)
    adata_combined.raw = adata_combined
    
    adata_combined.obs.shape
    
    target_n_cells = 46000
    current_n_cells = adata_combined.n_obs
    
    if current_n_cells > target_n_cells:
        # Randomly select indices to keep
        np.random.seed(42)  # for reproducibility
        keep_indices = np.random.choice(current_n_cells, size=target_n_cells, replace=False)
        
        # Subset the AnnData object
        adata_combined = adata_combined[keep_indices].copy()
        print(f"Downsampled from {current_n_cells} to {adata_combined.n_obs} cells")
    else:
        adata_combined = adata_combined.copy()
        print(f"No downsampling needed, {current_n_cells} cells <= {target_n_cells}")
    
    adata_combined.obs.shape
    
    ###################################
    # Compute Shannon entropy
    ###################################  
    
    # Get proportions of each group
    group_counts = adata_combined.obs["Group"].value_counts(normalize=True)
    p = group_counts.values
    
    # Compute Shannon entropy
    H = -(p * np.log(p)).sum()
    
    # Normalize by number of groups
    H_norm = H / np.log(len(p)) if len(p) > 1 else 0.0  # handle single-group edge case
    
    # Convert group proportions into a dictionary (GroupName: Proportion)
    p_dict = group_counts.to_dict()
    
    # Append results
    row = {
        'celltype': par_keep_cell_type,
        'input': "Diagnosis",
        'norm_shannon_entropy': H_norm,
    }
    
    # Merge the proportions into the row
    row.update(p_dict)
    
    data.append(row)

pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/class_imbalance_ALS_FTLD_control_python_BA4.csv"
results_df.to_csv(out_path, index=False)