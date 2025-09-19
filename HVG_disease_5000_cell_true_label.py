salloc -A def-tdurcan --time=0-2 -c 1 --mem=40g

module load StdEnv/2020 
module load python/3.8.10
module load r/4.3.1

export R_HOME=$(R RHOME) #2

PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  


if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("sva")

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ automated submission script

nano HVG_disease_5000_cell_true_label_Astro.sh

celltypes=(
"L3_L5"
"L2_L3"
"L4_L6"
"L4_L5"
"L5_L6"
"L6"
"PV"
"5HT3aR"
"Rosehip"
"SOM"
"Oligo"
"Astro"
"OPC"
"Micro"
"Mural"
"Endo"
"Fibro"
"L5"
)

SCRIPT_DIR="/home/fiorini9/scratch/machine_learning_ALS/scripts"

for ct in "${celltypes[@]}"; do
    echo "Creating job script and python script for ${ct}..."

    # Create the SLURM job script
    jobfile="${SCRIPT_DIR}/HVG_disease_5000_cell_true_label_${ct}.sh"
    cat > "$jobfile" <<EOF
#!/bin/bash
#SBATCH --account=def-grouleau
#SBATCH --time=00-05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=40g
#SBATCH --job-name=HVG_disease_5000_cell_true_label_${ct}
#SBATCH --error=${SCRIPT_DIR}/temp_error/job.%x-%j.err
#SBATCH --output=${SCRIPT_DIR}/temp_error/job.%x-%j.out

module load StdEnv/2020
module load python/3.8.10
module load r/4.3.1

export R_HOME=\$(R RHOME)

PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
source \$PYTHONENV0/bin/activate
export PYTHONPATH=\$PYTHONENV0/lib/python3.8/site-packages

cd ${SCRIPT_DIR}
python3.8 ${SCRIPT_DIR}/HVG_disease_5000_cell_true_label_${ct}.py
EOF
    
    # Create the SLURM job script
    pyfile="${SCRIPT_DIR}/HVG_disease_5000_cell_true_label_${ct}.py"
    cat > "$pyfile" <<EOF

## load libraries
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


sva = rpackages.importr('sva')
base = rpackages.importr('base')

#import pandas as pd
#from rpy2.robjects import pandas2ri
#pandas2ri.activate()

################################################################################################ 

###################################
# Parameters
###################################
par_prep = "CombatSeq"
#remove = ["SFTLD", "C9FTLD"]
#remove_li = ['C9-FTD']
par_keep_cell_type = "${ct}"


###################################
# Cell Specific parameters -- L4_L6
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
#lime_file_BA4 = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_All ALS_BA4_{par_keep_cell_type}_0_abs_case_control_narval_2.csv'
#lime_file_BA9 = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_All ALS_BA9_{par_keep_cell_type}_0_abs_case_control_narval_2.csv'

#features_BA4 = pd.read_csv(lime_file_BA4)["gene"].tolist()
#features_BA9 = pd.read_csv(lime_file_BA9)["gene"].tolist()

#features = list(set(features_BA4 + features_BA9))

#len(features)

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
mapping = {'C9ALS': 2, 'SALS': 1, 'SFTLD': 3, 'C9FTLD': 4,  'PN': 0}
adata_pineda_BA4.obs['Group'] = adata_pineda_BA4.obs['Group'].map(mapping)

## Add dataset column
adata_pineda_BA4.obs['Dataset'] = "Pineda"

## Add comprehensive donor column
adata_pineda_BA4.obs['Donor_comp'] = adata_pineda_BA4.obs["Donor"].astype(str) + "_" + adata_pineda_BA4.obs["Group"].astype(str)

set(adata_pineda_BA4.obs['Group'])
set(adata_pineda_BA4.obs['CellType'])
set(adata_pineda_BA4.obs['Region'])
set(adata_pineda_BA4.obs['Dataset'])
set(adata_pineda_BA4.obs['Donor_comp'])


###################################
# Pineda BA9
###################################
adata_pineda_BA9 = sc.read_h5ad(par_ann_data_Pineda_BA9)
adata_pineda_BA9 = adata_pineda_BA9[adata_pineda_BA9.obs['CellType'] == par_keep_cell_type]
adata_pineda_BA9 = adata_pineda_BA9[adata_pineda_BA9.obs['Region'] == 'BA9']
adata_pineda_BA9.obs_names = [f"Cell_{i:d}" for i in range(adata_pineda_BA9.n_obs)]
sc.pp.filter_cells(adata_pineda_BA9, min_genes=200) ## <------ moved here from line 272

set(adata_pineda_BA9.obs['Group'])

# Map disease status
mapping = {'C9ALS': 2, 'SALS': 1, 'SFTLD': 3, 'C9FTLD': 4,  'PN': 0}
adata_pineda_BA9.obs['Group'] = adata_pineda_BA9.obs['Group'].map(mapping)

## Add dataset column
adata_pineda_BA9.obs['Dataset'] = "Pineda"

## Add comprehensive donor column
adata_pineda_BA9.obs['Donor_comp'] = adata_pineda_BA9.obs["Donor"].astype(str) + "_"  + adata_pineda_BA9.obs["Group"].astype(str)

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
sc.pp.highly_variable_genes(adata_combined, flavor='seurat')

features = adata_combined.var_names[adata_combined.var['highly_variable']].tolist()
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
mapping = {'C9ALS': 2, 'SALS': 1, 'SFTLD': 3, 'C9FTLD': 4,  'PN': 0}
adata_pineda_BA4.obs['Group'] = adata_pineda_BA4.obs['Group'].map(mapping)

## Add dataset column
adata_pineda_BA4.obs['Dataset'] = "Pineda"

## Add comprehensive donor column
adata_pineda_BA4.obs['Donor_comp'] = adata_pineda_BA4.obs["Donor"].astype(str) + "_" + adata_pineda_BA4.obs["Group"].astype(str)

set(adata_pineda_BA4.obs['Group'])
set(adata_pineda_BA4.obs['CellType'])
set(adata_pineda_BA4.obs['Region'])
set(adata_pineda_BA4.obs['Dataset'])
set(adata_pineda_BA4.obs['Donor_comp'])


###################################
# Pineda BA9
###################################
adata_pineda_BA9 = sc.read_h5ad(par_ann_data_Pineda_BA9)
adata_pineda_BA9 = adata_pineda_BA9[adata_pineda_BA9.obs['CellType'] == par_keep_cell_type]
adata_pineda_BA9 = adata_pineda_BA9[adata_pineda_BA9.obs['Region'] == 'BA9']
adata_pineda_BA9.obs_names = [f"Cell_{i:d}" for i in range(adata_pineda_BA9.n_obs)]
sc.pp.filter_cells(adata_pineda_BA9, min_genes=200) ## <------ moved here from line 272

set(adata_pineda_BA9.obs['Group'])

# Map disease status
mapping = {'C9ALS': 2, 'SALS': 1, 'SFTLD': 3, 'C9FTLD': 4,  'PN': 0}
adata_pineda_BA9.obs['Group'] = adata_pineda_BA9.obs['Group'].map(mapping)

## Add dataset column
adata_pineda_BA9.obs['Dataset'] = "Pineda"

## Add comprehensive donor column
adata_pineda_BA9.obs['Donor_comp'] = adata_pineda_BA9.obs["Donor"].astype(str) + "_"  + adata_pineda_BA9.obs["Group"].astype(str)

set(adata_pineda_BA9.obs['Group'])
set(adata_pineda_BA9.obs['CellType'])
set(adata_pineda_BA9.obs['Region'])
set(adata_pineda_BA9.obs['Dataset'])
set(adata_pineda_BA9.obs['Donor_comp'])

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
mapping = {'C9-ALS': 2, 'C9-FTD': 4,  'Control': 0}
adata_li_BA4.obs['Group'] = adata_li_BA4.obs['Group'].map(mapping)

## Add dataset column
adata_li_BA4.obs['Dataset'] = "Li"

## Add proper region column
adata_li_BA4.obs['Region'] = "BA4"

## Add comprehensive donor column
adata_li_BA4.obs['Donor_comp'] = adata_li_BA4.obs["Donor"].astype(str) + "_" + adata_li_BA4.obs["Group"].astype(str)

set(adata_li_BA4.obs['Group'])
set(adata_li_BA4.obs['CellType'])
set(adata_li_BA4.obs['Region'])
set(adata_li_BA4.obs['Dataset'])
set(adata_li_BA4.obs['Donor_comp'])

###################################
# Li BA9
###################################
adata_li_BA9 = sc.read_h5ad(par_ann_data_Li_BA9)
adata_li_BA9 = adata_li_BA9[adata_li_BA9.obs['CellType'] == par_keep_cell_type]
adata_li_BA9 = adata_li_BA9[adata_li_BA9.obs['Region'] == 'medial frontal cortex']
adata_li_BA9.obs_names = [f"Cell_{i:d}" for i in range(adata_li_BA9.n_obs)]
sc.pp.filter_cells(adata_li_BA9, min_genes=200) ## <------ moved here from line 272

set(adata_li_BA9.obs['Group'])

# Map disease status
mapping = {'C9-ALS': 2, 'C9-FTD': 4,  'Control': 0}
adata_li_BA9.obs['Group'] = adata_li_BA9.obs['Group'].map(mapping)

## Add dataset column
adata_li_BA9.obs['Dataset'] = "Li"

## Add proper region column
adata_li_BA9.obs['Region'] = "BA9"

## Add comprehensive donor column
adata_li_BA9.obs['Donor_comp'] = adata_li_BA9.obs["Donor"].astype(str) + "_" + adata_li_BA9.obs["Group"].astype(str)

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
adata_limone_BA4.obs['Donor_comp'] = adata_limone_BA4.obs["Donor"].astype(str) + "_" + adata_limone_BA4.obs["Group"].astype(str)

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

keep_genes = common_genes.intersection(features)
len(keep_genes) ## USE THIS ONE
len(features)

###################################
# Create a unified object
###################################
keep_genes_1 = [g for g in adata_pineda_BA4.var_names if g in keep_genes]
adata_pineda_BA4 = adata_pineda_BA4[:, keep_genes_1]
len(adata_pineda_BA4.var_names)

keep_genes_1 = [g for g in adata_pineda_BA9.var_names if g in keep_genes]
adata_pineda_BA9 = adata_pineda_BA9[:, keep_genes_1]
len(adata_pineda_BA9.var_names)

keep_genes_1 = [g for g in adata_li_BA4.var_names if g in keep_genes]
adata_li_BA4 = adata_li_BA4[:, keep_genes_1]
len(adata_li_BA4.var_names)

keep_genes_1 = [g for g in adata_li_BA9.var_names if g in keep_genes]
adata_li_BA9 = adata_li_BA9[:, keep_genes_1]
len(adata_li_BA9.var_names)

keep_genes_1 = [g for g in adata_limone_BA4.var_names if g in keep_genes]
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

target_n_cells = 5000
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
# Print UMAP
################################### 
sc.pp.pca(adata_combined, n_comps=50)
sc.pp.neighbors(adata_combined)       
sc.tl.umap(adata_combined)  

sc.pl.umap(adata_combined, color="Group", show=False)
plt.savefig("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", bbox_inches="tight")

sc.pl.umap(adata_combined, color="Dataset", show=False)
plt.savefig("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", bbox_inches="tight")

sc.pl.umap(adata_combined, color="Region", show=False)
plt.savefig("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", bbox_inches="tight")

set(adata_combined.obs['Donor_comp'])


###################################
# Subset back out Pineda from the combine object
################################### 
#adata_combined_pineda = adata_combined[adata_combined.obs['Dataset'] == "Pineda"]

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


# Gradient Reversal Layer
class GradReverse(Function):
    @staticmethod
    def forward(ctx, x, lambd):
        ctx.lambd = lambd
        return x.view_as(x)
    
    @staticmethod
    def backward(ctx, grad_output):
        return grad_output.neg() * ctx.lambd, None

def grad_reverse(x, lambd=1.0):
    return GradReverse.apply(x, lambd)

def train_main_LOSO(main_model, domain_model, dataloader, epochs, device,
                    learning_rate=1e-3, lambda_domain=0.1):
    
    main_model.to(device)
    domain_model.to(device)
    
    optimizer_main = optim.Adam(main_model.parameters(), lr=learning_rate)
    optimizer_domain = optim.Adam(domain_model.parameters(), lr=learning_rate)
    
    loss_class = nn.CrossEntropyLoss()
    loss_domain = nn.CrossEntropyLoss()
    
    for epoch in range(epochs):
        main_model.train()
        domain_model.train()
        
        total_loss_c, total_loss_d = 0.0, 0.0
        correct, total = 0, 0
        num_batches = 0
        
        for X, y, d in dataloader:
            X, y, d = X.to(device), y.to(device), d.to(device)
            
            # Forward pass through main model
            y_out, shared = main_model(X)
            loss_c = loss_class(y_out, y)
            
            # Forward pass through domain model with GRL
            shared_rev = grad_reverse(shared, lambd=lambda_domain)
            d_out = domain_model(shared_rev)
            loss_d = loss_domain(d_out, d)
            
            # Total loss = classification + domain adversarial
            total_loss = loss_c + loss_d
            
            optimizer_main.zero_grad()
            optimizer_domain.zero_grad()
            total_loss.backward()
            optimizer_main.step()
            optimizer_domain.step()
            
            # Accumulate metrics
            total_loss_c += loss_c.item()
            total_loss_d += loss_d.item()
            num_batches += 1
            
            _, predicted = torch.max(y_out, 1)
            correct += (predicted == y).sum().item()
            total += y.size(0)
        
        avg_loss_c = total_loss_c / num_batches
        avg_loss_d = total_loss_d / num_batches
        accuracy = 100.0 * correct / total
        
        print(f"Epoch {epoch}, Avg Class Loss: {avg_loss_c:.4f}, "
              f"Avg Domain Loss: {avg_loss_d:.4f}, Accuracy: {accuracy:.2f}%")


def combat_foldwise_rpy2_adata(adata, holdout_donor, batch_key="Dataset", covariate_key="Group", donor_key="Donor_comp"):
    """
    Perform fold-wise ComBat correction using rpy2 + sva::ComBat.
    Fits ComBat on training donors only, then applies correction to both train and test donors.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object with .X (cells × genes) and .obs containing batch, donor, and covariates
    holdout_donor : str
        Donor ID to hold out for testing
    batch_key : str
        Column in adata.obs specifying batch (e.g., dataset/study)
    covariate_key : str
        Column in adata.obs specifying biological covariates to preserve (e.g., Group/disease status)
    donor_key : str
        Column in adata.obs specifying donor identity
    
    Returns
    -------
    combat_train : pd.DataFrame
        Batch-corrected expression matrix for training donors (genes × cells)
    combat_test : pd.DataFrame
        Batch-corrected expression matrix for test donor (genes × cells)
    """
    
    # --------------------------------------
    # 1. Split into training vs test indices
    # --------------------------------------
    train_idx = adata.obs[adata.obs[donor_key] != holdout_donor].index
    test_idx  = adata.obs[adata.obs[donor_key] == holdout_donor].index
    
    # Extract raw expression (genes × cells)
    expr = pd.DataFrame(
        adata.X.T if isinstance(adata.X, np.ndarray) else adata.X.toarray().T,
        index=adata.var_names,        # genes
        columns=adata.obs_names       # cells
    )
    
    expr_train = expr.loc[:, train_idx]
    expr_test  = expr.loc[:, test_idx]
    
    # --------------------------------------
    # 2. Build covariate matrices
    # --------------------------------------
    covars_train = pd.get_dummies(adata.obs.loc[train_idx, [covariate_key]], drop_first=True)
    covars_all   = pd.get_dummies(adata.obs.loc[train_idx.append(test_idx), [covariate_key]], drop_first=True)
    
    # --------------------------------------
    # 3. Convert training set to R and fit ComBat
    # --------------------------------------
    with localconverter(ro.default_converter + pandas2ri.converter):
        expr_train_r = ro.conversion.py2rpy(expr_train)
        batch_train_r = ro.conversion.py2rpy(adata.obs.loc[train_idx, batch_key])
        mod_train_r = ro.conversion.py2rpy(covars_train)
    
        combat_train_r = sva.ComBat(
            dat=expr_train_r,
            batch=batch_train_r,
            mod=mod_train_r,
            par_prior=True,
            prior_plots=False
        )
        combat_train = ro.conversion.rpy2py(combat_train_r)
    
    # Wrap back into DataFrame
    combat_train = pd.DataFrame(
        combat_train,
        index=expr_train.index,
        columns=expr_train.columns
    )
    
    # --------------------------------------
    # 4. Apply to full set (train+test) to extract test corrections
    # --------------------------------------
    expr_all = expr.loc[:, train_idx.append(test_idx)]
    
    with localconverter(ro.default_converter + pandas2ri.converter):
        expr_all_r = ro.conversion.py2rpy(expr_all)
        batch_all_r = ro.conversion.py2rpy(adata.obs.loc[train_idx.append(test_idx), batch_key])
        mod_all_r   = ro.conversion.py2rpy(covars_all)
    
        combat_all_r = sva.ComBat(
            dat=expr_all_r,
            batch=batch_all_r,
            mod=mod_all_r,
            par_prior=True,
            prior_plots=False
        )
        combat_all = ro.conversion.rpy2py(combat_all_r)
    
    # Wrap into DataFrame
    combat_all = pd.DataFrame(
        combat_all,
        index=expr_all.index,
        columns=expr_all.columns
    )
    
    # Extract corrected test cells
    combat_test = combat_all.loc[:, test_idx]
    
    return combat_train, combat_test



###################################
# Perform KNN LOSO with Pineda for model validation
###################################
num_genes = adata_combined.X.shape[1]
data = []

# LOSO loop
sample_IDs = set(adata_combined.obs['Donor_comp'])
#kNN_threshold_list = [.99, .95, .9, .85, .8]
kNN_threshold = 0.9
print(f"KNN Threshold: {kNN_threshold}")
training_epoch = 10
    

#donor = '235_BA4_2'
#train_corr, test_corr = combat_foldwise_rpy2_adata(adata_combined, donor, batch_key="Dataset", covariate_key="Group", donor_key="Donor_comp")

for donor in sample_IDs:
    
    train_corr, test_corr = combat_foldwise_rpy2_adata(adata_combined, donor, batch_key="Dataset", covariate_key="Group", donor_key="Donor_comp")
    
    for _ in range(5):
        print(f"Processing donor: {donor}")
        adata_train = adata_combined[adata_combined.obs['Donor_comp'] != donor]
        adata_test = adata_combined[adata_combined.obs['Donor_comp'] == donor]
        num_cells = adata_test.shape[0]
     
        #X_train = torch.FloatTensor(train_corr.toarray() if hasattr(train_corr, 'toarray') else train_corr)
        X_train = torch.FloatTensor(train_corr.values.T)
        y_train = torch.LongTensor(adata_train.obs['Group'].values)
        domains = pd.factorize(adata_train.obs['Dataset'])[0] ## Modified
        d_train = torch.LongTensor(domains)
     
        dataset = TensorDataset(X_train, y_train, d_train)
        input_size = adata_train.shape[1]
        num_classes = len(np.unique(y_train))
        num_domains = len(np.unique(domains))
     
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
     
        model_main = MainModel(input_size, num_classes)
        model_domain = DomainClassifier(25, num_domains)
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
     
        #train_main_LOSO(model_main, model_domain, loader, training_epoch, device) ## USED THIS ONE ORIGINALLY
        train_main_LOSO(model_main, model_domain, loader, epochs = 10, device=device, learning_rate=1e-3, lambda_domain=0.1)
        
        # Evaluation
        X_test  = torch.FloatTensor(test_corr.values.T)
        y_test = adata_test.obs['Group']
     
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
           
            print(y_all_pred)
            unique, counts = np.unique(y_all_pred, return_counts=True)
            count_dict = dict(zip(unique, counts))
            count_str = "; ".join([f'"{label}" = {count}' for label, count in count_dict.items()])
            print(count_str)
           
            data.append({
                'prep': par_prep, 'donor': donor,
                'group': 'ALS', 'celltype': par_keep_cell_type,
                'n_genes': num_genes, 'n_cells': len(y_true), 'train_epoch': training_epoch,
                'learning_rate': learning_rate, 'batch_size': batch_size,
                'test_accuracy_high': acc_high, 'test_accuracy_low': acc_low,
                'test_accuracy_all': acc_all, 'counts_pred':  count_str
            })

pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_disease_5000_cell_true_label_{par_keep_cell_type}_3.csv"
results_df.to_csv(out_path, index=False)
EOF

jobid=$(sbatch $jobfile | awk '{print $4}')
echo "Submitted job for ${CT} with JobID: ${jobid}"

done






