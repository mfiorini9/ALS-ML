## run this in Narval

salloc -A def-grouleau --time=0-4 -c 1 --mem=50g

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L3_L5 

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

nano optimal_LOO_0delta_absScaling_weighted_LIME_SALS_BA4_all_cells_narval_2.sh

#!/bin/bash  
#SBATCH --account=def-grouleau
#SBATCH --time=00-12:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=150g          # memory per cor
#SBATCH --job-name=optimal_LOO_0delta_absScaling_weighted_LIME_SALS_BA4_all_cells_narval_2
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/optimal_LOO_0delta_absScaling_weighted_LIME_SALS_BA4_all_cells_narval_2.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano optimal_LOO_0delta_absScaling_weighted_LIME_SALS_BA4_all_cells_narval_2.py



################################################################## clean version ------ USE THIS ONE

import pandas as pd
import numpy as np
import scanpy as sc
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
import seaborn as sns
import matplotlib.pyplot as plt

################################################################################################ 0

# Set parameters
#par_keep_cell_type = "L2_L3"
par_prep = "CombatSeq"
par_brain_region = "BA4"
par_status = "SALS"
remove = ["SFTLD", "C9FTLD", "C9ALS"]

#cell_types = [
#    "L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
#    "Rosehip", "SOM", "Oligo", "Astro", "OPC", "Micro",
#    "Mural", "Endo", "Fibro", "L5"
#]

# Loop through the list and print each item
#for par_keep_cell_type in cell_types:
#    print(par_keep_cell_type)


par_keep_cell_type = "L6"
    
# Load LIME-selected genes
lime_file = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_0_abs_case_control_narval_2.csv'
features = pd.read_csv(lime_file)["gene"].tolist()

# Load model report
report_file = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
report_df = pd.read_csv(report_file)
batch_size = int(report_df.at[0, 'batch_size'])
learning_rate = report_df.at[0, 'learning_rate']

# Load anndata
adata_file = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
adata = sc.read_h5ad(adata_file)

# Filter data
adata = adata[~adata.obs['Group'].isin(remove)]
adata = adata[adata.obs['CellType'] == par_keep_cell_type]
adata = adata[adata.obs['Region'] == par_brain_region]
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
sample_IDs = set(adata.obs['orig.ident'])

list(adata.obs.columns)

## Plt some QC variables
plt.figure(figsize=(12, 6))
sns.violinplot(x=adata.obs['orig.ident'], y=adata.obs['percent.ribo'])

plt.xticks(rotation=90)
plt.xlabel("orig.ident")
plt.ylabel("nFeature_RNA")
plt.title("Violin plot of nFeature_RNA across orig.ident")

# Save to PDF
plt.tight_layout()
plt.savefig("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf")  # <- This line saves the figure
plt.show()
