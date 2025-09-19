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

nano Cell_type_subtype_BA9.sh

#!/bin/bash  
#SBATCH --account=def-sfarhan
#SBATCH --time=00-05:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=100g          # memory per cor
#SBATCH --job-name=Cell_type_subtype_BA9
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

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/Cell_type_subtype_BA9.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano Cell_type_subtype_BA9.py


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

import numpy as np
import pandas as pd
from scipy.spatial.distance import euclidean, jensenshannon
from itertools import combinations
from scipy.stats import entropy, spearmanr, pearsonr

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
    # Cell Specific parameters 
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
    
    adata_combined = adata_pineda_BA9

    ###################################
    # Preprocess combined data
    ################################### 
    sc.pp.filter_genes(adata_combined, min_cells=3)
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    sc.pp.log1p(adata_combined)
    sc.pp.highly_variable_genes(adata_combined, flavor='seurat')
    
    features = adata_combined.var_names[adata_combined.var['highly_variable']].tolist()
    len(features)
    
    ###################################
    # Quantify cell subtype
    ###################################  
    
    # Expression matrix (cells × genes)
    expr = adata_combined.X  # dense or sparse
    obs = adata_combined.obs
    
    def shannon_entropy(counts):
        proportions = counts / counts.sum()
        return entropy(proportions, base=2)
    
    def compute_pseudobulk(expr, obs, groupby):
        """Aggregate expression per group (mean)."""
        df = pd.DataFrame(
            expr.toarray() if not isinstance(expr, np.ndarray) else expr,
            index=obs.index,
            columns=[f"gene{i}" for i in range(expr.shape[1])]
        )
        df[groupby] = obs[groupby].values
        return df.groupby(groupby).mean()
    
    # Count subtypes
    subtype_counts = obs["SubType"].value_counts()
    n_subtypes = len(subtype_counts)
    
    # Shannon entropy of subtype distribution
    H = shannon_entropy(subtype_counts)
    
    # Effective number of subtypes
    effective_subtypes = 2**H
    
    # Subtype divergence (pairwise pseudobulk distances)
    pseudobulk = compute_pseudobulk(expr, obs, "SubType")
    distances = []
    for s1, s2 in combinations(pseudobulk.index, 2):
        d = euclidean(pseudobulk.loc[s1], pseudobulk.loc[s2])
        # or use: d = jensenshannon(pseudobulk.loc[s1], pseudobulk.loc[s2])
        distances.append(d)
    avg_divergence = np.mean(distances) if distances else 0.0
    
    # Entropy-adjusted divergence
    adjusted_divergence = avg_divergence * effective_subtypes
    
    # Append results
    row = {
        "celltype": par_keep_cell_type,
        "NumSubtypes": n_subtypes,
        "ShannonEntropy": H,
        "AvgSubtypeDivergence": avg_divergence,
        "EffectiveSubtypes": effective_subtypes,
        "AdjustedDivergence": adjusted_divergence
    }
    
    data.append(row)

pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/Cell_type_subtype_BA9.csv"
results_df.to_csv(out_path, index=False)


