salloc -A def-sfarhan --time=0-8 -c 1 --mem=200g

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  

## generalizable model
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
import numpy as np
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix, vstack
from numpy import inf
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import NMF
import anndata as ad
from scETM import scETM, UnsupervisedTrainer, evaluate, prepare_for_transfer
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import TensorDataset
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
from scipy.sparse import csr_matrix

from sklearn.metrics import balanced_accuracy_score

import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report

###################################
# scanpy settings
###################################
sc.settings.verbosity = 3             
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


###################################
# Cell Specific parameters -- astro
###################################
par_ann_data = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B9.h5ad"
par_keep_cell_type = "astro"

###################################
# Load anndata
###################################
adata = sc.read_h5ad(par_ann_data)
adata.obs
adata = adata[adata.obs['Cell_Type'] == par_keep_cell_type]
adata.obs['Cell_Type']
print("genes: ", adata.var_names) 
print("cells: ", adata.obs_names) 
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]

###################################
# create new_annData object without LBD cell samples
###################################
adata.obs = adata.obs.reset_index() 
adata.obs.columns
set(adata.obs['Sample_ID'])

#adata.obs['Disease_Status']
#adata.obs
#adata.obs.shape
#adata.obs = adata.obs.reset_index() 

set(adata.obs['Disease_Status'])
set(adata.obs['Cell_Type'])


###################################
# Split dataset to keep certain subjects as source domain and others as target domain
###################################
# Get unique Sample_IDs for ALS and PN
als_ids = adata.obs[adata.obs['Sample_ID'].str.contains("ALS")]['Sample_ID'].unique()
pn_ids = adata.obs[adata.obs['Sample_ID'].str.contains("PN")]['Sample_ID'].unique()

# Function to split Sample_IDs into source and target with an 80-20 split
def split_sample_ids(ids, source_ratio=0.8):
    np.random.seed(42)  # For reproducibility
    np.random.shuffle(ids)  # Shuffle the IDs
    split_index = int(len(ids) * source_ratio)  # Calculate split index based on the source ratio
    return ids[:split_index], ids[split_index:]  # Source and target

# Split ALS and PN Sample_IDs
source_als, target_als = split_sample_ids(als_ids)
source_pn, target_pn = split_sample_ids(pn_ids)

# Combine source and target IDs
source_ids = np.concatenate([source_als, source_pn])
target_ids = np.concatenate([target_als, target_pn])

# Display the results
print(f"Source Sample_IDs: {source_ids}")
print(f"Target Sample_IDs: {target_ids}")

# Optionally, you can create new AnnData objects for the source and target datasets
source_dataset = adata[adata.obs['Sample_ID'].isin(source_ids)].copy()
target_dataset = adata[adata.obs['Sample_ID'].isin(target_ids)].copy()

# Check the shapes of the new datasets
print(f"Source dataset shape: {source_dataset.shape}")
print(f"Target dataset shape: {target_dataset.shape}")


# Preprocessing: Quality Control, Normalization, and Batch Correction
def preprocess_data(adata):
    sc.pp.filter_cells(adata, min_genes=200)  # Filter cells
    sc.pp.filter_genes(adata, min_cells=3)    # Filter genes
    sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize
    sc.pp.log1p(adata)                         # Log-transform
    sc.pp.highly_variable_genes(adata)        # Find HVGs
    adata.raw = adata  # Save raw counts for further analysis
    adata = adata[:, adata.var['highly_variable']]  # Subset to HVGs
    return adata

## Process train data
sc.pp.filter_cells(source_dataset, min_genes=200)  # Filter cells
sc.pp.filter_genes(source_dataset, min_cells=3)    # Filter genes
sc.pp.normalize_total(source_dataset, target_sum=1e4)  # Normalize
sc.pp.log1p(source_dataset)                         # Log-transform
sc.pp.highly_variable_genes(source_dataset)        # Find HVGs
source_dataset.raw = source_dataset  # Save raw counts for further analysis
adata_train = source_dataset[:, source_dataset.var['highly_variable']]  # Subset to HVGs


## Process test object
#sc.pp.filter_cells(target_dataset, min_genes=200)  # Filter cells
#sc.pp.filter_genes(target_dataset, min_cells=3)    # Filter genes
sc.pp.normalize_total(target_dataset, target_sum=1e4)  # Normalize
sc.pp.log1p(target_dataset)                         # Log-transform
target_dataset.raw = target_dataset  # Save raw counts for further analysis
adata_test = target_dataset[:, adata_train.var_names]  # Ensure the same genes are used


# Split training data into training and validation sets
X_train, X_val, y_train, y_val = train_test_split(
    adata_train.X, adata_train.obs['Disease_Status'], test_size=0.2, random_state=42
)

# Convert sparse matrices to dense format
X_train_dense = X_train.toarray() if hasattr(X_train, 'toarray') else X_train
X_val_dense = X_val.toarray() if hasattr(X_val, 'toarray') else X_val

# Scale features using the dense format
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train_dense)
X_val_scaled = scaler.transform(X_val_dense)

# Train MLP Classifier
mlp = MLPClassifier(hidden_layer_sizes=(100,50,3),max_iter=300, random_state=42)
mlp.fit(X_train_scaled, y_train)

# Validate on the validation set
y_pred_val = mlp.predict(X_val_scaled)
print("Validation Classification Report:\n", classification_report(y_val, y_pred_val))

# Transfer Learning: Fine-tune on a small subset of the test dataset if you have labels
# If you don't have labels, skip to classifying the test dataset.
#if 'Disease_Status' in adata_test.obs.columns:

# Scale the test data
X_test = adata_test.X
y_test = adata_test.obs['Disease_Status']

# Scale the test data
X_test_dense = X_test.toarray() if hasattr(X_test, 'toarray') else X_test
X_test_scaled = scaler.transform(X_test_dense)

# Randomly sample 10% of the test data for fine-tuning
n_fine_tune_samples = int(0.2 * len(y_test))  # 10% of the test set size
fine_tune_indices = np.random.choice(len(y_test), size=n_fine_tune_samples, replace=False)

# Create subsets for fine-tuning
X_fine_tune = X_test_scaled[fine_tune_indices]
y_fine_tune = y_test[fine_tune_indices]

# Create a mask to exclude fine-tuning samples for evaluation
mask = np.ones(len(y_test), dtype=bool)
mask[fine_tune_indices] = False

# Use the mask to create the evaluation set
X_eval = X_test_scaled[mask]
y_eval = y_test[mask]

# Fine-tune the model using the sampled data
mlp.fit(X_fine_tune, y_fine_tune)

# Evaluate on the remaining test dataset
y_pred_test = mlp.predict(X_eval)
print("Test Classification Report:\n", classification_report(y_eval, y_pred_test))



# Evaluate on the remaining test dataset
y_pred_test = mlp.predict(X_test_scaled)
print("Test Classification Report:\n", classification_report(y_test, y_pred_test))

