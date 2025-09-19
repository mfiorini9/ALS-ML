salloc -A def-sfarhan --time=0-8 -c 1 --mem=40g

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

## change Disease status labels
# Create a mapping dictionary
mapping = {'ALS': 1, 'ctrl': 0}

# Map the values
adata.obs['Disease_Status'] = adata.obs['Disease_Status'].map(mapping)

# Verify the change
print(adata.obs['Disease_Status'].value_counts())


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






















#######################################################
## Start here
#######################################################

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

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report

import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder
import torch

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

# Assuming the previous code has already set up the source_dataset and target_dataset
# Prepare data for the DataLoader

###################################
## Split the train set into train and validation
###################################
X_train, X_val, y_train, y_val = train_test_split(
    adata_train.X, adata_train.obs['Disease_Status'], test_size=0.2, random_state=42
)

X_train_dense = X_train.toarray() if hasattr(X_train, 'toarray') else X_train
X_val_dense = X_val.toarray() if hasattr(X_val, 'toarray') else X_val

###################################
## prepare tensors
###################################
###### Gene expression
X_train_tensor = torch.FloatTensor(X_train_dense)  # Gene expression matrix
X_val_tensor = torch.FloatTensor(X_val_dense)  # Gene expression matrix

###### Class labels
y_train_tensor = torch.LongTensor(y_train)

###### Sample ID
unique_samples = adata_train.obs['Sample_ID'].unique()
num_domains = len(unique_samples)

# Create a mapping from sample IDs to domain labels
sample_to_domain = {sample: idx for idx, sample in enumerate(unique_samples)}

# Create domain labels based on Sample_ID
adata_train.obs['Domain_Label'] = adata_train.obs['Sample_ID'].map(sample_to_domain)
y_train_indices = y_train.index
adata_train_subset = adata_train[y_train_indices]

domain_labels_tensor = torch.LongTensor(adata_train_subset.obs['Domain_Label'])


###################################
## Create data loaders
###################################

# Create dataset and DataLoader
train_dataset = TensorDataset(X_train_tensor, y_train_tensor, domain_labels_tensor)
train_dataloader = DataLoader(train_dataset, batch_size=32, shuffle=True)

###################################
## Define the models
###################################
# Define the main classification model
class MainModel(nn.Module):
    def __init__(self, input_size, num_classes):
        super(MainModel, self).__init__()
        self.shared = nn.Sequential(
            nn.Linear(input_size, 100), 
            nn.ReLU(),
            nn.Linear(100, 50),          # Second hidden layer
            nn.ReLU(),
            nn.Linear(50, 25),           # New third hidden layer
            nn.ReLU() 
        )
        self.classifier = nn.Linear(25, num_classes)
    
    def forward(self, x):
        shared_rep = self.shared(x.clone())  # Clone input to prevent modifications
        class_output = self.classifier(shared_rep)
        return class_output, shared_rep

# Define the domain classifier
class DomainClassifier(nn.Module):
    def __init__(self, input_size, num_domains):
        super(DomainClassifier, self).__init__()
        self.domain_classifier = nn.Sequential(
            nn.Linear(input_size, 25),
            nn.ReLU(),
            nn.Linear(25, num_domains)
        )
 
    def forward(self, x):
        return self.domain_classifier(x)

# Function to train the models
torch.autograd.set_detect_anomaly(True)

def train(main_model, domain_model, dataloader, num_epochs, device):
    main_model.to(device)
    domain_model.to(device)
 
    # Optimizers
    main_optimizer = optim.Adam(main_model.parameters(), lr=0.00001, weight_decay=0.1)
    domain_optimizer = optim.Adam(domain_model.parameters(), lr=0.001)
 
    # Loss functions
    class_loss_fn = nn.CrossEntropyLoss()
    domain_loss_fn = nn.CrossEntropyLoss()
 
    for epoch in range(num_epochs):
        main_model.train()
        domain_model.train()
 
        for X, y, domain_labels in dataloader:
            X, y, domain_labels = X.to(device), y.to(device), domain_labels.to(device)
 
            # Forward pass through the main model
            class_output, shared_rep = main_model(X)
            class_loss = class_loss_fn(class_output, y)
 
            # Use detach to avoid modifying shared_rep for the domain model
            domain_output = domain_model(shared_rep.detach())  
            domain_loss = domain_loss_fn(domain_output, domain_labels)
 
            # Combine losses with a negative sign for domain loss
            total_loss = class_loss - domain_loss
             
            # Update the main model
            main_optimizer.zero_grad()
            total_loss.backward(retain_graph=True)  # Do not retain graph here unless necessary
            main_optimizer.step()
 
            # Update the domain model
            domain_optimizer.zero_grad()
            domain_loss.backward()  # No need to retain graph here
            domain_optimizer.step()
 
        print(f"Epoch [{epoch+1}/{num_epochs}], Class Loss: {class_loss.item():.4f}, Domain Loss: {domain_loss.item():.4f}")

# Sample usage
# Define input size and number of classes
input_size = adata_train.shape[1]  # Number of highly variable genes
num_classes = len(np.unique(adata_train.obs['Disease_Status']))  # Number of unique disease classes
num_domains = len(np.unique(domain_labels_tensor.numpy()))  # Number of unique sample IDs or batches

# Initialize models
main_model = MainModel(input_size, num_classes)
domain_model = DomainClassifier(25, num_domains)  # 50 is the output size from shared layers

# Train the models
device = 'cuda' if torch.cuda.is_available() else 'cpu'
train(main_model, domain_model, train_dataloader, num_epochs=15, device=device)

#######################################
### Evaluation on the validation set
#######################################
# Prepare test data
X_val_tensor = torch.FloatTensor(X_val_dense)  # Gene expression matrix

main_model.eval()
with torch.no_grad():
    y_pred_val, _ = main_model(X_val_tensor)
    y_pred_val_labels = torch.argmax(y_pred_val, dim=1)
    print("Validation Classification Report:\n", classification_report(y_val, y_pred_val_labels.numpy()))

#######################################
### Evaluation on the test set
#######################################

X_test = adata_test.X
X_test_dense = X_test.toarray() if hasattr(X_test, 'toarray') else X_test
X_test_tensor = torch.FloatTensor(X_test_dense)

y_test = adata_test.obs['Disease_Status']

main_model.eval()
with torch.no_grad():
    y_pred_test, _ = main_model(X_test_tensor)
    y_pred_test_labels = torch.argmax(y_pred_test, dim=1)
    print("Test Classification Report:\n", classification_report(y_test, y_pred_test_labels.numpy()))

