salloc -A def-sfarhan --time=0-8 -c 1 --mem=40g

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  

## generalizable model
from bayes_opt import BayesianOptimization


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

import numpy as np
from sklearn.model_selection import KFold
from torch.utils.data import DataLoader, TensorDataset


from sklearn.metrics import classification_report, accuracy_score
import numpy as np
import torch
from torch.utils.data import DataLoader
from sklearn.model_selection import KFold

###################################
# scanpy settings
###################################
sc.settings.verbosity = 3             
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


###################################
# Cell Specific parameters -- astro
###################################
#par_ann_data = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B9.h5ad"
par_ann_data = "/home/fiorini9/scratch/machine_learning_ALS/combat_test_astro_int.h5ad"
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



#######################################################
## Perform the test train split
#######################################################

# Create an array of indices
indices = np.arange(adata.shape[0])

# Split indices into train and test
train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=42)

# Create training and testing AnnData objects
adata_train = adata[train_indices].copy()
adata_test = adata[test_indices].copy()

#X_train, X_test, y_train, y_test = train_test_split(
#    adata_train.X, adata_train.obs['Disease_Status'], test_size=0.2, random_state=42
#)

#X_train_dense = X_train.toarray() if hasattr(X_train, 'toarray') else X_train
#X_test_dense = X_test.toarray() if hasattr(X_test, 'toarray') else X_test


#######################################################
## Preprocess the train and test
#######################################################

## Process train data
sc.pp.filter_cells(adata_train, min_genes=200)  # Filter cells
sc.pp.filter_genes(adata_train, min_cells=3)    # Filter genes
sc.pp.normalize_total(adata_train, target_sum=1e4)  # Normalize
sc.pp.log1p(adata_train)                         # Log-transform
#sc.pp.highly_variable_genes(source_dataset, flavor='seurat_v3', n_top_genes = 7000)        # Find HVGs
sc.pp.highly_variable_genes(adata_train) 
adata_train.raw = adata_train  # Save raw counts for further analysis
adata_train = adata_train[:, adata_train.var['highly_variable']]  # Subset to HVGs


## Process test object
#sc.pp.filter_cells(target_dataset, min_genes=200)  # Filter cells
#sc.pp.filter_genes(target_dataset, min_cells=3)    # Filter genes
sc.pp.normalize_total(adata_test, target_sum=1e4)  # Normalize
sc.pp.log1p(adata_test)                         # Log-transform
adata_test.raw = adata_test  # Save raw counts for further analysis
adata_test = adata_test[:, adata_train.var_names]  # Ensure the same genes are used


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
    #main_optimizer = optim.Adam(main_model.parameters(), lr=0.0001, weight_decay=0.01)
    main_optimizer = optim.Adam(main_model.parameters(), lr=0.0001)
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


###################################
## 5 fold cross validation on the test set
###################################

# Number of folds
n_splits = 5
kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

########## Prepare your data as a TensorDataset ##########
## X
X_train = adata_train.X
X_train_dense = X_train.toarray() if hasattr(X_train, 'toarray') else X_train
X_train_tensor = torch.FloatTensor(X_train_dense)  # Gene expression matrix

## Y
y_train = adata_train.obs['Disease_Status']
y_train_tensor = torch.LongTensor(y_train)

## domain
unique_samples = adata_train.obs['Sample_ID'].unique()
num_domains = len(unique_samples)
# Create a mapping from sample IDs to domain labels
sample_to_domain = {sample: idx for idx, sample in enumerate(unique_samples)}
# Create domain labels based on Sample_ID
adata_train.obs['Domain_Label'] = adata_train.obs['Sample_ID'].map(sample_to_domain)
y_train_indices = y_train.index
adata_train_subset = adata_train[y_train_indices]
domain_labels_tensor = torch.LongTensor(adata_train_subset.obs['Domain_Label'])

## full dataset tensor
X_train_tensor.shape
y_train_tensor.shape
domain_labels_tensor.shape

full_dataset = TensorDataset(X_train_tensor, y_train_tensor, domain_labels_tensor)

num_classes = len(np.unique(adata_train.obs['Disease_Status']))  # Number of unique disease classes
num_domains = len(np.unique(domain_labels_tensor.numpy()))  # Number of unique sample IDs or batches
input_size = adata_train.shape[1]  # Number of highly variable genes


accuracies = []  # List to store accuracy for each fold

for fold, (train_idx, val_idx) in enumerate(kf.split(np.arange(len(full_dataset)))):
    print(f"Fold {fold + 1}/{n_splits}")
    
    # Create data loaders for the current fold
    train_subset = DataLoader(full_dataset, batch_size=32, sampler=torch.utils.data.SubsetRandomSampler(train_idx))
    val_subset = DataLoader(full_dataset, batch_size=32, sampler=torch.utils.data.SubsetRandomSampler(val_idx))
    
    # Initialize models for each fold
    main_model = MainModel(input_size, num_classes)
    domain_model = DomainClassifier(25, num_domains)  # 50 is the output size from shared layers
    
    # Train the models for the current fold
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    train(main_model, domain_model, train_subset, num_epochs=25, device=device)
    
    # Evaluation on the validation set
    main_model.eval()
    val_predictions = []
    val_labels = []
    
    with torch.no_grad():
        for X_val, y_val, domain_labels in val_subset:
            X_val, y_val, domain_labels = X_val.to(device), y_val.to(device), domain_labels.to(device)
            y_pred_val, _ = main_model(X_val)
            val_predictions.append(torch.argmax(y_pred_val, dim=1).cpu())
            val_labels.append(y_val.cpu())
    
    # Concatenate predictions and labels
    val_predictions = torch.cat(val_predictions).numpy()
    val_labels = torch.cat(val_labels).numpy()
    
    # Calculate accuracy for the current fold
    accuracy = accuracy_score(val_labels, val_predictions)
    accuracies.append(accuracy)
    print(f"Accuracy for Fold {fold + 1}: {accuracy:.4f}")

# Calculate and print the mean accuracy across all folds
mean_accuracy = np.mean(accuracies)
print(f"\nMean Accuracy across all folds: {mean_accuracy:.4f}")

###################################
## Train on the entire train set
###################################
# Sample usage
# Define input size and number of classes
input_size = adata_train.shape[1]  # Number of highly variable genes
num_classes = len(np.unique(adata_train.obs['Disease_Status']))  # Number of unique disease classes
num_domains = len(np.unique(domain_labels_tensor.numpy()))  # Number of unique sample IDs or batches

# Initialize models
main_model = MainModel(input_size, num_classes)
domain_model = DomainClassifier(25, num_domains)  # 50 is the output size from shared layers

# Create dataset and DataLoader
train_dataset = TensorDataset(X_train_tensor, y_train_tensor, domain_labels_tensor)
train_dataloader = DataLoader(train_dataset, batch_size=32, shuffle=True)

# Train the models
device = 'cuda' if torch.cuda.is_available() else 'cpu'
train(main_model, domain_model, train_dataloader, num_epochs=25, device=device)

###################################
## Evaluate on test set
###################################
X_test = adata_test.X
X_test_dense = X_test.toarray() if hasattr(X_test, 'toarray') else X_test
X_test_tensor = torch.FloatTensor(X_test_dense)

y_test = adata_test.obs['Disease_Status']

main_model.eval()
with torch.no_grad():
    y_pred_test, _ = main_model(X_test_tensor)
    y_pred_test_labels = torch.argmax(y_pred_test, dim=1)
    print("Test Classification Report:\n", classification_report(y_test, y_pred_test_labels.numpy()))


