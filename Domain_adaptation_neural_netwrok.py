salloc -A def-sfarhan --time=0-5 -c 1 --mem=40g

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
par_keep_cell_type = "olig"

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


# Preprocess the source data
sc.pp.normalize_total(source_dataset)
sc.pp.log1p(source_dataset)
sc.pp.highly_variable_genes(source_dataset, flavor='seurat_v3', n_top_genes = 10000)
source_dataset = source_dataset[:, source_dataset.var['highly_variable']]

# Preprocess the target data
sc.pp.normalize_total(target_dataset)
sc.pp.log1p(target_dataset)
target_dataset = target_dataset[:, source_dataset.var_names]  # Ensure the same genes are used

# Convert to tensors
X_source = torch.tensor(source_dataset.X.A, dtype=torch.float32)
source_dataset.obs['Disease_Status'] = source_dataset.obs['Disease_Status'].astype('category')
y_source = torch.tensor(source_dataset.obs['Disease_Status'].cat.codes.values, dtype=torch.long)
X_target = torch.tensor(target_dataset.X.A, dtype=torch.float32)

# Create DataLoader
source_dataset = TensorDataset(X_source, y_source)
target_dataset = TensorDataset(X_target)

source_loader = DataLoader(source_dataset, batch_size=256, shuffle=True) ## smaller batch size might help with generalizability
target_loader = DataLoader(target_dataset, batch_size=256, shuffle=False) ## smaller batch size might help with generalizability

# Define the Domain Adaptation Model with Batch Normalization
class DomainAdaptationModel(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(DomainAdaptationModel, self).__init__()
        self.feature_extractor = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),  # Added Batch Normalization
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),  # Added Batch Normalization
            nn.ReLU()
        )
        self.classifier = nn.Linear(hidden_dim, output_dim)
        
    def forward(self, x):
        features = self.feature_extractor(x)
        return features, self.classifier(features)  # Returns both features and class predictions

# Define the Adversarial Network with Batch Normalization
class AdversarialNetwork(nn.Module):
    def __init__(self, hidden_dim):
        super(AdversarialNetwork, self).__init__()
        self.discriminator = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),  # Added Batch Normalization
            nn.ReLU(),
            nn.Linear(hidden_dim, 1)
        )
        
    def forward(self, features):
        return self.discriminator(features)

# Initialize models
input_dim = X_source.shape[1]
hidden_dim = 512
output_dim = 2 #len(source.obs['disease_status'].cat.categories)

model = DomainAdaptationModel(input_dim, hidden_dim, output_dim)
adversarial_model = AdversarialNetwork(hidden_dim)

# Optimizers and Loss functions
optimizer = torch.optim.Adam(list(model.parameters()) + list(adversarial_model.parameters()), lr=0.001)
criterion_class = nn.CrossEntropyLoss()
criterion_adversarial = nn.BCEWithLogitsLoss()  # Binary Cross-Entropy for adversarial loss

# Training function with adversarial training
def train(model, adversarial_model, source_loader, target_loader, optimizer, criterion_class, criterion_adversarial, epochs=100):
    for epoch in range(epochs):
        model.train()
        for (source_data, source_labels) in source_loader:
            optimizer.zero_grad()
            # Forward pass through the main model
            source_features, source_outputs = model(source_data)
            loss_class = criterion_class(source_outputs, source_labels)
            # Adversarial training
            target_data = next(iter(target_loader))  # Sample target dataset, Only need one variable since there are no labels
            if isinstance(target_data, (list, tuple)):
                target_data = target_data[0]  # Unpack if it's a list or tuple
            target_features, _ = model(target_data)
            # Combine source and target features for adversarial loss
            all_features = torch.cat([source_features, target_features])
            domain_labels = torch.cat([torch.ones(source_features.size(0)), torch.zeros(target_features.size(0))])
            # Compute adversarial loss
            domain_predictions = adversarial_model(all_features)
            loss_adversarial = criterion_adversarial(domain_predictions.view(-1), domain_labels.float())
            # Total loss: Classification loss + Adversarial loss
            total_loss = loss_class + loss_adversarial
            total_loss.backward()
            optimizer.step()
        print(f'Epoch [{epoch + 1}/{epochs}], Loss: {total_loss.item():.4f}')

# Train the model with adversarial training
train(model, adversarial_model, source_loader, target_loader, optimizer, criterion_class, criterion_adversarial)


## Evaluate with balanced accuracy?
def evaluate(model, target_loader, true_labels):
    model.eval()
    predictions = []
    
    with torch.no_grad():
        for target_data in target_loader:
            # Use only the features for evaluation
            features = target_data[0]
            _, predicted = model(features)
            _, predicted_classes = torch.max(predicted, 1)
            predictions.extend(predicted_classes.numpy())
    
    # Calculate balanced accuracy
    balanced_acc = balanced_accuracy_score(true_labels, predictions)
    return balanced_acc

# Assuming you have the true labels for the target dataset
target_dataset = adata[adata.obs['Sample_ID'].isin(target_ids)].copy()
target_dataset.obs['Disease_Status'] = target_dataset.obs['Disease_Status'].astype('category')
true_labels = torch.tensor(target_dataset.obs['Disease_Status'].cat.codes.values, dtype=torch.long)

balanced_accuracy = evaluate(model, target_loader, true_labels)
print(f'Balanced Accuracy: {balanced_accuracy:.4f}')


