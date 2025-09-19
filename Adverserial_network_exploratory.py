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
par_test_anndata = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/test_pineda_astro_ALS_BA4.h5ad'
par_train_anndata = '/home/fiorini9/scratch/machine_learning/machine_learning_outs/train_pineda_astro_ALS_BA4.h5ad'

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

# Function to split Sample_IDs into source and target
def split_sample_ids(ids):
    np.random.seed(42)  # For reproducibility
    np.random.shuffle(ids)  # Shuffle the IDs
    split_index = len(ids) // 2  # Split in the middle
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

# EXPLORATORUY CODE WITH ANNDATA OBJECT ####################################################################################################################################################################################################################################
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset

# Assuming source_dataset and target_dataset are your source and target AnnData objects

# Extract features and labels for source and target datasets
# Replace 'disease_status' with the actual column name in your adata.obs
source_features = source_dataset.X  # Assuming it's in sparse format
source_labels = source_dataset.obs['Disease_Status'].astype('category').cat.codes

target_features = target_dataset.X  # Assuming it's in sparse format
target_labels = target_dataset.obs['Disease_Status'].astype('category').cat.codes

#class SyntheticDataset(Dataset):
#    def __init__(self, features, labels):
#        self.features = features
#        self.labels = labels
#
#    def __len__(self):
#        return len(self.labels)
#
#    def __getitem__(self, idx):
#        return self.features[idx].toarray().astype(np.float32), self.labels[idx]

class SyntheticDataset(Dataset):
    def __init__(self, features, labels):
        self.features = features
        self.labels = labels
            
    def __len__(self):
        return self.labels
    
    def __getitem__(self, idx):
        return self.features[idx].toarray().astype(np.float32), self.labels[idx].astype(np.float32)

#class SyntheticDataset(Dataset):
#    def __init__(self, features, labels):
#        self.features = torch.from_numpy(features.toarray()).float()  # Convert to tensor
#        self.labels = torch.tensor(labels, dtype=torch.float32)

#    def __len__(self):
#        return len(self.labels)

#    def __getitem__(self, idx):
#        return self.features[idx], self.labels[idx]


import torch
from torch.utils.data import Dataset
import scipy.sparse as sp

import torch
from torch.utils.data import Dataset
import scipy.sparse as sp

class SyntheticDataset(Dataset):
    def __init__(self, features, labels):
        # Convert sparse matrix to PyTorch sparse tensor
        if isinstance(features, sp.csr_matrix):
            # Convert CSR matrix to a PyTorch sparse tensor
            indices = torch.tensor(
                [features.nonzero()[0], features.nonzero()[1]], dtype=torch.long
            )
            values = torch.tensor(features.data, dtype=torch.float32)
            self.features = torch.sparse.FloatTensor(indices, values, torch.Size(features.shape))
        else:
            self.features = torch.FloatTensor(features)  # Convert to dense tensor if needed
        self.labels = torch.tensor(labels, dtype=torch.float32)
    def __len__(self):
        return len(self.labels)
    def __getitem__(self, idx):
        return self.features[idx], self.labels[idx]

# Assuming source_features and source_labels are defined appropriately

source_dataset = SyntheticDataset(source_features, source_labels)
target_dataset = SyntheticDataset(target_features, target_labels)

source_loader = DataLoader(source_dataset, batch_size=100, shuffle=True) ## smaller batch size might help with generalizability
target_loader = DataLoader(target_dataset, batch_size=100, shuffle=True)
       

# Define the model
class FeatureExtractor(nn.Module):
    def __init__(self, input_dim):
        super(FeatureExtractor, self).__init__()
        self.fc1 = nn.Linear(input_dim, 64)
        self.fc2 = nn.Linear(64, 32)
        
    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        return x

#If you have 5000 genes as input features, the choice of whether to use 5000 neurons in the input layer (or the first layer of your neural network) depends on your goals and the architecture of your model. Here are some considerations:
#
#1. Input Dimension
#Using 5000 Input Features: If your data consists of expression levels for 5000 genes, you should use 5000 input neurons in the input layer to represent each gene. This allows the model to learn from the full feature set.
#2. Model Complexity
#First Layer Size: The first hidden layer (e.g., self.fc1 = nn.Linear(input_dim, 64)) can be adjusted based on your specific needs:
#You might choose 64 neurons if you want to reduce dimensionality and capture the most important features.
#Alternatively, you could use more neurons (e.g., 128 or 256) in the first layer if you want to allow the model to learn more complex representations from the input data.
#3. Feature Selection
#Dimensionality Reduction Techniques: If you have a high-dimensional dataset with many genes, you might consider techniques like PCA (Principal Component Analysis) or t-SNE before feeding the data into the network. This can help reduce the dimensionality while retaining the most important variance in the data.
#4. Avoiding Overfitting
#Regularization: With a high number of features, especially in a small dataset, overfitting can be a concern. Consider using dropout layers, L2 regularization, or other techniques to mitigate this.
#5. Experimentation
#Trial and Error: Ultimately, the best approach often involves experimentation. Start with a model that has 5000 input features and a reasonable number of neurons in the first layer, then test different configurations to see what works best for your specific dataset and task.
#Summary
#If you have 5000 genes, you should use 5000 input features. For the first hidden layer, you can choose to reduce the number of neurons (e.g., to 64) based on your goals for dimensionality reduction and model complexity. Adjusting these parameters based on experimentation and validation performance will help you find the most effective architecture for your task.






class Classifier(nn.Module):
    def __init__(self, num_classes):
        super(Classifier, self).__init__()
        self.fc = nn.Linear(32, num_classes)  # Number of disease statuses
        
    def forward(self, x):
        return self.fc(x)

class DomainClassifier(nn.Module):
    def __init__(self):
        super(DomainClassifier, self).__init__()
        self.fc = nn.Linear(32, 2)  # Two domains (source, target)
        
    def forward(self, x):
        return self.fc(x)

# Training function
def train_dann(source_loader, target_loader, num_epochs=100):
    input_dim = source_features.shape[1]  # Number of features (genes)
    num_classes = len(np.unique(source_labels))  # Number of unique disease statuses
    
    feature_extractor = FeatureExtractor(input_dim)
    classifier = Classifier(num_classes)
    domain_classifier = DomainClassifier()
    
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(list(feature_extractor.parameters()) + 
                           list(classifier.parameters()) + 
                           list(domain_classifier.parameters()), lr=0.001)
    
    for epoch in range(num_epochs):
        feature_extractor.train()
        classifier.train()
        domain_classifier.train()
        
        for (source_data, source_labels), (target_data, _) in zip(source_loader, target_loader):
            optimizer.zero_grad()
            
            # Source domain
            features_s = feature_extractor(source_data)
            class_preds_s = classifier(features_s)
            domain_preds_s = domain_classifier(features_s)
            
            class_loss = criterion(class_preds_s, source_labels.long())
            domain_loss = criterion(domain_preds_s, torch.zeros(len(source_data)).long())  # Source domain label
            
            # Target domain (no labels)
            features_t = feature_extractor(target_data)
            domain_preds_t = domain_classifier(features_t)
            domain_loss_t = criterion(domain_preds_t, torch.ones(len(target_data)).long())  # Target domain label
            
            total_loss = class_loss + domain_loss + domain_loss_t
            total_loss.backward()
            optimizer.step()
        
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {total_loss.item():.4f}')
    
    return feature_extractor, classifier

# Train the model and retrieve the models
feature_extractor, classifier = train_dann(source_loader, target_loader)

# Evaluate the model on the target domain
def evaluate(target_loader, feature_extractor, classifier):
    feature_extractor.eval()
    classifier.eval()
    
    with torch.no_grad():
        for target_data, _ in target_loader:
            features_t = feature_extractor(target_data)
            class_preds_t = classifier(features_t)
            _, predicted = torch.max(class_preds_t, 1)
            print("Predictions on target domain:", predicted.numpy())

# Evaluate the model
evaluate(target_loader, feature_extractor, classifier)



# BASE CODE WITH SYNTHETIC DATA ####################################################################################################################################################################################################################################

# For simplicity, let's create two synthetic datasets: one for the source domain and another for the target domain.
class SyntheticDataset(Dataset):
    def __init__(self, n_samples, domain):
        self.n_samples = n_samples
        self.domain = domain
        
        # Generate data
        if domain == 'source':
            self.data = np.random.randn(n_samples, 2) + np.array([2, 2])
            self.labels = np.zeros(n_samples)  # Class 0
        else:
            self.data = np.random.randn(n_samples, 2) + np.array([-2, -2])
            self.labels = np.ones(n_samples)  # Class 1
            
    def __len__(self):
        return self.n_samples
    
    def __getitem__(self, idx):
        return self.data[idx].astype(np.float32), self.labels[idx].astype(np.float32)

source_dataset = SyntheticDataset(1000, domain='source')
target_dataset = SyntheticDataset(1000, domain='target')

source_loader = DataLoader(source_dataset, batch_size=32, shuffle=True)
target_loader = DataLoader(target_dataset, batch_size=32, shuffle=True)


class SyntheticDataset(Dataset):
    def __init__(self, n_samples, domain):
        self.n_samples = n_samples
        self.domain = domain
        
        # Generate data
        if domain == 'source':
            self.data = np.random.randn(n_samples, 2) + np.array([2, 2])
        else:
            self.data = np.random.randn(n_samples, 2) + np.array([-2, -2])
        
        # Randomly assign class labels 0 or 1
        self.labels = np.random.randint(0, 2, size=n_samples)  # Class label can be 0 or 1
            
    def __len__(self):
        return self.n_samples
    
    def __getitem__(self, idx):
        return self.data[idx].astype(np.float32), self.labels[idx].astype(np.float32)

# Create datasets with random class labels
source_dataset = SyntheticDataset(1000, domain='source')  # Random class labels 0 or 1
target_dataset = SyntheticDataset(1000, domain='target')  # Random class labels 0 or 1

source_loader = DataLoader(source_dataset, batch_size=32, shuffle=True)
target_loader = DataLoader(target_dataset, batch_size=32, shuffle=True)

# Inspecting the first few batches from the source_loader
print("Source Loader Batches:")
for i, (data, labels) in enumerate(source_loader):
    print(f"Batch {i}:")
    print("Data:\n", data.numpy())  # Convert to numpy for easier viewing
    print("Labels:\n", labels.numpy())
    if i == 2:  # Show only the first 3 batches
        break

# Inspecting the first few batches from the target_loader
print("\nTarget Loader Batches:")
for i, (data, labels) in enumerate(target_loader):
    print(f"Batch {i}:")
    print("Data:\n", data.numpy())  # Convert to numpy for easier viewing
    print("Labels:\n", labels.numpy())
    if i == 2:  # Show only the first 3 batches
        break


## Define the model
#We'll define a simple neural network that includes a shared feature extractor and domain classifier.
class FeatureExtractor(nn.Module):
    def __init__(self):
        super(FeatureExtractor, self).__init__()
        self.fc1 = nn.Linear(2, 64)
        self.fc2 = nn.Linear(64, 32)
        
    def forward(self, x):
        x = torch.relu(self.fc1(x))
        x = torch.relu(self.fc2(x))
        return x

class Classifier(nn.Module):
    def __init__(self):
        super(Classifier, self).__init__()
        self.fc = nn.Linear(32, 2)  # Two classes
        
    def forward(self, x):
        return self.fc(x)

class DomainClassifier(nn.Module):
    def __init__(self):
        super(DomainClassifier, self).__init__()
        self.fc = nn.Linear(32, 2)  # Two domains (source, target)
        
    def forward(self, x):
        return self.fc(x)

## Training function
# Define the training loop for the DANN model.

def train_dann(source_loader, target_loader, num_epochs=100):
    feature_extractor = FeatureExtractor()
    classifier = Classifier()
    domain_classifier = DomainClassifier()
    
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(list(feature_extractor.parameters()) + 
                           list(classifier.parameters()) + 
                           list(domain_classifier.parameters()), lr=0.001)
    
    for epoch in range(num_epochs):
        feature_extractor.train()
        classifier.train()
        domain_classifier.train()
        
        for (source_data, source_labels), (target_data, _) in zip(source_loader, target_loader):
            optimizer.zero_grad()
            
            # Source domain
            features_s = feature_extractor(source_data)
            class_preds_s = classifier(features_s)
            domain_preds_s = domain_classifier(features_s)
            
            class_loss = criterion(class_preds_s, source_labels.long())
            domain_loss = criterion(domain_preds_s, torch.zeros(len(source_data)).long())  # Source domain label
            
            # Target domain (no labels)
            features_t = feature_extractor(target_data)
            domain_preds_t = domain_classifier(features_t)
            domain_loss_t = criterion(domain_preds_t, torch.ones(len(target_data)).long())  # Target domain label
            
            total_loss = class_loss + domain_loss + domain_loss_t
            total_loss.backward()
            optimizer.step()
        
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {total_loss.item():.4f}')
    
    return feature_extractor, classifier  # Return models


# Train the model and retrieve the models
feature_extractor, classifier = train_dann(source_loader, target_loader)

## Evaludate the model on the tagret domain
def evaluate(target_loader, feature_extractor, classifier):
    feature_extractor.eval()
    classifier.eval()
    
    with torch.no_grad():
        for target_data, _ in target_loader:
            features_t = feature_extractor(target_data)
            class_preds_t = classifier(features_t)
            _, predicted = torch.max(class_preds_t, 1)
            print("Predictions on target domain:", predicted.numpy())

# Evaluate the model
evaluate(target_loader, feature_extractor, classifier)


####################################################################################################################################################################################################################################


def train_dann(source_loader, target_loader, num_epochs=100):
    feature_extractor = FeatureExtractor()
    classifier = Classifier()
    domain_classifier = DomainClassifier()
    
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(list(feature_extractor.parameters()) + 
                           list(classifier.parameters()) + 
                           list(domain_classifier.parameters()), lr=0.001)
    
    for epoch in range(num_epochs):
        feature_extractor.train()
        classifier.train()
        domain_classifier.train()
        
        for (source_data, source_labels), (target_data, _) in zip(source_loader, target_loader):
            optimizer.zero_grad()
            
            # Source domain
            features_s = feature_extractor(source_data)
            class_preds_s = classifier(features_s)
            domain_preds_s = domain_classifier(features_s)
            
            class_loss = criterion(class_preds_s, source_labels.long())
            domain_loss = criterion(domain_preds_s, torch.zeros(len(source_data)).long())  # Source domain label
            
            # Target domain (no labels)
            features_t = feature_extractor(target_data)
            domain_preds_t = domain_classifier(features_t)
            domain_loss_t = criterion(domain_preds_t, torch.ones(len(target_data)).long())  # Target domain label
            
            total_loss = class_loss + domain_loss + domain_loss_t
            total_loss.backward()
            optimizer.step()
        
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {total_loss.item():.4f}')

train_dann(source_loader, target_loader)


## Evaludate the model on the tagret domain
def evaluate(target_loader, feature_extractor, classifier):
    feature_extractor.eval()
    classifier.eval()
    
    with torch.no_grad():
        for target_data, _ in target_loader:
            features_t = feature_extractor(target_data)
            class_preds_t = classifier(features_t)
            _, predicted = torch.max(class_preds_t, 1)
            print("Predictions on target domain:", predicted.numpy())

# Evaluate the model
evaluate(target_loader, feature_extractor, classifier)


