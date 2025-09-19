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

# Map disease status
mapping = {'SALS': 1, 'PN': 0}
adata.obs['Group'] = adata.obs['Group'].map(mapping)

# Preprocess data
sc.pp.filter_cells(adata, min_genes=200)
adata = adata[:, features]
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
scaler = MaxAbsScaler()
adata.X = scaler.fit_transform(adata.X)

adata.raw = adata

# Subset to optimal features
adata = adata[:, features]
num_genes = adata.X.shape[1]
data = []

# Define models
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

class DomainClassifier(nn.Module):
    def __init__(self, input_size, num_domains):
        super().__init__()
        self.model = nn.Sequential(
            nn.Linear(input_size, 25), nn.ReLU(),
            nn.Linear(25, num_domains)
        )
    
    def forward(self, x):
        return self.model(x)

#def mixup_data_classwise(X, y, alpha=0.4, n_synthetic=1000):
#    N = X.shape[0]
#    y = y.numpy()
#    if y.ndim > 1 and y.shape[1] > 1:
#        class_indices = np.argmax(y, axis=1)
#    else:
#        class_indices = y.flatten()
    
#    synthetic_X = []
#    synthetic_y = []
#    unique_classes = np.unique(class_indices)
#    class_to_indices = {cls: np.where(class_indices == cls)[0] for cls in unique_classes}
    
#    for _ in range(n_synthetic):
#        cls = np.random.choice(unique_classes)
#        indices = class_to_indices[cls]
#        if len(indices) < 2:
#            continue
#        i1, i2 = np.random.choice(indices, size=2, replace=False)
#        lam = np.random.beta(alpha, alpha)
#        x_mix = lam * X[i1] + (1 - lam) * X[i2]
#        y_mix = lam * y[i1] + (1 - lam) * y[i2]
#        synthetic_X.append(x_mix)
#        synthetic_y.append(y_mix)
    
#    if len(synthetic_X) > 0:
#        synthetic_X = torch.FloatTensor(np.stack(synthetic_X))
#        synthetic_y = torch.LongTensor(np.round(np.stack(synthetic_y)).astype(int))
#        X_aug = torch.cat([X, synthetic_X], dim=0)
#        y_aug = torch.cat([torch.LongTensor(y), synthetic_y], dim=0)
#    else:
#        X_aug = torch.FloatTensor(X)
#        y_aug = torch.LongTensor(y)
    
#    return X_aug, y_aug


def mixup_data_classwise(X, y, domain, alpha=0.4, n_synthetic=1000):
    N = X.shape[0]
    
    # Ensure y and domain are numpy arrays
    if isinstance(y, torch.Tensor):
        y = y.numpy()
    if isinstance(domain, torch.Tensor):
        domain = domain.numpy()
    
    # Convert one-hot to class index if needed
    if y.ndim > 1 and y.shape[1] > 1:
        class_indices = np.argmax(y, axis=1)
    else:
        class_indices = y.flatten()
    
    synthetic_X = []
    synthetic_y = []
    synthetic_domain = []
    
    unique_classes = np.unique(class_indices)
    class_to_indices = {cls: np.where(class_indices == cls)[0] for cls in unique_classes}
    
    for _ in range(n_synthetic):
        cls = np.random.choice(unique_classes)
        indices = class_to_indices[cls]
        if len(indices) < 2:
            continue
        i1, i2 = np.random.choice(indices, size=2, replace=False)
        lam = np.random.beta(alpha, alpha)
        x_mix = lam * X[i1] + (1 - lam) * X[i2]
        y_mix = lam * y[i1] + (1 - lam) * y[i2]
        d_mix = domain[i1]  # <--- assign domain of i1
        synthetic_X.append(x_mix)
        synthetic_y.append(y_mix)
        synthetic_domain.append(d_mix)
    
    if len(synthetic_X) > 0:
        synthetic_X = torch.FloatTensor(np.stack(synthetic_X))
        synthetic_y = torch.LongTensor(np.round(np.stack(synthetic_y)).astype(int))
        synthetic_domain = torch.LongTensor(synthetic_domain)
     
        X_aug = torch.cat([torch.FloatTensor(X), synthetic_X], dim=0)
        y_aug = torch.cat([torch.LongTensor(y), synthetic_y], dim=0)
        domain_aug = torch.cat([torch.LongTensor(domain), synthetic_domain], dim=0)
    else:
        X_aug = torch.FloatTensor(X)
        y_aug = torch.LongTensor(y)
        domain_aug = torch.LongTensor(domain)
    
    return X_aug, y_aug, domain_aug


# Train function
def train(main_model, domain_model, dataloader, epochs, device):
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
        total_batches = 0
    
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
            total_batches += 1
    
        avg_loss_c = total_loss_c / total_batches
        avg_loss_d = total_loss_d / total_batches
        print(f"Epoch {epoch+1}/{epochs}: Classification Loss = {avg_loss_c:.4f}, Domain Loss = {avg_loss_d:.4f}")


# LOSO loop
for donor in sample_IDs:
    for _ in range(5):
        print(f"Processing donor: {donor}")
        adata_train = adata[adata.obs['orig.ident'] != donor]
        adata_test = adata[adata.obs['orig.ident'] == donor]
        num_cells = adata_test.shape[0]
        
        X_train = torch.FloatTensor(adata_train.X.toarray() if hasattr(adata_train.X, 'toarray') else adata_train.X)
        y_train = torch.LongTensor(adata_train.obs['Group'].values)
        domains = pd.factorize(adata_train.obs['Donor'])[0]
        d_train = torch.LongTensor(domains)
        
        # Mixup augmentation and append to real data
        X_synth, y_synth, d_synth = mixup_data_classwise(X_train, y_train, d_train,  alpha=0.4, n_synthetic=10000)
        
        dataset = TensorDataset(X_synth, y_synth, d_synth)
        input_size = adata_train.shape[1]
        num_classes = len(np.unique(y_synth))
        num_domains = len(np.unique(d_synth))
        
        loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
        
        model_main = MainModel(input_size, num_classes)
        model_domain = DomainClassifier(25, num_domains)
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        
        train(model_main, model_domain, loader, 10, device)
        
        # Evaluation
        X_test = torch.FloatTensor(adata_test.X.toarray() if hasattr(adata_test.X, 'toarray') else adata_test.X)
        y_test = adata_test.obs['Group']
        
        model_main.eval()
        with torch.no_grad():
            y_pred, _ = model_main(X_test)
            y_labels = torch.argmax(y_pred, dim=1)
            print(classification_report(y_test, y_labels.numpy()))
        
        acc = balanced_accuracy_score(y_test.values, y_labels.numpy())
        data.append({
            'prep': par_prep, 'donor': donor, 'region': par_brain_region,
            'group': par_status, 'celltype': par_keep_cell_type,
            'n_genes': num_genes, 'n_cells': num_cells,
            'learning_rate': learning_rate, 'batch_size': batch_size,
            'test_accuracy': acc, 'method': 'LIME'
        })

# Save results
pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_abs_scaling_narval_2_mixup_augmentation.csv"
results_df.to_csv(out_path, index=False)