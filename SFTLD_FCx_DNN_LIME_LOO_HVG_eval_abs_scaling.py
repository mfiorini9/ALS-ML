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

nano HVGs_LOO_absScaling_SFTLD_BA9_all_cells_narval_2.sh

#!/bin/bash  
#SBATCH --account=def-grouleau
#SBATCH --time=00-12:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=150g          # memory per cor
#SBATCH --job-name=HVGs_LOO_absScaling_SFTLD_BA9_all_cells_narval_2
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/HVGs_LOO_absScaling_SFTLD_BA9_all_cells_narval_2.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano HVGs_LOO_absScaling_SFTLD_BA9_all_cells_narval_2.py


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
par_brain_region = "BA9"
par_status = "SFTLD"
remove = ["C9ALS", "C9FTLD", "SALS"]

cell_types = [
    "L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Oligo", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5"
]

# Loop through the list and print each item
for par_keep_cell_type in cell_types:
    print(par_keep_cell_type)
    
    # Load model report
    report_file = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
    report_df = pd.read_csv(report_file)
    batch_size = int(report_df.at[0, 'batch_size'])
    learning_rate = report_df.at[0, 'learning_rate']
    
    # Set anndata path
    adata_file = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
    
    ##################################################################################################################################################################################################################
    # Get original HVGs
    ##################################################################################################################################################################################################################
    adata = sc.read_h5ad(adata_file)
    adata = adata[~adata.obs['Group'].isin(remove)]
    adata = adata[adata.obs['CellType'] == par_keep_cell_type]
    adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
    print("genes: ", adata.var_names) 
    print("cells: ", adata.obs_names) 
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    
    ## create a list of sample IDs
    sample_ID = set(adata.obs['orig.ident'])
    
    # clean AnnData object
    adata.obs = adata.obs.reset_index() 
    adata.obs.columns
    set(adata.obs['Group'])
    set(adata.obs['CellType'])
    set(adata.obs['Region'])
    
    ## change Disease status labels
    mapping = {'SFTLD': 1, 'PN': 0}
    
    # Map the values
    adata.obs['Group'] = adata.obs['Group'].map(mapping)
    
    # Verify the change
    print(adata.obs['Group'].value_counts())
    
    # Create an array of indices
    indices = np.arange(adata.shape[0])
    
    # Split indices into train and test
    train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=42)
    
    # Create training and testing AnnData objects
    adata_train = adata[train_indices].copy()
    adata_test = adata[test_indices].copy()
    
    ## Process train data
    sc.pp.filter_cells(adata_train, min_genes=200)  # Filter cells
    sc.pp.filter_genes(adata_train, min_cells=3)    # Filter genes
    sc.pp.normalize_total(adata_train, target_sum=1e4)  # Normalize
    sc.pp.log1p(adata_train)                         # Log-transform
    sc.pp.highly_variable_genes(adata_train) 
    adata_train.raw = adata_train  # Save raw counts for further analysis
    adata_train = adata_train[:, adata_train.var['highly_variable']]  # Subset to HVGs
    num_HVGs = adata_train.X.shape[1]
    
    HVGs_old = adata_train.var['features'].tolist()
    
    ## print to CSV for downstream use
    HVGs_old_df = pd.DataFrame(HVGs_old)
    HVGs_old_df.columns = ['HVGs']
    out_path_HVGs = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
    HVGs_old_df.to_csv(out_path_HVGs, index=False)
    
    ##################################################################################################################################################################################################################
    # Regular workflow
    ##################################################################################################################################################################################################################
    
    ## Read in adata
    adata = sc.read_h5ad(adata_file)
    
    # Filter data
    adata = adata[~adata.obs['Group'].isin(remove)]
    adata = adata[adata.obs['CellType'] == par_keep_cell_type]
    adata = adata[adata.obs['Region'] == par_brain_region]
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
    sample_IDs = set(adata.obs['orig.ident'])
    
    # Map disease status
    mapping = {'SFTLD': 1, 'PN': 0}
    adata.obs['Group'] = adata.obs['Group'].map(mapping)
    
    # Preprocess data
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[:, HVGs_old]                                       #### TRYING WITH THIS HERE, we will see how it affects the performance. 
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    scaler = MaxAbsScaler()
    adata.X = scaler.fit_transform(adata.X)
    
    adata.raw = adata
    
    # Subset to HVGs 
    adata = adata[:, HVGs_old]
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
        
            dataset = TensorDataset(X_train, y_train, d_train)
            input_size = adata_train.shape[1]
            num_classes = len(np.unique(y_train))
            num_domains = len(np.unique(domains))
        
            loader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
        
            model_main = MainModel(input_size, num_classes)
            model_domain = DomainClassifier(25, num_domains)
            device = 'cuda' if torch.cuda.is_available() else 'cpu'
        
            train(model_main, model_domain, loader, 10, device)      ##### NUM EPOCH IS HERE
        
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
    out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_LOO_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_abs_scaling_narval_2.csv"
    results_df.to_csv(out_path, index=False)