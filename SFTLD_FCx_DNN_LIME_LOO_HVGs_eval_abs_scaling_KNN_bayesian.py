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


nano HVGs_LOO_0delta_absScaling_weighted_Knn_threshold_sweep_LIME_SFTLD_BA9_all_cells_narval_2_1.sh

#!/bin/bash  
#SBATCH --account=def-grouleau
#SBATCH --time=03-00:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=150g          # memory per cor
#SBATCH --job-name=HVGs_LOO_0delta_absScaling_weighted_Knn_threshold_sweep_LIME_SFTLD_BA9_all_cells_narval_2_1
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/HVGs_LOO_0delta_absScaling_weighted_Knn_threshold_sweep_LIME_SFTLD_BA9_all_cells_narval_2_1.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano HVGs_LOO_0delta_absScaling_weighted_Knn_threshold_sweep_LIME_SFTLD_BA9_all_cells_narval_2_1.py



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
from sklearn.neighbors import KNeighborsClassifier


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
    
    # Load anndata
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
    
    # Subset to HVG features
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
    sample_IDs = set(adata.obs['orig.ident'])
    kNN_threshold_list = [.99, .95, .9, .85, .8]
    
    for kNN_threshold in kNN_threshold_list:
        print(f"KNN Threshold: {kNN_threshold}")
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
    
                train(model_main, model_domain, loader, 10, device)
    
                # Evaluation
                X_test = torch.FloatTensor(adata_test.X.toarray() if hasattr(adata_test.X, 'toarray') else adata_test.X)
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
    
                    data.append({
                        'prep': par_prep, 'donor': donor, 'region': par_brain_region,
                        'group': par_status, 'celltype': par_keep_cell_type,
                        'n_genes': num_genes, 'n_cells': len(y_true),
                        'learning_rate': learning_rate, 'batch_size': batch_size,
                        'test_accuracy_high': acc_high, 'test_accuracy_low': acc_low,
                        'test_accuracy_all': acc_all, 'method': 'LIME + kNN-filtered + prior-adjusted', 'kNN_thresh': kNN_threshold
                    })
    
    # Save results
    pd.set_option('display.max_rows', None)
    results_df = pd.DataFrame(data)
    out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVGs_gene_set_weighted_LOO_LIME_0_abs_case_control_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_abs_scaling_sweep_confidence_threshold_narval_2.csv"
    results_df.to_csv(out_path, index=False)