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


nano optimal_LOO_0delta_absScaling_weighted_Knn_threshold_sweep_LIME_SALS_BA9_all_cells_narval_2_1.sh

#!/bin/bash  
#SBATCH --account=def-grouleau
#SBATCH --time=03-00:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=150g          # memory per cor
#SBATCH --job-name=optimal_LOO_0delta_absScaling_weighted_Knn_threshold_sweep_LIME_SALS_BA9_all_cells_narval_2_1
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/optimal_LOO_0delta_absScaling_weighted_Knn_threshold_sweep_LIME_SALS_BA9_all_cells_narval_2_1.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano optimal_LOO_0delta_absScaling_weighted_Knn_threshold_sweep_LIME_SALS_BA9_all_cells_narval_2_1.py



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


################################################################################################ 

###################################
# Parameters
###################################
par_prep = "CombatSeq"
par_brain_region = "BA9"
par_status = "C9ALS"
remove = ["SFTLD", "C9FTLD", "SALS"]
par_brain_region_Li = "frontal cortex" ## Need to fix this. 
remove_li = ['C9-FTD']

#cell_types = [
#    "L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
#    "Rosehip", "SOM", "Oligo", "Astro", "OPC", "Micro",
#    "Mural", "Endo", "Fibro", "L5"
#]

# Loop through the list and print each item
#for par_keep_cell_type in cell_types:

par_keep_cell_type =  "L2_L3"


###################################
# Cell Specific parameters -- Mural
###################################
par_ann_data_Pineda = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
par_ann_data_Li = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
par_ann_data_Limone = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"


###################################
# Load information
###################################
print(par_keep_cell_type)  

# Load LIME-selected genes
lime_file = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_0_abs_case_control_narval_2.csv'
features = pd.read_csv(lime_file)["gene"].tolist()

# Load model report
report_file = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
report_df = pd.read_csv(report_file)
batch_size = int(report_df.at[0, 'batch_size'])
learning_rate = report_df.at[0, 'learning_rate']

###################################
# Load anndata -- Pineda
###################################
adata = sc.read_h5ad(par_ann_data_Pineda)
num_cells = adata.X.shape[0]
adata = adata[~adata.obs['Group'].isin(remove)]
adata = adata[adata.obs['CellType'] == par_keep_cell_type]
adata = adata[adata.obs['Region'] == par_brain_region]
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
sc.pp.filter_cells(adata, min_genes=200) ## <------ moved here from line 272

set(adata.obs['Group'])
set(adata.obs['CellType'])
set(adata.obs['Region'])

# Map disease status
mapping = {'C9ALS': 1, 'PN': 0}
adata.obs['Group'] = adata.obs['Group'].map(mapping)
set(adata.obs['Group'])

#sc.pp.normalize_total(adata, target_sum=1e4)
#sc.pp.log1p(adata)

###################################
# Load anndata -- Li
###################################
adata_li = sc.read_h5ad(par_ann_data_Li)
adata_li = adata_li[~adata_li.obs['Group'].isin(remove_li)]
adata_li = adata_li[adata_li.obs['CellType'] == par_keep_cell_type]
adata_li.obs_names = [f"Cell_{i:d}" for i in range(adata_li.n_obs)]
sc.pp.filter_cells(adata_li, min_genes=200) ## <------ moved here from line 272

set(adata_li.obs['Group'])
set(adata_li.obs['CellType'])
set(adata_li.obs['Region'])

# Map disease status
mapping = {'C9-ALS': 1, 'Control': 0}
adata_li.obs['Group'] = adata_li.obs['Group'].map(mapping)
set(adata_li.obs['Group'])

#sc.pp.normalize_total(adata_li, target_sum=1e4)
#sc.pp.log1p(adata_li)

###################################
# Overlapping genes across datasets
###################################
common_genes = adata.var_names.intersection(adata_li.var_names)
keep_genes = common_genes.intersection(features)
len(keep_genes) ## USE THIS ONE
len(features)


###################################
# Preprocess pineda
###################################
adata = adata[:, keep_genes]  
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
scaler = MaxAbsScaler()
adata.X = scaler.fit_transform(adata.X)

adata.raw = adata

keep_genes = adata.var_names ## redefine after applying min cell expression cutoff
len(keep_genes)

num_genes = adata.X.shape[1]
data = []

###################################
# Preprocess Li
###################################
adata_li = adata_li[:, keep_genes]  
sc.pp.normalize_total(adata_li, target_sum=1e4)
sc.pp.log1p(adata_li)
sc.pp.highly_variable_genes(adata_li)
adata_li.X = scaler.transform(adata_li.X)

adata_li.raw = adata_li

###################################
# Preliminary check
###################################
len(adata.var_names) == len(adata_li.var_names)

# Define models
#class MainModel(nn.Module):
#    def __init__(self, input_size, num_classes):
#        super().__init__()
 #       self.shared = nn.Sequential(
 #           nn.Linear(input_size, 100), nn.ReLU(),
 #           nn.Linear(100, 50), nn.ReLU(),
 #           nn.Linear(50, 25), nn.ReLU()
 #       )
 #       self.classifier = nn.Linear(25, num_classes)
    
 #   def forward(self, x):
 #       shared = self.shared(x)
 #       return self.classifier(shared), shared

#class DomainClassifier(nn.Module):
#    def __init__(self, input_size, num_domains):
#        super().__init__()
#        self.model = nn.Sequential(
#            nn.Linear(input_size, 25), nn.ReLU(),
#            nn.Linear(25, num_domains)
#        )
#    
#    def forward(self, x):
#        return self.model(x)

# Old Train function
#def train(main_model, domain_model, dataloader, epochs, device):
#    main_model.to(device)
#    domain_model.to(device)
#    optimizer_main = optim.Adam(main_model.parameters(), lr=learning_rate)
#    optimizer_domain = optim.Adam(domain_model.parameters(), lr=0.001)
#    loss_class = nn.CrossEntropyLoss()
#    loss_domain = nn.CrossEntropyLoss()
#   
#    for epoch in range(epochs):
#        main_model.train()
#        domain_model.train()
#        for X, y, d in dataloader:
#            X, y, d = X.to(device), y.to(device), d.to(device)
#            y_out, shared = main_model(X)
#            loss_c = loss_class(y_out, y)
#            d_out = domain_model(shared.detach())
#            loss_d = loss_domain(d_out, d)
#            loss = loss_c - loss_d
#    
#            optimizer_main.zero_grad()
#            loss.backward(retain_graph=True)
#           optimizer_main.step()
#    
#            optimizer_domain.zero_grad()
#            loss_d.backward()
#            optimizer_domain.step()



# New Train function
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
        num_batches = 0
        correct = 0
        total = 0
        
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
            num_batches += 1
        
            # Accuracy computation
            _, predicted = torch.max(y_out, 1)
            correct += (predicted == y).sum().item()
            total += y.size(0)
        
        avg_loss_c = total_loss_c / num_batches
        avg_loss_d = total_loss_d / num_batches
        accuracy = 100.0 * correct / total
        
        print(f"Epoch {epoch}, Avg Class Loss: {avg_loss_c:.4f}, Avg Domain Loss: {avg_loss_d:.4f}, Accuracy: {accuracy:.2f}%")
            
######################################################
## Using stuff below here. 
######################################################

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


## Train function with early stopping
def train_earlystop(main_model, domain_model, dataloader, epochs, device, patience=5):
    main_model.to(device)
    domain_model.to(device)
    optimizer_main = optim.Adam(main_model.parameters(), lr=learning_rate)
    optimizer_domain = optim.Adam(domain_model.parameters(), lr=0.001)
    loss_class = nn.CrossEntropyLoss()
    loss_domain = nn.CrossEntropyLoss()
    
    best_train_acc = 0.0
    patience_counter = 0
    best_main_state = None
    best_domain_state = None
    
    for epoch in range(epochs):
        main_model.train()
        domain_model.train()
        total_loss_c = 0.0
        total_loss_d = 0.0
        num_batches = 0
        correct = 0
        total = 0
    
        for X, y, d in dataloader:
            X, y, d = X.to(device), y.to(device), d.to(device)
    
            # Forward pass
            y_out, shared = main_model(X)
            loss_c = loss_class(y_out, y)
            d_out = domain_model(shared.detach())
            loss_d = loss_domain(d_out, d)
            loss = loss_c - loss_d
    
            # Update main model
            optimizer_main.zero_grad()
            loss.backward(retain_graph=True)
            optimizer_main.step()
    
            # Update domain model
            optimizer_domain.zero_grad()
            loss_d.backward()
            optimizer_domain.step()
    
            total_loss_c += loss_c.item()
            total_loss_d += loss_d.item()
            num_batches += 1
    
            # Accuracy
            _, predicted = torch.max(y_out, 1)
            correct += (predicted == y).sum().item()
            total += y.size(0)
    
        avg_loss_c = total_loss_c / num_batches
        avg_loss_d = total_loss_d / num_batches
        train_acc = 100.0 * correct / total
    
        print(f"Epoch {epoch}, "
              f"Avg Class Loss: {avg_loss_c:.4f}, "
              f"Avg Domain Loss: {avg_loss_d:.4f}, "
              f"Accuracy: {train_acc:.2f}%")
    
        # Early stopping check (based on training accuracy)
        if train_acc > best_train_acc:
            best_train_acc = train_acc
            best_main_state = main_model.state_dict()
            best_domain_state = domain_model.state_dict()
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch} "
                      f"(Best Train Acc: {best_train_acc:.2f}%)")
                main_model.load_state_dict(best_main_state)
                domain_model.load_state_dict(best_domain_state)
                return main_model, domain_model
                break

# NEED TO EVALUATE DIFFERENT TRAINING EPOCH VALUES AT SOME POINT. 


######################################################
## TEMPORARY TEST WITH GRL
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
######################################################
import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Function

# Gradient Reversal Layer definition
class GradientReversalFunction(Function):
    @staticmethod
    def forward(ctx, x, lambda_):
        ctx.lambda_ = lambda_
        return x.view_as(x)  # Identity in forward pass
     
    @staticmethod
    def backward(ctx, grad_output):
        # Reverse gradient in backward pass scaled by lambda
        return grad_output.neg() * ctx.lambda_, None

class GradientReversalLayer(nn.Module):
    def __init__(self, lambda_=1.0):
        super().__init__()
        self.lambda_ = lambda_
    
    def forward(self, x):
        return GradientReversalFunction.apply(x, self.lambda_)

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


# Train function with GRL and early stopping
def train_earlystop(main_model, domain_model, dataloader, epochs, device, learning_rate, patience=5,  grl_lambda=1.0):
    main_model.to(device)
    domain_model.to(device)
    grl = GradientReversalLayer(lambda_=grl_lambda).to(device)
    
    optimizer_main = optim.Adam(main_model.parameters(), lr=learning_rate)
    optimizer_domain = optim.Adam(domain_model.parameters(), lr=0.001)
    loss_class = nn.CrossEntropyLoss()
    loss_domain = nn.CrossEntropyLoss()
    
    best_train_acc = 0.0
    patience_counter = 0
    best_main_state = None
    best_domain_state = None
    
    for epoch in range(epochs):
        main_model.train()
        domain_model.train()
        total_loss_c = 0.0
        total_loss_d = 0.0
        num_batches = 0
        correct = 0
        total = 0
    
        for X, y, d in dataloader:
            X, y, d = X.to(device), y.to(device), d.to(device)
    
            # Forward pass
            y_out, shared = main_model(X)
    
            # Pass through GRL before domain classifier
            reversed_shared = grl(shared)
            d_out = domain_model(reversed_shared)
    
            loss_c = loss_class(y_out, y)
            loss_d = loss_domain(d_out, d)
    
            # Combine losses (add, since GRL reverses gradients)
            loss = loss_c + loss_d
        
            # Update both models jointly
            optimizer_main.zero_grad()
            optimizer_domain.zero_grad()
            loss.backward()
            optimizer_main.step()
            optimizer_domain.step()
        
            total_loss_c += loss_c.item()
            total_loss_d += loss_d.item()
            num_batches += 1
    
            # Accuracy on main task
            _, predicted = torch.max(y_out, 1)
            correct += (predicted == y).sum().item()
            total += y.size(0)
    
        avg_loss_c = total_loss_c / num_batches
        avg_loss_d = total_loss_d / num_batches
        train_acc = 100.0 * correct / total
    
        print(f"Epoch {epoch}, "
              f"Avg Class Loss: {avg_loss_c:.4f}, "
              f"Avg Domain Loss: {avg_loss_d:.4f}, "
              f"Accuracy: {train_acc:.2f}%")
    
        # Early stopping check (based on training accuracy)
        if train_acc > best_train_acc:
            best_train_acc = train_acc
            best_main_state = main_model.state_dict()
            best_domain_state = domain_model.state_dict()
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"Early stopping at epoch {epoch} "
                      f"(Best Train Acc: {best_train_acc:.2f}%)")
                main_model.load_state_dict(best_main_state)
                domain_model.load_state_dict(best_domain_state)
                return main_model, domain_model
    
    # If loop finishes without early stop, return best models anyway
    main_model.load_state_dict(best_main_state)
    domain_model.load_state_dict(best_domain_state)
    return main_model, domain_model



# LOSO loop
sample_IDs = set(adata.obs['orig.ident'])
#kNN_threshold_list = [.99, .95, .9, .85, .8]
kNN_threshold = 0.9
print(f"KNN Threshold: {kNN_threshold}")
training_epoch = 10

data2 = []

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
     
        model_main, model_domain = train_earlystop(model_main, model_domain, loader, 10, device, learning_rate, patience=5, grl_lambda=0.5)
     
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
     
            data2.append({
                'prep': par_prep, 'donor': donor, 'region': par_brain_region,
                'group': par_status, 'celltype': par_keep_cell_type,
                'n_genes': num_genes, 'n_cells': len(y_true), 'train_epoch': training_epoch,
                'learning_rate': learning_rate, 'batch_size': batch_size,
                'test_accuracy_high': acc_high, 'test_accuracy_low': acc_low,
                'test_accuracy_all': acc_all, 'method': 'LIME + kNN-filtered + prior-adjusted', 'kNN_thresh': kNN_threshold
            })






######################################################
######################################################
######################################################
######################################################
######################################################
######################################################


###################################
# Perform KNN LOSO with Pineda for model validation
###################################

# LOSO loop
sample_IDs = set(adata.obs['orig.ident'])
#kNN_threshold_list = [.99, .95, .9, .85, .8]
kNN_threshold = 0.9
print(f"KNN Threshold: {kNN_threshold}")
training_epoch = 10
    
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
     
        train(model_main, model_domain, loader, training_epoch, device)
     
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
                'n_genes': num_genes, 'n_cells': len(y_true), 'train_epoch': training_epoch,
                'learning_rate': learning_rate, 'batch_size': batch_size,
                'test_accuracy_high': acc_high, 'test_accuracy_low': acc_low,
                'test_accuracy_all': acc_all, 'method': 'LIME + kNN-filtered + prior-adjusted', 'kNN_thresh': kNN_threshold
            })

###################################
# Train model on entire Pineda dataset
# NEED TO EVALUATE DIFFERENT TRAINING EPOCH VALUES AT SOME POINT. 
###################################
adata_train = adata

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

#train(model_main, model_domain, loader, training_epoch, device)       

main_model, domain_model = train_earlystop(model_main, model_domain, loader, 100, device, learning_rate, patience=5, grl_lambda=0.5)



#all_loader = DataLoader(
#        TensorDataset(torch.tensor(X, dtype=torch.float32),
#                      torch.tensor(y, dtype=torch.long),
#                      torch.tensor(d, dtype=torch.long)),
#        batch_size=batch_size, shuffle=True
#    )
#    final_main_model = MainModel(input_size=X.shape[1], num_classes=len(np.unique(y)))
#    final_domain_model = DomainClassifier(input_size=25, num_domains=len(np.unique(d)))
    
#    train(final_main_model, final_domain_model, all_loader, epochs=epochs, device=device)
    
    # Save the final trained model
#    torch.save(final_main_model.state_dict(), "final_main_model.pth")
#    torch.save(final_domain_model.state_dict(), "final_domain_model.pth")
#    print("Final models saved.")

#    return fold_accuracies, final_main_model, final_domain_model



###################################
# Apply transfer learning on all Li subjects except that which is held out for model evaluation
###################################

###################################
# Evaluate model on left out Li subject with KNN bayesina approach
###################################


# Save results
pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_abs_scaling_sweep_confidence_threshold_narval_2.csv"
results_df.to_csv(out_path, index=False)