salloc -A def-tdurcan --time=0-4 -c 1 --mem=40g

## For Narval
module load StdEnv/2020 
module load python/3.8.10

python  


## For Beluga
module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Endo   ##### JOB ID: 43265174

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

nano test_SALS_BA4_regNMF_Endo_KNN_bayesian.sh

#!/bin/bash 
#SBATCH --account=def-sfarhan
#SBATCH --time=00-12:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=200g          # memory per cor
#SBATCH --job-name=test_SALS_BA4_regNMF_Endo_KNN_bayesian
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/test_SALS_BA4_regNMF_Endo_KNN_bayesian.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano test_SALS_BA4_regNMF_Endo_KNN_bayesian.py


# -----------> Job ID: 43321809

import sys
import os

## Add the specific site-packages path to sys.path
pythonenv_path = '/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA/lib/python3.8/site-packages'
sys.path.insert(0, pythonenv_path)

import scanpy as sc
import torch
import torch.nn.functional as F
from torch import nn
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset

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

from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.preprocessing import MaxAbsScaler

import scanpy as sc
import torch
import torch.nn.functional as F
from torch import nn
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset

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

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np
import scipy.linalg

import scanpy.external as sce

from cnmf import cNMF, Preprocess

from sklearn.preprocessing import MaxAbsScaler

import pandas as pd
import numpy as np
import scanpy as sc
import umap
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.backends.backend_pdf import PdfPages
import scipy.optimize
from scipy.spatial import procrustes
from scipy.linalg import orthogonal_procrustes

import torch.nn as nn
import torch.nn.functional as F

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Function
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.model_selection import GroupKFold
from sklearn.decomposition import PCA
import umap
import scanpy as sc
from sklearn.preprocessing import OneHotEncoder
from scipy.linalg import orthogonal_procrustes

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

from torch.utils.data import DataLoader, TensorDataset
#from torch_geometric.data import Data as GeoData
#from torch_geometric.nn import GCNConv
from sklearn.preprocessing import MaxAbsScaler
from sklearn.neighbors import KNeighborsClassifier

import scanpy as sc
import torch
import torch.nn.functional as F
from torch import nn
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset

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

from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.preprocessing import MaxAbsScaler

import scanpy as sc
import torch
import torch.nn.functional as F
from torch import nn
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset

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

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import numpy as np
import scipy.linalg

import scanpy.external as sce

from cnmf import cNMF, Preprocess

from sklearn.preprocessing import MaxAbsScaler

import pandas as pd
import numpy as np
import scanpy as sc
import umap
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.backends.backend_pdf import PdfPages
import scipy.optimize
from scipy.spatial import procrustes
from scipy.linalg import orthogonal_procrustes

import torch.nn as nn
import torch.nn.functional as F

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Function
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.model_selection import GroupKFold
from sklearn.decomposition import PCA
import umap
import scanpy as sc
from sklearn.preprocessing import OneHotEncoder
from scipy.linalg import orthogonal_procrustes


###################################
# Parameters
###################################
par_keep_cell_type = "Endo"
par_prep = "CombatSeq"
par_brain_region = "BA4"
par_status = "SALS"
remove = ["C9ALS", "C9FTLD", "SFTLD"]
par_brain_region_Li = "frontal cortex"
remove_li = ['C9-FTD']


###################################
# Cell Specific parameters -- Endo
###################################
par_ann_data_Pineda = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
par_ann_data_Li = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
par_ann_data_Limone = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"

###################################
# Load anndata -- Pineda
###################################
adata = sc.read_h5ad(par_ann_data_Pineda)
num_cells = adata.X.shape[0]
adata = adata[~adata.obs['Group'].isin(remove)]
adata = adata[adata.obs['CellType'] == par_keep_cell_type]
adata = adata[adata.obs['Region'] == par_brain_region]
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]

# Map disease status
mapping = {'SALS': 1, 'PN': 0}
adata.obs['Group'] = adata.obs['Group'].map(mapping)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

###################################
# Load anndata -- Li
###################################

adata_li = sc.read_h5ad(par_ann_data_Li)
adata_li = adata_li[~adata_li.obs['Group'].isin(remove_li)]
adata_li = adata_li[adata_li.obs['CellType'] == par_keep_cell_type]
adata_li.obs_names = [f"Cell_{i:d}" for i in range(adata_li.n_obs)]

sc.pp.normalize_total(adata_li, target_sum=1e4)
sc.pp.log1p(adata_li)

###################################
# Only keep overlapping genes and Apply MaxAbsScaler
###################################

# Ensure same genes before scaling
common_genes = adata.var_names.intersection(adata_li.var_names)
adata = adata[:, common_genes].copy()
adata.X.shape
adata_li = adata_li[:, common_genes].copy()
adata_li.X.shape #18344

scaler = MaxAbsScaler()
adata.X = scaler.fit_transform(adata.X)
#adata_li.X = scaler.transform(adata_li.X)

sc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor='seurat')
adata = adata[:, adata.var['highly_variable']].copy()
adata.X.shape
#adata_li = adata_li[:, adata.var_names].copy()

###################################
# Filter and QC
###################################
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.X.shape

###################################
# LOSO loop
###################################
    
sample_IDs = set(adata.obs['orig.ident'])
data = []

for donor in sample_IDs:
    print(f"Processing donor: {donor}")
    adata_train = adata[adata.obs['orig.ident'] != donor]
    ## Subset train to have equal number of disease and healthy cells
    group_counts = adata_train.obs['Group'].value_counts()
    min_cells = group_counts.min()
    balanced_indices = (
        adata_train.obs.groupby('Group')
        .sample(n=min_cells, random_state=42)
        .index
    )
    adata_train = adata_train[balanced_indices].copy()
    adata_test = adata[adata.obs['orig.ident'] == donor]
     
    adata_train.X.shape
    adata_test.X.shape
     
    ############################
    ## Apply NMF model to train
    ############################
    ## NMF model
    model = NMF(n_components=100, init='random', random_state=42, max_iter = 10000) 
     
    ## Fit the model
    model.fit(adata_train.X)
     
    ## W = transformed data matrix, V = original feature matrix
    W = model.transform(adata_train.X)
    H = model.components_
    W.shape
    H.shape
     
    H_df = pd.DataFrame(H, columns=adata_train.var_names)
    H_df
     
    W_df = pd.DataFrame(W)
    W_df
     
    ## Check class
    type(W), type(adata_train.X)
     
    ## Create dataframe
    new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
     
    ## Create new anndata object for method 3
    adata_train_NMF = sc.AnnData(X = W.copy(),
    obs = adata_train.obs.copy(),
    var = new_var,
    uns = adata_train.uns.copy(),
    obsm = adata_train.obsm.copy(),
    varm = adata_train.varm.copy(),
    layers = adata_train.layers.copy(),
    raw = adata_train.raw.copy(),
    dtype = "float32",
    shape = None,
    obsp = adata_train.obsp.copy(),
    varp = adata_train.varp
    )
    adata_train_NMF.__dict__['_raw'].__dict__['_var'] = adata_train.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
          
    ## Save method 3 data object
    #c_train.write(par_m3_train)
     
    ## Check the reconstruction error -- it is okay ish
    #model.reconstruction_err_
     
    #top_genes_dict = {}
    #for i in range(H_df.shape[0]):  # Iterate through each factor
    #    top_genes = H_df.iloc[i].nlargest(H_df.shape[1])
    #    top_genes_dict[f'feature_{i}'] = top_genes
    #top_genes_df = pd.DataFrame(top_genes_dict)
    #top_genes_df.to_csv(par_top_genes_factor)
     
    ############################
    ## Apply NMF model to test
    ############################
    ## NMF model
    W = model.transform(adata_test.X)
    H = model.components_
    W.shape
    H.shape
     
    H_df = pd.DataFrame(H, columns=adata_test.var_names)
    H_df
     
    W_df = pd.DataFrame(W)
    W_df
     
    ## Check class
    type(W), type(adata_test.X)
     
    ## Create dataframe
    new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
     
    ## Create new anndata object for method 3
    adata_test_NMF = sc.AnnData(X = W.copy(),
    obs = adata_test.obs.copy(),
    var = new_var,
    uns = adata_test.uns.copy(),
    obsm = adata_test.obsm.copy(),
    varm = adata_test.varm.copy(),
    layers = adata_test.layers.copy(),
    raw = adata_test.raw.copy(),
    dtype = "float32",
    shape = None,
    obsp = adata_test.obsp.copy(),
    varp = adata_test.varp
    )
    adata_test_NMF.__dict__['_raw'].__dict__['_var'] = adata_test.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
     
    ## Lets try to visualize the embeddings
    X_train = adata_train_NMF.X
    X_test = adata_test_NMF.X
     
    X_train.shape   
    X_test.shape  
     
    X_all = np.concatenate([X_train, X_test], axis=0)
     
    cell_metadata_train = adata_train.obs[['orig.ident', 'Group']].copy()
    cell_metadata_test = adata_test.obs[['orig.ident', 'Group']].copy()
     
    # Map Group numbers to names
    group_mapping = {0: 'PN', 1: 'SALS'}
    cell_metadata_train['Group'] = cell_metadata_train['Group'].map(group_mapping)
    cell_metadata_test['Group'] = cell_metadata_test['Group'].map(group_mapping)
     
    cell_metadata_train['Set'] = 'Train'
    cell_metadata_test['Set'] = 'Test'
     
    cell_metadata = pd.concat([cell_metadata_train, cell_metadata_test], axis=0).reset_index(drop=True)
     
    # Define a new label: "PlotGroup"
    cell_metadata['PlotGroup'] = np.where(
        cell_metadata['Set'] == 'Test', 'Test', cell_metadata['Group']
    )
                 
    ###################################
    # Implementation into the LOSO DNN
    ###################################
     
    num_genes = adata_train_NMF.X.shape[1]
                 
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
            
    # Bayesian loop
    learning_rate = 0.0008663099696291
    kNN_threshold_list = [.99, .95, .9, .85, .8]
    for kNN_threshold in kNN_threshold_list:
        for _ in range(5):
            print(f"Processing donor: {donor}")
            num_cells = adata_test_NMF.shape[0]
        
            X_train = torch.FloatTensor(adata_train_NMF.X.toarray() if hasattr(adata_train_NMF.X, 'toarray') else adata_train_NMF.X)
            y_train = torch.LongTensor(adata_train_NMF.obs['Group'].values)
            domains = pd.factorize(adata_train_NMF.obs['Donor'])[0]
            d_train = torch.LongTensor(domains)
        
            dataset = TensorDataset(X_train, y_train, d_train)
            input_size = X_train.shape[1]
            num_classes = len(np.unique(y_train))
            num_domains = len(np.unique(domains))
        
            #adata_train = adata[adata.obs['orig.ident'] != donor]
            #adata_test = adata[adata.obs['orig.ident'] == donor]
            #num_cells = adata_test.shape[0]
            
            #X_train = torch.FloatTensor(adata_train.X.toarray() if hasattr(adata_train.X, 'toarray') else adata_train.X)
            #y_train = torch.LongTensor(adata_train.obs['Group'].values)
            #domains = pd.factorize(adata_train.obs['Donor'])[0]
            #d_train = torch.LongTensor(domains)
            
            #dataset = TensorDataset(X_train, y_train, d_train)
            #input_size = adata_train.shape[1]
            #num_classes = len(np.unique(y_train))
            #num_domains = len(np.unique(domains))
            
            loader = DataLoader(dataset, batch_size=64, shuffle=True)
            
            model_main = MainModel(input_size, num_classes)
            model_domain = DomainClassifier(25, num_domains)
            device = 'cuda' if torch.cuda.is_available() else 'cpu'
            
            train(model_main, model_domain, loader, 25, device)
            
            # Evaluation
            #X_test = torch.FloatTensor(adata_test.X.toarray() if hasattr(adata_test.X, 'toarray') else adata_test.X)
            #y_test = adata_test.obs['Group']
            X_test = torch.FloatTensor(adata_test_NMF.X.toarray() if hasattr(adata_test_NMF.X, 'toarray') else adata_test_NMF.X)
            y_test = adata_test_NMF.obs['Group']
            
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
                    'learning_rate': learning_rate, 'batch_size': 64,
                    'test_accuracy_high': acc_high, 'test_accuracy_low': acc_low,
                    'test_accuracy_all': acc_all, 'method': 'LIME + kNN-filtered + prior-adjusted', 'kNN_thresh': kNN_threshold
                })
        
    # Delete obects to conserve memory
    del adata_train, adata_test
    del adata_train_NMF, adata_test_NMF
    del X_train, X_test, X_all, cell_metadata_train, cell_metadata_test, cell_metadata
    del model_main, model_domain

# Save results
pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data)

out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/regNMF_model_report2_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_{par_n_topics}_{par_n_epochs}_narval_2_KNN_Bayesian.csv"
results_df.to_csv(out_path, index=False)
        


