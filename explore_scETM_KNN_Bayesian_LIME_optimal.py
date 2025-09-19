salloc -A def-grouleau --time=0-8 -c 1 --mem=40g

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

nano test_SALS_BA4_regETM_Endo_epochs_sweep_200topics_narval_4.sh

#!/bin/bash 
#SBATCH --account=def-grouleau
#SBATCH --time=00-10:00           # time (DD-HH:MM)
#SBATCH --ntasks=1
#SBATCH --mem=300g          # memory per cor
#SBATCH --job-name=test_SALS_BA4_regETM_Endo_epochs_sweep_200topics_narval_4
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/test_SALS_BA4_regETM_Endo_epochs_sweep_200topics_narval_4.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano test_SALS_BA4_regETM_Endo_epochs_sweep_200topics_narval_4.py


# -----------> Job ID: 43321809

import sys
import os

# Add the specific site-packages path to sys.path
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

import lime
import lime.lime_tabular

from scipy.stats import zscore
from sklearn.neighbors import KNeighborsClassifier

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
import lime
import lime.lime_tabular
from joblib import Parallel, delayed
from sklearn.metrics import balanced_accuracy_score


###################################
# Parameters
###################################
par_keep_cell_type = "Endo"
par_prep = "CombatSeq"
par_brain_region = "BA4"
par_status = "SALS"
remove = ["C9ALS", "C9FTLD", "SFTLD"]
#par_keep_cell_type = "Endo"
par_brain_region_Li = "frontal cortex"
remove_li = ['C9-FTD']

###################################
# Cell topic-epoch parameters
###################################
## read in the optimal file
info_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_optimal_CombatSeq_{par_status}_{par_brain_region}_narval_2.csv"
df_info = pd.read_csv(info_path)

## subset to only retain cell type of interest
df_info = df_info[df_info["celltype"] == par_keep_cell_type]

par_n_topics = df_info["n_genes"].iloc[0]
par_n_epochs = df_info["n_epochs"].iloc[0]


## Start here
print(par_keep_cell_type)     
print(par_n_epochs)
print(par_n_topics)
     

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


#######################################################
## Create LIME dataframe to be populated for all donors
#######################################################
## modify to concatonate to an empty list.
LIME_df = pd.DataFrame()
LIME_df[0] = []
LIME_df[1] = []
LIME_df['cell_index'] = []
LIME_df['test_donor'] = []
LIME_df['Disease_status'] = []

# Initialize df to an empty DataFrame before the loop
LIME_df = pd.DataFrame(columns=['feature', 'importance', 'cell_index', 'test_donor', 'Disease_status'])

#######################################################
## Create optimal topic x gene dataframe to be populated for all donors
#######################################################
topic_df = pd.DataFrame(columns=['variable', 'value', 'donor', 'cell_type'])

#######################################################
## Create lime feature importance dataframe to be populated for all donors
#######################################################
LIME_feature_df = pd.DataFrame(columns=['feature', 'abs_importance', 'zscore', 'donor', 'cell_type'])

###################################
# LOSO loop
###################################

sample_IDs = set(adata.obs['orig.ident'])
data = []

for donor in sample_IDs:
    
    #donor = '191112_ALS_118_snRNA-C4' ## TEMP
     
    ## Create a donor specific LIME frame
    LIME_df_donor = pd.DataFrame()
    LIME_df_donor[0] = []
    LIME_df_donor[1] = []
    LIME_df_donor['cell_index'] = []
    LIME_df_donor['test_donor'] = []
    LIME_df_donor['Disease_status'] = []
    
    # Initialize df to an empty DataFrame before the loop
    LIME_df_donor = pd.DataFrame(columns=['feature', 'importance', 'cell_index', 'test_donor', 'Disease_status'])
    
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
    disease_status_test = set(adata_test.obs['Group'])
    
    adata_train.X.shape
    adata_test.X.shape
    
    ## Set the necesarry column
    adata_train.obs["batch_indices"] = adata_train.obs["orig.ident"]
    adata_train.obs["assigned_cluster"] = adata_train.obs["Group"]
    
    ## scETM model 
    model = scETM(adata_train.n_vars, adata_train.obs.batch_indices.nunique(), n_topics = par_n_topics)
    
    ## train the model 
    trainer = UnsupervisedTrainer(model, adata_train, test_ratio=0.1)
    trainer.train(n_epochs = par_n_epochs, eval_every = 5, eval_kwargs = dict(cell_type_col = 'assigned_cluster'), save_model_ckpt = False) ## increase 12000
    
    ## Get cell embeddings
    model.get_all_embeddings_and_nll(adata_train)
    delta, alpha, rho = map(pd.DataFrame, [adata_train.obsm['delta'], adata_train.uns['alpha'], adata_train.varm['rho']])
    delta.index = adata_train.obs_names
    rho.index = adata_train.var_names
    delta.shape, alpha.shape, rho.shape
    
    ## Get top 200 genes per topic
    beta = rho @ alpha.T  # (gene, topic)
    top_words = pd.DataFrame(adata_train.var_names.values[np.argsort(beta.values, axis=0)[:-201:-1]])  # (n_top, topic)
    #top_words.to_csv(par_top_genes_topic)
    
    ## save AnnData object
    adata_train_ETM = sc.AnnData(X=adata_train.X.copy(),
                    obs=adata_train.obs.copy(),
                    var=adata_train.var.copy(),
                    uns=adata_train.uns.copy(),
                    obsm=adata_train.obsm.copy(),
                    varm=adata_train.varm.copy(),
                    layers=adata_train.layers.copy(),
                    raw=adata_train.raw.copy(),
                    dtype="float32",
                    shape=None,
                    obsp=adata_train.obsp.copy(),
                    varp=adata_train.varp
                    )
    adata_train_ETM.__dict__['_raw'].__dict__['_var'] = adata_train.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
    
    ## Apply topic model to test set ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Set the necesarry column
    adata_test.obs["batch_indices"] = adata_test.obs["orig.ident"]
    adata_test.obs["assigned_cluster"] = adata_test.obs["Group"]
    
    ## Get cell embeddings
    model.get_all_embeddings_and_nll(adata_test)
    delta, alpha, rho = map(pd.DataFrame, [adata_test.obsm['delta'], adata_test.uns['alpha'], adata_test.varm['rho']])
    delta.index = adata_test.obs_names
    rho.index = adata_test.var_names
    delta.shape, alpha.shape, rho.shape
    
    ## save AnnData object
    adata_test_ETM = sc.AnnData(X=adata_test.X.copy(),
                    obs=adata_test.obs.copy(),
                    var=adata_test.var.copy(),
                    uns=adata_test.uns.copy(),
                    obsm=adata_test.obsm.copy(),
                    varm=adata_test.varm.copy(),
                    layers=adata_test.layers.copy(),
                    raw=adata_test.raw.copy(),
                    dtype="float32",
                    shape=None,
                    obsp=adata_test.obsp.copy(),
                    varp=adata_test.varp
                    )
    adata_test_ETM.__dict__['_raw'].__dict__['_var'] = adata_test.__dict__['_raw'].__dict__['_var'].rename(
        columns={'_index': 'features'})
    
    ## Lets try to visualize the embeddings
    X_train = adata_train_ETM.obsm["delta"] 
    X_test = adata_test_ETM.obsm["delta"]
    
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
    num_genes = adata_train_ETM.obsm["delta"].shape[1]
        
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
            
    # LOSO loop
    learning_rate = 0.0008663099696291
    
    print(f"Processing donor: {donor}")
    num_cells = adata_test_ETM.shape[0]
    
    X_train = torch.FloatTensor(adata_train_ETM.obsm["delta"].toarray() if hasattr(adata_train_ETM.obsm["delta"], 'toarray') else adata_train_ETM.obsm["delta"])
    y_train = torch.LongTensor(adata_train_ETM.obs['Group'].values)
    domains = pd.factorize(adata_train_ETM.obs['Donor'])[0]
    d_train = torch.LongTensor(domains)
    
    dataset = TensorDataset(X_train, y_train, d_train)
    input_size = X_train.shape[1]
    num_classes = len(np.unique(y_train))
    num_domains = len(np.unique(domains))
    
    loader = DataLoader(dataset, batch_size=64, shuffle=True)
    
    model_main = MainModel(input_size, num_classes)
    model_domain = DomainClassifier(25, num_domains)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    train(model_main, model_domain, loader, 25, device)      ##### NUM EPOCH IS HERE
    
    # Evaluation
    X_test = torch.FloatTensor(adata_test_ETM.obsm["delta"].toarray() if hasattr(adata_test_ETM.obsm["delta"], 'toarray') else adata_test_ETM.obsm["delta"])
    y_test = adata_test_ETM.obs['Group']
    
    model_main.eval()
    with torch.no_grad():
        y_pred, _ = model_main(X_test)
        y_labels = torch.argmax(y_pred, dim=1)
        print(classification_report(y_test, y_labels.numpy()))
    
    acc = balanced_accuracy_score(y_test.values, y_labels.numpy())
    
    if acc == 0:
        continue  # Skip to the next donor in the loop
    
    #data.append({
    #    'prep': par_prep, 'donor': donor, 'region': par_brain_region,
    #    'group': par_status, 'celltype': par_keep_cell_type,
    #    'n_genes': num_genes, 'n_cells': num_cells, 'n_epochs': par_n_epochs,
    #    'learning_rate': learning_rate, 'batch_size': 64,
    #    'test_accuracy': acc, 'method': 'LIME'
    #})
    
    ###################################
    ## LIME implementation
    ###################################
    #Get the cell index
    df_index = pd.DataFrame(y_test)
    df_index['predicted_label'] = y_labels.numpy()
    df_index['cell_index'] = df_index.index
    
    ## correct indices
    correct_indices = df_index[df_index['Group'] == df_index['predicted_label']].index
        
    ## convert array to dataframe?
    df=pd.DataFrame(X_train,) 
    df_test=pd.DataFrame(X_test,)
    df_test.index = df_index['cell_index']
            
    ## subset df index to only include correct indices
    df_test = df_test.loc[correct_indices]
    num_cells = np.shape(df_test.values)[0]  
    print(num_cells)
            
    ## run LIME for each cell
    #explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=['PN', par_status], feature_names=adata_test_ETM.var_names.tolist(), verbose=True, mode='classification')
    explainer = lime.lime_tabular.LimeTabularExplainer(df.values, class_names=['PN', par_status], feature_names = np.arange(0, adata_test_ETM.obsm["delta"].shape[1]).tolist(), verbose=True, mode='classification')
    
    #feature_names = np.arange(0, adata_test_ETM.obsm["delta"].shape[1])
    
    for i in range(num_cells):
        # Convert the NumPy array to a PyTorch tensor
        input_data = torch.tensor(df_test.values[i], dtype=torch.float32).unsqueeze(0)  # Add a batch dimension
        
        # Convert back to NumPy array for LIME
        input_data_np = input_data.numpy()  # Convert to NumPy array
        
        # Define a predict_proba function that works with NumPy arrays
        def predict_proba(input_data):
            # Convert input data to a PyTorch tensor
            input_tensor = torch.tensor(input_data, dtype=torch.float32)
            # Get predictions from the model
            with torch.no_grad():
                logits = model_main(input_tensor)
                
                # Check if logits is a tuple and extract the first element if necessary
                if isinstance(logits, tuple):
                    logits = logits[0]  # Adjust this based on your model's output structure
                
                probabilities = torch.softmax(logits, dim=1)  # Apply softmax
            return probabilities.numpy()  # Return as NumPy array
        
        # Use the model to get predictions
        exp = explainer.explain_instance(input_data_np[0], predict_proba, num_features=num_genes)  # Pass only the first item
        
        # Get explanation as list
        dd = exp.as_list()
        
        # Create a DataFrame from the explanation
        df2 = pd.DataFrame(dd, columns=['feature', 'importance'])
        df2['cell_index'] = df_test.index[i]  # Use the correct index
        df2['test_donor'] = donor  # Use the correct index
        df2['Disease_status'] = str(list(disease_status_test)[0])  # Use the correct index
        
        # Concatenate with the main DataFrame
        LIME_df_donor = pd.concat([LIME_df_donor, df2], ignore_index=True)  # Add ignore_index=True to avoid index conflicts
    
    ## Ensure that we are only taking the donor of interest
    LIME_df_donor = LIME_df_donor[LIME_df_donor["test_donor"] == donor]
    
    ###################################
    ## Identify the most important ETM topics according to LIME
    ###################################
    ## place dataframe in usable format
    LIME_df_donor['test'] = LIME_df_donor['feature'].str.count(' ')
    
    mask_4 = LIME_df_donor['test'] == 4
    LIME_df_donor.loc[mask_4, 'feature'] = LIME_df_donor.loc[mask_4, 'feature'].str.replace(
        r'^.*? .*? (.*?) .*.*$', r'\1', regex=True
    )
    
    mask_2 = LIME_df_donor['test'] == 2
    LIME_df_donor.loc[mask_2, 'feature'] = LIME_df_donor.loc[mask_2, 'feature'].str.replace(
        r'^(.*?) .* .*$', r'\1', regex=True
    )
    
    ## Take the LIME ABS value. 
    LIME_df_donor['abs_importance'] = LIME_df_donor['importance'].abs()
    
    ## Compute the mean importance of features across all cells
    mean_abs_importance = LIME_df_donor.groupby('feature')['abs_importance'].mean().reset_index()
    
    mean_abs_importance_sorted = mean_abs_importance.sort_values('abs_importance', ascending=False)
    
    plt.figure(figsize=(12, 6))
    sns.barplot(data=mean_abs_importance_sorted, x='feature', y='abs_importance', palette='viridis')
    
    plt.xticks(rotation=90)  # Rotate x-axis labels if many features
    plt.xlabel('Feature')
    plt.ylabel('Mean Absolute Importance')
    plt.title('Mean Absolute LIME Importance per Feature')
    plt.tight_layout()
    plt.show()
    
    plt.savefig("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", format='pdf')
    
    ## Z-score option
    z_threshold = 0
    
    mean_abs_importance['zscore'] = zscore(mean_abs_importance['abs_importance'])
    
    mean_abs_importance_sorted = mean_abs_importance.sort_values('zscore', ascending=False)
    mean_abs_importance_sorted['feature'] = mean_abs_importance_sorted['feature'].astype(str)
    
    plt.figure(figsize=(14, 6))
    plt.bar(mean_abs_importance_sorted['feature'], mean_abs_importance_sorted['zscore'], color='skyblue')
    plt.axhline(z_threshold, color='red', linestyle='--', label=f'Z = {z_threshold}')
    
    plt.xlabel('Feature')
    plt.ylabel('Z-score of abs_importance')
    plt.title('Z-score of Mean Absolute Importance per Feature')
    plt.xticks(rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    plt.savefig("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf", format='pdf')
    
    ## Get topic number with Z > 0
    unique_features = mean_abs_importance_sorted['feature'].unique()
    positive_zscore_df = mean_abs_importance_sorted[mean_abs_importance_sorted['zscore'] > z_threshold]
    unique_features = positive_zscore_df['feature'].unique()
    len(unique_features)
    set(unique_features)
    
    ## append
    positive_zscore_df = positive_zscore_df.copy()
    positive_zscore_df['donor'] = donor
    positive_zscore_df['cell_type'] = par_keep_cell_type
    
    LIME_feature_df = pd.concat([LIME_feature_df, positive_zscore_df], ignore_index=True)  # Add ignore_index=True to avoid index conflicts
    
    ###################################
    ## Retrain and evaluate using optimal ETM topics only. 
    ###################################
    ## Lets try to visualize the embeddings
    selected_columns = sorted([int(i) for i in unique_features])
    
    ## Top genes for important topics
    top_words_optimal = top_words.iloc[:, selected_columns]
    top_words_long = top_words_optimal.melt(ignore_index=False)
    top_words_long["donor"] = donor
    topic_df = pd.concat([topic_df, top_words_long], ignore_index=True)  # Add ignore_index=True to avoid index conflicts
    
    X_train = adata_train_ETM.obsm["delta"] 
    X_train = X_train[:, selected_columns]
    
    X_test = adata_test_ETM.obsm["delta"]
    X_test = X_test[:, selected_columns]
    
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
     
    num_genes = X_train.shape[1]
            
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
            
    # LOSO loop
    learning_rate = 0.0008663099696291
    kNN_threshold_list = [.99, .95, .9, .85, .8]
    for kNN_threshold in kNN_threshold_list:
        for _ in range(5):
            print(f"Processing donor: {donor}")
            num_cells = adata_test_ETM.shape[0]
        
            X_train = torch.FloatTensor(adata_train_ETM.obsm["delta"].toarray() if hasattr(adata_train_ETM.obsm["delta"], 'toarray') else adata_train_ETM.obsm["delta"])
            X_train = X_train[:, selected_columns]
            X_train.shape
            
            y_train = torch.LongTensor(adata_train_ETM.obs['Group'].values)
            domains = pd.factorize(adata_train_ETM.obs['Donor'])[0]
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
            X_test = torch.FloatTensor(adata_test_ETM.obsm["delta"].toarray() if hasattr(adata_test_ETM.obsm["delta"], 'toarray') else adata_test_ETM.obsm["delta"])
            X_test = X_test[:, selected_columns]
            
            y_test = adata_test_ETM.obs['Group']
            
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
    del model, trainer, adata_train_ETM, adata_test_ETM
    del X_train, X_test, X_all, cell_metadata_train, cell_metadata_test, cell_metadata
    del model_main, model_domain


###################################
## Print stuff
###################################

## All of these must be printed. 
results_df = pd.DataFrame(data)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_{par_n_topics}_{par_n_epochs}_narval_2_KNN_Bayesian_LIME_optimal.csv"
results_df.to_csv(out_path, index=False)


topic_df_df = pd.DataFrame(topic_df)
topic_df_df["cell_type"] = par_keep_cell_type
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_topic_gene_info_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_{par_n_topics}_{par_n_epochs}_narval_2_KNN_Bayesian_LIME_optimal.csv"
topic_df_df.to_csv(out_path, index=False)

LIME_feature_df_df = pd.DataFrame(LIME_feature_df)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_LIME_feature_importance_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_{par_n_topics}_{par_n_epochs}_narval_2_KNN_Bayesian_LIME_optimal.csv"
LIME_feature_df_df.to_csv(out_path, index=False)


## Top genes for important topics
#top_words_optimal = top_words.iloc[:, selected_columns]
#top_words_long = top_words_optimal.melt(ignore_index=False)
#top_words_long["donor"] = donor
#top_words_long.head()

# Save results
#pd.set_option('display.max_rows', None)
#results_df = pd.DataFrame(data)
    
# Save results
#pd.set_option('display.max_rows', None)
#results_df = pd.DataFrame(data)

#out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_{par_n_topics}_{par_n_epochs}_narval_2.csv"
#results_df.to_csv(out_path, index=False)
        

