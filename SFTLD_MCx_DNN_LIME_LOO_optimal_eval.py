## run this in beluga

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

nano optimal_LOO_0.9995_weighted_LIME_SFTLD_BA4_all_cells_narval_2.sh

#!/bin/bash  
#SBATCH --account=def-grouleau
#SBATCH --time=07-00:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=150g          # memory per cor
#SBATCH --job-name=optimal_LOO_0.9995_weighted_LIME_SFTLD_BA4_all_cells_narval_2
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/optimal_LOO_0.9995_weighted_LIME_SFTLD_BA4_all_cells_narval_2.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano optimal_LOO_0.9995_weighted_LIME_SFTLD_BA4_all_cells_narval_2.py



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


################################################################## 0.9995

# Set parameters
#par_keep_cell_type = "L2_L3"
par_prep = "CombatSeq"
par_brain_region = "BA4"
par_status = "SFTLD"
remove = ["SALS", "C9FTLD", "C9ALS"]

cell_types = [
    "L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Oligo", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5"
]

# Loop through the list and print each item
for par_keep_cell_type in cell_types:
    print(par_keep_cell_type)
    
    # Load LIME-selected genes
    lime_file = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_.9995_abs_case_control_narval_2.csv'
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
    mapping = {'SFTLD': 1, 'PN': 0}
    adata.obs['Group'] = adata.obs['Group'].map(mapping)
    
    # Preprocess data
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[:, features]                                       #### TRYING WITH THIS HERE, we will see how it affects the performance. 
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    adata.raw = adata
    
    # Subset to optimal features
    #                                                            adata = adata[:, features]
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
    out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
    results_df.to_csv(out_path, index=False)



################################################################## 0

# Set parameters
#par_keep_cell_type = "L2_L3"
par_prep = "CombatSeq"
par_brain_region = "BA4"
par_status = "SFTLD"
remove = ["SALS", "C9FTLD", "C9ALS"]

cell_types = [
    "L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Oligo", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5"
]

# Loop through the list and print each item
for par_keep_cell_type in cell_types:
    print(par_keep_cell_type)
    
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
    mapping = {'SFTLD': 1, 'PN': 0}
    adata.obs['Group'] = adata.obs['Group'].map(mapping)
    
    # Preprocess data
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[:, features]                                       #### TRYING WITH THIS HERE, we will see how it affects the performance. 
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    adata.raw = adata
    
    # Subset to optimal features
    #                                                            adata = adata[:, features]
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
    out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
    results_df.to_csv(out_path, index=False)

















































################################################################################################################ OLD not clean code

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
import lime
import lime.lime_tabular
from joblib import Parallel, delayed
from sklearn.metrics import balanced_accuracy_score

###################################
# scanpy settings
###################################
sc.settings.verbosity = 3             
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')


cell_types = [
    "L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Oligo", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5"
]

#par_keep_cell_type = "L2_L3"

# Loop through the list and print each item
for par_keep_cell_type in cell_types:
    print(par_keep_cell_type)
     
    
par_keep_cell_type = "L2_L3"
    
###################################
# Parameters
###################################
## Here we are keeping SFTLD and PN
par_prep = "CombatSeq"
par_brain_region = "BA4"
par_status = "SFTLD" ### IMPLEMENT THIS
remove = ["SALS", "C9FTLD", "C9ALS"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 
    
###################################
# Read in optimal gene set file. 
###################################
# File path for the CSV (adjust parameters accordingly)
file_path = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_.9995_narval_2.csv' ## Produced from LIME_outputs_explore3.R
    
#'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_.9995_narval_2.csv'
    
# Step 1: Read the CSV file
df = pd.read_csv(file_path)
features = df["gene"].tolist()
    
#################################
## Read in the previous report
#################################
file_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
model_report_df = pd.read_csv(file_path)
    
batch_size = int(model_report_df.at[0, 'batch_size'])
learning_rate = model_report_df.at[0, 'learning_rate']
    
###################################
# Cell Specific parameters -- L3_L5
###################################
data = []
par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
print(par_ann_data)
    
###################################
# Load anndata
###################################
adata = sc.read_h5ad(par_ann_data)
#num_cells = adata.X.shape[0] --> We are doing this at the test level now. 
adata.obs['Group']
adata = adata[~adata.obs['Group'].isin(remove)]
adata = adata[adata.obs['CellType'] == par_keep_cell_type]
adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
print("genes: ", adata.var_names) 
print("cells: ", adata.obs_names) 
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]

## create a list of sample IDs
sample_ID = set(adata.obs['orig.ident'])
    
###################################
# clean AnnData object
###################################
adata.obs = adata.obs.reset_index() 
adata.obs.columns
set(adata.obs['Group'])
set(adata.obs['CellType'])
set(adata.obs['Region'])
    
## change Disease status labels
# Create a mapping dictionary
mapping = {'SFTLD': 1, 'PN': 0}
    
# Map the values
adata.obs['Group'] = adata.obs['Group'].map(mapping)
    
# Verify the change
print(adata.obs['Group'].value_counts())
    
#######################################################
## Perform the test train split
#######################################################
# Create an array of indices
#indices = np.arange(adata.shape[0])
# Split indices into train and test
#train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=42)
# Create training and testing AnnData objects
#adata_train = adata[train_indices].copy()
#adata_test = adata[test_indices].copy()
        
#######################################################
## Preprocess the train and test objects bring this before the 
#######################################################

## Process the whole dataset together ------------------------------------------------ MIGHT BE USEFUL TO NORMALIZE COUNTS AFTER ONLY RETAINING THE GENES OF INTEREST?
sc.pp.filter_cells(adata, min_genes=200)  # Filter cells
sc.pp.filter_genes(adata, min_cells=3)    # Filter genes
sc.pp.normalize_total(adata, target_sum=1e4)  # Normalize
sc.pp.log1p(adata)                         # Log-transform
sc.pp.highly_variable_genes(adata) 
adata.raw = adata  # Save raw counts for further analysis
    
## Subset to only include the features of interest
adata_lim = adata[:, features]  # Subset to optimal
num_genes = adata_lim.X.shape[1]
    
#######################################################
## With optimal values
#######################################################  
    
for donor in sample_ID: 
    iter = 1
    while iter <= 3:
        print(f"Processing data for donor: {donor}")
        
        #######################################################
        ## Set leave one out test-train dataset
        #######################################################
        ## train object
        adata_train = adata_lim[adata_lim.obs['orig.ident'] != donor]
        set(adata_train.obs['orig.ident'])
        len(set(adata_train.obs['orig.ident']))
        set(adata_train.obs['Group'])
        adata_train.shape
        
        ## test object 
        adata_test = adata_lim[adata_lim.obs['orig.ident'] == donor]
        set(adata_test.obs['orig.ident'])
        len(set(adata_test.obs['orig.ident']))
        set(adata_test.obs['Group'])
        disease_status_test = set(adata_test.obs['Group'])
        adata_test.shape
        num_cells = adata_test.X.shape[0]
        
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
            
            def predict_proba(self, x):
                self.eval()
                with torch.no_grad():
                    logits = self(x)
                    return F.softmax(logits, dim=1)
                    
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
            
        # Function to train and evaluate the models
        torch.autograd.set_detect_anomaly(True)
            
        #def train_and_evaluate(learning_rate, weight_decay, batch_size):
        def train_and_evaluate(learning_rate, batch_size):
            # Split the dataset into training and validation sets (80-20 split)
            train_indices, val_indices = train_test_split(
                np.arange(len(full_dataset)), 
                test_size=0.2, 
                random_state=42
            )
            
            # Create DataLoaders for training and validation sets
            train_loader = DataLoader(full_dataset, batch_size=int(batch_size), sampler=torch.utils.data.SubsetRandomSampler(train_indices))
            val_subset = DataLoader(full_dataset, batch_size=int(batch_size), sampler=torch.utils.data.SubsetRandomSampler(val_indices))
            
            
            # Initialize models
            main_model = MainModel(input_size, num_classes)
            domain_model = DomainClassifier(25, num_domains)
            
            # Optimizers
            #main_optimizer = optim.Adam(main_model.parameters(), lr=float(learning_rate), weight_decay=float(weight_decay))
            main_optimizer = optim.Adam(main_model.parameters(), lr=float(learning_rate))
            domain_optimizer = optim.Adam(domain_model.parameters(), lr=0.001)
            
            # Loss functions
            class_loss_fn = nn.CrossEntropyLoss()
            domain_loss_fn = nn.CrossEntropyLoss()
            
            # Train the model
            for epoch in range(25):  # Fixed epochs for simplicity
                main_model.train()
                domain_model.train()
                
                for X, y, domain_labels in train_loader:
                    X, y, domain_labels = X.to(device), y.to(device), domain_labels.to(device)
                
                    class_output, shared_rep = main_model(X)
                    class_loss = class_loss_fn(class_output, y)
                
                    domain_output = domain_model(shared_rep.detach())
                    domain_loss = domain_loss_fn(domain_output, domain_labels)
                
                    total_loss = class_loss - domain_loss
                        
                    # Update the main model
                    main_optimizer.zero_grad()
                    total_loss.backward(retain_graph=True)  # Do not retain graph here unless necessary
                    main_optimizer.step()
            
                    # Update the domain model
                    domain_optimizer.zero_grad()
                    domain_loss.backward()  # No need to retain graph here
                    domain_optimizer.step()
                
            # Evaluation
            val_predictions = []
            val_labels = []
                
            with torch.no_grad():
                for X_val, y_val, domain_labels in val_subset:  # Use the last fold's validation set
                    X_val, y_val, domain_labels = X_val.to(device), y_val.to(device), domain_labels.to(device)
                    y_pred_val, _ = main_model(X_val)
                    val_predictions.append(torch.argmax(y_pred_val, dim=1).cpu())
                    val_labels.append(y_val.cpu())
                
            val_predictions = torch.cat(val_predictions).numpy()
            val_labels = torch.cat(val_labels).numpy()
                
            # Calculate accuracy
            balanced_accuracy = balanced_accuracy_score(val_labels, val_predictions)
            return balanced_accuracy
        
        #################################
        ## Bayesian approach to optimize parameters 
        #################################
        ########## Prepare your data as a TensorDataset ##########
        ## X
        X_train = adata_train.X
        X_train_dense = X_train.toarray() if hasattr(X_train, 'toarray') else X_train
        X_train_tensor = torch.FloatTensor(X_train_dense)  # Gene expression matrix
        
        ## Y
        y_train = adata_train.obs['Group']
        y_train_tensor = torch.LongTensor(y_train)
        
        ## domain
        unique_samples = adata_train.obs['Donor'].unique()
        num_domains = len(unique_samples)
        # Create a mapping from sample IDs to domain labels
        sample_to_domain = {sample: idx for idx, sample in enumerate(unique_samples)}
        # Create domain labels based on Sample_ID
        adata_train.obs['Domain_Label'] = adata_train.obs['Donor'].map(sample_to_domain)
        y_train_indices = y_train.index
        adata_train_subset = adata_train[y_train_indices]
        domain_labels_tensor = torch.LongTensor(adata_train_subset.obs['Domain_Label'])
        
        ## full dataset tensor
        X_train_tensor.shape
        y_train_tensor.shape
        domain_labels_tensor.shape
        full_dataset = TensorDataset(X_train_tensor, y_train_tensor, domain_labels_tensor)
        num_classes = len(np.unique(adata_train.obs['Group']))  # Number of unique disease classes
        num_domains = len(np.unique(domain_labels_tensor.numpy()))  # Number of unique sample IDs or batches
        input_size = adata_train.shape[1]  # Number of highly variable genes
        
        ###################################
        ## define train function
        ###################################
                
        def train(main_model, domain_model, dataloader, num_epochs, device):
            main_model.to(device)
            domain_model.to(device)
            
            # Optimizers
            #main_optimizer = optim.Adam(main_model.parameters(), lr=0.0001, weight_decay=0.01)
            main_optimizer = optim.Adam(main_model.parameters(), lr=learning_rate)
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
        ## Train on the entire train set
        ###################################
        # Sample usage
        # Define input size and number of classes
        input_size = adata_train.shape[1]  # Number of highly variable genes
        num_classes = len(np.unique(adata_train.obs['Group']))  # Number of unique disease classes
        num_domains = len(np.unique(domain_labels_tensor.numpy()))  # Number of unique sample IDs or batches
        
        # Initialize models
        main_model = MainModel(input_size, num_classes)
        domain_model = DomainClassifier(25, num_domains)  # 50 is the output size from shared layers
        
        # Create dataset and DataLoader
        train_dataset = TensorDataset(X_train_tensor, y_train_tensor, domain_labels_tensor)
        train_dataloader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        
        # Train the models
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        train(main_model, domain_model, train_dataloader, num_epochs=10, device=device)                                                                       ################################## CHANGE NUM EPOCHS HERE
        
        ###################################
        ## Evaluate on the left out sample
        ###################################
        X_test = adata_test.X
        X_test_dense = X_test.toarray() if hasattr(X_test, 'toarray') else X_test
        X_test_tensor = torch.FloatTensor(X_test_dense)
        
        y_test = adata_test.obs['Group']
        
        main_model.eval()
        with torch.no_grad():
            y_pred_test, _ = main_model(X_test_tensor)
            y_pred_test_labels = torch.argmax(y_pred_test, dim=1)
            probabilities = F.softmax(y_pred_test, dim=1)
            print("Test Classification Report:\n", classification_report(y_test, y_pred_test_labels.numpy()))
        
        # Convert predictions and true labels to numpy arrays
        y_pred_test_labels_np = y_pred_test_labels.numpy()
        y_test_np = y_test.values  # Ensure y_test is also a numpy array
                
        # Calculate accuracy
        test_balanced_accuracy = balanced_accuracy_score(y_test_np, y_pred_test_labels_np)
                
        ###################################
        ## Populate the dataframe
        ###################################
                
        # Append a dictionary with celltype and accuracy
        data.append({'prep': par_prep,
        'donor': donor,
        'region': par_brain_region,
        'group': par_status,
        'celltype': par_keep_cell_type, 
        'n_genes': num_genes,
        'n_cells': num_cells,
        'learning_rate': learning_rate,
        'batch_size': batch_size,
        'test_accuracy': test_balanced_accuracy,
        'method': 'LIME'
        })
        iter += 1
    
# Create a DataFrame from the list
df = pd.DataFrame(data)
    
## write csv file
file_name = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
df.to_csv(file_name, index=False)




