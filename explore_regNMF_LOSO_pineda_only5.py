salloc -A def-grouleau --time=0-8 -c 1 --mem=200g

module load StdEnv/2020 
module load python/3.8.10

python  

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L3_L5 

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

nano test_C9FTLD_BA9_regNMF_narval_2.sh

#!/bin/bash  
#SBATCH --account=def-grouleau
#SBATCH --time=04-00:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=150g          # memory per cor
#SBATCH --job-name=test_C9FTLD_BA9_regNMF_narval_2
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/test_C9FTLD_BA9_regNMF_narval_2.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano test_C9FTLD_BA9_regNMF_narval_2.py


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

###################################
# Parameters
###################################
par_prep = "CombatSeq"
par_brain_region = "BA9"
par_status = "C9FTLD"
remove = ["C9ALS", "SFTLD", "SALS"]
#par_keep_cell_type = "L2_L3"
par_brain_region_Li = "frontal cortex"
remove_li = ['C9-ALS']

cell_types = [
    "L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Oligo", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5"
]

for par_keep_cell_type in cell_types:
    print(par_keep_cell_type)
    
    ###################################
    # Cell Specific parameters -- Astro
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
    mapping = {'C9FTLD': 1, 'PN': 0}
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
    
    sc.pp.highly_variable_genes(adata, n_top_genes=10000, flavor='seurat')
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
    # Normalize and Log Transform
    ###################################
    #sc.pp.normalize_total(adata, target_sum=1e4)
    #sc.pp.log1p(adata)
    
    ###################################
    # Supervised Regularized NMF (PyTorch)
    ###################################
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    
    class SupervisedNMF(nn.Module):
        def __init__(self, n_cells, n_genes, n_components, lambda_subj=0.1, lambda_disease=1.0):
            super().__init__()
            self.W = nn.Parameter(torch.rand(n_cells, n_components, device=device))
            self.H = nn.Parameter(torch.rand(n_components, n_genes, device=device))
            self.lambda_subj = lambda_subj
            self.lambda_disease = lambda_disease
        
        def forward(self, X, subject_onehot, disease_labels):
            reconstruction = self.W @ self.H
            recon_loss = F.mse_loss(reconstruction, X)
        
            W_norm = self.W - self.W.mean(dim=0, keepdim=True)
            subj_norm = subject_onehot - subject_onehot.mean(dim=0, keepdim=True)
            disease_norm = disease_labels - disease_labels.mean()
        
            # Subject correlation penalty
            subj_corr = torch.matmul(W_norm.T, subj_norm) / (
                torch.norm(W_norm, dim=0).unsqueeze(1) * torch.norm(subj_norm, dim=0)
            )
            subject_penalty = torch.mean(torch.abs(subj_corr))
        
            # Disease correlation reward
            disease_corrs = torch.matmul(W_norm.T, disease_norm.unsqueeze(1)).squeeze() / (
                torch.norm(W_norm, dim=0) * torch.norm(disease_norm)
            )
            disease_reward = torch.mean(torch.abs(disease_corrs))
        
            loss = recon_loss + self.lambda_subj * subject_penalty - self.lambda_disease * disease_reward
            return loss, reconstruction
    
    def plot_latent_components_to_pdf(W_tensor, labels, num_components=10, output_path="latent_components.pdf"):
        W_np = W_tensor.detach().cpu().numpy()
        labels_np = labels.detach().cpu().numpy()
        
        with PdfPages(output_path) as pdf:
            for i in range(min(num_components, W_np.shape[1])):
                plt.figure(figsize=(6, 4))
                sns.violinplot(x=labels_np, y=W_np[:, i])
                plt.title(f"Component {i+1} by Disease Status")
                plt.xlabel("Disease Status")
                plt.ylabel(f"Component {i+1} Activation")
                plt.tight_layout()
                pdf.savefig()  # Save the current figure into the PDF
                plt.close()    # Close the figure to avoid displaying it
        
        print(f"Plots saved to {output_path}")
    
    ###################################
    # LOSO loop 
    ###################################
    
    sample_IDs = set(adata.obs['orig.ident'])
    data = []
    
    for donor in sample_IDs:
        # Split adata object
        print(f"Processing donor: {donor}")
        adata_train = adata[adata.obs['orig.ident'] != donor]
        adata_test = adata[adata.obs['orig.ident'] == donor]
        
        # Prepare tensors
        X = torch.tensor(adata_train.X.toarray(), dtype=torch.float32, device=device)
        subjects = adata_train.obs['orig.ident'].astype('category').cat.codes.values
        encoder = OneHotEncoder(sparse_output=False)
        subject_onehot = torch.tensor(encoder.fit_transform(subjects.reshape(-1, 1)), dtype=torch.float32, device=device)
        disease_labels = torch.tensor(adata_train.obs['Group'].values, dtype=torch.float32, device=device)
        
        # Train model
        model = SupervisedNMF(n_cells=X.shape[0], n_genes=X.shape[1], n_components=10, lambda_subj=1.0, lambda_disease=1.0).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)
        
        # Early stopping parameters
        patience = 10
        best_loss = float('inf')
        epochs_without_improvement = 0
        
        for epoch in range(200):
            optimizer.zero_grad()
            loss, reconstruction = model(X, subject_onehot, disease_labels)
            loss.backward()
            optimizer.step()
            with torch.no_grad():
                model.W.clamp_(min=0)
                model.H.clamp_(min=0)
            
            current_loss = loss.item()
            
            # Check for improvement
            if current_loss < best_loss - 1e-4:  # small threshold to avoid noise
                best_loss = current_loss
                epochs_without_improvement = 0
            else:
                epochs_without_improvement += 1
            
            # Print loss every epoch
            print(f"Epoch {epoch}: Loss = {current_loss:.4f}")
            
            # Early stopping check
            if epochs_without_improvement >= patience:
                print(f"Stopping early at epoch {epoch} (best loss = {best_loss:.4f})")
                break
        
        #✅ What to look for:
        #Is the loss decreasing consistently over epochs? That's a good sign of optimization progress.
        #Does the model generalize well? If you're holding out subjects (as in LOSO), see if W captures disease signal in the held-out subject.
        #Does the final reconstruction look good? You could check MSE separately.
        
        X_test = torch.tensor(adata_test.X.toarray(), dtype=torch.float32, device=device)
        W_test = torch.rand(X_test.shape[0], model.H.shape[0], device=device, requires_grad=True)
        
        optimizer_test = torch.optim.Adam([W_test], lr=1e-2)
        
        ###################################
        # Fit model to unseen subject ---------- NEED TO ADD THE STOP FUNCTION. 
        ###################################
        
        for epoch in range(2000):
            optimizer_test.zero_grad()
            reconstruction = W_test @ model.H  # Only H from trained model is used
            recon_loss = F.mse_loss(reconstruction, X_test)
            recon_loss.backward()
            optimizer_test.step()
            with torch.no_grad():
                W_test.clamp_(min=0)
            
            print(f"Test W fitting epoch {epoch}, loss = {recon_loss.item():.4f}")
        
        ###################################
        # Add to W components to the adata objects
        ###################################
        ## Train object
        W = model.W
        H = model.H
        
        plot_latent_components_to_pdf(model.W, disease_labels, num_components=10, output_path='/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf')
        
        W_np = model.W.detach().cpu().numpy()
        new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])
        
        adata_train_NMF = sc.AnnData(X = W_np.copy(),
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
        
        adata_train_NMF.__dict__['_raw'].__dict__['_var'] = adata_train_NMF.__dict__['_raw'].__dict__['_var'].rename(
            columns={'_index': 'features'})
        
        ## Test object
        W_test_np = W_test.detach().cpu().numpy()
        W_test_np.shape
        
        adata_test_NMF = sc.AnnData(X = W_test_np.copy(),
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
        
        adata_test_NMF.__dict__['_raw'].__dict__['_var'] = adata_test_NMF.__dict__['_raw'].__dict__['_var'].rename(
            columns={'_index': 'features'})
        
        ###################################
        # Implementation into the LOSO DNN
        ###################################
        
        num_genes = adata_test_NMF.X.shape[1]
            
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
        learning_rate = 0.0008663099696291
        
        for _ in range(5):
            print(f"Processing donor: {donor}")
            num_cells = adata_test_NMF.shape[0]
            
            X_train = torch.FloatTensor(adata_train_NMF.X.toarray() if hasattr(adata_train_NMF.X, 'toarray') else adata_train_NMF.X)
            y_train = torch.LongTensor(adata_train_NMF.obs['Group'].values)
            domains = pd.factorize(adata_train_NMF.obs['Donor'])[0]
            d_train = torch.LongTensor(domains)
            
            dataset = TensorDataset(X_train, y_train, d_train)
            input_size = adata_train_NMF.shape[1]
            num_classes = len(np.unique(y_train))
            num_domains = len(np.unique(domains))
            
            loader = DataLoader(dataset, batch_size=64, shuffle=True)
            
            model_main = MainModel(input_size, num_classes)
            model_domain = DomainClassifier(25, num_domains)
            device = 'cuda' if torch.cuda.is_available() else 'cpu'
            
            train(model_main, model_domain, loader, 10, device)      ##### NUM EPOCH IS HERE
            
            # Evaluation
            X_test = torch.FloatTensor(adata_test_NMF.X.toarray() if hasattr(adata_test_NMF.X, 'toarray') else adata_test_NMF.X)
            y_test = adata_test_NMF.obs['Group']
            
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
                'learning_rate': learning_rate, 'batch_size': 64,
                'test_accuracy': acc, 'method': 'LIME'
            })
    
    # Save results
    pd.set_option('display.max_rows', None)
    results_df = pd.DataFrame(data)
    
    out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/regNMF_model_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
    results_df.to_csv(out_path, index=False)






/home/fiorini9/scratch/machine_learning_ALS/model_outs/regNMF_model_report_CombatSeq_C9FTLD_BA9_L3_L5_narval_2.csv