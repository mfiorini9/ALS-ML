salloc -A def-grouleau --time=0-8 -c 1 --mem=200g

module load StdEnv/2020 
module load python/3.8.10

python  

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L3_L5   ##### JOB ID: 43265174

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

nano test_SFTLD_BA9_regNMF10_narval_2.sh

#!/bin/bash  
#SBATCH --account=def-grouleau
#SBATCH --time=01-00:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=150g          # memory per cor
#SBATCH --job-name=test_SFTLD_BA9_regNMF10_narval_2
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/test_SFTLD_BA9_regNMF10_narval_2.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano test_SFTLD_BA9_regNMF10_narval_2.py


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

#/home/fiorini9/scratch/machine_learning_ALS/model_outs/regNMF_model_report10_CombatSeq_SFTLD_BA9_L2_L3_narval_2.csv

###################################
# Parameters
###################################
par_prep = "CombatSeq"
par_brain_region = "BA9"
par_status = "SFTLD"
remove = ["C9ALS", "C9FTLD", "SALS"]
#par_keep_cell_type = "L2_L3"
par_brain_region_Li = "frontal cortex"
remove_li = ['C9-ALS']

cell_types = [ "L2_L3"]

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
    mapping = {'SFTLD': 1, 'PN': 0}
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
    # Supervised Regularized NMF with Domain Adaptation (PyTorch)
    ###################################
     
    # --- Device
     
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
     
    # --- Gradient Reversal Layer ---
    class GradReverse(Function):
        @staticmethod
        def forward(ctx, x, lambd):
            ctx.lambd = lambd
            return x.view_as(x)
     
        @staticmethod
        def backward(ctx, grad_output):
            return grad_output.neg() * ctx.lambd, None
     
    def grad_reverse(x, lambd=1.0):
        return GradReverse.apply(x, lambd)
     
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
                pdf.savefig()
                plt.close()
     
        print(f"Plots saved to {output_path}")
     
    class SupervisedNMF(nn.Module):
        def __init__(self, n_cells, n_genes, n_components, n_domains, lambda_subj=1, lambda_disease=20, lambda_domain=1):
            super().__init__()
            self.W = nn.Parameter(torch.rand(n_cells, n_components, device=device))
            self.H = nn.Parameter(torch.rand(n_components, n_genes, device=device))
            self.lambda_subj = lambda_subj
            self.lambda_disease = lambda_disease
            self.lambda_domain = lambda_domain
            self.n_domains = n_domains
     
            self.disease_predictor = nn.Linear(n_components, 1)
            self.domain_classifier = nn.Linear(n_components, n_domains)  # <-- output n_domains logits now
            self.bce_loss = nn.BCEWithLogitsLoss()
            self.ce_loss = nn.CrossEntropyLoss()  # <-- add CrossEntropyLoss for multi-domain
     
        def forward(self, X, subject_onehot, disease_labels):
            # Small Gaussian noise
            noise = 0.01 * torch.randn_like(X)
            X_noisy = X + noise
     
            reconstruction = self.W @ self.H
            recon_loss = F.mse_loss(reconstruction, X_noisy)
     
            # Subject decorrelation
            W_norm = self.W - self.W.mean(dim=0, keepdim=True)
            subj_norm = subject_onehot - subject_onehot.mean(dim=0, keepdim=True)
            subj_corr = torch.matmul(W_norm.T, subj_norm) / (
                torch.norm(W_norm, dim=0).unsqueeze(1) * torch.norm(subj_norm, dim=0)
            )
            subject_penalty = torch.mean(torch.abs(subj_corr))
     
            # Disease classification loss
            disease_logits = self.disease_predictor(self.W).squeeze()
            disease_classification_loss = self.bce_loss(disease_logits, disease_labels)
     
            # Domain adversarial loss
            rev_features = grad_reverse(self.W, lambd=1.0)
            domain_logits = self.domain_classifier(rev_features)
            domain_labels = subject_onehot.argmax(dim=1)  # <-- no modulo 2
            domain_classification_loss = self.ce_loss(domain_logits, domain_labels)  # <-- use CrossEntropyLoss
     
            # H penalty
            H_penalty = 1e-2 * torch.norm(self.H, p=2)
     
            # Total loss
            loss = (
                recon_loss
                + self.lambda_subj * subject_penalty
                + self.lambda_disease * disease_classification_loss
                + self.lambda_domain * domain_classification_loss
                + H_penalty
            )
            return loss, reconstruction, disease_logits
     
    ###################################
    # LOSO loop
    ###################################
     
    sample_IDs = [
        '210526_FTLD_225_snRNA-C6', '210526_PN_319_snRNA-E5', '210526_PN_303_snRNA-D9',
        '210526_FTLD_218_snRNA-C3', '210526_FTLD_226_snRNA-C7', '210526_FTLD_206_snRNA-B7',
        '210526_FTLD_204_snRNA-B6', '210526_FTLD_209_snRNA-B9', '210526_FTLD_228_snRNA-C8',
        '210526_FTLD_202_snRNA-B5', '210526_FTLD_217_snRNA-C2', '210526_FTLD_222_snRNA-C4',
        '210526_FTLD_207_snRNA-B8', '210526_FTLD_216_snRNA-C1'
    ]
    sample_IDs = set(adata.obs['orig.ident'])
    data = []
     
    for donor in sample_IDs:
        print(f"Processing donor: {donor}")
        adata_train = adata[adata.obs['orig.ident'] != donor]
        adata_test = adata[adata.obs['orig.ident'] == donor]
     
        X = torch.tensor(adata_train.X.toarray(), dtype=torch.float32, device=device)
        subjects = adata_train.obs['orig.ident'].astype('category').cat.codes.values
        encoder = OneHotEncoder(sparse_output=False)
        subject_onehot = torch.tensor(encoder.fit_transform(subjects.reshape(-1, 1)), dtype=torch.float32, device=device)
        disease_labels = torch.tensor(adata_train.obs['Group'].values, dtype=torch.float32, device=device)
     
        model = SupervisedNMF(
            n_cells=X.shape[0], 
            n_genes=X.shape[1], 
            n_components=75, 
            n_domains=subject_onehot.shape[1],  # NEW: number of domains = one-hot vector size
            lambda_subj=1, 
            lambda_disease=20, 
            lambda_domain=1
        ).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)
     
        patience = 10
        best_loss = float('inf')
        epochs_without_improvement = 0
     
        for epoch in range(100):
            optimizer.zero_grad()
            loss, reconstruction, disease_logits = model(X, subject_onehot, disease_labels)
            loss.backward()
            optimizer.step()
            with torch.no_grad():
                model.W.clamp_(min=0)
                model.H.clamp_(min=0)
     
            current_loss = loss.item()
            if current_loss < best_loss - 1e-4:
                best_loss = current_loss
                epochs_without_improvement = 0
            else:
                epochs_without_improvement += 1
     
            print(f"Epoch {epoch}: Loss = {current_loss:.4f}")
     
            if epochs_without_improvement >= patience:
                print(f"Stopping early at epoch {epoch} (best loss = {best_loss:.4f})")
                break
           
        # Compute subject prototypes
        #all_subject_ids = np.unique(subjects)
        W_train = model.W.detach()
        #train_subjects = torch.tensor(subjects, dtype=torch.long, device=device)
        #subject_prototypes = torch.stack([W_train[train_subjects == sid].mean(0) for sid in all_subject_ids])
     
        # --- Solve for W_test using NNLS ---
        H_np = model.H.detach().cpu().numpy()
        X_test_np = adata_test.X.toarray()
     
        W_test_list = []
        for i in range(X_test_np.shape[0]):
            x_i = X_test_np[i]
            w_i, _ = scipy.optimize.nnls(H_np.T, x_i)
            W_test_list.append(w_i)
     
        W_test_np = np.vstack(W_test_list)
        #W_test = torch.tensor(W_test_np, dtype=torch.float32, device=device)
       
        # --- Procrustes Alignment ---
        print("Running Procrustes alignment...")
       
        # Important: only apply to the feature space (not metadata
        W_train = model.W.detach().cpu().numpy()
        W_test_np = np.vstack(W_test_list)  # <-- still numpy!
        
        # Align
        W_train_mean = W_train.mean(axis=0)
        W_train_std = W_train.std(axis=0)
        
        W_test_mean = W_test_np.mean(axis=0)
        W_test_std = W_test_np.std(axis=0)
        
        # Align mean and std
        W_test_np_aligned = (W_test_np - W_test_mean) / (W_test_std + 1e-6)
        W_test_np_aligned = W_test_np_aligned * (W_train_std + 1e-6) + W_train_mean
        
        # Now *new tensor* for W_test
        W_test = torch.tensor(W_test_np_aligned, dtype=torch.float32, device=device)
       
        # --- Save latent components into new AnnData objects ---
        W_np = model.W.detach().cpu().numpy()
        new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W_np.shape[1])])
       
        adata_train_NMF = sc.AnnData(
            X=W_np.copy(),
            obs=adata_train.obs.copy(),
            var=new_var,
            uns=adata_train.uns.copy(),
            obsm=adata_train.obsm.copy(),
            varm=adata_train.varm.copy(),
            layers=adata_train.layers.copy(),
            raw=adata_train.raw.copy(),
            dtype="float32",
            obsp=adata_train.obsp.copy(),
            varp=adata_train.varp.copy(),
        )
        adata_train_NMF.raw.var.rename(columns={'_index': 'features'}, inplace=True)
       
        adata_test_NMF = sc.AnnData(
            X=W_test.detach().cpu().numpy().copy(),
            obs=adata_test.obs.copy(),
            var=new_var,
            uns=adata_test.uns.copy(),
            obsm=adata_test.obsm.copy(),
            varm=adata_test.varm.copy(),
            layers=adata_test.layers.copy(),
            raw=adata_test.raw.copy(),
            dtype="float32",
            obsp=adata_test.obsp.copy(),
            varp=adata_test.varp.copy(),
        )
        adata_test_NMF.raw.var.rename(columns={'_index': 'features'}, inplace=True)
       
        # --- Plot components violin (optional) ---
        plot_latent_components_to_pdf(
            model.W, disease_labels, num_components=10,
            output_path='/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp_violin.pdf'
        )
       
        # --- Collect all cell embeddings and metadata ---
        W_train = torch.tensor(W_train, dtype=torch.float32, device=W_test.device)
        W_all = torch.cat([W_train, W_test], dim=0).detach().cpu().numpy()
       
        cell_metadata_train = adata_train.obs[['orig.ident', 'Group']].copy()
        cell_metadata_test = adata_test.obs[['orig.ident', 'Group']].copy()
       
        # Map Group numbers to names
        group_mapping = {0: 'PN', 1: 'SFTLD'}
        cell_metadata_train['Group'] = cell_metadata_train['Group'].map(group_mapping)
        cell_metadata_test['Group'] = cell_metadata_test['Group'].map(group_mapping)
       
        cell_metadata_train['Set'] = 'Train'
        cell_metadata_test['Set'] = 'Test'
       
        cell_metadata = pd.concat([cell_metadata_train, cell_metadata_test], axis=0).reset_index(drop=True)
       
        # Define a new label: "PlotGroup"
        cell_metadata['PlotGroup'] = np.where(
            cell_metadata['Set'] == 'Test', 'Test', cell_metadata['Group']
        )
       
        # --- Open a PDF file to save PCA and UMAP plots ---
        with PdfPages('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/nmf_latent_space_visualizations.pdf') as pdf:
            # PCA plot
            pca = PCA(n_components=2)
            W_pca = pca.fit_transform(W_all)
            
            plt.figure(figsize=(8,6))
            sns.scatterplot(
                x=W_pca[:,0], y=W_pca[:,1], 
                hue=cell_metadata['PlotGroup'],
                s=20,
                palette={'SFTLD': 'tab:blue', 'PN': 'tab:orange', 'Test': 'tab:red'}
            )
            plt.title('PCA of NMF Latent Space')
            plt.xlabel('PC1')
            plt.ylabel('PC2')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            pdf.savefig()
            plt.close()
            
            # UMAP plot
            reducer = umap.UMAP(n_components=2, random_state=42)
            W_umap = reducer.fit_transform(W_all)
            
            plt.figure(figsize=(8,6))
            sns.scatterplot(
                x=W_umap[:,0], y=W_umap[:,1], 
                hue=cell_metadata['PlotGroup'],
                s=20,
                palette={'SFTLD': 'tab:blue', 'PN': 'tab:orange', 'Test': 'tab:red'}
            )
            plt.title('UMAP of NMF Latent Space')
            plt.xlabel('UMAP1')
            plt.ylabel('UMAP2')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.tight_layout()
            pdf.savefig()
            plt.close()
     
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
        
        for _ in range(3):
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
            
            train(model_main, model_domain, loader, 25, device)      ##### NUM EPOCH IS HERE
            
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
    
    out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/regNMF_model_report10_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
    results_df.to_csv(out_path, index=False)




