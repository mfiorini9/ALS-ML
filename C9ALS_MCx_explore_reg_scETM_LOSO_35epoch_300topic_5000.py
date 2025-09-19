salloc -A def-grouleau --time=0-8 -c 1 --mem=200g

module load StdEnv/2020 
module load python/3.8.10

python  

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
L3_L5   ##### JOB ID: 43265174

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

nano C9ALS_MCx_explore_reg_scETM_LOSO_35epoch_300topic_5000.sh

#!/bin/bash 
#SBATCH --account=def-grouleau
#SBATCH --time=07-00:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=300g          # memory per cor
#SBATCH --job-name=C9ALS_MCx_explore_reg_scETM_LOSO_35epoch_300topic_5000
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/C9ALS_MCx_explore_reg_scETM_LOSO_35epoch_300topic_5000.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano C9ALS_MCx_explore_reg_scETM_LOSO_35epoch_300topic_5000.py


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

###################################
# Parameters
###################################
par_prep = "CombatSeq"
par_brain_region = "BA4"
par_status = "C9ALS"
remove = ["SALS", "C9FTLD", "SFTLD"]
#par_keep_cell_type = "L2_L3"
par_brain_region_Li = "frontal cortex"
remove_li = ['C9-FTD']

cell_types = [
    "L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Oligo", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5"
]

n_topics_list = [300]
par_n_epochs = 35

for par_keep_cell_type in cell_types:
    print(par_keep_cell_type)   
     
    for par_n_topics in n_topics_list:
        print(par_n_topics)
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
        mapping = {'C9ALS': 1, 'PN': 0}
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
            adata_test = adata[adata.obs['orig.ident'] == donor]
        
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
            group_mapping = {0: 'PN', 1: 'C9ALS'}
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
        
            for _ in range(3):
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
                data.append({
                    'prep': par_prep, 'donor': donor, 'region': par_brain_region,
                    'group': par_status, 'celltype': par_keep_cell_type,
                    'n_genes': num_genes, 'n_cells': num_cells,
                    'learning_rate': learning_rate, 'batch_size': 64,
                    'test_accuracy': acc, 'method': 'LIME'
                })
      
            # Delete obects to conserve memory
            del adata_train, adata_test
            del model, trainer, adata_train_ETM, adata_test_ETM
            del X_train, X_test, X_all, cell_metadata_train, cell_metadata_test, cell_metadata
            del model_main, model_domain
      
        # Save results
        pd.set_option('display.max_rows', None)
        results_df = pd.DataFrame(data)
        
        out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/regETM_model_report2_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_{par_n_topics}_{par_n_epochs}_narval_2_5000.csv"
        results_df.to_csv(out_path, index=False)
        




