
salloc -A def-grouleau --time=0-1 -c 1 --mem=100g

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  

################################################################## clean version ------ Focusing on cells most similar to train
import pandas as pd
import numpy as np
import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.model_selection import train_test_split
from sklearn.metrics import balanced_accuracy_score, classification_report
from sklearn.neighbors import KNeighborsClassifier
from torch.utils.data import DataLoader, TensorDataset

# Set parameters
par_prep = "CombatSeq"
par_brain_region = "BA4"
par_status = "SFTLD"
remove = ["SALS", "C9FTLD", "C9ALS"]

cell_types = [
    "L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6", "PV", "5HT3aR", 
    "Rosehip", "SOM", "Oligo", "Astro", "OPC", "Micro",
    "Mural", "Endo", "Fibro", "L5"
]

par_keep_cell_type = "L2_L3"

lime_file = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_.9995_abs_case_control_narval_2.csv'
features = pd.read_csv(lime_file)["gene"].tolist()

report_file = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
report_df = pd.read_csv(report_file)
batch_size = int(report_df.at[0, 'batch_size'])
learning_rate = report_df.at[0, 'learning_rate']

adata_file = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
adata = sc.read_h5ad(adata_file)

adata = adata[~adata.obs['Group'].isin(remove)]
adata = adata[adata.obs['CellType'] == par_keep_cell_type]
adata = adata[adata.obs['Region'] == par_brain_region]
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
sample_IDs = set(adata.obs['orig.ident'])

mapping = {'SFTLD': 1, 'PN': 0}
adata.obs['Group'] = adata.obs['Group'].map(mapping)

adata = adata[:, features]
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata.raw = adata

num_genes = adata.X.shape[1]
data = []

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
    
            confidence_threshold = 0.9
            keep_indices = np.where(test_confidence >= confidence_threshold)[0]
    
            if len(keep_indices) > 0:
                y_logits = y_logits[keep_indices]
                y_labels = torch.argmax(y_logits, dim=1)
                y_true = y_test.iloc[keep_indices].values
            else:
                print("No test cells passed the confidence threshold.")
                y_labels = torch.argmax(y_logits, dim=1)
                y_true = y_test.values
    
            print(classification_report(y_true, y_labels.numpy()))
            acc = balanced_accuracy_score(y_true, y_labels.numpy())
    
            data.append({
                'prep': par_prep, 'donor': donor, 'region': par_brain_region,
                'group': par_status, 'celltype': par_keep_cell_type,
                'n_genes': num_genes, 'n_cells': len(y_true),
                'learning_rate': learning_rate, 'batch_size': batch_size,
                'test_accuracy': acc, 'method': 'LIME + kNN-filtered'
            })

pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
results_df.to_csv(out_path, index=False)











################################################################## clean version ------ TEST WITH GNN

import pandas as pd
import numpy as np
import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from sklearn.metrics import balanced_accuracy_score, classification_report
from torch.utils.data import DataLoader, TensorDataset
from torch_geometric.data import Data as GeoData
from torch_geometric.nn import GCNConv

# Set parameters
par_keep_cell_type = "L2_L3"
par_prep = "CombatSeq"
par_brain_region = "BA4"
par_status = "SFTLD"
remove = ["SALS", "C9FTLD", "C9ALS"]

# Load LIME-selected genes
lime_file = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_.9995_narval_2.csv'
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
adata = adata[:, features]
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata.raw = adata

num_genes = adata.X.shape[1]
data = []

# Define models
class GCNEncoder(nn.Module):
    def __init__(self, in_channels, hidden_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden_channels)
        self.conv2 = GCNConv(hidden_channels, in_channels)
    
    def forward(self, x, edge_index):
        x = F.relu(self.conv1(x, edge_index))
        x = self.conv2(x, edge_index)
        return x

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
def train(main_model, domain_model, gcn_encoder, graph_data, epochs, device):
    main_model.to(device)
    domain_model.to(device)
    gcn_encoder.to(device)
    optimizer_main = optim.Adam(main_model.parameters(), lr=learning_rate)
    optimizer_domain = optim.Adam(domain_model.parameters(), lr=0.001)
    optimizer_gcn = optim.Adam(gcn_encoder.parameters(), lr=0.001)
    loss_class = nn.CrossEntropyLoss()
    loss_domain = nn.CrossEntropyLoss()
    
    for epoch in range(epochs):
        main_model.train()
        domain_model.train()
        gcn_encoder.train()
    
        x = graph_data.x.to(device)
        y = graph_data.y.to(device)
        d = graph_data.d.to(device)
        edge_index = graph_data.edge_index.to(device)
    
        x = gcn_encoder(x, edge_index)
        y_out, shared = main_model(x)
        loss_c = loss_class(y_out, y)
        d_out = domain_model(shared.detach())
        loss_d = loss_domain(d_out, d)
        loss = loss_c - loss_d
    
        optimizer_main.zero_grad()
        optimizer_gcn.zero_grad()
        loss.backward(retain_graph=True)
        optimizer_main.step()
        optimizer_gcn.step()
    
        optimizer_domain.zero_grad()
        loss_d.backward()
        optimizer_domain.step()

# LOSO loop
for donor in sample_IDs:
    for _ in range(5):
        print(f"Processing donor: {donor}")
        adata_train = adata[adata.obs['orig.ident'] != donor].copy()
        adata_test = adata[adata.obs['orig.ident'] == donor].copy()
        num_cells = adata_test.shape[0]
    
        # Build donor-specific neighbor graphs
        donor_graphs = []
        X_train_list, y_train_list, d_train_list = [], [], []
        donor_to_index = {}
        start_idx = 0
    
        for donor_id in adata_train.obs['orig.ident'].unique():
            donor_adata = adata_train[adata_train.obs['orig.ident'] == donor_id].copy()
            sc.pp.neighbors(donor_adata, n_neighbors=3, use_rep='X')
            conn = donor_adata.obsp['connectivities']
            rows, cols = conn.nonzero()
            rows += start_idx
            cols += start_idx
            donor_graphs.append((rows, cols))
    
            X_donor = donor_adata.X.toarray() if hasattr(donor_adata.X, 'toarray') else donor_adata.X
            X_train_list.append(torch.FloatTensor(X_donor))
            y_train_list.append(torch.LongTensor(donor_adata.obs['Group'].values))
            d_train_list.append(torch.LongTensor([donor_to_index.setdefault(donor_id, len(donor_to_index))] * donor_adata.n_obs))
            start_idx += donor_adata.n_obs
    
        X_train = torch.cat(X_train_list)
        y_train = torch.cat(y_train_list)
        d_train = torch.cat(d_train_list)
    
        rows_all, cols_all = zip(*donor_graphs)
        edge_index = torch.tensor(np.vstack((np.concatenate(rows_all), np.concatenate(cols_all))), dtype=torch.long)
    
        graph_data = GeoData(x=X_train, y=y_train, edge_index=edge_index)
        graph_data.d = d_train
    
        input_size = X_train.shape[1]
        num_classes = len(torch.unique(y_train))
        num_domains = len(torch.unique(d_train))
    
        model_main = MainModel(input_size, num_classes)
        model_domain = DomainClassifier(25, num_domains)
        gcn_encoder = GCNEncoder(input_size, 64)
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
        train(model_main, model_domain, gcn_encoder, graph_data, 10, device)
    
        # Evaluation
        X_test = torch.FloatTensor(adata_test.X.toarray() if hasattr(adata_test.X, 'toarray') else adata_test.X)
        y_test = adata_test.obs['Group']
    
        model_main.eval()
        with torch.no_grad():
            X_test = X_test.to(device)
            y_pred, _ = model_main(X_test)
            y_labels = torch.argmax(y_pred, dim=1)
            print(classification_report(y_test, y_labels.cpu().numpy()))
    
        acc = balanced_accuracy_score(y_test.values, y_labels.cpu().numpy())
        data.append({
            'prep': par_prep, 'donor': donor, 'region': par_brain_region,
            'group': par_status, 'celltype': par_keep_cell_type,
            'n_genes': num_genes, 'n_cells': num_cells,
            'learning_rate': learning_rate, 'batch_size': batch_size,
            'test_accuracy': acc, 'method': 'LIME+GCN'
        })

# Save results
pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data)
out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0.9995_abs_case_control_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
results_df.to_csv(out_path, index=False)

















################################################ need to test this one

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, balanced_accuracy_score
import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
import pandas as pd

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

for par_keep_cell_type in cell_types:
    print(par_keep_cell_type)

par_keep_cell_type = "L2_L3"
par_prep = "CombatSeq"
par_brain_region = "BA9"
par_status = "SFTLD"
remove = ["SALS", "C9FTLD", "C9ALS"]

file_path = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_.9995_narval_2.csv'
df = pd.read_csv(file_path)
features = df["gene"].tolist()

file_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
model_report_df = pd.read_csv(file_path)
batch_size = int(model_report_df.at[0, 'batch_size'])
learning_rate = model_report_df.at[0, 'learning_rate']

adata = sc.read_h5ad(f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad")
adata = adata[~adata.obs['Group'].isin(remove)]
adata = adata[adata.obs['CellType'] == par_keep_cell_type]
adata = adata[adata.obs['Region'] == par_brain_region]
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]

sample_ID = set(adata.obs['orig.ident'])
adata.obs = adata.obs.reset_index() 
mapping = {'SFTLD': 1, 'PN': 0}
adata.obs['Group'] = adata.obs['Group'].map(mapping)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata.raw = adata
adata_lim = adata[:, features]
num_genes = adata_lim.X.shape[1]

data = []

import torch.nn.functional as F
class AttentionMILPool(nn.Module):
    def __init__(self, input_dim, hidden_dim=64):
        super().__init__()
        self.attention_a = nn.Linear(input_dim, hidden_dim)
        self.attention_b = nn.Linear(hidden_dim, 1)

    def forward(self, H):
        A = torch.tanh(self.attention_a(H))
        A = self.attention_b(A)
        A = torch.softmax(A, dim=0)
        M = torch.sum(A * H, dim=0)
        return M.unsqueeze(0)

class MainModel(nn.Module):
    def __init__(self, input_size, num_classes):
        super(MainModel, self).__init__()
        self.shared = nn.Sequential(
            nn.Linear(input_size, 100), nn.ReLU(),
            nn.Linear(100, 50), nn.ReLU(),
            nn.Linear(50, 25), nn.ReLU()
        )
        self.att_pool = AttentionMILPool(25)
        self.classifier = nn.Linear(25, num_classes)

    def forward(self, x):
        instance_features = self.shared(x)
        bag_rep = self.att_pool(instance_features)
        output = self.classifier(bag_rep)
        return output, instance_features

    def predict_subject(self, x):
        with torch.no_grad():
            instance_features = self.shared(x)
            bag_rep = self.att_pool(instance_features)
            output = self.classifier(bag_rep)
        return output

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

for donor in sample_ID:
    iter = 1
    while iter <= 3:
        adata_train = adata_lim[adata_lim.obs['orig.ident'] != donor]
        adata_test = adata_lim[adata_lim.obs['orig.ident'] == donor]
        disease_status_test = set(adata_test.obs['Group'])
        num_cells = adata_test.X.shape[0]

        X_train_tensor = torch.FloatTensor(adata_train.X.toarray())
        y_train_tensor = torch.LongTensor(adata_train.obs['Group'].values)
        unique_samples = adata_train.obs['Donor'].unique()
        sample_to_domain = {sample: idx for idx, sample in enumerate(unique_samples)}
        adata_train.obs['Domain_Label'] = adata_train.obs['Donor'].map(sample_to_domain)
        domain_labels_tensor = torch.LongTensor(adata_train.obs['Domain_Label'].values)

        full_dataset = TensorDataset(X_train_tensor, y_train_tensor, domain_labels_tensor)
        input_size = adata_train.shape[1]
        num_classes = len(np.unique(adata_train.obs['Group']))
        num_domains = len(np.unique(domain_labels_tensor.numpy()))

        def train(main_model, domain_model, dataloader, num_epochs, device):
            main_model.to(device)
            domain_model.to(device)
            main_optimizer = optim.Adam(main_model.parameters(), lr=learning_rate)
            domain_optimizer = optim.Adam(domain_model.parameters(), lr=0.001)
            class_loss_fn = nn.CrossEntropyLoss()
            domain_loss_fn = nn.CrossEntropyLoss()
            for epoch in range(num_epochs):
                for X, y, domain_labels in dataloader:
                    X, y, domain_labels = X.to(device), y.to(device), domain_labels.to(device)
                    output, shared_rep = main_model(X)
                    class_loss = class_loss_fn(output, y)
                    domain_output = domain_model(shared_rep.detach())
                    domain_loss = domain_loss_fn(domain_output, domain_labels)
                    total_loss = class_loss - domain_loss
                    main_optimizer.zero_grad()
                    total_loss.backward(retain_graph=True)
                    main_optimizer.step()
                    domain_optimizer.zero_grad()
                    domain_loss.backward()
                    domain_optimizer.step()

        main_model = MainModel(input_size, num_classes)
        domain_model = DomainClassifier(25, num_domains)
        train_dataloader = DataLoader(full_dataset, batch_size=batch_size, shuffle=True)
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        train(main_model, domain_model, train_dataloader, num_epochs=25, device=device)

        X_test_tensor = torch.FloatTensor(adata_test.X.toarray())
        y_test = adata_test.obs['Group'].values

        main_model.eval()
        with torch.no_grad():
            logits = main_model.predict_subject(X_test_tensor)
            pred_labels = torch.argmax(logits, dim=1).cpu().numpy()
            print("Test Classification Report:\n", classification_report(y_test, pred_labels))

        test_balanced_accuracy = balanced_accuracy_score(y_test, pred_labels)
        data.append({
            'prep': par_prep, 'donor': donor, 'region': par_brain_region,
            'group': par_status, 'celltype': par_keep_cell_type, 'n_genes': num_genes,
            'n_cells': num_cells, 'learning_rate': learning_rate,
            'batch_size': batch_size, 'test_accuracy': test_balanced_accuracy,
            'method': 'LIME_ATT_MIL'
        })
        iter += 1

df = pd.DataFrame(data)
file_name = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_ATT_MIL_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
df.to_csv(file_name, index=False)






###################################
## Full LOSO Script with Simple MIL (no GNN) + Attention Pooling + Domain Classifier + Per-Subject Bagging -- CURRENTLY USING THIS ONE.  
###################################

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, balanced_accuracy_score
import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
import pandas as pd

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
# Data Preprocessing
###################################
par_keep_cell_type = "L2_L3"
par_prep = "CombatSeq"
par_brain_region = "BA9"
par_status = "SFTLD"
par_domain = "orig.ident"
remove = ["SALS", "C9FTLD", "C9ALS"]

# Read optimal gene set
file_path = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_.9995_narval_2.csv'
df = pd.read_csv(file_path)
features = df["gene"].tolist()

# Read previous model report
file_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
model_report_df = pd.read_csv(file_path)
batch_size = int(model_report_df.at[0, 'batch_size'])
learning_rate = model_report_df.at[0, 'learning_rate']

# Load and preprocess data
par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
adata = sc.read_h5ad(par_ann_data)
adata = adata[~adata.obs['Group'].isin(remove)]
adata = adata[adata.obs['CellType'] == par_keep_cell_type]
adata = adata[adata.obs['Region'] == par_brain_region]
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
sample_ID = set(adata.obs['orig.ident'])
adata.obs = adata.obs.reset_index()
mapping = {'SFTLD': 1, 'PN': 0}
adata.obs['Group'] = adata.obs['Group'].map(mapping)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata.raw = adata
#adata_lim = adata[:, features]
adata_lim = adata

###################################
# Define Models
###################################
class FeatureExtractor(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(FeatureExtractor, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU()
        )
    
    def forward(self, x):
        return self.net(x)
     
class AttentionMIL(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(AttentionMIL, self).__init__()
        self.attention = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1)
        )
         
    def forward(self, x):
        attn_weights = self.attention(x)
        attn_weights = torch.softmax(attn_weights, dim=0)
        bag_rep = torch.sum(attn_weights * x, dim=0)
        return bag_rep, attn_weights

class MainModel(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_classes):
        super(MainModel, self).__init__()
        self.feature_extractor = FeatureExtractor(input_dim, hidden_dim)
        self.mil_pooling = AttentionMIL(hidden_dim, hidden_dim)
        self.classifier = nn.Linear(hidden_dim, num_classes)
         
    def forward(self, x):
        features = self.feature_extractor(x)
        bag_rep, _ = self.mil_pooling(features)
        logits = self.classifier(bag_rep.unsqueeze(0))
        return logits, bag_rep

class DomainClassifier(nn.Module):
    def __init__(self, input_size, num_domains):
        super(DomainClassifier, self).__init__()
        self.classifier = nn.Sequential(
            nn.Linear(input_size, 25),
            nn.ReLU(),
            nn.Linear(25, num_domains)
        )
     
    def forward(self, x):
        return self.classifier(x)

###################################
# LOSO MIL Training (Subject-Bagged)
###################################
results = []

for donor in sample_ID:
    for iter in range(1):
        print(f"Processing {par_domain.lower()}: {donor}, iteration {iter+1}")
        adata_train = adata_lim[adata_lim.obs['orig.ident'] != donor]
        adata_test = adata_lim[adata_lim.obs['orig.ident'] == donor]
     
        train_subjects = adata_train.obs['orig.ident'].unique()
        domain_map = {s: i for i, s in enumerate(train_subjects)}
     
        model = MainModel(input_dim=adata_train.shape[1], hidden_dim=64, num_classes=2)
        domain_model = DomainClassifier(input_size=64, num_domains=len(domain_map))
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model.to(device)
        domain_model.to(device)
     
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        domain_optimizer = optim.Adam(domain_model.parameters(), lr=0.001)
        criterion = nn.CrossEntropyLoss()
     
        for epoch in range(10):
            model.train()
            domain_model.train()
            total_class_loss = 0
            total_domain_loss = 0
            for subj in train_subjects:
                adata_subj = adata_train[adata_train.obs['orig.ident'] == subj]
                X = adata_subj.X.toarray() if hasattr(adata_subj.X, 'toarray') else adata_subj.X
                y = int(adata_subj.obs['Group'].values[0])
                domain = domain_map[subj]
     
                X_tensor = torch.FloatTensor(X).to(device)
                logits, rep = model(X_tensor)
     
                y_tensor = torch.tensor([y], dtype=torch.long).to(device)
                domain_tensor = torch.tensor([domain], dtype=torch.long).to(device)
     
                class_loss = criterion(logits, y_tensor)
                domain_pred = domain_model(rep.detach())
                domain_loss = criterion(domain_pred.unsqueeze(0), domain_tensor)
     
                loss = class_loss - domain_loss
                optimizer.zero_grad()
                loss.backward(retain_graph=True)
                optimizer.step()
     
                domain_optimizer.zero_grad()
                domain_loss.backward()
                domain_optimizer.step()
     
                total_class_loss += class_loss.item()
                total_domain_loss += domain_loss.item()
     
            print(f"Epoch {epoch+1}, Class Loss: {total_class_loss:.4f}, Domain Loss: {total_domain_loss:.4f}")
     
        model.eval()
        X_test = adata_test.X.toarray() if hasattr(adata_test.X, 'toarray') else adata_test.X
        y_test = adata_test.obs['Group'].astype(int).values
        X_test_tensor = torch.FloatTensor(X_test).to(device)
        with torch.no_grad():
            logits, _ = model(X_test_tensor)
            print("Logits:", logits)
            pred_label = torch.argmax(logits, dim=1).item()
     
        test_pred = [pred_label] * len(y_test)
        test_balanced_accuracy = balanced_accuracy_score(y_test, test_pred)
        print("Donor-Level Classification Report:\n", classification_report(y_test, test_pred))
        results.append({
            'donor': donor,
            'n_cells': X_test.shape[0],
            'test_accuracy': test_balanced_accuracy,
            'predicted_label': pred_label,
            'true_labels': list(y_test),
            'method': 'MIL'
        })
     
test_temp = pd.DataFrame(results)


                       donor  n_cells  test_accuracy  predicted_label                                        true_labels method
0   210526_FTLD_202_snRNA-B5       21            0.0                0  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...    MIL
1    210526_PN_318_snRNA-H12      231            1.0                0  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...    MIL
2     210526_PN_322_snRNA-E6      784            1.0                0  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...    MIL
3     210526_PN_309_snRNA-E2     1110            1.0                0  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...    MIL
4   210526_FTLD_217_snRNA-C2      293            0.0                0  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...    MIL
5   210526_FTLD_216_snRNA-C1      740            0.0                0  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...    MIL
6     210526_PN_301_snRNA-D7      713            1.0                0  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...    MIL
7   210526_FTLD_209_snRNA-B9      675            0.0                0  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...    MIL
8   210526_FTLD_204_snRNA-B6      259            0.0                0  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...    MIL
9    210526_PN_304_snRNA-D10      151            1.0                0  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...    MIL
10    210526_PN_303_snRNA-D9      761            1.0                0  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...    MIL
11    210526_PN_302_snRNA-D8      549            1.0                0  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...    MIL
12    210526_PN_325_snRNA-F5      735            1.0                0  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...    MIL
13    210526_PN_317_snRNA-E4     1488            1.0                0  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...    MIL
14  210526_FTLD_218_snRNA-C3      884            0.0                0  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...    MIL
15  210526_FTLD_225_snRNA-C6      517            0.0                0  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...    MIL
16  210526_FTLD_222_snRNA-C4      557            0.0                0  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...    MIL
17  210526_FTLD_206_snRNA-B7      250            0.0                0  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...    MIL




###################################
## Full LOSO Script with GNN-enhanced MIL + Attention + Domain Classifier + Per-Subject Bagging
###################################

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, balanced_accuracy_score
import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
import pandas as pd
from torch_geometric.nn import GCNConv
from torch_geometric.data import Data

###################################
# Data Preprocessing
###################################
par_keep_cell_type = "L2_L3"
par_prep = "CombatSeq"
par_brain_region = "BA9"
par_status = "SFTLD"
par_domain = "orig.ident"
remove = ["SALS", "C9FTLD", "C9ALS"]

file_path = f'/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_.9995_narval_2.csv'
df = pd.read_csv(file_path)
features = df["gene"].tolist()

file_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
model_report_df = pd.read_csv(file_path)
batch_size = int(model_report_df.at[0, 'batch_size'])
learning_rate = model_report_df.at[0, 'learning_rate']

par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
adata = sc.read_h5ad(par_ann_data)
adata = adata[~adata.obs['Group'].isin(remove)]
adata = adata[adata.obs['CellType'] == par_keep_cell_type]
adata = adata[adata.obs['Region'] == par_brain_region]
adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
sample_ID = set(adata.obs['orig.ident'])
adata.obs = adata.obs.reset_index()
mapping = {'SFTLD': 1, 'PN': 0}
adata.obs['Group'] = adata.obs['Group'].map(mapping)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
adata.raw = adata
adata_lim = adata[:, features]

###################################
# Define Models
###################################
class DenoisingAutoencoder(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(DenoisingAutoencoder, self).__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU()
        )
        self.decoder = nn.Sequential(
            nn.Linear(hidden_dim, input_dim),
            nn.Sigmoid()
        )

    def forward(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded, encoded

class GNNFeatureExtractor(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(GNNFeatureExtractor, self).__init__()
        self.gcn1 = GCNConv(input_dim, hidden_dim)
        self.gcn2 = GCNConv(hidden_dim, hidden_dim)

    def forward(self, x, edge_index):
        x = F.relu(self.gcn1(x, edge_index))
        x = F.relu(self.gcn2(x, edge_index))
        return x

class AttentionMIL(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(AttentionMIL, self).__init__()
        self.attention = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, x):
        attn_weights = self.attention(x)
        attn_weights = torch.softmax(attn_weights, dim=0)
        bag_rep = torch.sum(attn_weights * x, dim=0)
        return bag_rep, attn_weights

class MainModel(nn.Module):
    def __init__(self, input_dim, hidden_dim, num_classes):
        super(MainModel, self).__init__()
        self.feature_extractor = GNNFeatureExtractor(input_dim, hidden_dim)
        self.mil_pooling = AttentionMIL(hidden_dim, hidden_dim)
        self.classifier = nn.Linear(hidden_dim, num_classes)
        self.instance_classifier = nn.Linear(hidden_dim, num_classes)

    def forward(self, x, edge_index):
        features = self.feature_extractor(x, edge_index)
        bag_rep, _ = self.mil_pooling(features)
        logits = self.classifier(bag_rep.unsqueeze(0))
        cell_logits = self.instance_classifier(features)
        return logits, bag_rep, features, cell_logits

class DomainClassifier(nn.Module):
    def __init__(self, input_size, num_domains):
        super(DomainClassifier, self).__init__()
        self.classifier = nn.Sequential(
            nn.Linear(input_size, 25),
            nn.ReLU(),
            nn.Linear(25, num_domains)
        )

    def forward(self, x):
        return self.classifier(x)

###################################
# Graph Utility
###################################
def build_knn_graph(X, k=10):
    nbrs = NearestNeighbors(n_neighbors=k + 1).fit(X)
    knn_graph = nbrs.kneighbors_graph(X).tocoo()
    edge_index = torch.tensor(np.vstack((knn_graph.row, knn_graph.col)), dtype=torch.long)
    return edge_index

###################################
# LOSO Training with GNN-MIL
###################################
results = []
autoencoder = DenoisingAutoencoder(input_dim=adata_lim.shape[1], hidden_dim=64)
autoencoder_optimizer = optim.Adam(autoencoder.parameters(), lr=0.001)
autoencoder.train()

for epoch in range(10):
    for donor in sample_ID:
        adata_sub = adata_lim[adata_lim.obs['orig.ident'] == donor]
        X = adata_sub.X.toarray() if hasattr(adata_sub.X, 'toarray') else adata_sub.X
        X_tensor = torch.FloatTensor(X)
        noisy_X = X_tensor + 0.1 * torch.randn_like(X_tensor)
        decoded, _ = autoencoder(noisy_X)
        loss = F.mse_loss(decoded, X_tensor)
        autoencoder_optimizer.zero_grad()
        loss.backward()
        autoencoder_optimizer.step()

autoencoder.eval()

for donor in sample_ID:
    for iter in range(3):
        print(f"Processing {par_domain}: {donor}, Iteration {iter+1}")
        adata_train = adata_lim[adata_lim.obs['orig.ident'] != donor]
        adata_test = adata_lim[adata_lim.obs['orig.ident'] == donor]

        train_subjects = adata_train.obs['orig.ident'].unique()
        domain_map = {s: i for i, s in enumerate(train_subjects)}

        model = MainModel(input_dim=64, hidden_dim=64, num_classes=2)
        domain_model = DomainClassifier(input_size=64, num_domains=len(domain_map))
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        model.to(device)
        domain_model.to(device)
        autoencoder.to(device)

        optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        domain_optimizer = optim.Adam(domain_model.parameters(), lr=0.001)
        criterion = nn.CrossEntropyLoss()

        for epoch in range(25):
            model.train()
            total_loss = 0
            for subj in train_subjects:
                adata_subj = adata_train[adata_train.obs['orig.ident'] == subj]
                X = adata_subj.X.toarray() if hasattr(adata_subj.X, 'toarray') else adata_subj.X
                y = int(adata_subj.obs['Group'].values[0])
                domain = domain_map[subj]

                X_tensor = torch.FloatTensor(X).to(device)
                _, X_denoised = autoencoder(X_tensor)
                edge_index = build_knn_graph(X_denoised.cpu().numpy()).to(device)

                logits, rep, features, cell_logits = model(X_denoised, edge_index)
                y_tensor = torch.tensor([y], dtype=torch.long).to(device)
                domain_tensor = torch.tensor([domain], dtype=torch.long).to(device)
                instance_labels = torch.tensor([y] * X.shape[0], dtype=torch.long).to(device)

                class_loss = criterion(logits, y_tensor)
                domain_pred = domain_model(rep.detach())
                domain_loss = criterion(domain_pred.unsqueeze(0), domain_tensor)
                instance_loss = criterion(cell_logits, instance_labels)

                loss = class_loss + 0.2 * instance_loss - domain_loss
                optimizer.zero_grad()
                loss.backward(retain_graph=True)
                optimizer.step()

                domain_optimizer.zero_grad()
                domain_loss.backward()
                domain_optimizer.step()

        model.eval()
        X_test = adata_test.X.toarray() if hasattr(adata_test.X, 'toarray') else adata_test.X
        y_test = adata_test.obs['Group'].astype(int).values
        X_test_tensor = torch.FloatTensor(X_test).to(device)
        with torch.no_grad():
            _, X_test_denoised = autoencoder(X_test_tensor)
            edge_index_test = build_knn_graph(X_test_denoised.cpu().numpy()).to(device)
            logits, _, _, cell_logits = model(X_test_denoised, edge_index_test)
            pred_label = torch.argmax(logits, dim=1).item()
            cell_preds = torch.argmax(cell_logits, dim=1).cpu().numpy().tolist()

        test_pred = [pred_label] * len(y_test)
        test_balanced_accuracy = balanced_accuracy_score(y_test, test_pred)
        print("Donor-Level Classification Report:\n", classification_report(y_test, test_pred))
        results.append({
            'donor': donor,
            'n_cells': X_test.shape[0],
            'test_accuracy': test_balanced_accuracy,
            'predicted_label': pred_label,
            'true_labels': list(y_test),
            'cell_predictions': cell_preds,
            'method': 'GNN-MIL'
        })






















###################################
## Full LOSO Script with Simple MIL (no GNN) + Attention Pooling + Domain Classifier + Per-Subject Bagging
###################################


import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import DataLoader, TensorDataset
from sklearn.metrics import balanced_accuracy_score, classification_report
from sklearn.model_selection import train_test_split
import numpy as np

# Feature extractor
class FeatureExtractor(nn.Module):
    def __init__(self, input_size):
        super(FeatureExtractor, self).__init__()
        self.shared = nn.Sequential(
            nn.Linear(input_size, 100),
            nn.ReLU(),
            nn.Linear(100, 50),
            nn.ReLU(),
            nn.Linear(50, 25),
            nn.ReLU()
        )

    def forward(self, x):
        return self.shared(x)

# Attention-based pooling
class MILAttention(nn.Module):
    def __init__(self, input_dim, hidden_dim=128):
        super(MILAttention, self).__init__()
        self.attention = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, H):
        A = self.attention(H)  # [N, 1]
        A = torch.transpose(A, 1, 0)  # [1, N]
        A = F.softmax(A, dim=1)       # attention weights
        M = torch.mm(A, H)            # bag representation [1, input_dim]
        return M, A

# Bag-level classifier
class MILClassifier(nn.Module):
    def __init__(self, input_dim, num_classes):
        super(MILClassifier, self).__init__()
        self.fc = nn.Linear(input_dim, num_classes)

    def forward(self, M):
        return self.fc(M)

# Domain classifier
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


for donor in sample_ID:
    iter = 1
    while iter <= 3:
        print(f"\nProcessing data for donor: {donor} (iter {iter})")

        # Leave-one-subject-out split
        adata_train = adata_lim[adata_lim.obs['orig.ident'] != donor]
        adata_test = adata_lim[adata_lim.obs['orig.ident'] == donor]
        disease_status_test = set(adata_test.obs['Group'])
        num_cells = adata_test.shape[0]

        # Convert training data to tensors
        X_train = adata_train.X.toarray() if hasattr(adata_train.X, 'toarray') else adata_train.X
        X_train_tensor = torch.FloatTensor(X_train)
        y_train = adata_train.obs['Group']
        y_train_tensor = torch.LongTensor(y_train)

        unique_samples = adata_train.obs['Donor'].unique()
        sample_to_domain = {sample: idx for idx, sample in enumerate(unique_samples)}
        adata_train.obs['Domain_Label'] = adata_train.obs['Donor'].map(sample_to_domain)
        domain_labels_tensor = torch.LongTensor(adata_train.obs['Domain_Label'])

        full_dataset = TensorDataset(X_train_tensor, y_train_tensor, domain_labels_tensor)

        input_size = adata_train.shape[1]
        num_classes = len(np.unique(y_train_tensor.numpy()))
        num_domains = len(np.unique(domain_labels_tensor.numpy()))
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # Models
        feature_extractor = FeatureExtractor(input_size).to(device)
        domain_model = DomainClassifier(25, num_domains).to(device)
        mil_classifier = MILClassifier(25, num_classes).to(device)

        # Optimizers
        optimizer_feat = torch.optim.Adam(feature_extractor.parameters(), lr=learning_rate)
        optimizer_domain = torch.optim.Adam(domain_model.parameters(), lr=0.001)
        optimizer_mil = torch.optim.Adam(mil_classifier.parameters(), lr=learning_rate)

        # Losses
        class_loss_fn = nn.CrossEntropyLoss()
        domain_loss_fn = nn.CrossEntropyLoss()

        # Train loop
        dataloader = DataLoader(full_dataset, batch_size=batch_size, shuffle=True)
        for epoch in range(25):
            feature_extractor.train()
            domain_model.train()
            mil_classifier.train()

            for X, y, domain_labels in dataloader:
                X, y, domain_labels = X.to(device), y.to(device), domain_labels.to(device)

                H = feature_extractor(X)
                class_output = mil_classifier(H)
                domain_output = domain_model(H.detach())

                class_loss = class_loss_fn(class_output, y)
                domain_loss = domain_loss_fn(domain_output, domain_labels)
                total_loss = class_loss - domain_loss  # Adversarial

                optimizer_feat.zero_grad()
                optimizer_mil.zero_grad()
                total_loss.backward(retain_graph=True)
                optimizer_feat.step()
                optimizer_mil.step()

                optimizer_domain.zero_grad()
                domain_loss.backward()
                optimizer_domain.step()

        # Evaluation on left-out donor with MIL pooling
        X_test = adata_test.X.toarray() if hasattr(adata_test.X, 'toarray') else adata_test.X
        X_test_tensor = torch.FloatTensor(X_test).to(device)
        y_test = adata_test.obs['Group'].values

        feature_extractor.eval()
        mil_classifier.eval()

        with torch.no_grad():
            H_test = feature_extractor(X_test_tensor)
            attention = MILAttention(input_dim=H_test.shape[1]).to(device)
            bag_rep, attn_weights = attention(H_test)
            output = mil_classifier(bag_rep)
            pred_label = torch.argmax(output, dim=1).item()
            prob = F.softmax(output, dim=1).cpu().numpy()

        # Logging
        test_pred = [pred_label] * len(y_test)  # replicate for per-cell metrics
        test_balanced_accuracy = balanced_accuracy_score(y_test, test_pred)
        print("Donor-Level Classification Report:\n", classification_report(y_test, test_pred))

        data.append({
            'prep': par_prep,
            'donor': donor,
            'region': par_brain_region,
            'group': par_status,
            'celltype': par_keep_cell_type,
            'n_genes': num_genes,
            'n_cells': num_cells,
            'learning_rate': learning_rate,
            'batch_size': batch_size,
            'test_accuracy': test_balanced_accuracy,
            'method': 'MIL'
        })

        iter += 1


###################################
## Full LOSO Script with GNN + Attention-MIL + Domain Classifier + Real kNN Graph
###################################

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, balanced_accuracy_score
import numpy as np
import scanpy as sc
from sklearn.neighbors import NearestNeighbors
from torch_geometric.nn import GATConv
from torch_geometric.data import Data as PyGData, Batch as PyGBatch

# Define the models
class GNNEncoder(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(GNNEncoder, self).__init__()
        self.gnn1 = GATConv(input_dim, hidden_dim)
        self.gnn2 = GATConv(hidden_dim, hidden_dim)

    def forward(self, x, edge_index):
        x = F.relu(self.gnn1(x, edge_index))
        x = F.relu(self.gnn2(x, edge_index))
        return x

class AttentionMIL(nn.Module):
    def __init__(self, input_dim, hidden_dim):
        super(AttentionMIL, self).__init__()
        self.attention = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1)
        )

    def forward(self, x):
        weights = self.attention(x)
        weights = torch.softmax(weights, dim=0)
        weighted_avg = torch.sum(weights * x, dim=0)
        return weighted_avg, weights

class MainModel(nn.Module):
    def __init__(self, input_size, hidden_dim, num_classes):
        super(MainModel, self).__init__()
        self.gnn_encoder = GNNEncoder(input_size, hidden_dim)
        self.mil_pooling = AttentionMIL(hidden_dim, hidden_dim)
        self.classifier = nn.Linear(hidden_dim, num_classes)

    def forward(self, x, edge_index, batch):
        x = self.gnn_encoder(x, edge_index)
        pooled_outputs = []
        for i in torch.unique(batch):
            mask = batch == i
            x_i = x[mask]
            rep, _ = self.mil_pooling(x_i)
            pooled_outputs.append(rep)
        subject_reps = torch.stack(pooled_outputs)
        logits = self.classifier(subject_reps)
        return logits, subject_reps

class DomainClassifier(nn.Module):
    def __init__(self, input_size, num_domains):
        super(DomainClassifier, self).__init__()
        self.classifier = nn.Sequential(
            nn.Linear(input_size, 25),
            nn.ReLU(),
            nn.Linear(25, num_domains)
        )

    def forward(self, x):
        return self.classifier(x)

# Utility to build kNN graph
def build_knn_graph(X, k=10):
    nbrs = NearestNeighbors(n_neighbors=k+1, algorithm='ball_tree').fit(X)
    distances, indices = nbrs.kneighbors(X)
    edge_index = []
    for i in range(X.shape[0]):
        for j in indices[i][1:]:
            edge_index.append([i, j])
            edge_index.append([j, i])
    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    return edge_index

# Training function
def train(main_model, domain_model, dataloader, num_epochs, device):
    main_model.to(device)
    domain_model.to(device)

    main_optimizer = optim.Adam(main_model.parameters(), lr=0.0001)
    domain_optimizer = optim.Adam(domain_model.parameters(), lr=0.001)

    class_loss_fn = nn.CrossEntropyLoss()
    domain_loss_fn = nn.CrossEntropyLoss()

    for epoch in range(num_epochs):
        main_model.train()
        domain_model.train()

        for batch in dataloader:
            batch = batch.to(device)
            x, edge_index, batch_idx, y, domain_labels = batch.x, batch.edge_index, batch.batch, batch.y, batch.domain

            logits, reps = main_model(x, edge_index, batch_idx)
            class_loss = class_loss_fn(logits, y)

            domain_preds = domain_model(reps.detach())
            domain_loss = domain_loss_fn(domain_preds, domain_labels)

            total_loss = class_loss - domain_loss

            main_optimizer.zero_grad()
            total_loss.backward(retain_graph=True)
            main_optimizer.step()

            domain_optimizer.zero_grad()
            domain_loss.backward()
            domain_optimizer.step()

        print(f"Epoch {epoch+1}, Class Loss: {class_loss.item():.4f}, Domain Loss: {domain_loss.item():.4f}")

# Leave-One-Subject-Out Training Loop

data = []

for donor in sample_ID:
    for iter in range(3):
        print(f"Processing donor: {donor}, iteration {iter+1}")

        adata_train = adata_lim[adata_lim.obs['orig.ident'] != donor]
        adata_test = adata_lim[adata_lim.obs['orig.ident'] == donor]

        # Prepare training tensors
        X_train = adata_train.X.toarray() if hasattr(adata_train.X, 'toarray') else adata_train.X
        y_train = adata_train.obs['Group'].astype(int).values
        donor_labels = adata_train.obs['Donor'].values

        domain_map = {s: i for i, s in enumerate(np.unique(donor_labels))}
        domain_labels = np.array([domain_map[d] for d in donor_labels])

        X_tensor = torch.FloatTensor(X_train)
        y_tensor = torch.LongTensor(y_train)
        domain_tensor = torch.LongTensor(domain_labels)

        edge_index = build_knn_graph(X_tensor.numpy(), k=10)
        batch_tensor = torch.zeros(X_tensor.size(0), dtype=torch.long)

        data_obj = PyGData(x=X_tensor, edge_index=edge_index, y=y_tensor[0].unsqueeze(0), domain=domain_tensor[0].unsqueeze(0))
        train_loader = DataLoader(PyGBatch.from_data_list([data_obj]), batch_size=1, shuffle=True)

        input_size = X_tensor.shape[1]
        hidden_dim = 64
        num_classes = len(np.unique(y_tensor))
        num_domains = len(np.unique(domain_tensor))

        main_model = MainModel(input_size, hidden_dim, num_classes)
        domain_model = DomainClassifier(hidden_dim, num_domains)

        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        train(main_model, domain_model, train_loader, num_epochs=25, device=device)

        # Evaluation on test set
        X_test = adata_test.X.toarray() if hasattr(adata_test.X, 'toarray') else adata_test.X
        y_test = adata_test.obs['Group'].astype(int).values
        donor_name = adata_test.obs['Donor'].values[0]

        X_test_tensor = torch.FloatTensor(X_test).to(device)
        edge_index_test = build_knn_graph(X_test_tensor.cpu().numpy(), k=10).to(device)

        main_model.eval()
        with torch.no_grad():
            encoded = main_model.gnn_encoder(X_test_tensor, edge_index_test)
            pooled, _ = main_model.mil_pooling(encoded)
            pred = main_model.classifier(pooled.unsqueeze(0))
            pred_label = torch.argmax(pred, dim=1).item()

        test_pred = [pred_label] * len(y_test)
        test_balanced_accuracy = balanced_accuracy_score(y_test, test_pred)
        print("Donor-Level Classification Report:\n", classification_report(y_test, test_pred))

        data.append({
            'donor': donor,
            'n_cells': X_test.shape[0],
            'test_accuracy': test_balanced_accuracy,
            'predicted_label': pred_label,
            'true_labels': list(y_test),
            'method': 'GNN-MIL'
        })
