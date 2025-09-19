import scanpy as sc
import torch
import torch.nn.functional as F
from torch import nn
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns

###################################
# Parameters
###################################
par_prep = "CombatSeq"
par_brain_region = "BA9"
par_status = "SFTLD"
remove = ["C9ALS", "C9FTLD", "SALS"]
par_keep_cell_type = "L2_L3"

###################################
# Cell Specific parameters -- Astro
###################################
par_ann_data_Pineda = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"

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

###################################
# Filter and QC
###################################
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, min_counts=200)
sc.pp.filter_genes(adata, min_cells=3)

###################################
# Normalize and Log Transform
###################################
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

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

# Prepare tensors
X = torch.tensor(adata.X.toarray(), dtype=torch.float32, device=device)
subjects = adata.obs['orig.ident'].astype('category').cat.codes.values
encoder = OneHotEncoder(sparse_output=False)
subject_onehot = torch.tensor(encoder.fit_transform(subjects.reshape(-1, 1)), dtype=torch.float32, device=device)
disease_labels = torch.tensor(adata.obs['Group'].values, dtype=torch.float32, device=device)

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





###################################
# Visualize Latent Components
###################################
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import matplotlib.pyplot as plt

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

# Call the function to plot first 10 components
plot_latent_components_to_pdf(model.W, disease_labels, num_components=10, output_path='/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf')


###################################
# Add to W components to the addata object
###################################
## W = transformed data matrix, V = original feature matrix
W = model.W
H = model.H
W.shape
H.shape

W_np = model.W.detach().cpu().numpy()

new_var = pd.DataFrame(index=[f'feature_{i}' for i in range(W.shape[1])])

## Create new anndata object for method 3
adata_NMF = sc.AnnData(X = W_np.copy(),
  obs = adata.obs.copy(),
  var = new_var,
  uns = adata.uns.copy(),
  obsm = adata.obsm.copy(),
  varm = adata.varm.copy(),
  layers = adata.layers.copy(),
  raw = adata.raw.copy(),
  dtype = "float32",
  shape = None,
  obsp = adata.obsp.copy(),
  varp = adata.varp
  )

adata_NMF.__dict__['_raw'].__dict__['_var'] = adata_NMF.__dict__['_raw'].__dict__['_var'].rename(
    columns={'_index': 'features'})



###################################
# Implementation into the LOSO DNN
###################################

num_genes = adata_NMF.X.shape[1]
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
learning_rate = 0.0008663099696291
sample_IDs = set(adata_NMF.obs['orig.ident'])

for donor in sample_IDs:
    for _ in range(3):
        print(f"Processing donor: {donor}")
        adata_train = adata_NMF[adata_NMF.obs['orig.ident'] != donor]
        adata_test = adata_NMF[adata_NMF.obs['orig.ident'] == donor]
        num_cells = adata_test.shape[0]
    
        X_train = torch.FloatTensor(adata_train.X.toarray() if hasattr(adata_train.X, 'toarray') else adata_train.X)
        y_train = torch.LongTensor(adata_train.obs['Group'].values)
        domains = pd.factorize(adata_train.obs['Donor'])[0]
        d_train = torch.LongTensor(domains)
    
        dataset = TensorDataset(X_train, y_train, d_train)
        input_size = adata_train.shape[1]
        num_classes = len(np.unique(y_train))
        num_domains = len(np.unique(domains))
    
        loader = DataLoader(dataset, batch_size=64, shuffle=True)
    
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
            'learning_rate': learning_rate, 'batch_size': 64,
            'test_accuracy': acc, 'method': 'LIME'
        })

# Save results
pd.set_option('display.max_rows', None)
results_df = pd.DataFrame(data)


#out_path = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
#results_df.to_csv(out_path, index=False)


210526_PN_322_snRNA-E6

###################################
# Fit model to unseen subject
###################################
def fit_new_subject(adata_new, trained_model, n_iter=1000, lr=1e-2):
    sc.pp.filter_cells(adata_new, min_genes=200)
    sc.pp.filter_cells(adata_new, min_counts=200)
    sc.pp.filter_genes(adata_new, min_cells=3)
    sc.pp.normalize_total(adata_new, target_sum=1e4)
    sc.pp.log1p(adata_new)

    X_new = torch.tensor(adata_new.X.toarray(), dtype=torch.float32, device=device)
    n_cells_new = X_new.shape[0]
    n_components = trained_model.H.shape[0]

    W_new = nn.Parameter(torch.rand(n_cells_new, n_components, device=device))
    H_fixed = trained_model.H.detach()

    optimizer = torch.optim.Adam([W_new], lr=lr)

    for i in range(n_iter):
        optimizer.zero_grad()
        recon = W_new @ H_fixed
        loss = F.mse_loss(recon, X_new)
        loss.backward()
        optimizer.step()
        with torch.no_grad():
            W_new.clamp_(min=0)
        if i % 200 == 0:
            print(f"[New Subject] Iter {i}: Loss = {loss.item():.4f}")

    return W_new

# Example usage for unseen subject
# adata_new = sc.read_h5ad("/path/to/new_subject.h5ad")
# W_new = fit_new_subject(adata_new, model)





###################################
# Clustering based on Latent Components
###################################
def cluster_latent_components(W_tensor, n_clusters=5):
    W_np = W_tensor.detach().cpu().numpy()
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(W_np)
    return clusters

# Cluster cells and add to adata.obs
adata.obs['latent_cluster'] = cluster_latent_components(model.W, n_clusters=5)

# Visualize cluster composition by disease status
plt.figure(figsize=(6, 4))
sns.countplot(x='latent_cluster', hue='Group', data=adata.obs)
plt.title("Cluster Composition by Disease Status")
plt.tight_layout()
plt.show()

