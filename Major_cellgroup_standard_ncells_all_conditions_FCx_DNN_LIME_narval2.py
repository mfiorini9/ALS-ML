## run this in beluga

salloc -A def-grouleau --time=0-4 -c 4 --mem=150g

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

python  

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

nano LIME_all_conditions_BA9_all_major_narval_standard_n.sh

#!/bin/bash  
#SBATCH --account=def-sfarhan
#SBATCH --time=00-10:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=150g          # memory per cor
#SBATCH --job-name=LIME_all_conditions_BA9_all_major_narval_standard_n
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

python3.8 /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_all_conditions_BA9_all_major_narval_standard_n.py

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nano LIME_all_conditions_BA9_all_major_narval_standard_n.py

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


###################################
# Create loop for all major cell types
###################################

par_keep_celltypes = ["Ex", "In", "Glia", "Vasc"]

# Loop over the different par_keep_cell_type values
for keep_cell_type in par_keep_celltypes:
    print(f"Processing for cell type: {keep_cell_type}")
       
    ###################################
    # Parameters
    ###################################
    ## Here we are keeping C9FTLD and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA9"
    par_status = "all_conditions" ### IMPLEMENT THIS
    #remove = ["SALS", "SFTLD", "C9ALS"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 
    par_keep_cell_type = keep_cell_type
       
    ###################################
    # Cell Specific parameters -- Ex
    ###################################
    data = []
    par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
    print(par_ann_data)
       
    ###################################
    # Load anndata
    ###################################
    adata = sc.read_h5ad(par_ann_data)
    num_cells = adata.X.shape[0]
    adata.obs['Group']
    #adata = adata[~adata.obs['Group'].isin(remove)]
    adata = adata[adata.obs['Major_CellType'] == par_keep_cell_type]
    adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
    print("genes: ", adata.var_names) 
    print("cells: ", adata.obs_names) 
    adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
       
    ###################################
    # clean AnnData object
    ###################################
    #adata.obs = adata.obs.reset_index() 
    adata.obs.columns
    set(adata.obs['Group'])
    set(adata.obs['Major_CellType'])
    set(adata.obs['Region'])
       
    ## change Disease status labels
    # Create a mapping dictionary
    mapping = {'SALS': 1, 'C9ALS': 1, 'SFTLD': 1, 'C9FTLD': 1, 'PN': 0}
       
    # Map the values
    adata.obs['Group'] = adata.obs['Group'].map(mapping)
       
    # Verify the change
    print(adata.obs['Group'].value_counts())
       
    ###################################
    # Subset to 881 cells
    ###################################
    ## filter cells before subset
    sc.pp.filter_cells(adata, min_genes=200)  # Filter cells
       
    # Check the number of cells
    num_cells = adata.shape[0]
       
    # Subset to 881 cells if there are more than 881 cells
    if num_cells > 881:
        # Randomly select 881 cells
        adata = adata[np.random.choice(adata.obs_names, size=881, replace=False)]
        
    # View the new shape of the AnnData object
    print(adata.shape)
       
    adata.obs = adata.obs.reset_index() 
       
    #######################################################
    ## Perform the test train split
    #######################################################
    # Create an array of indices
    indices = np.arange(adata.shape[0])
       
    # Split indices into train and test
    train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=42)
       
    # Create training and testing AnnData objects
    adata_train = adata[train_indices].copy()
    adata_test = adata[test_indices].copy()
       
    #######################################################
    ## Preprocess the train and test objects
    #######################################################
       
    ## Process train data
    sc.pp.filter_genes(adata_train, min_cells=3)    # Filter genes
    sc.pp.normalize_total(adata_train, target_sum=1e4)  # Normalize
    sc.pp.log1p(adata_train)                         # Log-transform
    #sc.pp.highly_variable_genes(source_dataset, flavor='seurat_v3', n_top_genes = 7000)        # Find HVGs
    sc.pp.highly_variable_genes(adata_train) 
    adata_train.raw = adata_train  # Save raw counts for further analysis
    adata_train = adata_train[:, adata_train.var['highly_variable']]  # Subset to HVGs
       
    ## Process test object
    sc.pp.normalize_total(adata_test, target_sum=1e4)  # Normalize
    sc.pp.log1p(adata_test)                         # Log-transform
    adata_test.raw = adata_test  # Save raw counts for further analysis
    adata_test = adata_test[:, adata_train.var_names]  # Ensure the same genes are used
        
    num_HVGs = adata_train.X.shape[1]
        
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
       
    # Parameter bounds for Bayesian Optimization
    pbounds = {
        'learning_rate': (1e-6, 1e-3), 
        #'weight_decay': (1e-5, 1e-1),
        'batch_size': (16, 64)
    }
       
    # Initialize Bayesian Optimization
    optimizer = BayesianOptimization(
        f=train_and_evaluate,
        pbounds=pbounds,
        random_state=42,
    )
       
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
       
    # Perform Bayesian Optimization
    optimizer.maximize(init_points=5, n_iter=20)
       
    # Output the best parameters found
    print("Best parameters:", optimizer.max)
       
    ## set variables
    batch_size = round(optimizer.max['params']['batch_size'])
    learning_rate = optimizer.max['params']['learning_rate']
       
    ###################################
    ## define 5 fold cross validation function
    ###################################
       
    # Function to train the models
    torch.autograd.set_detect_anomaly(True)
        
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
    ## 5 fold cross validation on the train set
    ###################################
       
    # Number of folds
    n_splits = 5
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
       
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
    # Create domain labels based on Donor
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
       
    accuracies = []  # List to store accuracy for each fold
       
    for fold, (train_idx, val_idx) in enumerate(kf.split(np.arange(len(full_dataset)))):
        print(f"Fold {fold + 1}/{n_splits}")
            
        # Create data loaders for the current fold
        train_subset = DataLoader(full_dataset, batch_size=batch_size, sampler=torch.utils.data.SubsetRandomSampler(train_idx))
        val_subset = DataLoader(full_dataset, batch_size=32, sampler=torch.utils.data.SubsetRandomSampler(val_idx))
            
        # Initialize models for each fold
        main_model = MainModel(input_size, num_classes)
        domain_model = DomainClassifier(25, num_domains)  # 50 is the output size from shared layers
            
        # Train the models for the current fold
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        train(main_model, domain_model, train_subset, num_epochs=25, device=device)
            
        # Evaluation on the validation set
        main_model.eval()
        val_predictions = []
        val_labels = []
            
        with torch.no_grad():
            for X_val, y_val, domain_labels in val_subset:
                X_val, y_val, domain_labels = X_val.to(device), y_val.to(device), domain_labels.to(device)
                y_pred_val, _ = main_model(X_val)
                val_predictions.append(torch.argmax(y_pred_val, dim=1).cpu())
                val_labels.append(y_val.cpu())
        
        # Concatenate predictions and labels
        val_predictions = torch.cat(val_predictions).numpy()
        val_labels = torch.cat(val_labels).numpy()
        
        # Calculate accuracy for the current fold
        balanced_accuracy = balanced_accuracy_score(val_labels, val_predictions)
        accuracies.append(balanced_accuracy)
        print(f"Accuracy for Fold {fold + 1}: {balanced_accuracy:.4f}")
       
    # Calculate and print the mean accuracy across all folds
    mean_balanced_accuracy = np.mean(accuracies)
    print(f"\nMean Accuracy across all folds: {mean_balanced_accuracy:.4f}")
       
    # Convert each float to a string and join with semicolons
    balanced_accuracies_str = "; ".join(map(str, accuracies))
       
    # Print the result
    print(balanced_accuracies_str)
       
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
    train(main_model, domain_model, train_dataloader, num_epochs=25, device=device)
       
    ###################################
    ## Evaluate on test set
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
    'region': par_brain_region,
    'group': par_status,
    'celltype': par_keep_cell_type, 
    'n_HVGs': num_HVGs,
    'n_cells': num_cells,
    'learning_rate': learning_rate,
    'batch_size': batch_size,
    'five_fold_accuracies': balanced_accuracies_str,
    'mean_five_fold_accuracies': mean_balanced_accuracy,
    'test_accuracy': test_balanced_accuracy,
    })
       
    ###################################
    ## Export summary csv file 
    ###################################
    # Create a DataFrame from the list
    df = pd.DataFrame(data)
       
    ## write csv file
    file_name = f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_{par_prep}_{par_status}_{par_brain_region}_{par_keep_cell_type}_narval_2.csv"
    df.to_csv(file_name, index=False)


