#run in Narval
salloc -A def-grouleau --time=0-8 -c 1 --mem=100g

module load StdEnv/2023
module load r/4.4.0
R

############################################################################ Print metadata for Pineda
############################################################################
############################################################################
############################################################################
## Code BA4 SALS
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA4"
    par_status = "SALS" ### IMPLEMENT THIS
    remove = ["C9ALS", "C9FTLD", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                    "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'SALS': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA4 C9ALS
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA4"
    par_status = "C9ALS" ### IMPLEMENT THIS
    remove = ["SALS", "C9FTLD", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                    "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'C9ALS': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA4 SFTLD
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python 

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA4"
    par_status = "SFTLD" ### IMPLEMENT THIS
    remove = ["SALS", "C9FTLD", "C9ALS"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                    "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
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
        indices = np.arange(adata.shape[0])
        
        # Split indices into train and test
        train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=42)
        
        # Create training and testing AnnData objects
        adata_train = adata[train_indices].copy()
        adata_test = adata[test_indices].copy()
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA4 C9FTLD
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python 

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA4"
    par_status = "C9FTLD" ### IMPLEMENT THIS
    remove = ["SALS", "SFTLD", "C9ALS"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                    "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'C9FTLD': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA9 SALS
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA9"
    par_status = "SALS" ### IMPLEMENT THIS
    remove = ["C9ALS", "C9FTLD", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                    "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'SALS': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA9 C9ALS
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA9"
    par_status = "C9ALS" ### IMPLEMENT THIS
    remove = ["SALS", "C9FTLD", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                    "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'C9ALS': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA9 SFTLD
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python 

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA9"
    par_status = "SFTLD" ### IMPLEMENT THIS
    remove = ["SALS", "C9FTLD", "C9ALS"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                    "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
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
        indices = np.arange(adata.shape[0])
        
        # Split indices into train and test
        train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=42)
        
        # Create training and testing AnnData objects
        adata_train = adata[train_indices].copy()
        adata_test = adata[test_indices].copy()
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA9 C9FTLD
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python 

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA9"
    par_status = "C9FTLD" ### IMPLEMENT THIS
    remove = ["SALS", "SFTLD", "C9ALS"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                    "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'C9FTLD': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##

## Code BA4 all conditions 
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA4"
    par_status = "all_conditions"
    #remove = ["C9ALS", "C9FTLD", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                    "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        #adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['CellType'])
        set(adata.obs['Region'])
        
        ## Add a disease category column to metadata
        adata.obs['Disease_Group'] = adata.obs['Group']
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'SALS': 1, 'C9ALS': 1, 'SFTLD': 1, 'C9FTLD': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA9 all conditions 
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA9"
    par_status = "all_conditions"
    #remove = ["C9ALS", "C9FTLD", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                    "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                    "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        #adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['CellType'])
        set(adata.obs['Region'])
        
        ## Add a disease category column to metadata
        adata.obs['Disease_Group'] = adata.obs['Group']
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'SALS': 1, 'C9ALS': 1, 'SFTLD': 1, 'C9FTLD': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


####################################################################################################
#################################################################################################### major cell groups


## Code BA4 SALS
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA4"
    par_status = "SALS" ### IMPLEMENT THIS
    remove = ["C9ALS", "C9FTLD", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["Ex", "In", "Glia", "Vasc"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['Major_CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['Major_CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'SALS': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA4 C9ALS
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA4"
    par_status = "C9ALS" ### IMPLEMENT THIS
    remove = ["SALS", "C9FTLD", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["Ex", "In", "Glia", "Vasc"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['Major_CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['Major_CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'C9ALS': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA4 SFTLD
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA4"
    par_status = "SFTLD" ### IMPLEMENT THIS
    remove = ["SALS", "C9FTLD", "C9ALS"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["Ex", "In", "Glia", "Vasc"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['Major_CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['Major_CellType'])
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
        indices = np.arange(adata.shape[0])
        
        # Split indices into train and test
        train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=42)
        
        # Create training and testing AnnData objects
        adata_train = adata[train_indices].copy()
        adata_test = adata[test_indices].copy()
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA4 C9FTLD
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA4"
    par_status = "C9FTLD" ### IMPLEMENT THIS
    remove = ["SALS", "C9ALS", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["Ex", "In", "Glia", "Vasc"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['Major_CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['Major_CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'C9FTLD': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA9 SALS
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA9"
    par_status = "SALS" ### IMPLEMENT THIS
    remove = ["C9ALS", "C9FTLD", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["Ex", "In", "Glia", "Vasc"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['Major_CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['Major_CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'SALS': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA9 C9ALS
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA9"
    par_status = "C9ALS" ### IMPLEMENT THIS
    remove = ["SALS", "C9FTLD", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["Ex", "In", "Glia", "Vasc"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['Major_CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['Major_CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'C9ALS': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA9 SFTLD
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA9"
    par_status = "SFTLD" ### IMPLEMENT THIS
    remove = ["SALS", "C9FTLD", "C9ALS"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["Ex", "In", "Glia", "Vasc"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['Major_CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['Major_CellType'])
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
        indices = np.arange(adata.shape[0])
        
        # Split indices into train and test
        train_indices, test_indices = train_test_split(indices, test_size=0.2, random_state=42)
        
        # Create training and testing AnnData objects
        adata_train = adata[train_indices].copy()
        adata_test = adata[test_indices].copy()
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##


## Code BA9 C9FTLD
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    ###################################
    # Libraries
    ###################################
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from sklearn.model_selection import train_test_split

    ###################################
    # Parameters
    ###################################
    ## Here we are keeping sALS and PN
    par_prep = "CombatSeq"
    par_brain_region = "BA9"
    par_status = "C9FTLD" ### IMPLEMENT THIS
    remove = ["SALS", "C9ALS", "SFTLD"] ## C9ALS, C9FTLD, SALS, SFTLD; choose 3 of the 4. 

    ###################################
    # Set cell type list
    ###################################
    celltype_list = ["Ex", "In", "Glia", "Vasc"]

    ###################################
    # Loop through cell types
    ###################################
    for par_keep_cell_type in celltype_list:  
        print(par_keep_cell_type)
        
        ###################################
        # Cell Specific parameters -- astro
        ###################################
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[~adata.obs['Group'].isin(remove)]
        adata = adata[adata.obs['Major_CellType'] == par_keep_cell_type]
        adata = adata[adata.obs['Region'] == par_brain_region] ## check to see if this work
        print("genes: ", adata.var_names) 
        print("cells: ", adata.obs_names) 
        adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
        
        ###################################
        # create new_annData object with only sALS and not c9ALS
        ###################################
        adata.obs = adata.obs.reset_index() 
        adata.obs.columns
        set(adata.obs['Group'])
        set(adata.obs['Major_CellType'])
        set(adata.obs['Region'])
        
        ## change Disease status labels
        # Create a mapping dictionary
        mapping = {'C9FTLD': 1, 'PN': 0}
        
        # Map the values
        adata.obs['Group'] = adata.obs['Group'].map(mapping)
        
        # Verify the change
        print(adata.obs['Group'].value_counts())
        
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
        
        # Print
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda_narval.csv")
##



############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ SALS BA4
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_SALS_BA4.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_SALS_BA4
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_SALS_BA4.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_SALS_BA4.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA4"
    par_status = "SALS"
    par_remove_group = c("C9ALS", "SFTLD", "C9FTLD")
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")


    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)
        
        ## we are going to merge the meta datas to be able to perform sample specific LIME analyses
        # Test annData object meta data
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #7839476

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME)

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID)
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average PD
        ######################################
        #### Compute sample specific Z scores
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
            
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        z_scores_donor= "fill"
        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
        
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
            LIME_avg$Sample_ID <- i

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor)
            z_scores_donor <- (data-mean(data))/sd(data)
            z_scores_donor

            LIME_avg$z_scores_donor <- z_scores_donor

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
        
        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

        Complete_df <- merge(filler, filler_1, by = "feature")
        
        Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
        
        # Initialize an empty list to store results
        express_fill_list <- list()

        # Loop through unique feature names in LIME_merge_PD$feature
        unique_features <- unique(LIME_merge_PD$feature)

        for (i in unique_features) {
            # Extract gene expression for the current feature i
            expression <- seu@assays$RNA@data[i, ]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the feature
            percent_expressed <- length(cells_expressing) / ncol(seu)
            
            # Store the results in a temporary data frame
            express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
            
            # Append the temporary data frame to the list
            express_fill_list[[i]] <- express_fill_temp
        }

        # Combine all data frames in the list into a single data frame
        express_fill <- do.call(rbind, express_fill_list)

        ## normalize so that they all add to 1
        express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        percent_expressed= "fill"
        normalized_percent_expressed= "fill"
        Mean_feature_importance_donor_weighted= "fill"
        z_scores_weighted_donor= "fill"

        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

        #i = "191112_ALS_110_snRNA-B9.RDS"
        for (i in unique(LIME_merge_PD$Sample_ID)) {
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

            ## weigh by percent expression
            LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
            LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
            z_scores <- (data-mean(data))/sd(data)
            z_scores

            LIME_avg$z_scores_weighted_donor <- z_scores

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

        Complete_df_weighted <- merge(filler, filler_1, by = "feature")
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")


    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##

############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ C9ALS BA4
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_C9ALS_BA4.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_C9ALS_BA4
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_C9ALS_BA4.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_C9ALS_BA4.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA4"
    par_status = "C9ALS"
    par_remove_group = c("SALS", "SFTLD", "C9FTLD")
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")


    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)
        
        ## we are going to merge the meta datas to be able to perform sample specific LIME analyses
        # Test annData object meta data
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #7839476

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME)

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID)
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average PD
        ######################################
        #### Compute sample specific Z scores
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
            
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        z_scores_donor= "fill"
        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
        
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
            LIME_avg$Sample_ID <- i

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor)
            z_scores_donor <- (data-mean(data))/sd(data)
            z_scores_donor

            LIME_avg$z_scores_donor <- z_scores_donor

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
        
        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

        Complete_df <- merge(filler, filler_1, by = "feature")
        
        Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
        
        # Initialize an empty list to store results
        express_fill_list <- list()

        # Loop through unique feature names in LIME_merge_PD$feature
        unique_features <- unique(LIME_merge_PD$feature)

        for (i in unique_features) {
            # Extract gene expression for the current feature i
            expression <- seu@assays$RNA@data[i, ]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the feature
            percent_expressed <- length(cells_expressing) / ncol(seu)
            
            # Store the results in a temporary data frame
            express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
            
            # Append the temporary data frame to the list
            express_fill_list[[i]] <- express_fill_temp
        }

        # Combine all data frames in the list into a single data frame
        express_fill <- do.call(rbind, express_fill_list)

        ## normalize so that they all add to 1
        express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        percent_expressed= "fill"
        normalized_percent_expressed= "fill"
        Mean_feature_importance_donor_weighted= "fill"
        z_scores_weighted_donor= "fill"

        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

        #i = "191112_ALS_110_snRNA-B9.RDS"
        for (i in unique(LIME_merge_PD$Sample_ID)) {
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

            ## weigh by percent expression
            LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
            LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
            z_scores <- (data-mean(data))/sd(data)
            z_scores

            LIME_avg$z_scores_weighted_donor <- z_scores

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

        Complete_df_weighted <- merge(filler, filler_1, by = "feature")
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")


    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##

############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ SFTLD BA4
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_SFTLD_BA4.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_SFTLD_BA4
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_SFTLD_BA4.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_SFTLD_BA4.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA4"
    par_status = "SFTLD"
    par_remove_group = c("C9ALS", "SALS", "C9FTLD")
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")


    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)
        
        ## we are going to merge the meta datas to be able to perform sample specific LIME analyses
        # Test annData object meta data
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #7839476

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME)

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID)
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average PD
        ######################################
        #### Compute sample specific Z scores
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
            
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        z_scores_donor= "fill"
        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
        
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
            LIME_avg$Sample_ID <- i

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor)
            z_scores_donor <- (data-mean(data))/sd(data)
            z_scores_donor

            LIME_avg$z_scores_donor <- z_scores_donor

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
        
        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

        Complete_df <- merge(filler, filler_1, by = "feature")
        
        Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
        
        # Initialize an empty list to store results
        express_fill_list <- list()

        # Loop through unique feature names in LIME_merge_PD$feature
        unique_features <- unique(LIME_merge_PD$feature)

        for (i in unique_features) {
            # Extract gene expression for the current feature i
            expression <- seu@assays$RNA@data[i, ]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the feature
            percent_expressed <- length(cells_expressing) / ncol(seu)
            
            # Store the results in a temporary data frame
            express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
            
            # Append the temporary data frame to the list
            express_fill_list[[i]] <- express_fill_temp
        }

        # Combine all data frames in the list into a single data frame
        express_fill <- do.call(rbind, express_fill_list)

        ## normalize so that they all add to 1
        express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        percent_expressed= "fill"
        normalized_percent_expressed= "fill"
        Mean_feature_importance_donor_weighted= "fill"
        z_scores_weighted_donor= "fill"

        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

        #i = "191112_ALS_110_snRNA-B9.RDS"
        for (i in unique(LIME_merge_PD$Sample_ID)) {
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

            ## weigh by percent expression
            LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
            LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
            z_scores <- (data-mean(data))/sd(data)
            z_scores

            LIME_avg$z_scores_weighted_donor <- z_scores

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

        Complete_df_weighted <- merge(filler, filler_1, by = "feature")
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")


    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##

############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ C9FTLD BA4
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_C9FTLD_BA4.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_C9FTLD_BA4
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_C9FTLD_BA4.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_C9FTLD_BA4.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA4"
    par_status = "C9FTLD"
    par_remove_group = c("C9ALS", "SFTLD", "SALS")
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")


    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)
        
        ## we are going to merge the meta datas to be able to perform sample specific LIME analyses
        # Test annData object meta data
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #7839476

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME)

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID)
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average PD
        ######################################
        #### Compute sample specific Z scores
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
            
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        z_scores_donor= "fill"
        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
        
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
            LIME_avg$Sample_ID <- i

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor)
            z_scores_donor <- (data-mean(data))/sd(data)
            z_scores_donor

            LIME_avg$z_scores_donor <- z_scores_donor

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
        
        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

        Complete_df <- merge(filler, filler_1, by = "feature")
        
        Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
        
        # Initialize an empty list to store results
        express_fill_list <- list()

        # Loop through unique feature names in LIME_merge_PD$feature
        unique_features <- unique(LIME_merge_PD$feature)

        for (i in unique_features) {
            # Extract gene expression for the current feature i
            expression <- seu@assays$RNA@data[i, ]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the feature
            percent_expressed <- length(cells_expressing) / ncol(seu)
            
            # Store the results in a temporary data frame
            express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
            
            # Append the temporary data frame to the list
            express_fill_list[[i]] <- express_fill_temp
        }

        # Combine all data frames in the list into a single data frame
        express_fill <- do.call(rbind, express_fill_list)

        ## normalize so that they all add to 1
        express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        percent_expressed= "fill"
        normalized_percent_expressed= "fill"
        Mean_feature_importance_donor_weighted= "fill"
        z_scores_weighted_donor= "fill"

        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

        #i = "191112_ALS_110_snRNA-B9.RDS"
        for (i in unique(LIME_merge_PD$Sample_ID)) {
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

            ## weigh by percent expression
            LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
            LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
            z_scores <- (data-mean(data))/sd(data)
            z_scores

            LIME_avg$z_scores_weighted_donor <- z_scores

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

        Complete_df_weighted <- merge(filler, filler_1, by = "feature")
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")


    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##


############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ SALS BA9
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_SALS_BA9.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_SALS_BA9
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_SALS_BA9.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_SALS_BA9.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA9"
    par_status = "SALS"
    par_remove_group = c("C9ALS", "SFTLD", "C9FTLD")
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")


    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)
        
        ## we are going to merge the meta datas to be able to perform sample specific LIME analyses
        # Test annData object meta data
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #7839476

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME)

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID)
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average PD
        ######################################
        #### Compute sample specific Z scores
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
            
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        z_scores_donor= "fill"
        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
        
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
            LIME_avg$Sample_ID <- i

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor)
            z_scores_donor <- (data-mean(data))/sd(data)
            z_scores_donor

            LIME_avg$z_scores_donor <- z_scores_donor

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
        
        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

        Complete_df <- merge(filler, filler_1, by = "feature")
        
        Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
        
        # Initialize an empty list to store results
        express_fill_list <- list()

        # Loop through unique feature names in LIME_merge_PD$feature
        unique_features <- unique(LIME_merge_PD$feature)

        for (i in unique_features) {
            # Extract gene expression for the current feature i
            expression <- seu@assays$RNA@data[i, ]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the feature
            percent_expressed <- length(cells_expressing) / ncol(seu)
            
            # Store the results in a temporary data frame
            express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
            
            # Append the temporary data frame to the list
            express_fill_list[[i]] <- express_fill_temp
        }

        # Combine all data frames in the list into a single data frame
        express_fill <- do.call(rbind, express_fill_list)

        ## normalize so that they all add to 1
        express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        percent_expressed= "fill"
        normalized_percent_expressed= "fill"
        Mean_feature_importance_donor_weighted= "fill"
        z_scores_weighted_donor= "fill"

        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

        #i = "191112_ALS_110_snRNA-B9.RDS"
        for (i in unique(LIME_merge_PD$Sample_ID)) {
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

            ## weigh by percent expression
            LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
            LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
            z_scores <- (data-mean(data))/sd(data)
            z_scores

            LIME_avg$z_scores_weighted_donor <- z_scores

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

        Complete_df_weighted <- merge(filler, filler_1, by = "feature")
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")


    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##

############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ C9ALS BA9
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_C9ALS_BA9.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_C9ALS_BA9
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_C9ALS_BA9.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_C9ALS_BA9.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA9"
    par_status = "C9ALS"
    par_remove_group = c("SALS", "SFTLD", "C9FTLD")
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")


    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)
        
        ## we are going to merge the meta datas to be able to perform sample specific LIME analyses
        # Test annData object meta data
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #7839476

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME)

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID)
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average PD
        ######################################
        #### Compute sample specific Z scores
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
            
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        z_scores_donor= "fill"
        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
        
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
            LIME_avg$Sample_ID <- i

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor)
            z_scores_donor <- (data-mean(data))/sd(data)
            z_scores_donor

            LIME_avg$z_scores_donor <- z_scores_donor

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
        
        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

        Complete_df <- merge(filler, filler_1, by = "feature")
        
        Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
        
        # Initialize an empty list to store results
        express_fill_list <- list()

        # Loop through unique feature names in LIME_merge_PD$feature
        unique_features <- unique(LIME_merge_PD$feature)

        for (i in unique_features) {
            # Extract gene expression for the current feature i
            expression <- seu@assays$RNA@data[i, ]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the feature
            percent_expressed <- length(cells_expressing) / ncol(seu)
            
            # Store the results in a temporary data frame
            express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
            
            # Append the temporary data frame to the list
            express_fill_list[[i]] <- express_fill_temp
        }

        # Combine all data frames in the list into a single data frame
        express_fill <- do.call(rbind, express_fill_list)

        ## normalize so that they all add to 1
        express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        percent_expressed= "fill"
        normalized_percent_expressed= "fill"
        Mean_feature_importance_donor_weighted= "fill"
        z_scores_weighted_donor= "fill"

        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

        #i = "191112_ALS_110_snRNA-B9.RDS"
        for (i in unique(LIME_merge_PD$Sample_ID)) {
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

            ## weigh by percent expression
            LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
            LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
            z_scores <- (data-mean(data))/sd(data)
            z_scores

            LIME_avg$z_scores_weighted_donor <- z_scores

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

        Complete_df_weighted <- merge(filler, filler_1, by = "feature")
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")


    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##

############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ SFTLD BA9
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_SFTLD_BA9.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_SFTLD_BA9
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_SFTLD_BA9.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_SFTLD_BA9.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA9"
    par_status = "SFTLD"
    par_remove_group = c("C9ALS", "SALS", "C9FTLD")
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")


    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)
        
        ## we are going to merge the meta datas to be able to perform sample specific LIME analyses
        # Test annData object meta data
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #7839476

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME)

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID)
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average PD
        ######################################
        #### Compute sample specific Z scores
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
            
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        z_scores_donor= "fill"
        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
        
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
            LIME_avg$Sample_ID <- i

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor)
            z_scores_donor <- (data-mean(data))/sd(data)
            z_scores_donor

            LIME_avg$z_scores_donor <- z_scores_donor

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
        
        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

        Complete_df <- merge(filler, filler_1, by = "feature")
        
        Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
        
        # Initialize an empty list to store results
        express_fill_list <- list()

        # Loop through unique feature names in LIME_merge_PD$feature
        unique_features <- unique(LIME_merge_PD$feature)

        for (i in unique_features) {
            # Extract gene expression for the current feature i
            expression <- seu@assays$RNA@data[i, ]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the feature
            percent_expressed <- length(cells_expressing) / ncol(seu)
            
            # Store the results in a temporary data frame
            express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
            
            # Append the temporary data frame to the list
            express_fill_list[[i]] <- express_fill_temp
        }

        # Combine all data frames in the list into a single data frame
        express_fill <- do.call(rbind, express_fill_list)

        ## normalize so that they all add to 1
        express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        percent_expressed= "fill"
        normalized_percent_expressed= "fill"
        Mean_feature_importance_donor_weighted= "fill"
        z_scores_weighted_donor= "fill"

        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

        #i = "191112_ALS_110_snRNA-B9.RDS"
        for (i in unique(LIME_merge_PD$Sample_ID)) {
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

            ## weigh by percent expression
            LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
            LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
            z_scores <- (data-mean(data))/sd(data)
            z_scores

            LIME_avg$z_scores_weighted_donor <- z_scores

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

        Complete_df_weighted <- merge(filler, filler_1, by = "feature")
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")


    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##

############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ C9FTLD BA9
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_C9FTLD_BA9.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_C9FTLD_BA9
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_C9FTLD_BA9.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_C9FTLD_BA9.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA9"
    par_status = "C9FTLD"
    par_remove_group = c("C9ALS", "SFTLD", "SALS")
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")


    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)
        
        ## we are going to merge the meta datas to be able to perform sample specific LIME analyses
        # Test annData object meta data
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #7839476

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME)

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID)
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average PD
        ######################################
        #### Compute sample specific Z scores
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
            
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        z_scores_donor= "fill"
        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
        
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
            LIME_avg$Sample_ID <- i

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor)
            z_scores_donor <- (data-mean(data))/sd(data)
            z_scores_donor

            LIME_avg$z_scores_donor <- z_scores_donor

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
        
        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

        Complete_df <- merge(filler, filler_1, by = "feature")
        
        Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
        
        # Initialize an empty list to store results
        express_fill_list <- list()

        # Loop through unique feature names in LIME_merge_PD$feature
        unique_features <- unique(LIME_merge_PD$feature)

        for (i in unique_features) {
            # Extract gene expression for the current feature i
            expression <- seu@assays$RNA@data[i, ]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the feature
            percent_expressed <- length(cells_expressing) / ncol(seu)
            
            # Store the results in a temporary data frame
            express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
            
            # Append the temporary data frame to the list
            express_fill_list[[i]] <- express_fill_temp
        }

        # Combine all data frames in the list into a single data frame
        express_fill <- do.call(rbind, express_fill_list)

        ## normalize so that they all add to 1
        express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        percent_expressed= "fill"
        normalized_percent_expressed= "fill"
        Mean_feature_importance_donor_weighted= "fill"
        z_scores_weighted_donor= "fill"

        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

        #i = "191112_ALS_110_snRNA-B9.RDS"
        for (i in unique(LIME_merge_PD$Sample_ID)) {
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

            ## weigh by percent expression
            LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
            LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
            z_scores <- (data-mean(data))/sd(data)
            z_scores

            LIME_avg$z_scores_weighted_donor <- z_scores

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

        Complete_df_weighted <- merge(filler, filler_1, by = "feature")
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")


    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##


###########################################################################
## All conditions models
###########################################################################


############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ All conditions model, LIME computation  all together BA4 
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_all_conditions_BA4.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_all_conditions_BA4
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_all_conditions_BA4.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_all_conditions_BA4.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA4"
    par_status = "all_conditions"
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

    celltype2 = "L3_L5"

    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #29196100

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME) #28199800

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID) ## CHANGE THIS TO INCLUDE GROUP
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average 
        ######################################
        #### Compute sample specific Z scores
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
            
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        z_scores_donor= "fill"
        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
        
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
            LIME_avg$Sample_ID <- i

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor)
            z_scores_donor <- (data-mean(data))/sd(data)
            z_scores_donor

            LIME_avg$z_scores_donor <- z_scores_donor

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
        
        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

        Complete_df <- merge(filler, filler_1, by = "feature")
        
        Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
        
        # Initialize an empty list to store results
        express_fill_list <- list()

        # Loop through unique feature names in LIME_merge_PD$feature
        unique_features <- unique(LIME_merge_PD$feature)

        for (i in unique_features) {
            # Extract gene expression for the current feature i
            expression <- seu@assays$RNA@data[i, ]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the feature
            percent_expressed <- length(cells_expressing) / ncol(seu)
            
            # Store the results in a temporary data frame
            express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
            
            # Append the temporary data frame to the list
            express_fill_list[[i]] <- express_fill_temp
        }

        # Combine all data frames in the list into a single data frame
        express_fill <- do.call(rbind, express_fill_list)

        ## normalize so that they all add to 1
        express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        percent_expressed= "fill"
        normalized_percent_expressed= "fill"
        Mean_feature_importance_donor_weighted= "fill"
        z_scores_weighted_donor= "fill"

        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

            ## weigh by percent expression
            LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
            LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
            z_scores <- (data-mean(data))/sd(data)
            z_scores

            LIME_avg$z_scores_weighted_donor <- z_scores

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

        Complete_df_weighted <- merge(filler, filler_1, by = "feature")
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")


    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##


############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ All conditions model, LIME computation all together BA9
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_all_conditions_BA9.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_all_conditions_BA9
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_all_conditions_BA9.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_all_conditions_BA9.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA9"
    par_status = "all_conditions"
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

    celltype2 = "L3_L5"

    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #29196100

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME) #28199800

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID) ## CHANGE THIS TO INCLUDE GROUP
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average 
        ######################################
        #### Compute sample specific Z scores
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
            
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        z_scores_donor= "fill"
        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
        
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
            LIME_avg$Sample_ID <- i

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor)
            z_scores_donor <- (data-mean(data))/sd(data)
            z_scores_donor

            LIME_avg$z_scores_donor <- z_scores_donor

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
        
        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

        Complete_df <- merge(filler, filler_1, by = "feature")
        
        Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        LIME_merge_PD <- subset(LIME_merge, Group == 1)
        
        # Initialize an empty list to store results
        express_fill_list <- list()

        # Loop through unique feature names in LIME_merge_PD$feature
        unique_features <- unique(LIME_merge_PD$feature)

        for (i in unique_features) {
            # Extract gene expression for the current feature i
            expression <- seu@assays$RNA@data[i, ]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the feature
            percent_expressed <- length(cells_expressing) / ncol(seu)
            
            # Store the results in a temporary data frame
            express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
            
            # Append the temporary data frame to the list
            express_fill_list[[i]] <- express_fill_temp
        }

        # Combine all data frames in the list into a single data frame
        express_fill <- do.call(rbind, express_fill_list)

        ## normalize so that they all add to 1
        express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                
        ## Create filler frame
        feature = "fill"
        Mean_feature_importance_donor= "fill"
        Sample_ID= "fill"
        percent_expressed= "fill"
        normalized_percent_expressed= "fill"
        Mean_feature_importance_donor_weighted= "fill"
        z_scores_weighted_donor= "fill"

        filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

        for (i in unique(LIME_merge_PD$Sample_ID)) {
            temper <- subset(LIME_merge_PD, Sample_ID == i)
            
            ## Compute average feature importance across cells
            LIME_avg <- temper %>%
                group_by(feature) %>%
                summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

            ## weigh by percent expression
            LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
            LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

            ## compute Z scores
            data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
            z_scores <- (data-mean(data))/sd(data)
            z_scores

            LIME_avg$z_scores_weighted_donor <- z_scores

            filler <- rbind(filler, LIME_avg) 
        }     
        
        filler <- subset(filler, feature != "fill")
        
        ## take the average Z score for gene
        filler_1 <- filler %>%
            group_by(feature) %>%
            summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

        filler_1 <- data.frame(filler_1)

        filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

        Complete_df_weighted <- merge(filler, filler_1, by = "feature")
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")


    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##

############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ All conditions model, LIME computation separate for each disease group  BA4 
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_all_conditions_disease_group_sep_BA4.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=02-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_all_conditions_disease_group_sep_BA4
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_all_conditions_disease_group_sep_BA4.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_all_conditions_disease_group_sep_BA4.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA4"
    par_status = "all_conditions"
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

    celltype2 = "L3_L5"

    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")
        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #29196100

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME) #28199800

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID, Disease_Group) ## CHANGE THIS TO INCLUDE GROUP
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average 
        ######################################
        #### Compute sample specific Z scores
        disease_groups <- c("SALS", "C9ALS",  "SFTLD",  "C9FTLD")
        
        for(disease_group in unique(disease_groups)){
            ## calculate the average importance for each gene  
            LIME_merge_PD <- subset(LIME_merge, Group == 1)
            LIME_merge_PD <- subset(LIME_merge_PD, Disease_Group == disease_group)
                
            ## Create filler frame
            feature = "fill"
            Mean_feature_importance_donor= "fill"
            Sample_ID= "fill"
            z_scores_donor= "fill"
            filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

            for (i in unique(LIME_merge_PD$Sample_ID)) {
            
                temper <- subset(LIME_merge_PD, Sample_ID == i)
                
                ## Compute average feature importance across cells
                LIME_avg <- temper %>%
                    group_by(feature) %>%
                    summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

                ## compute Z scores
                data <- abs(LIME_avg$Mean_feature_importance_donor)
                z_scores_donor <- (data-mean(data))/sd(data)
                z_scores_donor

                LIME_avg$z_scores_donor <- z_scores_donor

                filler <- rbind(filler, LIME_avg) 
            }     
            
            filler <- subset(filler, feature != "fill")
            
            ## take the average Z score for gene
            filler_1 <- filler %>%
                group_by(feature) %>%
                summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
            
            filler_1 <- data.frame(filler_1)

            filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

            Complete_df <- merge(filler, filler_1, by = "feature")
            
            Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
            write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',disease_group, '_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        }

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        #### Compute sample specific Z scores
        disease_groups <- c("SALS", "C9ALS",  "SFTLD",  "C9FTLD")
        
        for(disease_group in unique(disease_groups)){
            print(disease_group)
            ## calculate the average importance for each gene  
            LIME_merge_PD <- subset(LIME_merge, Group == 1)
            LIME_merge_PD <- subset(LIME_merge_PD, Disease_Group == disease_group)
        
            # Initialize an empty list to store results
            express_fill_list <- list()

            # Loop through unique feature names in LIME_merge_PD$feature
            unique_features <- unique(LIME_merge_PD$feature)

            for (i in unique_features) {
                # Extract gene expression for the current feature i
                expression <- seu@assays$RNA@data[i, ]
                
                # Count cells where expression > 0
                cells_expressing <- expression[expression > 0]
                
                # Calculate percent of cells expressing the feature
                percent_expressed <- length(cells_expressing) / ncol(seu)
                
                # Store the results in a temporary data frame
                express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
                
                # Append the temporary data frame to the list
                express_fill_list[[i]] <- express_fill_temp
            }

            # Combine all data frames in the list into a single data frame
            express_fill <- do.call(rbind, express_fill_list)

            ## normalize so that they all add to 1
            express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                    
            ## Create filler frame
            feature = "fill"
            Mean_feature_importance_donor= "fill"
            Sample_ID= "fill"
            percent_expressed= "fill"
            normalized_percent_expressed= "fill"
            Mean_feature_importance_donor_weighted= "fill"
            z_scores_weighted_donor= "fill"

            filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

            for (i in unique(LIME_merge_PD$Sample_ID)) {
                temper <- subset(LIME_merge_PD, Sample_ID == i)
                
                ## Compute average feature importance across cells
                LIME_avg <- temper %>%
                    group_by(feature) %>%
                    summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                    LIME_avg$Sample_ID <- i

                ## weigh by percent expression
                LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
                LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

                ## compute Z scores
                data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
                z_scores <- (data-mean(data))/sd(data)
                z_scores

                LIME_avg$z_scores_weighted_donor <- z_scores

                filler <- rbind(filler, LIME_avg) 
            }     
            
            filler <- subset(filler, feature != "fill")
            
            ## take the average Z score for gene
            filler_1 <- filler %>%
                group_by(feature) %>%
                summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

            filler_1 <- data.frame(filler_1)

            filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

            Complete_df_weighted <- merge(filler, filler_1, by = "feature")
            write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',disease_group, '_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        }    

    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    disease_groups <- c("SALS", "C9ALS",  "SFTLD",  "C9FTLD")
        
    
    
    prep <- function(celltype2){
        for(disease_group in unique(disease_groups)){
            ## read in csv file
            print(disease_group)
            print(celltype2)
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',disease_group, '_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered',par_prep,'_',par_status,'_',disease_group, '_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        }
    }    

    for(j in unique(celltype_list)){
        prep(j)
    }
##


############################################################################ Compute LIME Z score weighted and unweighted -- broad cell types
############################################################################
############################################################################
############################################################################ All conditions model, LIME computation separate for each disease group  BA9 
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score_all_conditions_disease_group_sep_BA9.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=02-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score_all_conditions_disease_group_sep_BA9
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score_all_conditions_disease_group_sep_BA9.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score_all_conditions_disease_group_sep_BA9.R
   

    ## Load libraries
    library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(withr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
    library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(backports, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(stringr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(rstatix, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
    library(remotes, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") 
    
    #library(SeuratData, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    #library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")
    #/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2
    #install.packages("ggplot2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("withr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("backports", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("ggpubr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("SeuratData", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("remotes", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #install.packages("hdf5r", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )
    #remotes::install_github("mojaveazure/seurat-disk", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA9"
    par_status = "all_conditions"
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

    celltype2 = "L3_L5"

    prep <- function(celltype2){
        ## Code
        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda_narval.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")
        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'_narval.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) 
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) #29196100

        ## check overlap
        test$X %in% LIME$cell_index
        
        ## only keep correctly classified instances
        LIME <- LIME[LIME$Group == LIME$predicted_label,]
        nrow(LIME) #28199800

        ## subset the metadata object
        test_2 <- subset(test, X %in% LIME$cell_index)
        test_2 <- test_2 %>% dplyr::select(X, Sample_ID, Disease_Group) ## CHANGE THIS TO INCLUDE GROUP
        nrow(test_2) == length(unique(LIME$cell_index))
        test_2 <- test_2[!duplicated(test_2$X),]
        nrow(test_2)== length(unique(LIME$cell_index))

        LIME_merge <- merge(LIME, test_2, by.x = "cell_index", by.y = "X", all.X = T, all.y = F )
        nrow(LIME_merge) == nrow(LIME)

        table(LIME_merge$Group, LIME_merge$Sample_ID)
        length(unique(LIME_merge$Sample_ID))
        
        ######################################
        ## Z score average 
        ######################################
        #### Compute sample specific Z scores
        disease_groups <- c("SALS", "C9ALS",  "SFTLD",  "C9FTLD")
        
        for(disease_group in unique(disease_groups)){
            ## calculate the average importance for each gene  
            LIME_merge_PD <- subset(LIME_merge, Group == 1)
            LIME_merge_PD <- subset(LIME_merge_PD, Disease_Group == disease_group)
                
            ## Create filler frame
            feature = "fill"
            Mean_feature_importance_donor= "fill"
            Sample_ID= "fill"
            z_scores_donor= "fill"
            filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, z_scores_donor)

            for (i in unique(LIME_merge_PD$Sample_ID)) {
            
                temper <- subset(LIME_merge_PD, Sample_ID == i)
                
                ## Compute average feature importance across cells
                LIME_avg <- temper %>%
                    group_by(feature) %>%
                    summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                LIME_avg$Sample_ID <- i

                ## compute Z scores
                data <- abs(LIME_avg$Mean_feature_importance_donor)
                z_scores_donor <- (data-mean(data))/sd(data)
                z_scores_donor

                LIME_avg$z_scores_donor <- z_scores_donor

                filler <- rbind(filler, LIME_avg) 
            }     
            
            filler <- subset(filler, feature != "fill")
            
            ## take the average Z score for gene
            filler_1 <- filler %>%
                group_by(feature) %>%
                summarize(Mean_z_score_across_donors = mean(as.numeric(z_scores_donor), na.rm=TRUE))
            
            filler_1 <- data.frame(filler_1)

            filler_1 <- filler_1[order( -(filler_1$Mean_z_score_across_donors)),]

            Complete_df <- merge(filler, filler_1, by = "feature")
            
            Complete_df <- Complete_df[order( -(Complete_df$Mean_z_score_across_donors)),]
            write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_unweighted_',par_prep,'_',par_status,'_',disease_group, '_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        }

        ######################################
        ## Weighted by percent expression
        ##################################### 
        ## calculate the average importance for each gene  
        #### Compute sample specific Z scores
        disease_groups <- c("SALS", "C9ALS",  "SFTLD",  "C9FTLD")
        
        for(disease_group in unique(disease_groups)){
            print(disease_group)
            ## calculate the average importance for each gene  
            LIME_merge_PD <- subset(LIME_merge, Group == 1)
            LIME_merge_PD <- subset(LIME_merge_PD, Disease_Group == disease_group)
        
            # Initialize an empty list to store results
            express_fill_list <- list()

            # Loop through unique feature names in LIME_merge_PD$feature
            unique_features <- unique(LIME_merge_PD$feature)

            for (i in unique_features) {
                # Extract gene expression for the current feature i
                expression <- seu@assays$RNA@data[i, ]
                
                # Count cells where expression > 0
                cells_expressing <- expression[expression > 0]
                
                # Calculate percent of cells expressing the feature
                percent_expressed <- length(cells_expressing) / ncol(seu)
                
                # Store the results in a temporary data frame
                express_fill_temp <- data.frame(percent_expressed = percent_expressed, feature = i)
                
                # Append the temporary data frame to the list
                express_fill_list[[i]] <- express_fill_temp
            }

            # Combine all data frames in the list into a single data frame
            express_fill <- do.call(rbind, express_fill_list)

            ## normalize so that they all add to 1
            express_fill$normalized_percent_expressed <- express_fill$percent_expressed/sum(express_fill$percent_expressed)
                    
            ## Create filler frame
            feature = "fill"
            Mean_feature_importance_donor= "fill"
            Sample_ID= "fill"
            percent_expressed= "fill"
            normalized_percent_expressed= "fill"
            Mean_feature_importance_donor_weighted= "fill"
            z_scores_weighted_donor= "fill"

            filler <- data.frame(feature, Mean_feature_importance_donor, Sample_ID, percent_expressed, normalized_percent_expressed, Mean_feature_importance_donor_weighted, z_scores_weighted_donor)

            for (i in unique(LIME_merge_PD$Sample_ID)) {
                temper <- subset(LIME_merge_PD, Sample_ID == i)
                
                ## Compute average feature importance across cells
                LIME_avg <- temper %>%
                    group_by(feature) %>%
                    summarize(Mean_feature_importance_donor = mean(importance, na.rm=TRUE))
                    LIME_avg$Sample_ID <- i

                ## weigh by percent expression
                LIME_avg <- merge(LIME_avg,express_fill, by = "feature"  )
                LIME_avg$Mean_feature_importance_donor_weighted <- LIME_avg$Mean_feature_importance_donor*LIME_avg$normalized_percent_expressed

                ## compute Z scores
                data <- abs(LIME_avg$Mean_feature_importance_donor_weighted)
                z_scores <- (data-mean(data))/sd(data)
                z_scores

                LIME_avg$z_scores_weighted_donor <- z_scores

                filler <- rbind(filler, LIME_avg) 
            }     
            
            filler <- subset(filler, feature != "fill")
            
            ## take the average Z score for gene
            filler_1 <- filler %>%
                group_by(feature) %>%
                summarize(Mean_z_score_weighted_across_donor = mean(as.numeric(z_scores_weighted_donor), na.rm=TRUE))

            filler_1 <- data.frame(filler_1)

            filler_1 <- filler_1[order( -(filler_1$Mean_z_score_weighted_across_donor)),]

            Complete_df_weighted <- merge(filler, filler_1, by = "feature")
            write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',disease_group, '_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        }    

    }

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################

    disease_groups <- c("SALS", "C9ALS",  "SFTLD",  "C9FTLD")
        
    
    
    prep <- function(celltype2){
        for(disease_group in unique(disease_groups)){
            ## read in csv file
            print(disease_group)
            print(celltype2)
            df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_',par_prep,'_',par_status,'_',disease_group, '_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
            df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
            df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

            df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
            write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/pineda_LIME_weighted_ordered',par_prep,'_',par_status,'_',disease_group, '_',par_brain_region,'_',celltype2,'_narval.csv'), sep = ",")
        }
    }    

    for(j in unique(celltype_list)){
        prep(j)
    }
##


