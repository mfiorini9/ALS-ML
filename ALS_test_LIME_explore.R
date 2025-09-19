salloc -A def-sfarhan --time=0-8 -c 6 --mem=100g

module load StdEnv/2020 
module load r/4.2.2 
R


############################################################################ Print metadata for Pineda
############################################################################
############################################################################
############################################################################
## Code BA4 SALS
    module load StdEnv/2020 

    module load python/3.8.10
    PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
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
    remove = "C9ALS"

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
        par_ann_data = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_ALS_{par_brain_region}_{par_keep_cell_type}_int.h5ad"
        print(par_ann_data)
        
        ###################################
        # Load anndata
        ###################################
        adata = sc.read_h5ad(par_ann_data)
        num_cells = adata.X.shape[0]
        adata.obs['Group']
        adata = adata[adata.obs['Group'] != remove]
        adata = adata[adata.obs['CellType'] == par_keep_cell_type]
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
        
        adata_train.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda.csv")
        
        adata_test.obs.to_csv(f"/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_{par_status}_{par_brain_region}_{par_keep_cell_type}_pineda.csv")

##

############################################################################ Compute LIME Z score weighted and unweighted -- broad cell tpes
############################################################################
############################################################################
############################################################################
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score.R
   
    ## Load libraries
    library(Seurat)
    library(ggplot2)
    library(Seurat)
    library(ggplot2)
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(tidyr)
    library(stringr)
    library(dplyr)
    library(ggrepel)
    library(phateR)
    library(ComplexHeatmap)
    library(rstatix)
    library(inflection, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(SeuratData, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(SeuratDisk, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(SeuratDisk, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/R/x86_64-pc-linux-gnu-library/4.2")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA4"
    par_status = "SALS"
    par_remove_group = "C9ALS"
    par_prep = "CombatSeq"
    
    prep <- function(celltype2){
        ## Code
        seu <-LoadH5Seurat(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_ALS_",par_brain_region,"_",celltype2,"_int.h5Seurat"))    
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
        
        meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda.csv')
        meta_data <- read.delim(meta_data, header = T, sep = ",")

        test <- meta_data
        
        ## Fold 1
        DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'.csv')
        DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'.csv')
        LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
        
        colnames(LIME) <- c("X", "X0", "importance", "cell_index")

        LIME$test <- str_count(LIME$X0, ' ')

        LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
        LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
        
        ## merge with cell index and remove incorrectly classified cells
        index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
        LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
        nrow(LIME) #4249410
        
        index <- index %>% dplyr::select(Group, predicted_label, cell_index)
        
        LIME <- merge(LIME, index, by = "cell_index", all = T)
        nrow(LIME) 

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
        write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'.csv'), sep = ",")

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
        write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'.csv'), sep = ",")


    }

    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################
    par_brain_region = "BA4"
    par_status = "SALS"
    par_remove_group = "C9ALS"
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##


############################################################################ Compute LIME Z score weighted and unweighted -- refined cell tpes
############################################################################
############################################################################
############################################################################
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_Z_score.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=01-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=50g          # memory per cor
    #SBATCH --job-name=LIME_Z_score
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_Z_score.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_Z_score.R
   
    ## Load libraries
    library(Seurat)
    library(ggplot2)
    library(Seurat)
    library(ggplot2)
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(tidyr)
    library(stringr)
    library(dplyr)
    library(ggrepel)
    library(phateR)
    library(ComplexHeatmap)
    library(rstatix)
    library(inflection, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(SeuratData, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(SeuratDisk, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    ######################################################################
    ## this function assigns a percentile rank to mean Z scores and plots summary figures: both weighted by percent expressed and not weighted.
    ######################################################################
    par_brain_region = "BA4"
    par_status = "SALS"
    par_remove_group = "C9ALS"
    par_prep = "CombatSeq"
    
    celltype2 = "L3_L5"
    
    prep <- function(celltype2){
        ## Code
        seu1 <-LoadH5Seurat(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_ALS_",par_brain_region,"_",celltype2,"_int.h5Seurat"))    
        ncol(seu1) 
        DefaultAssay(seu1)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu1@meta.data$Group) ##
        full_label <- unique(seu1@meta.data$full_label) ##

        for(t in unique(full_label)){
            seu <- seu1
            ## Subset to remove group that we do not want
            unique(seu@meta.data$full_label)
                
                xx <- unique(seu@meta.data$full_label)
                xx <- xx[xx %in% c(t)]

                Idents(seu) <- "full_label"
                seu=subset(seu,idents=xx)

            unique(seu@meta.data$full_label) 
            nrow(seu@meta.data)
        
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
            
            meta_data <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_',par_status,'_',par_brain_region,'_',celltype2,'_pineda.csv')
            meta_data <- read.delim(meta_data, header = T, sep = ",")

            test <- meta_data
            
            ## Fold 1
            DNN_LIME_m4_fold1 <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'.csv')
            DNN_cell_index_m4_fold1 = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_index_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'.csv')
            LIME <- read.delim(DNN_LIME_m4_fold1, header = T, sep = ",")
            
            colnames(LIME) <- c("X", "X0", "importance", "cell_index")

            LIME$test <- str_count(LIME$X0, ' ')

            LIME$feature[LIME$test == 4] <- sub(".* .* (.*) .* .*", "\\1", LIME$X0[LIME$test == 4])
            LIME$feature[LIME$test == 2] <- sub("(.*) .* .*", "\\1", LIME$X0[LIME$test == 2])
            
            ## merge with cell index and remove incorrectly classified cells
            index <- read.delim(DNN_cell_index_m4_fold1, header = T, sep = ",")
            LIME <- LIME %>% dplyr::select(X0, importance, cell_index, feature )
            nrow(LIME) #4249410
            
            index <- index %>% dplyr::select(Group, predicted_label, cell_index)
            
            LIME <- merge(LIME, index, by = "cell_index", all = T)
            nrow(LIME) 

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
            write.csv(Complete_df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_unweighted_',par_prep,'_',par_status,'_',par_brain_region,'_',t,'.csv'), sep = ",")

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
            write.csv(Complete_df_weighted, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',t,'.csv'), sep = ",")
        }

    }

    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

    for(j in unique(celltype_list)){
        prep(j)
    }

    ######################################################################
    ## this function orders the weighted percent expression Z scores and re-prints the dataframe
    ######################################################################
    par_brain_region = "BA4"
    par_status = "SALS"
    par_remove_group = "C9ALS"
    par_prep = "CombatSeq"
    
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",   "T_Cell",  "Mural",   "Endo",    "Fibro", "L5")

    prep <- function(celltype2){
        ## read in csv file
        df <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_weighted_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'.csv'), sep = ",")
        df <- aggregate(. ~ feature, data = df, FUN = function(x) x[1])
        df$Mean_z_score_weighted_across_donor <- as.numeric(df$Mean_z_score_weighted_across_donor)

        df <- df[order( -(df$Mean_z_score_weighted_across_donor)),]
        write.csv(df, file = paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_weighted_ordered_',par_prep,'_',par_status,'_',par_brain_region,'_',celltype2,'.csv'), sep = ",")
    }

    for(j in unique(celltype_list)){
        prep(j)
    }
##
