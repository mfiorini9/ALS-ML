salloc -A def-grouleau --time=0-8 -c 1 --mem=300g

module load StdEnv/2023
module load r/4.4.0
R

## load libraries
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
library(labeling, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(scCustomize, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

#install.packages("scCustomize", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )

################################################################## ## TXd test
##################################################################
##################################################################
##################################################################
##################################################################
################################################################## Pineda

## TXd clean
    
    ## NOTES: ########################################################################################
    
    #"The TxD score is a quantification of transcriptional dysregulation. 
    #It represents the change in transcriptome-wide gene expression of each 
    #subtype in disease from its respective PN expression profile. 
    #The divergence score is the Euclidean distance between the median disease 
    #and corresponding PN covariate-corrected, pseudo-bulk expression profiles for each cell type."
    
    ## Functions: ######################################################################################## 

    # Create a function to calculate the pseudo-bulk expression profile (median of gene expression per cell type)
    calculate_pseudo_bulk <- function(seurat_obj, celltype) {
        # Subset the Seurat object by celltype
        celltype_data <- subset(seurat_obj, CellType == celltype)
        
        # Calculate pseudo-bulk expression
        pseudo_bulk_profile <- AggregateExpression(celltype_data, return.seurat = TRUE, verbose = FALSE, group.by = "Sample_ID")
        return(pseudo_bulk_profile)
    }

    # Function to calculate the Euclidean distance
    calculate_divergence_score <- function(disease_expr, pn_expr) {
        distance <- sqrt(sum((disease_expr - pn_expr)^2))
        return(distance)
    }

    complete_workflow_TXd <- function(brain_region, cell_class){
        seu <- readRDS(paste0('/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_',brain_region,'_',cell_class,'_int.rds'))
        ncol(seu)

        ## List of cell types
        Groups <- unique(seu$Group)

        for (group in Groups){
        
            ## Subset the Seurat object for disease and PN groups
            seurat_obj_disease <- subset(seu, Group == group)
            seurat_obj_pn <- subset(seu, Group == "PN")

            ## List of cell types
            celltypes <- unique(seu$CellType)

            ## Create empty lists to store pseudo-bulk profiles
            pseudo_bulk_disease <- list()
            pseudo_bulk_pn <- list()

            ## Loop through cell types and calculate pseudo-bulk profiles for both disease and PN
            for (celltype in celltypes) {
                pseudo_bulk_disease[[celltype]] <- calculate_pseudo_bulk(seurat_obj_disease, celltype)
                pseudo_bulk_pn[[celltype]] <- calculate_pseudo_bulk(seurat_obj_pn, celltype)
            }

            ## Calculate the divergence score for each cell type
            divergence_scores <- list()


            for (celltype in celltypes) {
                disease_expr <- pseudo_bulk_disease[[celltype]]@assays$RNA@layers$scale.data
                rownames(disease_expr) <- rownames(seurat_obj_disease)
                
                pn_expr <- pseudo_bulk_pn[[celltype]]@assays$RNA@layers$scale.data
                rownames(pn_expr) <- rownames(seurat_obj_pn)
                
                disease_expr <- as.matrix(disease_expr)
                pn_expr <- as.matrix(pn_expr)

                disease_expr <- apply(disease_expr, 1, median)
                pn_expr <- apply(pn_expr, 1, median)

                divergence_scores[[celltype]] <- calculate_divergence_score(disease_expr, pn_expr)
            }

            ## Dataframe
            score_results <- data.frame(
            celltype = celltypes,
            divergence_score = unlist(divergence_scores)
            )
            score_results$Region <- brain_region
            score_results$Group <- group

            fill <- rbind(fill, score_results)


        }
        return(fill)
    }
        
    ###########################
    ## Ex BA4
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA4"
    cell_class = "Ex"
    BA4_Ex_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## In BA4
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA4"
    cell_class = "In"
    BA4_In_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Glia BA4
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA4"
    cell_class = "Glia"
    BA4_Glia_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Vasc BA4
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA4"
    cell_class = "Vasc"
    BA4_Vasc_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Ex BA9
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA9"
    cell_class = "Ex"
    BA9_Ex_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## In BA9
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA9"
    cell_class = "In"
    BA9_In_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Glia BA9
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA9"
    cell_class = "Glia"
    BA9_Glia_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Vasc BA9
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA9"
    cell_class = "Vasc"
    BA9_Vasc_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Bind all data frames
    ###########################
    all_txd <- rbind(BA4_Ex_df, BA4_In_df, BA4_Glia_df, BA4_Vasc_df, BA9_Ex_df, BA9_In_df, BA9_Glia_df, BA9_Vasc_df)
    all_txd <- subset(all_txd, celltype != "fill")

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ DNN classification accuracies (Need to compute this for the accuracy correlation)

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9ALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SFTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
        #cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9ALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SFTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## USE THIS ONE: merge all: only keep TEST accuracies
    ###########################
    ## code
        
        ## Bind and process BA4
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
                ## TEMP: remove all conditions
                BA4_bind_all <- subset(BA4_bind_all, condition != "all_conditions")
        
        ## Bind and process BA9
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
                ## TEMP: remove all conditions
                BA9_bind_all <- subset(BA9_bind_all, condition != "all_conditions")


                ## TEMP: we are going to start with nly the TEST values. --> if not 
                BA4_bind_all_test <- subset(BA4_bind_all, test == "test")
                BA9_bind_all_test <- subset(BA9_bind_all, test == "test")
        
        ## Bind and process BA4 and BA9
        bind_all_accuracy_test <- rbind(BA4_bind_all_test, BA9_bind_all_test)

        
        
    ##

    ###########################
    ## DONT USE THIS ONE: merge all: MEAN across test and validation acccuracies
    ###########################
    ## code 
        ## Bind and process BA4
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
                ## TEMP: remove all conditions
                BA4_bind_all <- subset(BA4_bind_all, condition != "all_conditions")
        
        summary_df_BA4 <- BA4_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(mean_accuracy = mean(accuracy, na.rm = TRUE)) %>% as.data.frame()

        ## Bind and process BA9
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
                ## TEMP: remove all conditions
                BA9_bind_all <- subset(BA9_bind_all, condition != "all_conditions")

        summary_df_BA9 <- BA9_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(mean_accuracy = mean(accuracy, na.rm = TRUE)) %>% as.data.frame()

        
        ## Bind and process BA4 and BA9
        bind_all_accuracy_test <- rbind(summary_df_BA4, summary_df_BA9)   
    ##

    ###########################
    ## DONT USE THIS ONE: merge all: MEDIAN across test and validation acccuracies
    ###########################
    ## code 
        ## Bind and process BA4
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
                ## TEMP: remove all conditions
                BA4_bind_all <- subset(BA4_bind_all, condition != "all_conditions")
        
        summary_df_BA4 <- BA4_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(median_accuracy = median(accuracy, na.rm = TRUE)) %>% as.data.frame()

        ## Bind and process BA9
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
                ## TEMP: remove all conditions
                BA9_bind_all <- subset(BA9_bind_all, condition != "all_conditions")

        summary_df_BA9 <- BA9_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(median_accuracy = median(accuracy, na.rm = TRUE)) %>% as.data.frame()

        ## Bind and process BA4 and BA9
        bind_all_accuracy_test <- rbind(summary_df_BA4, summary_df_BA9)   
    ##



##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Corrrelation and plot
## code
    ## USE THIS ONE: code for test accuracies only
        ## merge TXd with accuracies (test)

        colnames(bind_all_accuracy_test) <- c('accuracy', 'test', 'celltype', 'Group', 'Region')

        all_Txd_accuracy <- merge(all_txd, bind_all_accuracy_test, by = c('celltype', 'Region', 'Group'))

        ## normalize TXd score
        all_Txd_accuracy$normalized_divergence <- (all_Txd_accuracy$divergence_score - min(all_Txd_accuracy$divergence_score)) / 
                                (max(all_Txd_accuracy$divergence_score) - min(all_Txd_accuracy$divergence_score))

        ## plot
        ggplot(all_Txd_accuracy, aes(x = normalized_divergence, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("DNN accuracy") + 
        xlab("Normalized transcriptional\ndivergence from control profiles")  
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)
    ##

    ## DONT USE THIS ONE: code for test and validation accuracies MEAN.
        ## merge TXd with accuracies (test)

        colnames(bind_all_accuracy_test) <- c('celltype', 'Group', 'Region', 'accuracy')

        all_Txd_accuracy <- merge(all_txd, bind_all_accuracy_test, by = c('celltype', 'Region', 'Group'))

        ## normalize TXd score
        all_Txd_accuracy$normalized_divergence <- (all_Txd_accuracy$divergence_score - min(all_Txd_accuracy$divergence_score)) / 
                                (max(all_Txd_accuracy$divergence_score) - min(all_Txd_accuracy$divergence_score))

        ## plot
        ggplot(all_Txd_accuracy, aes(x = normalized_divergence, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("DNN accuracy") + 
        xlab("Normalized transcriptional\ndivergence from control profiles") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)
    ##

    ## DONT USE THIS ONE: code for test and validation accuracies MEDIAN.
        ## merge TXd with accuracies (test)

        colnames(bind_all_accuracy_test) <- c('celltype', 'Group', 'Region', 'accuracy')

        all_Txd_accuracy <- merge(all_txd, bind_all_accuracy_test, by = c('celltype', 'Region', 'Group'))

        ## normalize TXd score
        all_Txd_accuracy$normalized_divergence <- (all_Txd_accuracy$divergence_score - min(all_Txd_accuracy$divergence_score)) / 
                                (max(all_Txd_accuracy$divergence_score) - min(all_Txd_accuracy$divergence_score))

        ## plot
        ggplot(all_Txd_accuracy, aes(x = normalized_divergence, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("DNN accuracy") + 
        xlab("Normalized transcriptional\ndivergence from control profiles") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)
    ##
##


################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################ Major groups


################################################################## ## TXd test
##################################################################
##################################################################
##################################################################
##################################################################
################################################################## Pineda

## TXd clean
    
    ## NOTES: ########################################################################################
    
    #"The TxD score is a quantification of transcriptional dysregulation. 
    #It represents the change in transcriptome-wide gene expression of each 
    #subtype in disease from its respective PN expression profile. 
    #The divergence score is the Euclidean distance between the median disease 
    #and corresponding PN covariate-corrected, pseudo-bulk expression profiles for each cell type."
    
    ## Functions: ######################################################################################## 

    # Create a function to calculate the pseudo-bulk expression profile (median of gene expression per cell type)
    calculate_pseudo_bulk <- function(seurat_obj, celltype) {
        # Subset the Seurat object by celltype
        celltype_data <- seurat_obj
        
        # Calculate pseudo-bulk expression
        pseudo_bulk_profile <- AggregateExpression(celltype_data, return.seurat = TRUE, verbose = FALSE, group.by = "Sample_ID")
        return(pseudo_bulk_profile)
    }

    # Function to calculate the Euclidean distance
    calculate_divergence_score <- function(disease_expr, pn_expr) {
        distance <- sqrt(sum((disease_expr - pn_expr)^2))
        return(distance)
    }

    complete_workflow_TXd <- function(brain_region, cell_class){
        seu <- readRDS(paste0('/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_',brain_region,'_',cell_class,'_int.rds'))
        ncol(seu)

        ## List of cell types
        Groups <- unique(seu$Group)

        for (group in Groups){
        
            ## Subset the Seurat object for disease and PN groups
            seurat_obj_disease <- subset(seu, Group == group)
            seurat_obj_pn <- subset(seu, Group == "PN")

            ## List of cell types
            #celltypes <- unique(seu$CellType)

            ## Create empty lists to store pseudo-bulk profiles
            pseudo_bulk_disease <- list()
            pseudo_bulk_pn <- list()

            ## Loop through cell types and calculate pseudo-bulk profiles for both disease and PN
            #for (celltype in celltypes) {
            pseudo_bulk_disease[[cell_class]] <- calculate_pseudo_bulk(seurat_obj_disease, celltype)
            pseudo_bulk_pn[[cell_class]] <- calculate_pseudo_bulk(seurat_obj_pn, celltype)
            #}

            ## Calculate the divergence score for each cell type
            divergence_scores <- list()


            #for (celltype in celltypes) {
                disease_expr <- pseudo_bulk_disease[[cell_class]]@assays$RNA@layers$scale.data
                rownames(disease_expr) <- rownames(seurat_obj_disease)
                
                pn_expr <- pseudo_bulk_pn[[cell_class]]@assays$RNA@layers$scale.data
                rownames(pn_expr) <- rownames(seurat_obj_pn)
                
                disease_expr <- as.matrix(disease_expr)
                pn_expr <- as.matrix(pn_expr)

                disease_expr <- apply(disease_expr, 1, median)
                pn_expr <- apply(pn_expr, 1, median)

                divergence_scores[[cell_class]] <- calculate_divergence_score(disease_expr, pn_expr)
            #}

            ## Dataframe
            score_results <- data.frame(
            celltype = cell_class,
            divergence_score = unlist(divergence_scores)
            )
            score_results$Region <- brain_region
            score_results$Group <- group

            fill <- rbind(fill, score_results)


        }
        return(fill)
    }
        
    ###########################
    ## Ex BA4
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA4"
    cell_class = "Ex"
    BA4_Ex_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## In BA4
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA4"
    cell_class = "In"
    BA4_In_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Glia BA4
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA4"
    cell_class = "Glia"
    BA4_Glia_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Vasc BA4
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA4"
    cell_class = "Vasc"
    BA4_Vasc_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Ex BA9
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA9"
    cell_class = "Ex"
    BA9_Ex_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## In BA9
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA9"
    cell_class = "In"
    BA9_In_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Glia BA9
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA9"
    cell_class = "Glia"
    BA9_Glia_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Vasc BA9
    ###########################

    ## Reset the filler frame
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")

    ## parameters
    brain_region = "BA9"
    cell_class = "Vasc"
    BA9_Vasc_df <- complete_workflow_TXd(brain_region, cell_class)

    ###########################
    ## Bind all data frames
    ###########################
    all_txd <- rbind(BA4_Ex_df, BA4_In_df, BA4_Glia_df, BA4_Vasc_df, BA9_Ex_df, BA9_In_df, BA9_Glia_df, BA9_Vasc_df)
    all_txd <- subset(all_txd, celltype != "fill")

    ## save file


##



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ DNN classification accuracies (Need to compute this for the accuracy correlation)

## code
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()
        
        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9ALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SFTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4 -- skip for now. 
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()
        
        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9ALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')


        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SFTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9 -- skip for now. 
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_standard_n_CombatSeq_',condition,'_',region,'_',i,'_narval.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## USE THIS ONE: merge all: only keep TEST accuracies
    ###########################
    ## code
        
        ## Bind and process BA4
        #BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4)
                ## TEMP: remove all conditions
                BA4_bind_all <- subset(BA4_bind_all, condition != "all_conditions")
        
        ## Bind and process BA9
        #BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9)
                ## TEMP: remove all conditions
                BA9_bind_all <- subset(BA9_bind_all, condition != "all_conditions")


                ## TEMP: we are going to start with nly the TEST values. --> if not 
                BA4_bind_all_test <- subset(BA4_bind_all, test == "test")
                BA9_bind_all_test <- subset(BA9_bind_all, test == "test")
        
        ## Bind and process BA4 and BA9
        bind_all_accuracy_test <- rbind(BA4_bind_all_test, BA9_bind_all_test)

        
        
    ##


    ###########################
    ## DONT USE THIS ONE: merge all: MEAN across test and validation acccuracies
    ###########################
    ## code 
        ## Bind and process BA4
        #BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4)
                ## TEMP: remove all conditions
                BA4_bind_all <- subset(BA4_bind_all, condition != "all_conditions")
        
        summary_df_BA4 <- BA4_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(mean_accuracy = mean(accuracy, na.rm = TRUE)) %>% as.data.frame()

        ## Bind and process BA9
        #BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9)
                ## TEMP: remove all conditions
                BA9_bind_all <- subset(BA9_bind_all, condition != "all_conditions")

        summary_df_BA9 <- BA9_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(mean_accuracy = mean(accuracy, na.rm = TRUE)) %>% as.data.frame()

        
        ## Bind and process BA4 and BA9
        bind_all_accuracy_test <- rbind(summary_df_BA4, summary_df_BA9)   
    ##

    ###########################
    ## DONT USE THIS ONE: merge all: MEDIAN across test and validation acccuracies
    ###########################
    ## code 
        ## Bind and process BA4
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
                ## TEMP: remove all conditions
                BA4_bind_all <- subset(BA4_bind_all, condition != "all_conditions")
        
        summary_df_BA4 <- BA4_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(median_accuracy = median(accuracy, na.rm = TRUE)) %>% as.data.frame()

        ## Bind and process BA9
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
                ## TEMP: remove all conditions
                BA9_bind_all <- subset(BA9_bind_all, condition != "all_conditions")

        summary_df_BA9 <- BA9_bind_all %>%
            group_by(celltype, condition, region) %>%
            summarize(median_accuracy = median(accuracy, na.rm = TRUE)) %>% as.data.frame()

        ## Bind and process BA4 and BA9
        bind_all_accuracy_test <- rbind(summary_df_BA4, summary_df_BA9)   
    ##

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Corrrelation and plot
## code
    ## USE THIS ONE: code for test accuracies only
        ## merge TXd with accuracies (test)

        colnames(bind_all_accuracy_test) <- c('accuracy', 'test', 'celltype', 'Group', 'Region')

        all_Txd_accuracy <- merge(all_txd, bind_all_accuracy_test, by = c('celltype', 'Region', 'Group'))

        ## normalize TXd score
        all_Txd_accuracy$normalized_divergence <- (all_Txd_accuracy$divergence_score - min(all_Txd_accuracy$divergence_score)) / 
                                (max(all_Txd_accuracy$divergence_score) - min(all_Txd_accuracy$divergence_score))

        ## plot
        ggplot(all_Txd_accuracy, aes(x = normalized_divergence, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("DNN accuracy") + 
        xlab("Normalized transcriptional\ndivergence from control profiles")  
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)

        ## BA4 only
        all_Txd_accuracy_BA4 <- subset(all_Txd_accuracy, Region == "BA4")
        ggplot(all_Txd_accuracy_BA4, aes(x = normalized_divergence, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("DNN accuracy") + 
        xlab("Normalized transcriptional\ndivergence from control profiles")  
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)


        ## BA9 only
        all_Txd_accuracy_BA9 <- subset(all_Txd_accuracy, Region == "BA9")
        ggplot(all_Txd_accuracy_BA9, aes(x = normalized_divergence, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("DNN accuracy") + 
        xlab("Normalized transcriptional\ndivergence from control profiles")  
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)

    ##

    ## DONT USE THIS ONE: code for test and validation accuracies MEAN.
        ## merge TXd with accuracies (test)

        colnames(bind_all_accuracy_test) <- c('celltype', 'Group', 'Region', 'accuracy')

        all_Txd_accuracy <- merge(all_txd, bind_all_accuracy_test, by = c('celltype', 'Region', 'Group'))

        ## normalize TXd score
        all_Txd_accuracy$normalized_divergence <- (all_Txd_accuracy$divergence_score - min(all_Txd_accuracy$divergence_score)) / 
                                (max(all_Txd_accuracy$divergence_score) - min(all_Txd_accuracy$divergence_score))

        ## plot
        ggplot(all_Txd_accuracy, aes(x = normalized_divergence, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("DNN accuracy") + 
        xlab("Normalized transcriptional\ndivergence from control profiles") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)
    ##

    ## DONT USE THIS ONE: code for test and validation accuracies MEDIAN.
        ## merge TXd with accuracies (test)

        colnames(bind_all_accuracy_test) <- c('celltype', 'Group', 'Region', 'accuracy')

        all_Txd_accuracy <- merge(all_txd, bind_all_accuracy_test, by = c('celltype', 'Region', 'Group'))

        ## normalize TXd score
        all_Txd_accuracy$normalized_divergence <- (all_Txd_accuracy$divergence_score - min(all_Txd_accuracy$divergence_score)) / 
                                (max(all_Txd_accuracy$divergence_score) - min(all_Txd_accuracy$divergence_score))

        ## plot
        ggplot(all_Txd_accuracy, aes(x = normalized_divergence, y = accuracy)) + 
        theme_bw() + 
        geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "solid") + # Add line of best fit
        geom_point(colour = "darkgrey") +
        stat_cor(method = "pearson", colour = "black", label.y = 1,) +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.title.x = element_text(face = "bold"),
            axis.title.y = element_text(face = "bold")
        ) +
        ylab("DNN accuracy") + 
        xlab("Normalized transcriptional\ndivergence from control profiles") 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width = 3)
    ##
##








#### OLD CODE ################################################################################################

## TXd temp
    #"The TxD score is a quantification of transcriptional dysregulation. 
    #It represents the change in transcriptome-wide gene expression of each 
    #subtype in disease from its respective PN expression profile. 
    #The divergence score is the Euclidean distance between the median disease 
    #and corresponding PN covariate-corrected, pseudo-bulk expression profiles for each cell type."   
    
    
    ## we can maybe use these (covariat carrected) objects for TXd. 
    #'/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_BA9_',seu_celltype,'_int.rds'

    ##############
    ## Excitatory neurons BA4
    ##############
    
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_BA4_Ex_int.rds')
    ncol(seu)

    
    ## log normalize
    #seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    #seu <- ScaleData(seu)   

    # Assuming your Seurat object is called `seurat_obj`
    # Subset the Seurat object for disease and PN groups
    seurat_obj_disease <- subset(seu, Group == "C9ALS")
    seurat_obj_pn <- subset(seu, Group == "PN")

    # List of cell types (assuming 'celltype' is in the metadata)
    celltypes <- unique(seu$CellType)

    # Create a function to calculate the pseudo-bulk expression profile (median of gene expression per cell type)
    calculate_pseudo_bulk <- function(seurat_obj, celltype) {
        # Subset the Seurat object by celltype
        #celltype_data <- subset(seurat_obj, subset = celltype == celltype)
        celltype_data <- subset(seurat_obj, CellType == celltype)
        
        # Calculate median expression for each gene in the cell type
        
        #pseudo_bulk_profile <- AverageExpression(celltype_data, return.seurat = TRUE, verbose = FALSE, group.by = "Sample_ID")
        #return(pseudo_bulk_profile)

        pseudo_bulk_profile <- AggregateExpression(celltype_data, return.seurat = TRUE, verbose = FALSE, group.by = "Sample_ID")
        return(pseudo_bulk_profile)
    }

   
                ##################### test
                    #celltype = "Astro"
                    #celltype_data <- subset(seurat_obj_disease, CellType == celltype)
                    
                    # Calculate median expression for each gene in the cell type
                    #pseudo_bulk_profile <- AverageExpression(celltype_data, return.seurat = TRUE, verbose = FALSE, group.by = "Sample_ID")
                    #return(pseudo_bulk_profile)
                ##################### test
   

    # Create empty lists to store pseudo-bulk profiles
    pseudo_bulk_disease <- list()
    pseudo_bulk_pn <- list()

    # Loop through cell types and calculate pseudo-bulk profiles for both disease and PN
    for (celltype in celltypes) {
        pseudo_bulk_disease[[celltype]] <- calculate_pseudo_bulk(seurat_obj_disease, celltype)
        pseudo_bulk_pn[[celltype]] <- calculate_pseudo_bulk(seurat_obj_pn, celltype)
    }

                # Calculate TxD score (change in gene expression)
                #txd_scores <- list()

                #for (celltype in celltypes) {
                #    disease_expr <- pseudo_bulk_disease[[celltype]]@assays$RNA@layers$data ## maybe we use scaled data????????????????????
                #    pn_expr <- pseudo_bulk_pn[[celltype]]@assays$RNA@layers$data
                    
                #    # TxD score is the absolute difference in gene expression between disease and PN
                #    txd_scores[[celltype]] <- abs(disease_expr - pn_expr)
                #}

    # Function to calculate the Euclidean distance
    calculate_divergence_score <- function(disease_expr, pn_expr) {
        # Calculate the Euclidean distance between disease and PN profiles
        distance <- sqrt(sum((disease_expr - pn_expr)^2))
        return(distance)
    }

    # Calculate the divergence score for each cell type
    divergence_scores <- list()

    for (celltype in celltypes) {
        disease_expr <- pseudo_bulk_disease[[celltype]]@assays$RNA@layers$scale.data
        rownames(disease_expr) <- rownames(seurat_obj_disease)
        
        pn_expr <- pseudo_bulk_pn[[celltype]]@assays$RNA@layers$scale.data
        rownames(pn_expr) <- rownames(seurat_obj_pn)
        
        disease_expr <- as.matrix(disease_expr)
        pn_expr <- as.matrix(pn_expr)

        ## calculate the median for each gene
        disease_expr <- apply(disease_expr, 1, median)
        pn_expr <- apply(pn_expr, 1, median)

        # Calculate divergence score
        divergence_scores[[celltype]] <- calculate_divergence_score(disease_expr, pn_expr)
    }


    # Create a data frame with the scores
    score_results <- data.frame(
    celltype = celltypes,
    divergence_score = unlist(divergence_scores)
    )

    # View the results
    print(score_results)
##
