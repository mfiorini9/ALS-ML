salloc -A def-sfarhan --time=0-2 -c 1 --mem=50g

module load StdEnv/2023
module load r/4.4.0
R

library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(ggpubr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(tidyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(dplyr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

##########
## Load SALS_MCx
##########
## code ##
    # List of cell types
    cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
    status <- 'SALS'
    brain_region <- 'BA4'

    # Create an empty list to store dataframes
    df_list <- list()

    # Loop through each cell type
    for (i in cell_type) {
    # Construct the file path for the current cell type
    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

    # Check if the file exists
    if (file.exists(file_path)) {
        # Read the CSV file
        current_df <- read.csv(file_path)
        
        # Optionally, add the cell type as a new column to identify the source of each row
        current_df$dataset <- paste0(status, "_", brain_region)
        
        # Add the dataframe to the list
        df_list[[i]] <- current_df
    } else {
        warning(paste("File does not exist:", file_path))
    }
    }

    # Merge all dataframes into one
    final_SALS_MCx <- do.call(rbind, df_list)

##

##########
## Load C9ALS_MCx
##########
## code ##
    # List of cell types
    cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
    status <- 'C9ALS'
    brain_region <- 'BA4'

    # Create an empty list to store dataframes
    df_list <- list()

    # Loop through each cell type
    for (i in cell_type) {
    # Construct the file path for the current cell type
    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

    # Check if the file exists
    if (file.exists(file_path)) {
        # Read the CSV file
        current_df <- read.csv(file_path)
        
        # Optionally, add the cell type as a new column to identify the source of each row
        current_df$dataset <- paste0(status, "_", brain_region)
        
        # Add the dataframe to the list
        df_list[[i]] <- current_df
    } else {
        warning(paste("File does not exist:", file_path))
    }
    }

    # Merge all dataframes into one
    final_C9ALS_MCx <- do.call(rbind, df_list)

##

##########
## Load SALS_FCx
##########
## code ##
    # List of cell types
    cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
    status <- 'SALS'
    brain_region <- 'BA9'

    # Create an empty list to store dataframes
    df_list <- list()

    # Loop through each cell type
    for (i in cell_type) {
    # Construct the file path for the current cell type
    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

    # Check if the file exists
    if (file.exists(file_path)) {
        # Read the CSV file
        current_df <- read.csv(file_path)
        
        # Optionally, add the cell type as a new column to identify the source of each row
        current_df$dataset <- paste0(status, "_", brain_region)
        
        # Add the dataframe to the list
        df_list[[i]] <- current_df
    } else {
        warning(paste("File does not exist:", file_path))
    }
    }

    # Merge all dataframes into one
    final_SALS_FCx <- do.call(rbind, df_list)

##

##########
## Load C9ALS_FCx
##########
## code ##
    # List of cell types
    cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
    status <- 'C9ALS'
    brain_region <- 'BA9'

    # Create an empty list to store dataframes
    df_list <- list()

    # Loop through each cell type
    for (i in cell_type) {
    # Construct the file path for the current cell type
    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

    # Check if the file exists
    if (file.exists(file_path)) {
        # Read the CSV file
        current_df <- read.csv(file_path)
        
        # Optionally, add the cell type as a new column to identify the source of each row
        current_df$dataset <- paste0(status, "_", brain_region)
        
        # Add the dataframe to the list
        df_list[[i]] <- current_df
    } else {
        warning(paste("File does not exist:", file_path))
    }
    }

    # Merge all dataframes into one
    final_C9ALS_FCx <- do.call(rbind, df_list)

##

##########
## Load SFTLD_MCx
##########
## code ##
    # List of cell types
    cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
    status <- 'SFTLD'
    brain_region <- 'BA4'

    # Create an empty list to store dataframes
    df_list <- list()

    # Loop through each cell type
    for (i in cell_type) {
    # Construct the file path for the current cell type
    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

    # Check if the file exists
    if (file.exists(file_path)) {
        # Read the CSV file
        current_df <- read.csv(file_path)
        
        # Optionally, add the cell type as a new column to identify the source of each row
        current_df$dataset <- paste0(status, "_", brain_region)
        
        # Add the dataframe to the list
        df_list[[i]] <- current_df
    } else {
        warning(paste("File does not exist:", file_path))
    }
    }

    # Merge all dataframes into one
    final_SFTLD_MCx <- do.call(rbind, df_list)

##

##########
## Load C9FTLD_MCx
##########
## code ##
    # List of cell types
    cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
    status <- 'C9FTLD'
    brain_region <- 'BA4'

    # Create an empty list to store dataframes
    df_list <- list()

    # Loop through each cell type
    for (i in cell_type) {
    # Construct the file path for the current cell type
    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

    # Check if the file exists
    if (file.exists(file_path)) {
        # Read the CSV file
        current_df <- read.csv(file_path)
        
        # Optionally, add the cell type as a new column to identify the source of each row
        current_df$dataset <- paste0(status, "_", brain_region)
        
        # Add the dataframe to the list
        df_list[[i]] <- current_df
    } else {
        warning(paste("File does not exist:", file_path))
    }
    }

    # Merge all dataframes into one
    final_C9FTLD_MCx <- do.call(rbind, df_list)

##

##########
## Load SFTLD_FCx
##########
## code ##
    # List of cell types
    cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
    status <- 'SFTLD'
    brain_region <- 'BA9'

    # Create an empty list to store dataframes
    df_list <- list()

    # Loop through each cell type
    for (i in cell_type) {
    # Construct the file path for the current cell type
    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

    # Check if the file exists
    if (file.exists(file_path)) {
        # Read the CSV file
        current_df <- read.csv(file_path)
        
        # Optionally, add the cell type as a new column to identify the source of each row
        current_df$dataset <- paste0(status, "_", brain_region)
        
        # Add the dataframe to the list
        df_list[[i]] <- current_df
    } else {
        warning(paste("File does not exist:", file_path))
    }
    }

    # Merge all dataframes into one
    final_SFTLD_FCx <- do.call(rbind, df_list)

##

##########
## Load C9FTLD_FCx
##########
## code ##
    # List of cell types
    cell_type <- c("5HT3aR", 'Astro', 'Endo', 'Fibro', 'L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5_L6', 'L5', 'L6', 'Micro', 'Mural', 'Oligo', 'OPC', 'PV', 'Rosehip', 'SOM',  'T_Cell')
    status <- 'C9FTLD'
    brain_region <- 'BA9'

    # Create an empty list to store dataframes
    df_list <- list()

    # Loop through each cell type
    for (i in cell_type) {
    # Construct the file path for the current cell type
    file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',status,'_',brain_region,'_',i,'_narval_2_LOO.csv')

    # Check if the file exists
    if (file.exists(file_path)) {
        # Read the CSV file
        current_df <- read.csv(file_path)
        
        # Optionally, add the cell type as a new column to identify the source of each row
        current_df$dataset <- paste0(status, "_", brain_region)
        
        # Add the dataframe to the list
        df_list[[i]] <- current_df
    } else {
        warning(paste("File does not exist:", file_path))
    }
    }

    # Merge all dataframes into one
    final_C9FTLD_FCx <- do.call(rbind, df_list)

##





##########
## Bind all and plot MCx
##########
## code
    ########################
    ## Cell-wise accuracy
    ########################

    #total_bind <- rbind(final_df_kamath, final_df_wang, final_df_smajic)

    total_bind <- rbind(final_SALS_MCx, final_C9ALS_MCx, final_SFTLD_MCx, final_C9FTLD_MCx)

    #total_bind$dataset <- factor(total_bind$dataset, levels = c('kamath', 'wang', 'smajic'))
    #total_bind$method <- factor(total_bind$method, levels = c("optimal", "random_HVGs", "random_genes"))

    #total_bind <- subset(total_bind, dataset != 'smajic')
    #total_bind <- subset(total_bind, method != 'random_HVGs')

    total_bind <- total_bind %>%
    group_by(celltype, donor, dataset) %>%
    arrange(desc(test_accuracy)) %>%    # Arrange X0 in descending order
    slice_head(n = 1)  

    ## median bar plot
    summary_stats <- total_bind %>%
    group_by(celltype, dataset) %>%
    summarise(
        median_test_accuracy = median(test_accuracy), 
        mean_test_accuracy = mean(test_accuracy),           # Calculate the median of X0
        sd_test_accuracy = sd(test_accuracy),                  # Calculate the standard deviation of X0
        mad_test_accuracy = median(abs(test_accuracy - median(test_accuracy))),  # Calculate MAD: median of absolute deviations from the median
        .groups = "drop"                 # Drop the grouping after summarising
    )

    
    summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5', 'L5_L6', 'L6', '5HT3aR', 'PV', 'Rosehip', 'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC', 'T_Cell' ))
    summary_stats$dataset <- factor(summary_stats$dataset, levels = c("SALS_BA4", "C9ALS_BA4", "SFTLD_BA4", "C9FTLD_BA4"))

    ggplot(summary_stats, aes(x = dataset, y = median_test_accuracy, fill = dataset)) + 
    geom_bar(stat = "identity", position = "dodge", linewidth = 1) +  # Bar plot
    geom_errorbar(aes(ymin=median_test_accuracy - mad_test_accuracy, ymax=median_test_accuracy + mad_test_accuracy), width=.2,
            position=position_dodge(.9)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values = c("orange", "red", "blue", "purple")) +   
    theme_bw() +
    theme(panel.grid = element_blank(),
        #strip.background =element_rect(fill="white", colour = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = c("black")),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(size = 6, hjust = 0.5),
        legend.position = "none")+
    ylab("Accuracy") + 
    coord_cartesian(ylim=c(0.0,1.01)) +
    facet_wrap(~celltype, ncol =19) +
    scale_x_discrete(labels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))

    ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_test_leave_one_out_bar.pdf',sep=""), height = 3.5, width = 13) 

    ## compute average
    mean(summary_stats$median_test_accuracy)


##

##########
## Bind all and plot FCx
##########
## code
    ########################
    ## Cell-wise accuracy
    ########################

    #total_bind <- rbind(final_df_kamath, final_df_wang, final_df_smajic)

    total_bind <- rbind(final_SALS_FCx, final_C9ALS_FCx, final_SFTLD_FCx, final_C9FTLD_FCx)

    #total_bind$dataset <- factor(total_bind$dataset, levels = c('kamath', 'wang', 'smajic'))
    #total_bind$method <- factor(total_bind$method, levels = c("optimal", "random_HVGs", "random_genes"))

    #total_bind <- subset(total_bind, dataset != 'smajic')
    #total_bind <- subset(total_bind, method != 'random_HVGs')

    total_bind <- total_bind %>%
    group_by(celltype, donor, dataset) %>%
    arrange(desc(test_accuracy)) %>%    # Arrange X0 in descending order
    slice_head(n = 1)  

    ## median bar plot
    summary_stats <- total_bind %>%
    group_by(celltype, dataset) %>%
    summarise(
        median_test_accuracy = median(test_accuracy), 
        mean_test_accuracy = mean(test_accuracy),           # Calculate the median of X0
        sd_test_accuracy = sd(test_accuracy),                  # Calculate the standard deviation of X0
        mad_test_accuracy = median(abs(test_accuracy - median(test_accuracy))),  # Calculate MAD: median of absolute deviations from the median
        .groups = "drop"                 # Drop the grouping after summarising
    )

    
    summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3', 'L3_L5', 'L4_L5', 'L4_L6', 'L5', 'L5_L6', 'L6', '5HT3aR', 'PV', 'Rosehip', 'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC', 'T_Cell' ))
    summary_stats$dataset <- factor(summary_stats$dataset, levels = c("SALS_BA9", "C9ALS_BA9", "SFTLD_BA9", "C9FTLD_BA9"))

    ggplot(summary_stats, aes(x = dataset, y = median_test_accuracy, fill = dataset)) + 
    geom_bar(stat = "identity", position = "dodge", linewidth = 1) +  # Bar plot
    geom_errorbar(aes(ymin=median_test_accuracy - mad_test_accuracy, ymax=median_test_accuracy + mad_test_accuracy), width=.2,
            position=position_dodge(.9)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values = c("orange", "red", "blue", "purple")) +   
    theme_bw() +
    theme(panel.grid = element_blank(),
        #strip.background =element_rect(fill="white", colour = "white"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = c("black")),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(size = 6, hjust = 0.5),
        legend.position = "none")+
    ylab("Accuracy") + 
    coord_cartesian(ylim=c(0.0,1.01)) +
    facet_wrap(~celltype, ncol =19) +
    scale_x_discrete(labels = c("SALS", "C9ALS", "SFTLD", "C9FTLD"))

    ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_test_leave_one_out_bar.pdf',sep=""), height = 3.5, width = 13)

    ## compute average
    mean(summary_stats$median_test_accuracy) 

##