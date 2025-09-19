salloc -A def-sfarhan --time=0-4 -c 1 --mem=10g

module load StdEnv/2023
module load r/4.4.0
R

## Load libraries
#library(MAST)
#library(Seurat, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") ##good
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
library(labeling, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(tidyverse, lib="/home/fiorini9/scracth/R")
#install.packages("RColorBrewer", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(readr)
library(reshape2)

library(tidyverse, lib="/home/fiorini9/scracth/R")
#install.packages("RColorBrewer", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(readr)
library(reshape2)
library(ggh4x)
library(Hmisc)


###############################################################################################
## HVGs BA4 
###############################################################################################

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ ALS v FTLD v Control BA4 

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_All_HVGs_combat_LOSO_generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        #write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        #BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 5)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)


        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("ALS", "FTLD", "control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 3 class random BA4

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_3_class_randomized_BA4_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        #length(unique(BA4_bind_all$donor))

        #write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        #BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_min(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "1",
                str_ends(donor, "_2") ~ "2",
                str_ends(donor, "_0") ~ "0",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("1", "2", "0"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##



###############################################################################################
## HVGs BA9
###############################################################################################

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ ALS v FTLD v Control BA9

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_All_HVGs_combat_LOSO_generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        #write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        #BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 5)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)


        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("ALS", "FTLD", "control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 3 class random BA9

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_3_class_randomized_BA9_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        #length(unique(BA4_bind_all$donor))

        #write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        #BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_min(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "1",
                str_ends(donor, "_2") ~ "2",
                str_ends(donor, "_0") ~ "0",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("1", "2", "0"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##


FOCUSING ON THIS BELOW


###############################################################################################
## HVGs BA4 + BA9 together. 
###############################################################################################


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ ALS v FTLD v Control

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_combat_LOSO_generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.6, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)


        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("ALS", "FTLD", "control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Sporadic v C9orf72 v Control

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_combat_LOSO_generalizable_sporadic_C9_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
                str_ends(donor, "_0") ~ "Control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.6, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)


        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("Sporadic", "C9orf72", "Control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf', height = 3, width = 5)
        
    ##  

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ SALS v C9ALS v SFTLD v C9FTLD v Control

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_combat_LOSO_generalizable_subtypes_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_HVG_generalizable_subtype_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "SALS",
                str_ends(donor, "_2") ~ "C9ALS",
                str_ends(donor, "_3") ~ "SFTLD",
                str_ends(donor, "_4") ~ "C9FTLD",
                str_ends(donor, "_0") ~ "Control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        

        summary_stats <- data.frame(summary_stats)


        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.3, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("SALS", "C9ALS",  "SFTLD", "C9FTLD", "Control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)

    
        
    ##  

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Cell-level model performance summary plot

## code
    ###########################
    ## ALS, FTLD, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        group_summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- rbind(overall_summary_stats, group_summary_stats)

        summary_stats$input <- "Diagnosis"
        summary_stats_diagnosis <- summary_stats
        summary_stats_diagnosis$group[summary_stats_diagnosis$group == "control"] <- "Control"
    ##

    ###########################
    ## Sporadic, C9orf72, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        group_summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- rbind(overall_summary_stats, group_summary_stats)

        summary_stats$input <- "Genetic group"
        summary_stats_genetic <- summary_stats
        summary_stats_genetic$group[summary_stats_genetic$group == "control"] <- "Control"
    ##

    
    ###########################
    ## All subtypes
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
                    mutate(group = case_when(
                        str_ends(donor, "_1") ~ "SALS",
                        str_ends(donor, "_2") ~ "C9ALS",
                        str_ends(donor, "_3") ~ "SFTLD",
                        str_ends(donor, "_4") ~ "C9FTLD",
                        str_ends(donor, "_0") ~ "Control",
                        TRUE ~ NA_character_
                    ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        group_summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- rbind(overall_summary_stats, group_summary_stats)

        summary_stats$input <- "Disease group"
        summary_stats_group <- summary_stats
        summary_stats_group$group[summary_stats_group$group == "control"] <- "Control"

    ##


    ###########################
    ## Plot
    ###########################
    ## code
        total_accuracy_combin <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_group)
        
        ## set factor levels
        total_accuracy_combin$input <- factor(total_accuracy_combin$input, levels = c("Diagnosis", "Genetic group", "Disease group"))
        total_accuracy_combin$group <- factor(total_accuracy_combin$group, levels = rev(c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control")))
        total_accuracy_combin$celltype[total_accuracy_combin$celltype == "Rosehip"] <- "Rose"
        total_accuracy_combin$celltype <- factor(total_accuracy_combin$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rose',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        total_accuracy_combin <- total_accuracy_combin %>%
        group_by(input, group) %>%
        tidyr::complete(celltype) %>%
        ungroup()

                                ## plot diagnosis bar
                                total_accuracy_combin_diagnosis <- subset(total_accuracy_combin, input == "Diagnosis")
                                total_accuracy_combin_diagnosis$group <- factor(total_accuracy_combin_diagnosis$group, levels = c("Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control"))

                                ggplot(total_accuracy_combin_diagnosis, aes(x = group, y = mean_accuracy, fill = group, colour = group)) + 
                                geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
                                    position=position_dodge(.9), colour = "black") +
                                geom_bar(stat="identity", position=position_dodge()) +
                                geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                theme_bw() + 
                                theme(
                                    legend.position = "none",
                                    panel.grid = element_blank(),
                                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_text(face = "bold")
                                ) +
                                facet_grid(input~celltype, scales = "free_x") +
                                ylab("Mean LOSO\naccuracy") +
                                scale_fill_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966",
                                                        "Sporadic" = "darkyellow",
                                                        "C9orf72" = "darkpurple")) +
                                scale_colour_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966",
                                                        "Sporadic" = "darkyellow",
                                                        "C9orf72" = "darkpurple")) 
                                #coord_cartesian(ylim=c(0.5,1.0)) 
                                ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)

                                ## plot genetic group
                                total_accuracy_combin_genetic <- subset(total_accuracy_combin, input == "Genetic group")
                                total_accuracy_combin_genetic$group <- factor(total_accuracy_combin_genetic$group, levels = c("Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control"))

                                ggplot(total_accuracy_combin_genetic, aes(x = group, y = mean_accuracy, fill = group, colour = group)) + 
                                geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
                                    position=position_dodge(.9), colour = "black") +
                                geom_bar(stat="identity", position=position_dodge()) +
                                geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                theme_bw() + 
                                theme(
                                    legend.position = "none",
                                    panel.grid = element_blank(),
                                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_text(face = "bold")
                                ) +
                                facet_grid(input~celltype, scales = "free_x") +
                                ylab("Mean LOSO\naccuracy") +
                                scale_fill_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966",
                                                        "Sporadic" = "goldenrod2",
                                                        "C9orf72" = "purple4")) +
                                scale_colour_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966",
                                                        "Sporadic" = "goldenrod2",
                                                        "C9orf72" = "purple4")) 
                                #coord_cartesian(ylim=c(0.5,1.0)) 
                                ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)


                                ## plot disease group bar
                                total_accuracy_combin_group <- subset(total_accuracy_combin, input == "Disease group")
                                total_accuracy_combin_group$group <- factor(total_accuracy_combin_group$group, levels = c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Control"))

                                mean(total_accuracy_combin_group$mean_accuracy, na.rm = TRUE)

                                ggplot(total_accuracy_combin_group, aes(x = group, y = mean_accuracy, fill = group, colour = group)) + 
                                geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
                                    position=position_dodge(.9), colour = "black") +
                                geom_bar(stat="identity", position=position_dodge()) +
                                geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                theme_bw() + 
                                theme(
                                    legend.position = "none",
                                    panel.grid = element_blank(),
                                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_text(face = "bold")
                                ) +
                                facet_grid(input~celltype, scales = "free_x") +
                                ylab("Mean LOSO\naccuracy") +
                                scale_fill_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966")) +
                                scale_colour_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966")) 
                                #coord_cartesian(ylim=c(0.5,1.0)) 
                                ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)


        ## plot heatmap
        ggplot(total_accuracy_combin, aes(x = celltype, y = group, fill = mean_accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(
            colour = "black",
            size = 8,
            angle = 45,
            hjust = 1,   # aligns the end of text with tick
            vjust = 1    # pushes text closer to axis
            ),
            #axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour = "black", size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested(input ~ ., scales = "free", space = "free") +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3.55, width = 5)

        ## Save final file
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_final_performance_cell_level_whole.csv')
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_final_performance_cell_level_top.csv')
  
    ##


##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Sample-level model performance summary plot

## code
    ###########################
    ## ALS, FTLD, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        # overall donor-wise accuracy:
        sample_acc_overall <- donor_acc %>%
        group_by(celltype) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct), 
            n_donors = n()
        )

        sample_acc_overall$group <- "Overall"

        sample_acc_overall <- sample_acc_overall %>% dplyr::select(celltype, group,   accuracy, sd_accuracy, n_donors)

        # group donor-wise accuracy:
        group_acc_sample <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct),
            n_donors = n()
        )

        summary_stats_donor <- rbind(sample_acc_overall, group_acc_sample)

        summary_stats_donor$input <- "Diagnosis"
        summary_stats_diagnosis <- summary_stats_donor
        summary_stats_diagnosis$group[summary_stats_diagnosis$group == "control"] <- "Control"
    ##

    ###########################
    ## Sporadic, C9orf72, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        # overall donor-wise accuracy:
        sample_acc_overall <- donor_acc %>%
        group_by(celltype) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct), 
            n_donors = n()
        )

        sample_acc_overall$group <- "Overall"

        sample_acc_overall <- sample_acc_overall %>% dplyr::select(celltype, group,   accuracy, sd_accuracy, n_donors)

        # group donor-wise accuracy:
        group_acc_sample <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct),
            n_donors = n()
        )

        summary_stats_donor <- rbind(sample_acc_overall, group_acc_sample)

        summary_stats_donor$input <- "Genetic group"
        summary_stats_genetic <- summary_stats_donor
        summary_stats_genetic$group[summary_stats_genetic$group == "control"] <- "Control"
    ##

    ###########################
    ## All subtypes
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
                    mutate(group = case_when(
                        str_ends(donor, "_1") ~ "SALS",
                        str_ends(donor, "_2") ~ "C9ALS",
                        str_ends(donor, "_3") ~ "SFTLD",
                        str_ends(donor, "_4") ~ "C9FTLD",
                        str_ends(donor, "_0") ~ "Control",
                        TRUE ~ NA_character_
                    ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        
        ############ Sample level. 
        
        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.3, 1, 0),
            .groups = "drop"
        )

        # overall donor-wise accuracy:
        sample_acc_overall <- donor_acc %>%
        group_by(celltype) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct), 
            n_donors = n()
        )

        sample_acc_overall$group <- "Overall"

        sample_acc_overall <- sample_acc_overall %>% dplyr::select(celltype, group,   accuracy, sd_accuracy, n_donors)

        # group donor-wise accuracy:
        group_acc_sample <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct),
            n_donors = n()
        )

        summary_stats_donor <- rbind(sample_acc_overall, group_acc_sample)

        summary_stats_donor$input <- "Disease group"
        summary_stats_group <- summary_stats_donor
        summary_stats_group$group[summary_stats_group$group == "control"] <- "Control"

    ##
        
    ###########################
    ## Plot
    ###########################
    ## code
        total_accuracy_combin <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_group)
        
        ## set factor levels
        total_accuracy_combin$input <- factor(total_accuracy_combin$input, levels = c("Diagnosis", "Genetic group", "Disease group"))
        total_accuracy_combin$group <- factor(total_accuracy_combin$group, levels = rev(c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control")))
        total_accuracy_combin$celltype[total_accuracy_combin$celltype == "Rosehip"] <- "Rose"
        total_accuracy_combin$celltype <- factor(total_accuracy_combin$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rose',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        total_accuracy_combin <- total_accuracy_combin %>%
        group_by(input, group) %>%
        tidyr::complete(celltype) %>%
        ungroup()

                                    ## plot diagnosis bar
                                    total_accuracy_combin_diagnosis <- subset(total_accuracy_combin, input == "Diagnosis")
                                    total_accuracy_combin_diagnosis$group <- factor(total_accuracy_combin_diagnosis$group, levels = c("Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control"))

                                    ggplot(total_accuracy_combin_diagnosis, aes(x = group, y = accuracy, fill = group, colour = group)) + 
                                    geom_errorbar(aes(ymin=accuracy - sd_accuracy, ymax=accuracy + sd_accuracy), width=.2,
                                        position=position_dodge(.9), colour = "black") +
                                    geom_bar(stat="identity", position=position_dodge()) +
                                    geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                    theme_bw() + 
                                    theme(
                                        legend.position = "none",
                                        panel.grid = element_blank(),
                                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                        axis.title.x = element_blank(),
                                        axis.title.y = element_text(face = "bold")
                                    ) +
                                    facet_grid(input~celltype, scales = "free_x") +
                                    ylab("Mean LOSO\naccuracy") +
                                    scale_fill_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966",
                                                            "Sporadic" = "darkyellow",
                                                            "C9orf72" = "darkpurple")) +
                                    scale_colour_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966",
                                                            "Sporadic" = "darkyellow",
                                                            "C9orf72" = "darkpurple")) 
                                    #coord_cartesian(ylim=c(0.5,1.0)) 
                                    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)

                                    ## plot genetic group
                                    total_accuracy_combin_genetic <- subset(total_accuracy_combin, input == "Genetic group")
                                    total_accuracy_combin_genetic$group <- factor(total_accuracy_combin_genetic$group, levels = c("Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control"))

                                    ggplot(total_accuracy_combin_genetic, aes(x = group, y = accuracy, fill = group, colour = group)) + 
                                    geom_errorbar(aes(ymin=0, ymax=accuracy + sd_accuracy), width=.2,
                                        position=position_dodge(.9), colour = "black") +
                                    geom_bar(stat="identity", position=position_dodge()) +
                                    geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                    theme_bw() + 
                                    theme(
                                        legend.position = "none",
                                        panel.grid = element_blank(),
                                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                        axis.title.x = element_blank(),
                                        axis.title.y = element_text(face = "bold")
                                    ) +
                                    facet_grid(input~celltype, scales = "free_x") +
                                    ylab("Mean LOSO\naccuracy") +
                                    scale_fill_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966",
                                                            "Sporadic" = "goldenrod2",
                                                            "C9orf72" = "purple4")) +
                                    scale_colour_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966",
                                                            "Sporadic" = "goldenrod2",
                                                            "C9orf72" = "purple4")) 
                                    #coord_cartesian(ylim=c(0.5,1.0)) 
                                    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)


                                    ## plot disease group bar
                                    total_accuracy_combin_group <- subset(total_accuracy_combin, input == "Disease group")
                                    total_accuracy_combin_group$group <- factor(total_accuracy_combin_group$group, levels = c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Control"))

                                    mean(total_accuracy_combin_group$accuracy, na.rm = TRUE)

                                    ggplot(total_accuracy_combin_group, aes(x = group, y = accuracy, fill = group, colour = group)) + 
                                    geom_errorbar(aes(ymin=0, ymax=accuracy + sd_accuracy), width=.2,
                                        position=position_dodge(.9), colour = "black") +
                                    geom_bar(stat="identity", position=position_dodge()) +
                                    geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                    theme_bw() + 
                                    theme(
                                        legend.position = "none",
                                        panel.grid = element_blank(),
                                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                        axis.title.x = element_blank(),
                                        axis.title.y = element_text(face = "bold")
                                    ) +
                                    facet_grid(input~celltype, scales = "free_x") +
                                    ylab("Mean LOSO\naccuracy") +
                                    scale_fill_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966")) +
                                    scale_colour_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966")) 
                                    #coord_cartesian(ylim=c(0.5,1.0)) 
                                    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)

        
        ## plot heatmap
        ggplot(total_accuracy_combin, aes(x = celltype, y = group, fill = accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(accuracy, 2)), size = 3, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(
            colour = "black",
            size = 9.5,
            angle = 45,
            hjust = 1,   # aligns the end of text with tick
            vjust = 1    # pushes text closer to axis
            ),
            #axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour = "black", size = 9.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested(input ~ ., scales = "free", space = "free") +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3.55, width = 8)

        ## Save final file
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_final_performance_sample_level.csv')
        
        #total_accuracy_combin_diagnosis <- subset(total_accuracy_combin, input == "Diagnosis")
        #mean(total_accuracy_combin_diagnosis$accuracy, na.rm = T) #0.94
        #total_accuracy_combin_group <- subset(total_accuracy_combin, input == "Disease group")
        #mean(total_accuracy_combin_group$accuracy, na.rm = T) #0.94
  
    ##


##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Sample-level aggregate performance across cell types

## code
    ###########################
    ## ALS, FTLD, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')

        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(sample_mean_accuracy),        
                .groups = "drop"               
            )

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 

        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        donor_acc <- data.frame(donor_acc)

        ## Calculate donor aggregate across cell types
        donor_perf <- donor_acc %>%
        group_by(donor, group) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf <- donor_perf %>%
        group_by(group) %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf <- data.frame(group_perf)

        # overall donor-wise accuracy:
        donor_perf <- donor_acc %>%
        group_by(donor) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf_overall <- donor_perf %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf_overall <- data.frame(group_perf_overall)

        group_perf_overall$group <- "Overall"

        group_perf_overall <- group_perf_overall %>% dplyr::select(group, prop_donors_correct, n_donors)

        summary_stats_donor <- rbind(group_perf, group_perf_overall)

        summary_stats_donor$input <- "Diagnosis"
        summary_stats_diagnosis <- summary_stats_donor
        summary_stats_diagnosis$group[summary_stats_diagnosis$group == "control"] <- "Control"
    ##

    ###########################
    ## Sporadic, C9orf72, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(sample_mean_accuracy),        
                .groups = "drop"               
            )

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        donor_acc <- data.frame(donor_acc)

        ## Calculate donor aggregate across cell types
        donor_perf <- donor_acc %>%
        group_by(donor, group) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf <- donor_perf %>%
        group_by(group) %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf <- data.frame(group_perf)

        # overall donor-wise accuracy:
        donor_perf <- donor_acc %>%
        group_by(donor) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf_overall <- donor_perf %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf_overall <- data.frame(group_perf_overall)

        group_perf_overall$group <- "Overall"

        group_perf_overall <- group_perf_overall %>% dplyr::select(group, prop_donors_correct, n_donors)
            
        summary_stats_donor <- rbind(group_perf, group_perf_overall)

        summary_stats_donor$input <- "Genetic group"
        summary_stats_genetic <- summary_stats_donor
        summary_stats_genetic$group[summary_stats_genetic$group == "control"] <- "Control"
    ##

    ###########################
    ## All subtypes
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(sample_mean_accuracy),        
                .groups = "drop"               
            )

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "SALS",
                str_ends(donor, "_2") ~ "C9ALS",
                str_ends(donor, "_3") ~ "SFTLD",
                str_ends(donor, "_4") ~ "C9FTLD",
                str_ends(donor, "_0") ~ "Control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.3, 1, 0),
            .groups = "drop"
        )

        donor_acc <- data.frame(donor_acc)

        ## Calculate donor aggregate across cell types
        donor_perf <- donor_acc %>%
        group_by(donor, group) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf <- donor_perf %>%
        group_by(group) %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf <- data.frame(group_perf)

        # overall donor-wise accuracy:
        donor_perf <- donor_acc %>%
        group_by(donor) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf_overall <- donor_perf %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf_overall <- data.frame(group_perf_overall)

        group_perf_overall$group <- "Overall"

        group_perf_overall <- group_perf_overall %>% dplyr::select(group, prop_donors_correct, n_donors)
            
        summary_stats_donor <- rbind(group_perf, group_perf_overall)

        summary_stats_donor$input <- "Disease group"
        summary_stats_group <- summary_stats_donor
        summary_stats_group$group[summary_stats_group$group == "control"] <- "Control"

    ##
        
    ###########################
    ## Plot
    ###########################
    ## code
        total_accuracy_combin <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_group)
        
        ## set factor levels
        total_accuracy_combin$input <- factor(total_accuracy_combin$input, levels = c("Diagnosis", "Genetic group", "Disease group"))
        
        total_accuracy_combin$group_info <- paste0(total_accuracy_combin$group, " (n = ", total_accuracy_combin$n_donors, ")" )

        total_accuracy_combin$group_info <- factor(total_accuracy_combin$group_info, levels = rev(c( "Overall (n = 96)", "ALS (n = 42)", "SALS (n = 20)", "C9ALS (n = 22)", 
                                                                                                     "FTLD (n = 29)", "SFTLD (n = 13)", "C9FTLD (n = 16)", "Sporadic (n = 33)", "C9orf72 (n = 38)", "Control (n = 25)")))
        
        ## plot heatmap
        ggplot(total_accuracy_combin, aes(x = "Aggregate", y = group_info, fill = prop_donors_correct)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(prop_donors_correct, 2)), size = 3, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(
            colour = "black",
            size = 9.5,
            angle = 45,
            hjust = 1,   # aligns the end of text with tick
            vjust = 1    # pushes text closer to axis
            ),
            #axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour = "black", size = 9.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested(input ~ ., scales = "free", space = "free") +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3.55, width = 2)
        
        #total_accuracy_combin_diagnosis <- subset(total_accuracy_combin, input == "Diagnosis")
        #mean(total_accuracy_combin_diagnosis$accuracy, na.rm = T) #0.94
        #total_accuracy_combin_group <- subset(total_accuracy_combin, input == "Disease group")
        #mean(total_accuracy_combin_group$accuracy, na.rm = T) #0.94
  
    ##


##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ randomized 3 class

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_3_class_randomized_CombatSeq_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_3_class_randomized.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_3_class_randomized.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_min(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "1",
                str_ends(donor, "_2") ~ "2",
                str_ends(donor, "_0") ~ "0",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("1", "2", "0"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)

        mean(summary_stats$mean_accuracy)
        
    ##  

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ randomized 5 class

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5_class_randomized_CombatSeq_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_5_class_randomized.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_5_class_randomized.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_min(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "1",
                str_ends(donor, "_2") ~ "2",
                str_ends(donor, "_3") ~ "3",
                str_ends(donor, "_4") ~ "4",
                str_ends(donor, "_0") ~ "0",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)


        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("1", "2", "3", "4", "0"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)

        mean(summary_stats$mean_accuracy)



        
    ##  

##

###############################################################################################
## NMF BA4 + BA9 together. 
###############################################################################################

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ ALS v FTLD v Control

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_combat_LOSO_generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.6, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)


        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("ALS", "FTLD", "control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Sporadic v C9orf72 v Control

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_combat_LOSO_generalizable_sporadic_C9_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
                str_ends(donor, "_0") ~ "Control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.6, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)


        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("Sporadic", "C9orf72", "Control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf', height = 3, width = 5)
        
    ##  

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ SALS v C9ALS v SFTLD v C9FTLD v Control

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_combat_LOSO_generalizable_subtypes_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_HVG_generalizable_subtype_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "SALS",
                str_ends(donor, "_2") ~ "C9ALS",
                str_ends(donor, "_3") ~ "SFTLD",
                str_ends(donor, "_4") ~ "C9FTLD",
                str_ends(donor, "_0") ~ "Control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        

        summary_stats <- data.frame(summary_stats)


        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.3, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("SALS", "C9ALS",  "SFTLD", "C9FTLD", "Control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
 
    ##  

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Cell-level model performance summary plot

## code
    ###########################
    ## ALS, FTLD, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        group_summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- rbind(overall_summary_stats, group_summary_stats)

        summary_stats$input <- "Diagnosis"
        summary_stats_diagnosis <- summary_stats
        summary_stats_diagnosis$group[summary_stats_diagnosis$group == "control"] <- "Control"
    ##

    ###########################
    ## Sporadic, C9orf72, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        group_summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- rbind(overall_summary_stats, group_summary_stats)

        summary_stats$input <- "Genetic group"
        summary_stats_genetic <- summary_stats
        summary_stats_genetic$group[summary_stats_genetic$group == "control"] <- "Control"
    ##

    
    ###########################
    ## All subtypes
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
                    mutate(group = case_when(
                        str_ends(donor, "_1") ~ "SALS",
                        str_ends(donor, "_2") ~ "C9ALS",
                        str_ends(donor, "_3") ~ "SFTLD",
                        str_ends(donor, "_4") ~ "C9FTLD",
                        str_ends(donor, "_0") ~ "Control",
                        TRUE ~ NA_character_
                    ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        group_summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- rbind(overall_summary_stats, group_summary_stats)

        summary_stats$input <- "Disease group"
        summary_stats_group <- summary_stats
        summary_stats_group$group[summary_stats_group$group == "control"] <- "Control"

    ##

    ###########################
    ## Plot
    ###########################
    ## code
        total_accuracy_combin <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_group)
        
        ## set factor levels
        total_accuracy_combin$input <- factor(total_accuracy_combin$input, levels = c("Diagnosis", "Genetic group", "Disease group"))
        total_accuracy_combin$group <- factor(total_accuracy_combin$group, levels = rev(c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control")))
        total_accuracy_combin$celltype[total_accuracy_combin$celltype == "Rosehip"] <- "Rose"
        total_accuracy_combin$celltype <- factor(total_accuracy_combin$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rose',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        total_accuracy_combin <- total_accuracy_combin %>%
        group_by(input, group) %>%
        tidyr::complete(celltype) %>%
        ungroup()

                                ## plot diagnosis bar
                                total_accuracy_combin_diagnosis <- subset(total_accuracy_combin, input == "Diagnosis")
                                total_accuracy_combin_diagnosis$group <- factor(total_accuracy_combin_diagnosis$group, levels = c("Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control"))

                                ggplot(total_accuracy_combin_diagnosis, aes(x = group, y = mean_accuracy, fill = group, colour = group)) + 
                                geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
                                    position=position_dodge(.9), colour = "black") +
                                geom_bar(stat="identity", position=position_dodge()) +
                                geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                theme_bw() + 
                                theme(
                                    legend.position = "none",
                                    panel.grid = element_blank(),
                                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_text(face = "bold")
                                ) +
                                facet_grid(input~celltype, scales = "free_x") +
                                ylab("Mean LOSO\naccuracy") +
                                scale_fill_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966",
                                                        "Sporadic" = "darkyellow",
                                                        "C9orf72" = "darkpurple")) +
                                scale_colour_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966",
                                                        "Sporadic" = "darkyellow",
                                                        "C9orf72" = "darkpurple")) 
                                #coord_cartesian(ylim=c(0.5,1.0)) 
                                ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)

                                ## plot genetic group
                                total_accuracy_combin_genetic <- subset(total_accuracy_combin, input == "Genetic group")
                                total_accuracy_combin_genetic$group <- factor(total_accuracy_combin_genetic$group, levels = c("Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control"))

                                ggplot(total_accuracy_combin_genetic, aes(x = group, y = mean_accuracy, fill = group, colour = group)) + 
                                geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
                                    position=position_dodge(.9), colour = "black") +
                                geom_bar(stat="identity", position=position_dodge()) +
                                geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                theme_bw() + 
                                theme(
                                    legend.position = "none",
                                    panel.grid = element_blank(),
                                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_text(face = "bold")
                                ) +
                                facet_grid(input~celltype, scales = "free_x") +
                                ylab("Mean LOSO\naccuracy") +
                                scale_fill_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966",
                                                        "Sporadic" = "goldenrod2",
                                                        "C9orf72" = "purple4")) +
                                scale_colour_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966",
                                                        "Sporadic" = "goldenrod2",
                                                        "C9orf72" = "purple4")) 
                                #coord_cartesian(ylim=c(0.5,1.0)) 
                                ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)


                                ## plot disease group bar
                                total_accuracy_combin_group <- subset(total_accuracy_combin, input == "Disease group")
                                total_accuracy_combin_group$group <- factor(total_accuracy_combin_group$group, levels = c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Control"))

                                mean(total_accuracy_combin_group$mean_accuracy, na.rm = TRUE)

                                ggplot(total_accuracy_combin_group, aes(x = group, y = mean_accuracy, fill = group, colour = group)) + 
                                geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
                                    position=position_dodge(.9), colour = "black") +
                                geom_bar(stat="identity", position=position_dodge()) +
                                geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                theme_bw() + 
                                theme(
                                    legend.position = "none",
                                    panel.grid = element_blank(),
                                    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                    axis.title.x = element_blank(),
                                    axis.title.y = element_text(face = "bold")
                                ) +
                                facet_grid(input~celltype, scales = "free_x") +
                                ylab("Mean LOSO\naccuracy") +
                                scale_fill_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966")) +
                                scale_colour_manual(values = c("Overall" = "black",
                                                        "ALS" = "darkred",
                                                        "FTLD" = "darkblue",
                                                        "SALS" = "orange",
                                                        "C9ALS"  = "red",
                                                        "SFTLD" = "blue",
                                                        "C9FTLD"  = "purple",
                                                        "Control" = "#339966")) 
                                #coord_cartesian(ylim=c(0.5,1.0)) 
                                ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)


        ## plot heatmap
        ggplot(total_accuracy_combin, aes(x = celltype, y = group, fill = mean_accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(
            colour = "black",
            size = 8,
            angle = 45,
            hjust = 1,   # aligns the end of text with tick
            vjust = 1    # pushes text closer to axis
            ),
            #axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour = "black", size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested(input ~ ., scales = "free", space = "free") +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3.55, width = 5)

        ## Save final file
        #write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_final_performance_cell_level_whole.csv')
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_final_performance_cell_level_top.csv')
  
    ##



##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Sample-level model performance summary plot

## code
    ###########################
    ## ALS, FTLD, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        # overall donor-wise accuracy:
        sample_acc_overall <- donor_acc %>%
        group_by(celltype) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct), 
            n_donors = n()
        )

        sample_acc_overall$group <- "Overall"

        sample_acc_overall <- sample_acc_overall %>% dplyr::select(celltype, group,   accuracy, sd_accuracy, n_donors)

        # group donor-wise accuracy:
        group_acc_sample <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct),
            n_donors = n()
        )

        summary_stats_donor <- rbind(sample_acc_overall, group_acc_sample)

        summary_stats_donor$input <- "Diagnosis"
        summary_stats_diagnosis <- summary_stats_donor
        summary_stats_diagnosis$group[summary_stats_diagnosis$group == "control"] <- "Control"
    ##

    ###########################
    ## Sporadic, C9orf72, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        # overall donor-wise accuracy:
        sample_acc_overall <- donor_acc %>%
        group_by(celltype) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct), 
            n_donors = n()
        )

        sample_acc_overall$group <- "Overall"

        sample_acc_overall <- sample_acc_overall %>% dplyr::select(celltype, group,   accuracy, sd_accuracy, n_donors)

        # group donor-wise accuracy:
        group_acc_sample <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct),
            n_donors = n()
        )

        summary_stats_donor <- rbind(sample_acc_overall, group_acc_sample)

        summary_stats_donor$input <- "Genetic group"
        summary_stats_genetic <- summary_stats_donor
        summary_stats_genetic$group[summary_stats_genetic$group == "control"] <- "Control"
    ##

    ###########################
    ## All subtypes
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
                    mutate(group = case_when(
                        str_ends(donor, "_1") ~ "SALS",
                        str_ends(donor, "_2") ~ "C9ALS",
                        str_ends(donor, "_3") ~ "SFTLD",
                        str_ends(donor, "_4") ~ "C9FTLD",
                        str_ends(donor, "_0") ~ "Control",
                        TRUE ~ NA_character_
                    ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        
        ############ Sample level. 
        
        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.3, 1, 0),
            .groups = "drop"
        )

        # overall donor-wise accuracy:
        sample_acc_overall <- donor_acc %>%
        group_by(celltype) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct), 
            n_donors = n()
        )

        sample_acc_overall$group <- "Overall"

        sample_acc_overall <- sample_acc_overall %>% dplyr::select(celltype, group,   accuracy, sd_accuracy, n_donors)

        # group donor-wise accuracy:
        group_acc_sample <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            sd_accuracy = sd(donor_correct),
            n_donors = n()
        )

        summary_stats_donor <- rbind(sample_acc_overall, group_acc_sample)

        summary_stats_donor$input <- "Disease group"
        summary_stats_group <- summary_stats_donor
        summary_stats_group$group[summary_stats_group$group == "control"] <- "Control"

    ##

    ###########################
    ## Plot
    ###########################
    ## code
        total_accuracy_combin <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_group)
        
        ## set factor levels
        total_accuracy_combin$input <- factor(total_accuracy_combin$input, levels = c("Diagnosis", "Genetic group", "Disease group"))
        total_accuracy_combin$group <- factor(total_accuracy_combin$group, levels = rev(c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control")))
        total_accuracy_combin$celltype[total_accuracy_combin$celltype == "Rosehip"] <- "Rose"
        total_accuracy_combin$celltype <- factor(total_accuracy_combin$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rose',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        total_accuracy_combin <- total_accuracy_combin %>%
        group_by(input, group) %>%
        tidyr::complete(celltype) %>%
        ungroup()

                                    ## plot diagnosis bar
                                    total_accuracy_combin_diagnosis <- subset(total_accuracy_combin, input == "Diagnosis")
                                    total_accuracy_combin_diagnosis$group <- factor(total_accuracy_combin_diagnosis$group, levels = c("Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control"))

                                    ggplot(total_accuracy_combin_diagnosis, aes(x = group, y = accuracy, fill = group, colour = group)) + 
                                    geom_errorbar(aes(ymin=accuracy - sd_accuracy, ymax=accuracy + sd_accuracy), width=.2,
                                        position=position_dodge(.9), colour = "black") +
                                    geom_bar(stat="identity", position=position_dodge()) +
                                    geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                    theme_bw() + 
                                    theme(
                                        legend.position = "none",
                                        panel.grid = element_blank(),
                                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                        axis.title.x = element_blank(),
                                        axis.title.y = element_text(face = "bold")
                                    ) +
                                    facet_grid(input~celltype, scales = "free_x") +
                                    ylab("Mean LOSO\naccuracy") +
                                    scale_fill_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966",
                                                            "Sporadic" = "darkyellow",
                                                            "C9orf72" = "darkpurple")) +
                                    scale_colour_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966",
                                                            "Sporadic" = "darkyellow",
                                                            "C9orf72" = "darkpurple")) 
                                    #coord_cartesian(ylim=c(0.5,1.0)) 
                                    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)

                                    ## plot genetic group
                                    total_accuracy_combin_genetic <- subset(total_accuracy_combin, input == "Genetic group")
                                    total_accuracy_combin_genetic$group <- factor(total_accuracy_combin_genetic$group, levels = c("Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Sporadic", "C9orf72", "Control"))

                                    ggplot(total_accuracy_combin_genetic, aes(x = group, y = accuracy, fill = group, colour = group)) + 
                                    geom_errorbar(aes(ymin=0, ymax=accuracy + sd_accuracy), width=.2,
                                        position=position_dodge(.9), colour = "black") +
                                    geom_bar(stat="identity", position=position_dodge()) +
                                    geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                    theme_bw() + 
                                    theme(
                                        legend.position = "none",
                                        panel.grid = element_blank(),
                                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                        axis.title.x = element_blank(),
                                        axis.title.y = element_text(face = "bold")
                                    ) +
                                    facet_grid(input~celltype, scales = "free_x") +
                                    ylab("Mean LOSO\naccuracy") +
                                    scale_fill_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966",
                                                            "Sporadic" = "goldenrod2",
                                                            "C9orf72" = "purple4")) +
                                    scale_colour_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966",
                                                            "Sporadic" = "goldenrod2",
                                                            "C9orf72" = "purple4")) 
                                    #coord_cartesian(ylim=c(0.5,1.0)) 
                                    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)


                                    ## plot disease group bar
                                    total_accuracy_combin_group <- subset(total_accuracy_combin, input == "Disease group")
                                    total_accuracy_combin_group$group <- factor(total_accuracy_combin_group$group, levels = c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Control"))

                                    mean(total_accuracy_combin_group$accuracy, na.rm = TRUE)

                                    ggplot(total_accuracy_combin_group, aes(x = group, y = accuracy, fill = group, colour = group)) + 
                                    geom_errorbar(aes(ymin=0, ymax=accuracy + sd_accuracy), width=.2,
                                        position=position_dodge(.9), colour = "black") +
                                    geom_bar(stat="identity", position=position_dodge()) +
                                    geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
                                    geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
                                    theme_bw() + 
                                    theme(
                                        legend.position = "none",
                                        panel.grid = element_blank(),
                                        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
                                        axis.title.x = element_blank(),
                                        axis.title.y = element_text(face = "bold")
                                    ) +
                                    facet_grid(input~celltype, scales = "free_x") +
                                    ylab("Mean LOSO\naccuracy") +
                                    scale_fill_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966")) +
                                    scale_colour_manual(values = c("Overall" = "black",
                                                            "ALS" = "darkred",
                                                            "FTLD" = "darkblue",
                                                            "SALS" = "orange",
                                                            "C9ALS"  = "red",
                                                            "SFTLD" = "blue",
                                                            "C9FTLD"  = "purple",
                                                            "Control" = "#339966")) 
                                    #coord_cartesian(ylim=c(0.5,1.0)) 
                                    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)

        
        ## plot heatmap
        ggplot(total_accuracy_combin, aes(x = celltype, y = group, fill = accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(accuracy, 2)), size = 3, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(
            colour = "black",
            size = 9.5,
            angle = 45,
            hjust = 1,   # aligns the end of text with tick
            vjust = 1    # pushes text closer to axis
            ),
            #axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour = "black", size = 9.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested(input ~ ., scales = "free", space = "free") +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3.55, width = 8)

        ## Save final file
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_final_performance_sample_level.csv')
        
        #total_accuracy_combin_diagnosis <- subset(total_accuracy_combin, input == "Diagnosis")
        #mean(total_accuracy_combin_diagnosis$accuracy, na.rm = T) #0.94
        #total_accuracy_combin_group <- subset(total_accuracy_combin, input == "Disease group")
        #mean(total_accuracy_combin_group$accuracy, na.rm = T) #0.94
  
    ##
        
##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Sample-level aggregate performance across cell types

## code
    ###########################
    ## ALS, FTLD, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')

        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(sample_mean_accuracy),        
                .groups = "drop"               
            )

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 

        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        donor_acc <- data.frame(donor_acc)

        ## Calculate donor aggregate across cell types
        donor_perf <- donor_acc %>%
        group_by(donor, group) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf <- donor_perf %>%
        group_by(group) %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf <- data.frame(group_perf)

        # overall donor-wise accuracy:
        donor_perf <- donor_acc %>%
        group_by(donor) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf_overall <- donor_perf %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf_overall <- data.frame(group_perf_overall)

        group_perf_overall$group <- "Overall"

        group_perf_overall <- group_perf_overall %>% dplyr::select(group, prop_donors_correct, n_donors)

        summary_stats_donor <- rbind(group_perf, group_perf_overall)

        summary_stats_donor$input <- "Diagnosis"
        summary_stats_diagnosis <- summary_stats_donor
        summary_stats_diagnosis$group[summary_stats_diagnosis$group == "control"] <- "Control"
    ##

    ###########################
    ## Sporadic, C9orf72, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(sample_mean_accuracy),        
                .groups = "drop"               
            )

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        donor_acc <- data.frame(donor_acc)

        ## Calculate donor aggregate across cell types
        donor_perf <- donor_acc %>%
        group_by(donor, group) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf <- donor_perf %>%
        group_by(group) %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf <- data.frame(group_perf)

        # overall donor-wise accuracy:
        donor_perf <- donor_acc %>%
        group_by(donor) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf_overall <- donor_perf %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf_overall <- data.frame(group_perf_overall)

        group_perf_overall$group <- "Overall"

        group_perf_overall <- group_perf_overall %>% dplyr::select(group, prop_donors_correct, n_donors)
            
        summary_stats_donor <- rbind(group_perf, group_perf_overall)

        summary_stats_donor$input <- "Genetic group"
        summary_stats_genetic <- summary_stats_donor
        summary_stats_genetic$group[summary_stats_genetic$group == "control"] <- "Control"
    ##

    ###########################
    ## All subtypes
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(sample_mean_accuracy),        
                .groups = "drop"               
            )

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "SALS",
                str_ends(donor, "_2") ~ "C9ALS",
                str_ends(donor, "_3") ~ "SFTLD",
                str_ends(donor, "_4") ~ "C9FTLD",
                str_ends(donor, "_0") ~ "Control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.3, 1, 0),
            .groups = "drop"
        )

        donor_acc <- data.frame(donor_acc)

        ## Calculate donor aggregate across cell types
        donor_perf <- donor_acc %>%
        group_by(donor, group) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf <- donor_perf %>%
        group_by(group) %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf <- data.frame(group_perf)

        # overall donor-wise accuracy:
        donor_perf <- donor_acc %>%
        group_by(donor) %>%
        summarise(
            prop_correct = mean(donor_correct),
            donor_correct_agg = if_else(prop_correct > 0.6, 1, 0)
        )

        donor_perf <- data.frame(donor_perf)

        group_perf_overall <- donor_perf %>%
        summarise(
            prop_donors_correct = mean(donor_correct_agg),
            n_donors = n()
        )

        group_perf_overall <- data.frame(group_perf_overall)

        group_perf_overall$group <- "Overall"

        group_perf_overall <- group_perf_overall %>% dplyr::select(group, prop_donors_correct, n_donors)
            
        summary_stats_donor <- rbind(group_perf, group_perf_overall)

        summary_stats_donor$input <- "Disease group"
        summary_stats_group <- summary_stats_donor
        summary_stats_group$group[summary_stats_group$group == "control"] <- "Control"

    ##
        
    ###########################
    ## Plot
    ###########################
    ## code
        total_accuracy_combin <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_group)
        
        ## set factor levels
        total_accuracy_combin$input <- factor(total_accuracy_combin$input, levels = c("Diagnosis", "Genetic group", "Disease group"))
        
        total_accuracy_combin$group_info <- paste0(total_accuracy_combin$group, " (n = ", total_accuracy_combin$n_donors, ")" )

        total_accuracy_combin$group_info <- factor(total_accuracy_combin$group_info, levels = rev(c( "Overall (n = 96)", "ALS (n = 42)", "SALS (n = 20)", "C9ALS (n = 22)", 
                                                                                                     "FTLD (n = 29)", "SFTLD (n = 13)", "C9FTLD (n = 16)", "Sporadic (n = 33)", "C9orf72 (n = 38)", "Control (n = 25)")))
        
        ## plot heatmap
        ggplot(total_accuracy_combin, aes(x = "Aggregate", y = group_info, fill = prop_donors_correct)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(prop_donors_correct, 2)), size = 3, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(
            colour = "black",
            size = 9.5,
            angle = 45,
            hjust = 1,   # aligns the end of text with tick
            vjust = 1    # pushes text closer to axis
            ),
            #axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour = "black", size = 9.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested(input ~ ., scales = "free", space = "free") +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3.55, width = 2)
        
        #total_accuracy_combin_diagnosis <- subset(total_accuracy_combin, input == "Diagnosis")
        #mean(total_accuracy_combin_diagnosis$accuracy, na.rm = T) #0.94
        #total_accuracy_combin_group <- subset(total_accuracy_combin, input == "Disease group")
        #mean(total_accuracy_combin_group$accuracy, na.rm = T) #0.94
  
    ##


##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ randomized 3 class

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_3_class_randomized_CombatSeq_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_min(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "1",
                str_ends(donor, "_2") ~ "2",
                str_ends(donor, "_0") ~ "0",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.6, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)


        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("1", "2", "0"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ randomized 5 class

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_3_class_randomized_CombatSeq_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        length(unique(BA4_bind_all$donor))

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_min(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "1",
                str_ends(donor, "_2") ~ "2",
                str_ends(donor, "_0") ~ "0",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.6, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)


        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("1", "2", "0"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##


###############################################################################################
## HVG: BA4 + BA9 normalized n cells
###############################################################################################

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ ALS v FTLD v Control

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_diagnosis_5000_cell_true_label_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        #length(unique(BA4_bind_all$donor))

        #write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        #BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("ALS", "FTLD", "control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Sporadic v C9 v Control

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_genetic_5000_cell_true_label_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        #length(unique(BA4_bind_all$donor))

        #write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        #BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("ALS", "FTLD", "control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ SALS v C9ALS v SFTLD v C9FTLD v Control

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_disease_5000_cell_true_label_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        #length(unique(BA4_bind_all$donor))

        #write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        #BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "SALS",
                str_ends(donor, "_2") ~ "C9ALS",
                str_ends(donor, "_3") ~ "SFTLD",
                str_ends(donor, "_4") ~ "C9FTLD",
                str_ends(donor, "_0") ~ "Control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("SALS", "C9ALS", "SFTLD", "C9FTLD", "Control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
        
    ##  

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 3 class random

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_3_class_5000_cell_random_label_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation_1 <- Pineda_validation[,-1]
    ##

    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_3_class_randomized_CombatSeq_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation_2 <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        
        BA4_bind_all <- rbind(Pineda_validation_1, Pineda_validation_2)

        #length(unique(BA4_bind_all$donor))

        #write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        #BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_min(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "1",
                str_ends(donor, "_2") ~ "2",
                str_ends(donor, "_0") ~ "0",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("1", "2", "0"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)

        mean(summary_stats$mean_accuracy)
        
    ##  

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 5 class random

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5_class_5000_cell_random_label_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation_1 <- Pineda_validation[,-1]
    ##

    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_3_class_randomized_CombatSeq_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation_2 <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        
        #BA4_bind_all <- rbind(Pineda_validation_1, Pineda_validation_2)
        BA4_bind_all <- rbind(Pineda_validation_1)

        #length(unique(BA4_bind_all$donor))

        #write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        #BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_min(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "1",
                str_ends(donor, "_2") ~ "2",
                str_ends(donor, "_3") ~ "3",
                str_ends(donor, "_4") ~ "4",
                str_ends(donor, "_0") ~ "0",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("1", "2", "3", "4", "0"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)

        mean(summary_stats$mean_accuracy)
        
    ##  

##

###############################################################################################
## NMF: BA4 + BA9 normalized n cells
###############################################################################################


NOT USING BEYOND THIS POINT. 






################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 100 NMF fix: ALS vs FTLD vs Control with combat

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_combat_LOSO_generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/100_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/100_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.6, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)
        ## This looks really good once again. 



        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("ALS", "FTLD", "control"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC')))
        summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = group, y = celltype, fill = mean_accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(mean_accuracy, 2)), size = 2.5, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour = c("black")),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background =element_rect(fill="white", colour = "white")
        ) +
        facet_grid(. ~ run) +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        #scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp3.pdf', height = 4, width = 6)
        
    ##  

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 100 NMF fix: SALS vs C9ALS vs SFTLD vs C9FTLD vs Control with combat

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_combat_LOSO_generalizable_subtypes_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
        
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/100_NMF_fix_LOSO_combat_generalizable_subtype_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/100_NMF_fix_LOSO_combat_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "SALS",
                str_ends(donor, "_2") ~ "C9ALS",
                str_ends(donor, "_3") ~ "SFTLD",
                str_ends(donor, "_4") ~ "C9FTLD",
                str_ends(donor, "_0") ~ "Control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        

        summary_stats <- data.frame(summary_stats)


        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)

        ## THIS LOOKS REALLY GOOD. 

        ## Sample-wise accuracy when considering all cell types together??


        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = c("ALS", "FTLD", "control"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC')))
        summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = group, y = celltype, fill = mean_accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(mean_accuracy, 2)), size = 2.5, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
            axis.text.y = element_text(colour = c("black")),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background =element_rect(fill="white", colour = "white")
        ) +
        facet_grid(. ~ run) +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"))   # Use the "RdYlBu" palette for continuous data
        #scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        #scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp3.pdf', height = 4, width = 6)
        
    ##  

##



################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 100 NMF: Cell-level model performance summary plot

## code
    ###########################
    ## ALS, FTLD, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/100_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        group_summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- rbind(overall_summary_stats, group_summary_stats)

        summary_stats$input <- "Diagnosis"
        summary_stats_diagnosis <- summary_stats
        summary_stats_diagnosis$group[summary_stats_diagnosis$group == "control"] <- "Control"
    ##

    ###########################
    ## All subtypes
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/100_NMF_fix_LOSO_combat_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
                    mutate(group = case_when(
                        str_ends(donor, "_1") ~ "SALS",
                        str_ends(donor, "_2") ~ "C9ALS",
                        str_ends(donor, "_3") ~ "SFTLD",
                        str_ends(donor, "_4") ~ "C9FTLD",
                        str_ends(donor, "_0") ~ "Control",
                        TRUE ~ NA_character_
                    ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        group_summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- rbind(overall_summary_stats, group_summary_stats)

        summary_stats$input <- "Disease group"
        summary_stats_group <- summary_stats
        summary_stats_group$group[summary_stats_group$group == "control"] <- "Control"

    ##
        
    ###########################
    ## Plot
    ###########################
    ## code
        total_accuracy_combin <- rbind(summary_stats_diagnosis, summary_stats_group)
        
        ## set factor levels
        total_accuracy_combin$input <- factor(total_accuracy_combin$input, levels = c("Diagnosis", "Disease group"))
        total_accuracy_combin$group <- factor(total_accuracy_combin$group, levels = rev(c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Control")))
        total_accuracy_combin$celltype[total_accuracy_combin$celltype == "Rosehip"] <- "Rose"
        total_accuracy_combin$celltype <- factor(total_accuracy_combin$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rose',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        total_accuracy_combin <- total_accuracy_combin %>%
        group_by(input, group) %>%
        tidyr::complete(celltype) %>%
        ungroup()

        ## plot diagnosis bar
        total_accuracy_combin_diagnosis <- subset(total_accuracy_combin, input == "Diagnosis")
        total_accuracy_combin_diagnosis$group <- factor(total_accuracy_combin_diagnosis$group, levels = c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Control"))

        mean(total_accuracy_combin_diagnosis$mean_accuracy, na.rm = TRUE)

        ggplot(total_accuracy_combin_diagnosis, aes(x = group, y = mean_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(input~celltype, scales = "free_x") +
        ylab("Mean LOSO\naccuracy") +
        scale_fill_manual(values = c("Overall" = "black",
                                "ALS" = "darkred",
                                "FTLD" = "darkblue",
                                "SALS" = "orange",
                                "C9ALS"  = "red",
                                "SFTLD" = "blue",
                                "C9FTLD"  = "purple",
                                "Control" = "#339966")) +
        scale_colour_manual(values = c("Overall" = "black",
                                "ALS" = "darkred",
                                "FTLD" = "darkblue",
                                "SALS" = "orange",
                                "C9ALS"  = "red",
                                "SFTLD" = "blue",
                                "C9FTLD"  = "purple",
                                "Control" = "#339966")) 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)

        ## plot disease group bar
        total_accuracy_combin_group <- subset(total_accuracy_combin, input == "Disease group")
        total_accuracy_combin_group$group <- factor(total_accuracy_combin_group$group, levels = c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Control"))

        mean(total_accuracy_combin_group$mean_accuracy, na.rm = TRUE)

        ggplot(total_accuracy_combin_group, aes(x = group, y = mean_accuracy, fill = group, colour = group)) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey") +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(input~celltype, scales = "free_x") +
        ylab("Mean LOSO\naccuracy") +
        scale_fill_manual(values = c("Overall" = "black",
                                "ALS" = "darkred",
                                "FTLD" = "darkblue",
                                "SALS" = "orange",
                                "C9ALS"  = "red",
                                "SFTLD" = "blue",
                                "C9FTLD"  = "purple",
                                "Control" = "#339966")) +
        scale_colour_manual(values = c("Overall" = "black",
                                "ALS" = "darkred",
                                "FTLD" = "darkblue",
                                "SALS" = "orange",
                                "C9ALS"  = "red",
                                "SFTLD" = "blue",
                                "C9FTLD"  = "purple",
                                "Control" = "#339966")) 
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2, width = 13)


        ## plot heatmap
        ggplot(total_accuracy_combin, aes(x = celltype, y = group, fill = mean_accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(
            colour = "black",
            size = 8,
            angle = 45,
            hjust = 1,   # aligns the end of text with tick
            vjust = 1    # pushes text closer to axis
            ),
            #axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour = "black", size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested(input ~ ., scales = "free", space = "free") +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
  
    ##


##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ 100 NMF: Sample-level model performance summary plot

## code
    ###########################
    ## ALS, FTLD, Control
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/100_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        # overall donor-wise accuracy:
        sample_acc_overall <- donor_acc %>%
        group_by(celltype) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        sample_acc_overall$group <- "Overall"

        sample_acc_overall <- sample_acc_overall %>% dplyr::select(celltype, group,   accuracy, n_donors)

        # group donor-wise accuracy:
        group_acc_sample <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        summary_stats_donor <- rbind(sample_acc_overall, group_acc_sample)

        summary_stats_donor$input <- "Diagnosis"
        summary_stats_diagnosis <- summary_stats_donor
        summary_stats_diagnosis$group[summary_stats_diagnosis$group == "control"] <- "Control"
    ##

    ###########################
    ## All subtypes
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/100_NMF_fix_LOSO_combat_generalizable_subtype_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## Calculate cell-type specific overall accuracy
        overall_summary_stats <- BA4_bind_all %>%
            group_by(celltype) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        overall_summary_stats$group <- "Overall"

        overall_summary_stats <- overall_summary_stats %>% dplyr::select(celltype, group, median_accuracy, mean_accuracy, sd_accuracy, mad_accuracy)
                
        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
                    mutate(group = case_when(
                        str_ends(donor, "_1") ~ "SALS",
                        str_ends(donor, "_2") ~ "C9ALS",
                        str_ends(donor, "_3") ~ "SFTLD",
                        str_ends(donor, "_4") ~ "C9FTLD",
                        str_ends(donor, "_0") ~ "Control",
                        TRUE ~ NA_character_
                    ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        
        ############ Sample level. 
        
        ## Sample wise accuracy -- group
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        # overall donor-wise accuracy:
        sample_acc_overall <- donor_acc %>%
        group_by(celltype) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        sample_acc_overall$group <- "Overall"

        sample_acc_overall <- sample_acc_overall %>% dplyr::select(celltype, group,   accuracy, n_donors)

        # group donor-wise accuracy:
        group_acc_sample <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

       summary_stats_donor <- rbind(sample_acc_overall, group_acc_sample)

        summary_stats_donor$input <- "Disease group"
        summary_stats_group <- summary_stats_donor
        summary_stats_group$group[summary_stats_group$group == "control"] <- "Control"

    ##
        
    ###########################
    ## Plot
    ###########################
    ## code
        total_accuracy_combin <- rbind(summary_stats_diagnosis, summary_stats_group)
        
        ## set factor levels
        total_accuracy_combin$input <- factor(total_accuracy_combin$input, levels = c("Diagnosis", "Disease group"))
        total_accuracy_combin$group <- factor(total_accuracy_combin$group, levels = rev(c( "Overall", "ALS", "SALS", "C9ALS", "FTLD", "SFTLD", "C9FTLD", "Control")))
        total_accuracy_combin$celltype[total_accuracy_combin$celltype == "Rosehip"] <- "Rose"
        total_accuracy_combin$celltype <- factor(total_accuracy_combin$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rose',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

        total_accuracy_combin <- total_accuracy_combin %>%
        group_by(input, group) %>%
        tidyr::complete(celltype) %>%
        ungroup()

        ## plot heatmap
        ggplot(total_accuracy_combin, aes(x = celltype, y = group, fill = accuracy)) +
        theme_bw() + 
        geom_tile() +
        geom_text(aes(label = round(accuracy, 2)), size = 2.25, colour = "black") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(
            colour = "black",
            size = 8,
            angle = 45,
            hjust = 1,   # aligns the end of text with tick
            vjust = 1    # pushes text closer to axis
            ),
            #axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour = "black", size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested(input ~ ., scales = "free", space = "free") +
        ylab("Accuracy") +
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)

        total_accuracy_combin_diagnosis <- subset(total_accuracy_combin, input == "Diagnosis")
        mean(total_accuracy_combin_diagnosis$accuracy, na.rm = T) #0.94
        total_accuracy_combin_group <- subset(total_accuracy_combin, input == "Disease group")
        mean(total_accuracy_combin_group$accuracy, na.rm = T) #0.94
  
    ##


##

END FOCUS. 


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ BA4 All HVGs fix: ALS vs FTLD vs Control with combat

## code
    #########
    ## Validation on Pineda
    #########
    ## code
        #'/home/fiorini9/scratch/machine_learning_ALS/model_outs/generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'.csv'
        # List of cell types
        cell_type <- c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",   "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",  "Micro",  "Mural",   "Endo",    "Fibro", "L5")

        # Create an empty list to store dataframes
        df_list <- list()
        
        # Loop through each cell type
        for (i in cell_type) {
        # Construct the file path for the current cell type
        #file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_All_HVGs_combat_LOSO_generalizable_ALS_FTLD_control_Pineda_optimal_gene_set_weighted_LOO_LIME_0_abs_case_control_report_CombatSeq_ALS_',i,'_3.csv')
 
        # Check if the file exists
        if (file.exists(file_path)) {
            # Read the CSV file
            current_df <- read.csv(file_path)
            current_df$run <- "validation"
            
            # Add the dataframe to the list
            df_list[[i]] <- current_df
        } else {
            warning(paste("File does not exist:", file_path))
        }
        }

        # Merge all dataframes into one
        Pineda_validation <- do.call(rbind, df_list)
        Pineda_validation <- Pineda_validation[,-1]
    ##

    ###########################
    ## merge and plot
    ###########################
    ## Mean: code for barplot
        BA4_bind_all <- Pineda_validation

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        #BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region, n_cells) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )
         nrow(BA4_bind_all)

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)

        ## compute mean        
        summary_stats <- BA4_bind_all %>%
            group_by(celltype,group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

    
        summary_stats <- data.frame(summary_stats)

        ## Sample wise accuracy
        donor_acc <- BA4_bind_all %>%
        group_by(donor, celltype, group) %>%
        summarise(
            donor_correct = ifelse(mean(sample_mean_accuracy) > 0.5, 1, 0),
            .groups = "drop"
        )

        # If you want overall donor-wise accuracy per group:
        group_acc <- donor_acc %>%
        group_by(celltype, group) %>%
        summarise(
            accuracy = mean(donor_correct),
            n_donors = n()
        )

        group_acc <- data.frame(group_acc)

        ## set factor levels
        summary_stats$group <- factor(summary_stats$group, levels = (rev(c("ALS", "FTLD", "control"))))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))
        #summary_stats$run <- factor(summary_stats$run, levels = c('validation','evaluation_unbalanced','evaluation_1X','evaluation_2X'))

        ## Plot heatmap
        ggplot(summary_stats, aes(x = celltype, y = group, fill = mean_accuracy)) +
                theme_bw() + 
                geom_tile() +
                geom_text(aes(label = round(mean_accuracy, 2)), size = 2.25, colour = "black") +
                theme(
                    legend.position = "none",
                    panel.grid = element_blank(),
                    axis.text.x = element_text(
                    colour = "black",
                    size = 8,
                    angle = 45,
                    hjust = 1,   # aligns the end of text with tick
                    vjust = 1    # pushes text closer to axis
                    ),
                    #axis.ticks.x = element_blank(),
                    axis.text.y = element_text(colour = "black", size = 8),
                    axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    strip.background = element_rect(fill="lightgrey", colour = "white"),
                    strip.text = element_text(size= 8, face="bold", colour = "black"),
                    strip.placement = "outside"
                ) +
                #facet_nested(input ~ ., scales = "free", space = "free") +
                ylab("Accuracy") +
                scale_x_discrete(expand = c(0,0)) +
                scale_y_discrete(expand = c(0,0)) +
                scale_fill_gradientn(colors = brewer.pal(9, "RdYlGn"), limits = c(0, 1), na.value = "white" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width = 5)
  
        
    ##  

##