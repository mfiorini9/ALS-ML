salloc -A def-tdurcan --time=0-8 -c 1 --mem=10g

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

Performance with NORMALIZED CELL COUNTS


###############################################################################################
## HVG: BA4 + BA9 normalized n cells
###############################################################################################

##################################
##################################
##################################
##################################
##################################
################################## Diagnosis

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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_diagnosis_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Genetic

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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_genetic_5000_cell_true_label.csv')
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Disease

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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_disease_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Cell-level model performance summary plot

## code
    ###########################
    ## Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_diagnosis_5000_cell_true_label.csv')
        
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
    ## Genetic
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_genetic_5000_cell_true_label.csv')
        
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
    ## Disease
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_disease_5000_cell_true_label.csv')
        
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
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5000_final_performance_cell_level_top.csv')
  
    ##


##



###############################################################################################
## NMF: BA4 + BA9 normalized n cells
###############################################################################################

##################################
##################################
##################################
##################################
##################################
################################## Diagnosis

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_diagnosis_5000_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_diagnosis_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Genetic

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_genetic_5000_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_genetic_5000_cell_true_label.csv')
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Disease

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_disease_5000_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_disease_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Cell-level model performance summary plot

## code
    ###########################
    ## Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_diagnosis_5000_cell_true_label.csv')
        
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
    ## Genetic
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_genetic_5000_cell_true_label.csv')
        
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
    ## Disease
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_disease_5000_cell_true_label.csv')
        
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
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_5000_final_performance_cell_level_top.csv')
    ##


##


###############################################################################################
## HVG: BA4 normalized n cells
###############################################################################################

##################################
##################################
##################################
##################################
##################################
################################## Diagnosis

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5000_diagnosis_BA4_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_diagnosis_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Genetic

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5000_genetic_BA4_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_genetic_5000_cell_true_label.csv')
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Disease

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5000_disease_BA4_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_disease_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Cell-level model performance summary plot

## code
    ###########################
    ## Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_diagnosis_5000_cell_true_label.csv')
        
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
    ## Genetic
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_genetic_5000_cell_true_label.csv')
        
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
    ## Disease
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_disease_5000_cell_true_label.csv')
        
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
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_5000_final_performance_cell_level_top.csv')
    ##


##


###############################################################################################
## NMF: BA4 normalized n cells
###############################################################################################

##################################
##################################
##################################
##################################
##################################
################################## Diagnosis

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_diagnosis_5000_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_diagnosis_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Genetic

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_genetic_5000_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_genetic_5000_cell_true_label.csv')
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Disease

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_disease_5000_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_disease_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Cell-level model performance summary plot

## code
    ###########################
    ## Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_diagnosis_5000_cell_true_label.csv')
        
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
    ## Genetic
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_genetic_5000_cell_true_label.csv')
        
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
    ## Disease
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_disease_5000_cell_true_label.csv')
        
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
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_5000_final_performance_cell_level_top.csv')
    ##


##

###############################################################################################
## HVG: BA9 normalized n cells
###############################################################################################

##################################
##################################
##################################
##################################
##################################
################################## Diagnosis

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5000_diagnosis_BA9_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_diagnosis_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Genetic

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5000_genetic_BA9_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_genetic_5000_cell_true_label.csv')
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Disease

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5000_disease_BA9_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_disease_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Cell-level model performance summary plot

## code
    ###########################
    ## Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_diagnosis_5000_cell_true_label.csv')
        
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
    ## Genetic
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_genetic_5000_cell_true_label.csv')
        
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
    ## Disease
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_disease_5000_cell_true_label.csv')
        
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
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_5000_final_performance_cell_level_top.csv')
    ##


##

###############################################################################################
## NMF: BA9 normalized n cells
###############################################################################################

##################################
##################################
##################################
##################################
##################################
################################## Diagnosis

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_diagnosis_5000_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_diagnosis_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Genetic

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_genetic_5000_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_genetic_5000_cell_true_label.csv')
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Disease

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
        file_path <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_disease_5000_cell_true_label_',i,'_3.csv')
        
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

        write.csv(BA4_bind_all, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_disease_5000_cell_true_label.csv')        
    ##  

##

##################################
##################################
##################################
##################################
##################################
################################## Cell-level model performance summary plot

## code
    ###########################
    ## Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_diagnosis_5000_cell_true_label.csv')
        
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
    ## Genetic
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_genetic_5000_cell_true_label.csv')
        
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
    ## Disease
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_disease_5000_cell_true_label.csv')
        
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
        write.csv(total_accuracy_combin, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_5000_final_performance_cell_level_top.csv')
    ##
##