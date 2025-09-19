salloc -A def-grouleau --time=0-4 -c 1 --mem=200g
salloc -A def-sfarhan --time=0-4 -c 1 --mem=100g
salloc -A def-tdurcan --time=0-4 -c 1 --mem=200g

def-grouleau

module load StdEnv/2023
module load r/4.4.0
R

## Load libraries
library(MAST)
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
library(ggsankey)
library(ggalluvial)

#install.packages("remotes")
#remotes::install_github("davidsjoberg/ggsankey")

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Performance by broad cell type
## code
    ###########################
    ## All HVGs
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/All_HVG_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, run) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- data.frame(summary_stats)

        ## Set validation and evakuation columns.
        summary_stats$run_type[summary_stats$run == "validation"] <- "Validation"
        summary_stats$run_type[summary_stats$run == "evaluation_unbalanced"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_1X"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_2X"] <- "Evaluation"

        ## Create a validation dataset
        summary_stats_val <- subset(summary_stats, run_type == "Validation")

        ## Create a evaluation dataset
        summary_stats_eval <- subset(summary_stats, run_type == "Evaluation")

        summary_stats_eval_best <- summary_stats_eval %>%
            group_by(celltype) %>%
            slice_max(order_by = mean_accuracy, n = 1) %>%
            ungroup()

        ## Create an overall accuracy dataframe
        overall_accuracy <- rbind(summary_stats_val, summary_stats_eval_best)
        overall_accuracy <- data.frame(overall_accuracy)
        overall_accuracy$run_type[overall_accuracy$run == "validation"] <- "Validation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_unbalanced"] <- "Evaluation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_1X"] <- "Evaluation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_2X"] <- "Evaluation"
        overall_accuracy$key <- paste0(overall_accuracy$celltype, '_', overall_accuracy$run )
        overall_accuracy$group <- "Overall"

        table(overall_accuracy$celltype)

        overall_accuracy$major_group[overall_accuracy$celltype == "L2_L3"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L3_L5"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L4_L5"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L4_L6"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L5"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L5_L6"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L6"] <- "Excitatory"

        overall_accuracy$major_group[overall_accuracy$celltype == "5HT3aR"] <- "Inhibitory"
        overall_accuracy$major_group[overall_accuracy$celltype == "PV"] <- "Inhibitory"
        overall_accuracy$major_group[overall_accuracy$celltype == "Rosehip"] <- "Inhibitory"
        overall_accuracy$major_group[overall_accuracy$celltype == "SOM"] <- "Inhibitory"

        overall_accuracy$major_group[overall_accuracy$celltype == "Astro"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "Endo"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "Fibro"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "Micro"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "Mural"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "Oligo"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "OPC"] <- "non-neuronal"

        ## Compute disease group specific
        summary_stats_group <- overall_accuracy %>%
            group_by(major_group) %>%
            summarise(
                median_accuracy = median(mean_accuracy), 
                mean_accuracy = mean(mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(mean_accuracy - median(mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats_group <- data.frame(summary_stats_group)
        #major_group median_accuracy mean_accuracy sd_accuracy mad_accuracy
        #Excitatory       0.9534029     0.9341193          NA            0
        #Inhibitory       0.9321111     0.9162656          NA            0
        #non-neuronal       0.9525628     0.9052889          NA            0

        overall_accuracy_HVG <- overall_accuracy
        overall_accuracy_HVG$input <- "HVG"

    ## 

    ###########################
    ## NMF
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, run) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- data.frame(summary_stats)

        ## Set validation and evakuation columns.
        summary_stats$run_type[summary_stats$run == "validation"] <- "Validation"
        summary_stats$run_type[summary_stats$run == "evaluation_unbalanced"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_1X"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_2X"] <- "Evaluation"

        ## Create a validation dataset
        summary_stats_val <- subset(summary_stats, run_type == "Validation")

        ## Create a evaluation dataset
        summary_stats_eval <- subset(summary_stats, run_type == "Evaluation")

        summary_stats_eval_best <- summary_stats_eval %>%
            group_by(celltype) %>%
            slice_max(order_by = mean_accuracy, n = 1) %>%
            ungroup()

        ## Create an overall accuracy dataframe
        overall_accuracy <- rbind(summary_stats_val, summary_stats_eval_best)
        overall_accuracy <- data.frame(overall_accuracy)
        overall_accuracy$run_type[overall_accuracy$run == "validation"] <- "Validation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_unbalanced"] <- "Evaluation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_1X"] <- "Evaluation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_2X"] <- "Evaluation"
        overall_accuracy$key <- paste0(overall_accuracy$celltype, '_', overall_accuracy$run )
        overall_accuracy$group <- "Overall"

        table(overall_accuracy$celltype)

        overall_accuracy$major_group[overall_accuracy$celltype == "L2_L3"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L3_L5"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L4_L5"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L4_L6"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L5"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L5_L6"] <- "Excitatory"
        overall_accuracy$major_group[overall_accuracy$celltype == "L6"] <- "Excitatory"

        overall_accuracy$major_group[overall_accuracy$celltype == "5HT3aR"] <- "Inhibitory"
        overall_accuracy$major_group[overall_accuracy$celltype == "PV"] <- "Inhibitory"
        overall_accuracy$major_group[overall_accuracy$celltype == "Rosehip"] <- "Inhibitory"
        overall_accuracy$major_group[overall_accuracy$celltype == "SOM"] <- "Inhibitory"

        overall_accuracy$major_group[overall_accuracy$celltype == "Astro"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "Endo"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "Fibro"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "Micro"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "Mural"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "Oligo"] <- "non-neuronal"
        overall_accuracy$major_group[overall_accuracy$celltype == "OPC"] <- "non-neuronal"

        ## Compute disease group specific
        summary_stats_group <- overall_accuracy %>%
            group_by(major_group) %>%
            summarise(
                median_accuracy = median(mean_accuracy), 
                mean_accuracy = mean(mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(mean_accuracy - median(mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats_group <- data.frame(summary_stats_group)
        #major_group  median_accuracy mean_accuracy sd_accuracy mad_accuracy
        #<chr>                  <dbl>         <dbl>       <dbl>        <dbl>
        #1 Excitatory             0.903         0.870          NA            0
        #2 Inhibitory             0.835         0.787          NA            0
        #3 non-neuronal           0.751         0.783          NA            0

        overall_accuracy_NMF <- overall_accuracy
        overall_accuracy_NMF$input <- "NMF"

    ##

    ###########################
    ## Combine plot
    ###########################
    ## code
        overall_accuracy_total <- rbind(overall_accuracy_HVG, overall_accuracy_NMF)
        
    
        ggplot(overall_accuracy_total, aes(x = major_group, y = mean_accuracy, fill = major_group)) + 
        theme_bw() + 
        geom_boxplot(outlier.shape = NA) +
        geom_line(data = means_df, aes(x = major_group, y = mean_val, group = input),
            inherit.aes = FALSE, color = "red", linewidth = 0.4) +
        stat_summary(fun = mean, geom = "point", shape = 23, size = 1, fill = "red", colour = "red") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~input, scales = "free_x") +
        scale_fill_manual(values = c("Excitatory" = "#00A8A8",
                                "Inhibitory"  = "#9AD3DA",
                                "non-neuronal" = "#FF924D")) +
        ylab("Accuracy") +
        scale_x_discrete(labels = c("Excitatory", "Inhibitory", "Non-neuronal"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width =2.25)
    ##

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Performance vs n cells.
## code
    ###########################
    ## All HVGs
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/All_HVG_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, run) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
        summary_stats <- data.frame(summary_stats)

        ## Set validation and evakuation columns.
        summary_stats$run_type[summary_stats$run == "validation"] <- "Validation"
        summary_stats$run_type[summary_stats$run == "evaluation_unbalanced"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_1X"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_2X"] <- "Evaluation"

        ## Create a validation dataset
        summary_stats_val <- subset(summary_stats, run_type == "Validation")

        ## Create a evaluation dataset
        summary_stats_eval <- subset(summary_stats, run_type == "Evaluation")

        summary_stats_eval_best <- summary_stats_eval %>%
            group_by(celltype) %>%
            slice_max(order_by = mean_accuracy, n = 1) %>%
            ungroup()

        ## Create an overall accuracy dataframe
        overall_accuracy <- rbind(summary_stats_val, summary_stats_eval_best)
        overall_accuracy <- data.frame(overall_accuracy)
        overall_accuracy$run_type[overall_accuracy$run == "validation"] <- "Validation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_unbalanced"] <- "Evaluation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_1X"] <- "Evaluation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_2X"] <- "Evaluation"
        overall_accuracy$key <- paste0(overall_accuracy$celltype, '_', overall_accuracy$run )
        overall_accuracy$group <- "Overall"

        table(overall_accuracy$celltype)

        overall_accuracy_best <- overall_accuracy


        #############
        
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/All_HVG_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        BA4_bind_all$key <- paste0(BA4_bind_all$celltype, '_', BA4_bind_all$run )
        BA4_bind_all_best <- subset(BA4_bind_all, key %in% overall_accuracy_best$key)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all_best %>%
            group_by(celltype, run, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- data.frame(summary_stats)
        

        summary_stats$run_type[summary_stats$run == "validation"] <- "Validation"
        summary_stats$run_type[summary_stats$run == "evaluation_unbalanced"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_1X"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_2X"] <- "Evaluation"

        summary_stats$merge_key <- paste0(summary_stats$run_type, '_', summary_stats$group, '_', summary_stats$celltype )

        ## read in ncell counts
        n_cell <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/outs/cell_count_broad_disease_group.csv')

        n_cell <- subset(n_cell, Dataset %in% c("Pineda et al.", "Li et al. and Limone et al."))

        n_cell$run_type[n_cell$Dataset == "Pineda et al."] <- "Validation"
        n_cell$run_type[n_cell$Dataset == "Li et al. and Limone et al."] <- "Evaluation"
        n_cell$Group[n_cell$Group == "Control"] <- "control"

        n_cell$merge_key <- paste0(n_cell$run_type, '_', n_cell$Group, '_', n_cell$Cell_Type )

        n_cell <- n_cell %>% dplyr::select(Freq, merge_key)

        ## merge
        nrow(n_cell) == nrow(summary_stats)
        merge <- merge(summary_stats, n_cell, by = "merge_key")
        nrow(merge) == nrow(summary_stats)

        merge$input = "HVG"
        merge_HVG <- merge
    ##

    ###########################
    ## NMF
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, run) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
        summary_stats <- data.frame(summary_stats)

        ## Set validation and evakuation columns.
        summary_stats$run_type[summary_stats$run == "validation"] <- "Validation"
        summary_stats$run_type[summary_stats$run == "evaluation_unbalanced"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_1X"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_2X"] <- "Evaluation"

        ## Create a validation dataset
        summary_stats_val <- subset(summary_stats, run_type == "Validation")

        ## Create a evaluation dataset
        summary_stats_eval <- subset(summary_stats, run_type == "Evaluation")

        summary_stats_eval_best <- summary_stats_eval %>%
            group_by(celltype) %>%
            slice_max(order_by = mean_accuracy, n = 1) %>%
            ungroup()

        ## Create an overall accuracy dataframe
        overall_accuracy <- rbind(summary_stats_val, summary_stats_eval_best)
        overall_accuracy <- data.frame(overall_accuracy)
        overall_accuracy$run_type[overall_accuracy$run == "validation"] <- "Validation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_unbalanced"] <- "Evaluation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_1X"] <- "Evaluation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_2X"] <- "Evaluation"
        overall_accuracy$key <- paste0(overall_accuracy$celltype, '_', overall_accuracy$run )
        overall_accuracy$group <- "Overall"

        table(overall_accuracy$celltype)

        overall_accuracy_best <- overall_accuracy


        #############
        
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        BA4_bind_all$key <- paste0(BA4_bind_all$celltype, '_', BA4_bind_all$run )
        BA4_bind_all_best <- subset(BA4_bind_all, key %in% overall_accuracy_best$key)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all_best %>%
            group_by(celltype, run, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- data.frame(summary_stats)
        

        summary_stats$run_type[summary_stats$run == "validation"] <- "Validation"
        summary_stats$run_type[summary_stats$run == "evaluation_unbalanced"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_1X"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_2X"] <- "Evaluation"

        summary_stats$merge_key <- paste0(summary_stats$run_type, '_', summary_stats$group, '_', summary_stats$celltype )

        ## read in ncell counts
        n_cell <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/outs/cell_count_broad_disease_group.csv')

        n_cell <- subset(n_cell, Dataset %in% c("Pineda et al.", "Li et al. and Limone et al."))

        n_cell$run_type[n_cell$Dataset == "Pineda et al."] <- "Validation"
        n_cell$run_type[n_cell$Dataset == "Li et al. and Limone et al."] <- "Evaluation"
        n_cell$Group[n_cell$Group == "Control"] <- "control"

        n_cell$merge_key <- paste0(n_cell$run_type, '_', n_cell$Group, '_', n_cell$Cell_Type )

        n_cell <- n_cell %>% dplyr::select(Freq, merge_key)

        ## merge
        nrow(n_cell) == nrow(summary_stats)
        merge <- merge(summary_stats, n_cell, by = "merge_key")
        nrow(merge) == nrow(summary_stats)

        merge$input = "NMF"
        merge_NMF <- merge
    ##

    ###########################
    ## bind and plot
    ###########################
    ## code
        merge_all <- rbind(merge_HVG, merge_NMF)
       
        ## Plot scatter plot
        cell_type_levels <- c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC')
        colors <- c(
            "L3_L5" = "#f3c300",
            "L2_L3" = "#f38400",
            "L4_L6" = "#a1caf1",
            "L4_L5" = "#be0032",
            "L5_L6" = "#c2b280",
            "L5" = "#008856",
            "L6" = "#2b3d26",
            "PV" = "#e25822",
            "5HT3aR" = "#e68fac",
            "Rosehip" = "#0067a5",
            "SOM" = "#f99379",
            "Oligo" = "#604e97",
            "Astro" = "#875692",
            "OPC" = "#f6a600",
            "Micro" = "#b3446c",
            "T_Cell" = "#dcd300",
            "Mural" = "#882d17",
            "Endo" = "#8db600",
            "Fibro" = "#654522")

        merge_all$celltype <- factor(merge_all$celltype, levels =  cell_type_levels)
        merge_all$run_type <- factor(merge_all$run_type, levels =  c("Validation", "Evaluation"))

        ggplot(merge_all, aes(x = Freq, y = mean_accuracy)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_point(aes(colour = celltype, shape = run_type)) +
            theme(
                legend.position = "right",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black"),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_grid(.~input) +
            scale_colour_manual(values = colors) +
            scale_shape_manual(values = c(16,15)) +
            ylab("Accuracy") +
            xlab("N cells") +
            stat_cor(
                method = "pearson",
                label.x = 2.75,  # X-coordinate in data units
                label.y = 0.01, # Y-coordinate in data units
                size = 2.5) +
            scale_x_log10(
                breaks = c(1, 1e1, 1e2, 1e3, 1e4, 1e5),
                labels = c("1", "1e1", "1e2", "1e3", "1e4", "1e5")
            ) 
        
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 5)

        ggplot(merge_all, aes(x = Freq, y = mean_accuracy)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_point(aes(colour = celltype, shape = run_type)) +
            theme(
                legend.position = "right",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black"),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_grid(group~input) +
            scale_colour_manual(values = colors) +
            scale_shape_manual(values = c(16,15)) +
            ylab("Accuracy") +
            xlab("N cells") +
            stat_cor(
                method = "pearson",
                label.x = 2.75,  # X-coordinate in data units
                label.y = 0.01, # Y-coordinate in data units
                size = 2.5) +
            scale_x_log10(
                breaks = c(1, 1e1, 1e2, 1e3, 1e4, 1e5),
                labels = c("1", "1e1", "1e2", "1e3", "1e4", "1e5")
            ) 

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 5, width = 5)
    ##
##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Missclassification bubble plot
## code: HVGs ALS, FTLD, Control
    ## read in file
    BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
    
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

    ## Create group label
    BA4_bind_all <- BA4_bind_all %>%
        mutate(group = case_when(
            str_ends(donor, "_1") ~ "ALS",
            str_ends(donor, "_2") ~ "FTLD",
            str_ends(donor, "_0") ~ "control",
            TRUE ~ NA_character_
        ))

    ## Create true class label
    BA4_bind_all <- BA4_bind_all %>%
        mutate(true_class = case_when(
            str_ends(donor, "_1") ~ "1",
            str_ends(donor, "_2") ~ "2",
            str_ends(donor, "_0") ~ "0",
            TRUE ~ NA_character_
        ))

    ## reduce limit cells
    BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
    
    ## Lets reduce the number of columns
    BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, region, group, test_accuracy_all, true_class, counts_pred)

    # Example: df is your dataframe
    df_parsed <- BA4_bind_all %>%
    # Expand counts_pred into multiple rows
    mutate(counts_pred = strsplit(counts_pred, ";")) %>%
    unnest(counts_pred) %>%
    mutate(
        counts_pred = str_trim(counts_pred),
        class = str_extract(counts_pred, '"\\d+"'),
        class = str_replace_all(class, '"', ""),   # remove quotes
        count = as.numeric(str_extract(counts_pred, "\\d+$"))
    ) %>%
    select(-counts_pred)

    # Now compute proportions per donor + celltype
    df_props <- df_parsed %>%
    group_by(donor, celltype, region, group, true_class) %>%
    mutate(
        total = sum(count),
        proportion = count / total
    ) %>%
    ungroup()

    df_props <- data.frame(df_props)

    df_avg_props <- df_props %>%
    group_by(celltype, group, class) %>%
    summarise(
        mean_proportion = mean(proportion, na.rm = TRUE),
        sd_proportion   = sd(proportion, na.rm = TRUE),
        n_donors        = n_distinct(donor)
    ) %>%
    ungroup()

    df_avg_props <- data.frame(df_avg_props)

    ## prepare plot columns
    df_avg_props$true_class[df_avg_props$group == "ALS"] <- "1"
    df_avg_props$true_class[df_avg_props$group == "FTLD"] <- "2"
    df_avg_props$true_class[df_avg_props$group == "control"] <- "0"

    df_avg_props$predicted_class <- df_avg_props$class

    df_avg_props$celltype <- factor(df_avg_props$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

    temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
        theme_bw() + 
        #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
        geom_tile() +
        theme(
            legend.position = "right",
            legend.key.size = unit(0.2, "cm"),      
            legend.text = element_text(size = 6),   
            legend.title = element_text(size = 7),  
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title = element_text(face = "bold", size = 8),
            strip.text = element_text(size = 8, face = "bold")
        ) +
        facet_grid(.~celltype) +
        #scale_fill_viridis_c(direction = -1) +
        scale_fill_gradient(low = "white", high = "black") +   # <-- white = low, darkred = high
        #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
        #scale_colour_manual(values = "black") +
        #scale_shape_manual(values = c(16,15)) +
        ylab("Predicted label") +
        xlab("True label") +
        scale_x_discrete(labels = c("Control", "ALS", "FTLD")) +
        scale_y_discrete(labels = c("Control", "ALS", "FTLD"))

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.75, width = 13)

##


## code: HVGs ALS, FTLD, Control All  cells together swankey plot.
    ## read in file
    BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
    
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

    ## Create group label
    BA4_bind_all <- BA4_bind_all %>%
        mutate(group = case_when(
            str_ends(donor, "_1") ~ "ALS",
            str_ends(donor, "_2") ~ "FTLD",
            str_ends(donor, "_0") ~ "control",
            TRUE ~ NA_character_
        ))

    ## Create true class label
    BA4_bind_all <- BA4_bind_all %>%
        mutate(true_class = case_when(
            str_ends(donor, "_1") ~ "1",
            str_ends(donor, "_2") ~ "2",
            str_ends(donor, "_0") ~ "0",
            TRUE ~ NA_character_
        ))

    ## reduce limit cells
    BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
    
    ## Lets reduce the number of columns
    BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, region, group, test_accuracy_all, true_class, counts_pred)

    # Example: df is your dataframe
    df_parsed <- BA4_bind_all %>%
    # Expand counts_pred into multiple rows
    mutate(counts_pred = strsplit(counts_pred, ";")) %>%
    unnest(counts_pred) %>%
    mutate(
        counts_pred = str_trim(counts_pred),
        class = str_extract(counts_pred, '"\\d+"'),
        class = str_replace_all(class, '"', ""),   # remove quotes
        count = as.numeric(str_extract(counts_pred, "\\d+$"))
    ) %>%
    select(-counts_pred)

    # Now compute proportions per donor + celltype
    df_props <- df_parsed %>%
    group_by(donor, celltype, region, group, true_class) %>%
    mutate(
        total = sum(count),
        proportion = count / total
    ) %>%
    ungroup()

    df_props <- data.frame(df_props)

    df_props <- data.frame(df_props)

    df_avg_props <- df_props %>%
    group_by(group, class) %>%
    summarise(
        mean_proportion = mean(proportion, na.rm = TRUE),
        sd_proportion   = sd(proportion, na.rm = TRUE),
        n_donors        = n_distinct(donor)
    ) %>%
    ungroup()

    df_avg_props <- data.frame(df_avg_props)

    ## prepare plot columns
    df_avg_props$true_class[df_avg_props$group == "ALS"] <- "1"
    df_avg_props$true_class[df_avg_props$group == "FTLD"] <- "2"
    df_avg_props$true_class[df_avg_props$group == "control"] <- "0"

    df_avg_props$predicted_class <- df_avg_props$class

    df_avg_props$true_class[df_avg_props$true_class == "1"] <- "ALS"
    df_avg_props$true_class[df_avg_props$true_class == "2"] <- "FTLD"
    df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"

    df_avg_props$predicted_class[df_avg_props$predicted_class == "1"] <- "ALS"
    df_avg_props$predicted_class[df_avg_props$predicted_class == "2"] <- "FTLD"
    df_avg_props$predicted_class[df_avg_props$predicted_class == "0"] <- "Control"

    
    df_avg_props$predicted_class <- factor(df_avg_props$predicted_class, levels = c("Control", "ALS", "FTLD"))
    df_avg_props$true_class <- factor(df_avg_props$true_class, levels = c("Control", "ALS", "FTLD"))

    ## TEST SWANKEY PLOT all cell types together
    ggplot(data = df_avg_props, aes(axis1 = true_class, axis2 = predicted_class, y = mean_proportion)) +
        geom_alluvium(aes(fill = group)) +
        geom_stratum() +
        geom_text(stat = "stratum",
                    aes(label = after_stat(stratum))) +
        scale_x_discrete(limits = c("True", "Predicted"),
                        expand = c(0.15, 0.05)) +
        theme_void() + 
        #facet_grid(.~celltype)
        theme(legend.position = "none") +
        scale_fill_manual(values = c("ALS" = "darkred",
                        "FTLD" = "darkblue",
                        "control" = "#339966")) 
     ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 5, width = 5)

##


## code: HVGs SALS, C9ALS, SFTLD, C9FTLD, Control
    ## read in file
    BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/All_HVGs_fix_HVG_generalizable_subtype_combined.csv')
    
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

    ## Create group label
    BA4_bind_all <- BA4_bind_all %>%
        mutate(group = case_when(
            str_ends(donor, "_1") ~ "SALS",
            str_ends(donor, "_2") ~ "C9ALS",
            str_ends(donor, "_3") ~ "SFTLD",
            str_ends(donor, "_4") ~ "C9FTLD",
            str_ends(donor, "_0") ~ "control",
            TRUE ~ NA_character_
        ))

    ## Create true class label
    BA4_bind_all <- BA4_bind_all %>%
        mutate(true_class = case_when(
            str_ends(donor, "_1") ~ "1",
            str_ends(donor, "_2") ~ "2",
            str_ends(donor, "_3") ~ "3",
            str_ends(donor, "_4") ~ "4",
            str_ends(donor, "_0") ~ "0",
            TRUE ~ NA_character_
        ))

    ## reduce limit cells
    BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
    
    ## Lets reduce the number of columns
    BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, region, group, test_accuracy_all, true_class, counts_pred)

    # Example: df is your dataframe
    df_parsed <- BA4_bind_all %>%
    # Expand counts_pred into multiple rows
    mutate(counts_pred = strsplit(counts_pred, ";")) %>%
    unnest(counts_pred) %>%
    mutate(
        counts_pred = str_trim(counts_pred),
        class = str_extract(counts_pred, '"\\d+"'),
        class = str_replace_all(class, '"', ""),   # remove quotes
        count = as.numeric(str_extract(counts_pred, "\\d+$"))
    ) %>%
    select(-counts_pred)

    # Now compute proportions per donor + celltype
    df_props <- df_parsed %>%
    group_by(donor, celltype, region, group, true_class) %>%
    mutate(
        total = sum(count),
        proportion = count / total
    ) %>%
    ungroup()

    df_props <- data.frame(df_props)

    df_avg_props <- df_props %>%
    group_by(celltype, group, class) %>%
    summarise(
        mean_proportion = mean(proportion, na.rm = TRUE),
        sd_proportion   = sd(proportion, na.rm = TRUE),
        n_donors        = n_distinct(donor)
    ) %>%
    ungroup()

    df_avg_props <- data.frame(df_avg_props)

    ## prepare plot columns
    df_avg_props$true_class[df_avg_props$group == "SALS"] <- "1"
    df_avg_props$true_class[df_avg_props$group == "C9ALS"] <- "2"
    df_avg_props$true_class[df_avg_props$group == "SFTLD"] <- "3"
    df_avg_props$true_class[df_avg_props$group == "C9FTLD"] <- "4"
    df_avg_props$true_class[df_avg_props$group == "control"] <- "0"

    df_avg_props$predicted_class <- df_avg_props$class

    df_avg_props$celltype <- factor(df_avg_props$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

    
    temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
        theme_bw() + 
        #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
        geom_tile() +
        theme(
            legend.position = "right",
            legend.key.size = unit(0.2, "cm"),      
            legend.text = element_text(size = 6),   
            legend.title = element_text(size = 7),  
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title = element_text(face = "bold", size = 8),
            strip.text = element_text(size = 8, face = "bold")
        ) +
        facet_grid(.~celltype) +
        #scale_fill_viridis_c(direction = -1) +
        scale_fill_gradient(low = "white", high = "black") +   # <-- white = low, darkred = high
        #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
        #scale_colour_manual(values = "black") +
        #scale_shape_manual(values = c(16,15)) +
        ylab("Predicted label") +
        xlab("True label") +
        scale_x_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD")) +
        scale_y_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD"))

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.75, width = 13)

    ## Flow Arrow plot
    # Example dataset (df_avg_props)
    # It should have: true_class, predicted_class, mean_proportion, celltype

    #df_avg_props$true_class[df_avg_props$true_class == "1"] <- "SALS"
    #df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9ALS"
    #df_avg_props$true_class[df_avg_props$true_class == "3"] <- "SFTLD"
    #df_avg_props$true_class[df_avg_props$true_class == "4"] <- "C9FTLD"
    #df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"

    #df_avg_props$predicted_class[df_avg_props$predicted_class == "1"] <- "SALS"
    #df_avg_props$predicted_class[df_avg_props$predicted_class == "2"] <- "C9ALS"
    #df_avg_props$predicted_class[df_avg_props$predicted_class == "3"] <- "SFTLD"
    #df_avg_props$predicted_class[df_avg_props$predicted_class == "4"] <- "C9FTLD"
    #df_avg_props$predicted_class[df_avg_props$predicted_class == "0"] <- "Control"


    #df_avg_props$predicted_class <- factor(df_avg_props$predicted_class, levels = c("Control", "SALS", "C9ALS", "SFTLD", "C9FTLD"))
    #df_avg_props$true_class <- factor(df_avg_props$true_class, levels = c("Control", "SALS", "C9ALS", "SFTLD", "C9FTLD"))

    ## TEST SWANKEY PLOT all cell types together
    #ggplot(data = df_avg_props, aes(axis1 = true_class, axis2 = predicted_class, y = mean_proportion)) +
    #    geom_alluvium(aes(fill = group)) +
    #    geom_stratum() +
    #    geom_text(stat = "stratum",
    #                aes(label = after_stat(stratum))) +
    #    scale_x_discrete(limits = c("True", "Predicted"),
    #                    expand = c(0.15, 0.05)) +
    #    theme_void() + 
    #    #facet_grid(.~celltype)
    #    theme(legend.position = "none") +
    #    scale_fill_manual(values = c("SALS" = "orange",
    #                    "C9ALS"  = "red",
    #                    "SFTLD" = "blue",
    #                    "C9FTLD"  = "purple",
    #                    "control" = "#339966")) 


    #df_avg_props$major_group[df_avg_props$celltype == "L2_L3"] <- "Excitatory"
    #df_avg_props$major_group[df_avg_props$celltype == "L3_L5"] <- "Excitatory"
    #df_avg_props$major_group[df_avg_props$celltype == "L4_L5"] <- "Excitatory"
    #df_avg_props$major_group[df_avg_props$celltype == "L4_L6"] <- "Excitatory"
    #df_avg_props$major_group[df_avg_props$celltype == "L5"] <- "Excitatory"
    #df_avg_props$major_group[df_avg_props$celltype == "L5_L6"] <- "Excitatory"
    #df_avg_props$major_group[df_avg_props$celltype == "L6"] <- "Excitatory"

    #df_avg_props$major_group[df_avg_props$celltype == "5HT3aR"] <- "Inhibitory"
    #df_avg_props$major_group[df_avg_props$celltype == "PV"] <- "Inhibitory"
    #df_avg_props$major_group[df_avg_props$celltype == "Rosehip"] <- "Inhibitory"
    #df_avg_props$major_group[df_avg_props$celltype == "SOM"] <- "Inhibitory"

    #df_avg_props$major_group[df_avg_props$celltype == "Astro"] <- "non-neuronal"
    #df_avg_props$major_group[df_avg_props$celltype == "Endo"] <- "non-neuronal"
    #df_avg_props$major_group[df_avg_props$celltype == "Fibro"] <- "non-neuronal"
    #df_avg_props$major_group[df_avg_props$celltype == "Micro"] <- "non-neuronal"
    #df_avg_props$major_group[df_avg_props$celltype == "Mural"] <- "non-neuronal"
    #df_avg_props$major_group[df_avg_props$celltype == "Oligo"] <- "non-neuronal"
    #df_avg_props$major_group[df_avg_props$celltype == "OPC"] <- "non-neuronal"


    
    ## TEST SWANKEY PLOT all cell types together
    #ggplot(data = df_avg_props, aes(axis1 = true_class, axis2 = predicted_class, y = mean_proportion)) +
    #    geom_alluvium(aes(fill = group)) +
    #    geom_stratum() +
    #    geom_text(stat = "stratum",
    #                aes(label = after_stat(stratum))) +
    #    scale_x_discrete(limits = c("True", "Predicted"),
    #                    expand = c(0.15, 0.05)) +
    #    theme_void() + 
    #    #facet_grid(.~celltype)
    #    theme(legend.position = "none") +
    #    scale_fill_manual(values = c("SALS" = "orange",
    #                    "C9ALS"  = "red",
    #                    "SFTLD" = "blue",
    #                    "C9FTLD"  = "purple",
    #                    "control" = "#339966")) 

    ## Exc neurons
    #df_avg_props_exc <- subset(df_avg_props, major_group == "Excitatory")
    #ggplot(data = df_avg_props_exc, aes(axis1 = true_class, axis2 = predicted_class, y = mean_proportion)) +
    #    geom_alluvium(aes(fill = group)) +
    #    geom_stratum() +
    #    geom_text(stat = "stratum",
    #                aes(label = after_stat(stratum))) +
    #    scale_x_discrete(limits = c("True", "Predicted"),
    #                    expand = c(0.15, 0.05)) +
    #    theme_void() + 
    #    #facet_grid(.~celltype)
    #    theme(legend.position = "none") +
    #    scale_fill_manual(values = c("SALS" = "orange",
    #                    "C9ALS"  = "red",
    #                    "SFTLD" = "blue",
    #                    "C9FTLD"  = "purple",
    #                    "control" = "#339966")) 

    #  ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 5, width = 5)

    ## Non-neuron
    #df_avg_props_non <- subset(df_avg_props, major_group == "non-neuronal")
    #ggplot(data = df_avg_props_non, aes(axis1 = true_class, axis2 = predicted_class, y = mean_proportion)) +
    #    geom_alluvium(aes(fill = group)) +
    #    geom_stratum() +
    #    geom_text(stat = "stratum",
    #                aes(label = after_stat(stratum))) +
    #    scale_x_discrete(limits = c("True", "Predicted"),
    #                    expand = c(0.15, 0.05)) +
    #    theme_void() + 
    #    #facet_grid(.~celltype)
    #    theme(legend.position = "none") +
    #    scale_fill_manual(values = c("SALS" = "orange",
    #                    "C9ALS"  = "red",
    #                    "SFTLD" = "blue",
    #                    "C9FTLD"  = "purple",
    #                    "control" = "#339966")) 

    #  ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 5, width = 5)

##

## code: HVGs SALS, C9ALS, SFTLD, C9FTLD, Control all cells together swankey plot
    ## read in file
    BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/All_HVGs_fix_HVG_generalizable_subtype_combined.csv')
    
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

    ## Create group label
    BA4_bind_all <- BA4_bind_all %>%
        mutate(group = case_when(
            str_ends(donor, "_1") ~ "SALS",
            str_ends(donor, "_2") ~ "C9ALS",
            str_ends(donor, "_3") ~ "SFTLD",
            str_ends(donor, "_4") ~ "C9FTLD",
            str_ends(donor, "_0") ~ "control",
            TRUE ~ NA_character_
        ))

    ## Create true class label
    BA4_bind_all <- BA4_bind_all %>%
        mutate(true_class = case_when(
            str_ends(donor, "_1") ~ "1",
            str_ends(donor, "_2") ~ "2",
            str_ends(donor, "_3") ~ "3",
            str_ends(donor, "_4") ~ "4",
            str_ends(donor, "_0") ~ "0",
            TRUE ~ NA_character_
        ))

    ## reduce limit cells
    BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
    
    ## Lets reduce the number of columns
    BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, region, group, test_accuracy_all, true_class, counts_pred)

    # Example: df is your dataframe
    df_parsed <- BA4_bind_all %>%
    # Expand counts_pred into multiple rows
    mutate(counts_pred = strsplit(counts_pred, ";")) %>%
    unnest(counts_pred) %>%
    mutate(
        counts_pred = str_trim(counts_pred),
        class = str_extract(counts_pred, '"\\d+"'),
        class = str_replace_all(class, '"', ""),   # remove quotes
        count = as.numeric(str_extract(counts_pred, "\\d+$"))
    ) %>%
    select(-counts_pred)

    # Now compute proportions per donor + celltype
    df_props <- df_parsed %>%
    group_by(donor, celltype, region, group, true_class) %>%
    mutate(
        total = sum(count),
        proportion = count / total
    ) %>%
    ungroup()

    df_props <- data.frame(df_props)

    df_avg_props <- df_props %>%
        group_by(group, class) %>%
        summarise(
            mean_proportion = mean(proportion, na.rm = TRUE),
            sd_proportion   = sd(proportion, na.rm = TRUE),
            n_donors        = n_distinct(donor)
        ) %>%
        ungroup()

        df_avg_props <- data.frame(df_avg_props)

    ## prepare plot columns
    df_avg_props$true_class[df_avg_props$group == "SALS"] <- "1"
    df_avg_props$true_class[df_avg_props$group == "C9ALS"] <- "2"
    df_avg_props$true_class[df_avg_props$group == "SFTLD"] <- "3"
    df_avg_props$true_class[df_avg_props$group == "C9FTLD"] <- "4"
    df_avg_props$true_class[df_avg_props$group == "control"] <- "0"

    df_avg_props$predicted_class <- df_avg_props$class
    
    df_avg_props$true_class[df_avg_props$true_class == "1"] <- "SALS"
    df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9ALS"
    df_avg_props$true_class[df_avg_props$true_class == "3"] <- "SFTLD"
    df_avg_props$true_class[df_avg_props$true_class == "4"] <- "C9FTLD"
    df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"

    df_avg_props$predicted_class[df_avg_props$predicted_class == "1"] <- "SALS"
    df_avg_props$predicted_class[df_avg_props$predicted_class == "2"] <- "C9ALS"
    df_avg_props$predicted_class[df_avg_props$predicted_class == "3"] <- "SFTLD"
    df_avg_props$predicted_class[df_avg_props$predicted_class == "4"] <- "C9FTLD"
    df_avg_props$predicted_class[df_avg_props$predicted_class == "0"] <- "Control"


    df_avg_props$predicted_class <- factor(df_avg_props$predicted_class, levels = c("Control", "SALS", "C9ALS", "SFTLD", "C9FTLD"))
    df_avg_props$true_class <- factor(df_avg_props$true_class, levels = c("Control", "SALS", "C9ALS", "SFTLD", "C9FTLD"))

    ## TEST SWANKEY PLOT all cell types together
    ggplot(data = df_avg_props, aes(axis1 = true_class, axis2 = predicted_class, y = mean_proportion)) +
        geom_alluvium(aes(fill = group)) +
        geom_stratum() +
        geom_text(stat = "stratum",
                    aes(label = after_stat(stratum))) +
        scale_x_discrete(limits = c("True", "Predicted"),
                        expand = c(0.15, 0.05)) +
        theme_void() + 
        #facet_grid(.~celltype)
        theme(legend.position = "none") +
        scale_fill_manual(values = c("SALS" = "orange",
                        "C9ALS"  = "red",
                        "SFTLD" = "blue",
                        "C9FTLD"  = "purple",
                        "control" = "#339966")) 
             ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 5, width = 5)

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Performance vs TxD
 
## code
    
    ###########################
    ## load libraries
    ###########################
    ##
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
        #library(limma)
        library(sva)

    ##

    ###########################
    ## Prep merge object
    ###########################
    ## code
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_all_celltypes_lim_narval.rds')
        dim(seu)
        
        ## Subset to 100,000
        n_cells <- 100000
        set.seed(123)  # for reproducibility
        cells_to_keep <- sample(colnames(seu), n_cells)
        seu_sub <- subset(seu, cells = cells_to_keep)
        dim(seu_sub)
        rm(seu)

        ## Prep
        seu_sub <- NormalizeData(seu_sub)
        seu_sub <- FindVariableFeatures(seu_sub, selection.method = "vst", nfeatures = 10000)
        seu_sub <- ScaleData(seu_sub
            #features = rownames(seu_sub),
            #vars.to.regress = c("orig.ident")
            )
        seu_sub <- RunPCA(seu_sub, features = VariableFeatures(seu_sub))
        seu_sub <- FindNeighbors(seu_sub, dims = 1:20)
        seu_sub <- FindClusters(seu_sub, resolution = 0.5)
        seu_sub <- RunUMAP(seu_sub, dims = 1:20)
        
        ## Plot UMAPs
        colors <- c(
                "L3_L5" = "#f3c300",
                "L2_L3" = "#f38400",
                "L4_L6" = "#a1caf1",
                "L4_L5" = "#be0032",
                "L5_L6" = "#c2b280",
                "L5" = "#008856",
                "L6" = "#2b3d26",
                "PV" = "#e25822",
                "5HT3aR" = "#e68fac",
                "Rosehip" = "#0067a5",
                "SOM" = "#f99379",
                "Oligo" = "#604e97",
                "Astro" = "#875692",
                "OPC" = "#f6a600",
                "Micro" = "#b3446c",
                "T_Cell" = "#dcd300",
                "Mural" = "#882d17",
                "Endo" = "#8db600",
                "Fibro" = "#654522")
        
        str(seu_sub@meta.data)

        ## By Celltype
        DimPlot(seu_sub, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE, pt.size = 0.1) +
            theme_void()+
            theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            plot.title = element_blank()) + 
            xlab("UMAP1") + ylab("UMAP2") +
            scale_color_manual(values =  colors)

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 10, width = 10)

        ## By Region
        DimPlot(seu_sub, reduction = "umap", group.by = "Region", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE, pt.size = 0.1) +
            theme_void()+
            theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            plot.title = element_blank()) + 
            xlab("UMAP1") + ylab("UMAP2") #+
            #scale_color_manual(values =  colors)

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 10, width = 10)

        ## By Region
        DimPlot(seu_sub, reduction = "umap", group.by = "Condition", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE, pt.size = 0.1) +
            theme_void()+
            theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            plot.title = element_blank()) + 
            xlab("UMAP1") + ylab("UMAP2") #+
            #scale_color_manual(values =  colors)

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 10, width = 10)


        #WILL ALSO TRY WITH COMBAT CORRECTION.
        
        #TRY USING MINIMUM NUMBER OF CELLS -- doesnt do shit

        trying to co-variate correct the pseudobulk profiles -- will try to visualize with PCA. 
        If doesnt work can covariate correct with all cell types together, right now were correcting cell type specific. 
    ##

    ###########################
    ## Compute TxD
    ###########################
    ## code
    
        ## NOTES: ########################################################################################
        
        #"The TxD score is a quantification of transcriptional dysregulation. 
        #It represents the change in transcriptome-wide gene expression of each 
        #cell type in disease from its respective PN expression profile. 
        #The divergence score is the Euclidean distance between the median disease 
        #and corresponding PN covariate-corrected, pseudo-bulk expression profiles for each cell type."
        
            #Try with same number of cells
            #seu_sub <- subset(seu_sub, subset = CellType != "T_Cell")
            # Check how many cells per cell type
            #table(seu_sub$CellType)

            # How many cells to keep per cell type?
            #n_cells <- min(table(seu_sub$CellType))

            #set.seed(123)  # reproducibility
            #cells_keep <- seu_sub@meta.data %>%
            #tibble::rownames_to_column("cell") %>%  # keep cell names
            #group_by(CellType) %>%
            #sample_n(n_cells) %>%
            #pull(cell)   # extract the cell barcodes

            # Subset Seurat object
            #seu_down <- subset(seu_sub, cells = cells_keep)

            # Check balance
            #table(seu_down$CellType)




            ## MEthod 1
            # Create a function to calculate the pseudo-bulk expression profile (median of gene expression per cell type)
            calculate_pseudo_bulk <- function(seurat_obj) {
                pseudo_bulk_profile <- AggregateExpression(seurat_obj, 
                                                            return.seurat = TRUE, 
                                                            verbose = FALSE, 
                                                            group.by = c("Sample_ID", "Condition")#,
                                                            #assays = "RNA",
                                                            #slot = "scale.data"
                                                            )
                return(pseudo_bulk_profile)
            }

            # Function to calculate the Euclidean distance
            calculate_divergence_score <- function(disease_expr, pn_expr) {
                distance <- sqrt(sum((disease_expr - pn_expr)^2))
                return(distance)
            }

            #complete_workflow_TXd <- function(seu_sub){

            ## RUN FROM HERE    
                
                ## List of cell types
                Condition <- unique(seu_sub$Condition)
                Condition <- c("ALS", "FTLD")
                celltypes <- unique(seu_sub$CellType)


                fill <- data.frame(
                celltype = "fill",
                divergence_score = 0,
                Condition = "fill")
            
            
                for (condition in Condition){
                    for (celltype in celltypes) {

                    #condition <- "ALS" ##### TEMP
                    #celltype <- "L3_L5" ## TEMP

                    ## Subset the Seurat object for disease and PN groups
                    sue_obj_lim <- subset(seu_sub, CellType == celltype)
                    print(unique(sue_obj_lim@meta.data$CellType))

                    keep_condition <- c(condition, "PN")    

                    sue_obj_lim <- subset(sue_obj_lim, Condition %in% keep_condition)
                    print(unique(sue_obj_lim@meta.data$Condition))
    
                    ## Remove samples with only one sample
                    num_cells <- data.frame(table(sue_obj_lim@meta.data$orig.ident))
                    num_cells <- subset(num_cells, Freq > 1)
                    keep_donor <- unique(num_cells$Var1)

                    Idents(sue_obj_lim) <- "orig.ident"
                    sue_obj_lim=subset(sue_obj_lim,idents=keep_donor)
                    table(sue_obj_lim@meta.data$orig.ident)

                    m <- sue_obj_lim@assays$RNA@counts
                    adjusted <- ComBat_seq(as.matrix(m), batch=as.numeric(as.factor(sue_obj_lim@meta.data$orig.ident)))

                    min(m)
                    max(m)
                    min(adjusted)
                    max(adjusted)
                    
                    
                    ## Add feature and cell names to combat-corrected matrix
                    #test <- data.frame(seu_lim@assays$RNA@counts)
                    rownames(adjusted) <- rownames(sue_obj_lim)
                    colnames(adjusted) <- rownames(sue_obj_lim@meta.data)

                    ## Create Seurat assay with combat-generated batch corrected matrix
                    assay.v5 <- CreateAssay5Object(counts = as.matrix(adjusted))
                    sue_obj_lim_corr <- CreateSeuratObject(assay.v5)
                    sue_obj_lim_corr <- NormalizeData(sue_obj_lim_corr)
                    sue_obj_lim_corr <- AddMetaData(sue_obj_lim_corr, sue_obj_lim@meta.data)
                    #pbmc3k_slim <- ScaleData(pbmc3k_slim, verbose = FALSE)
                    #pbmc3k_slim <- FindVariableFeatures(pbmc3k_slim)
                    #pbmc3k_slim <- RunPCA(pbmc3k_slim, npcs = 25, verbose = FALSE)
                    #pbmc3k_slim <- RunUMAP(pbmc3k_slim, dims = 1:25, n.neighbors =45)
                    #pbmc3k_slim <- RunTSNE(pbmc3k_slim, dims = 1:25)
                    #dim(pbmc3k_slim)
                        
                    pseudo_bulk <- calculate_pseudo_bulk(sue_obj_lim)
                    print(dim(pseudo_bulk))
                    
                    pseudo_bulk_corr <- calculate_pseudo_bulk(sue_obj_lim_corr)
                    print(dim(pseudo_bulk_corr))

                    ## Covariate correction with LIMA test
                    # Example: remove batch and donor effects
                    #pb_meta <- data.frame(pseudo_bulk@meta.data)
                    #batch <- pb_meta$orig.ident
                    #covariates <- model.matrix(~ 0 + Condition, data = pb_meta)  # leave Condition effect intact

                    #expr_pb <- pseudo_bulk@assays$RNA@layers$counts
                    #print(dim(expr_pb))
                    #max(expr_pb)
                    #min(expr_pb)

                    #expr_pb_corrected <- removeBatchEffect(expr_pb,
                                                        #batch = batch,
                                                        #design = covariates)

                    #max(expr_pb_corrected)
                    #min(expr_pb_corrected)

                    ## Combat correction
                    #m <- pseudo_bulk@assays$RNA@layers$counts
                        


                    ######################################################## PRINT PCAs
                    expr_pb <- pseudo_bulk@assays$RNA@layers$counts
                    expr_pb_df <- data.frame(expr_pb)
                    expr_pb_filtered <- expr_pb_df[apply(expr_pb_df, 1, var) > 0, ]
                    pca_raw <- prcomp(t(expr_pb_filtered), scale. = TRUE)
                    
                    expr_pb_corr <- pseudo_bulk_corr@assays$RNA@layers$counts
                    expr_pb_corr_df <- data.frame(expr_pb_corr)
                    expr_pb_corrected_filtered <- expr_pb_corr_df[apply(expr_pb_corr_df, 1, var) > 0, ]
                    pca_corr <- prcomp(t(expr_pb_corrected_filtered), scale. = TRUE)

                    # Extract first 2 PCs
                    pca_raw_df <- data.frame(
                    Sample = colnames(expr_pb_filtered),
                    PC1 = pca_raw$x[,1],
                    PC2 = pca_raw$x[,2],
                    orig.ident = pseudo_bulk@meta.data$orig.ident,
                    Condition = pseudo_bulk@meta.data$Condition
                    )

                    pca_corr_df <- data.frame(
                    Sample = colnames(expr_pb_corrected_filtered),
                    PC1 = pca_corr$x[,1],
                    PC2 = pca_corr$x[,2],
                    orig.ident = pseudo_bulk_corr$orig.ident,
                    Condition = pseudo_bulk_corr$Condition
                    )

                    # Plot before correction
                    #p1 <- ggplot(pca_raw_df, aes(x = PC1, y = PC2, color = orig.ident)) +
                    #geom_point(size = 3, alpha = 0.8) +
                    #theme_bw() +
                    #theme(legend.position = "none") +
                    #labs(title = "PCA - Before Batch Correction")
                    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 10, width = 10)

                    #p2 <- ggplot(pca_raw_df, aes(x = PC1, y = PC2, color = Condition)) +
                    #geom_point(size = 3, alpha = 0.8) +
                    #theme_bw() +
                    #theme(legend.position = "bottom") +
                    #labs(title = "PCA - Before Batch Correction")
                    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 10, width = 10)

                    #p1 + p2
                    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 5, width = 10)


                    # Plot after correction
                    #p1 <- ggplot(pca_corr_df, aes(x = PC1, y = PC2, color = orig.ident)) +
                    #geom_point(size = 3, alpha = 0.8) +
                    #theme_bw() +
                    #theme(legend.position = "none") +
                    #labs(title = "PCA - Before Batch Correction")
                    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 10, width = 10)

                    #p2 <- ggplot(pca_corr_df, aes(x = PC1, y = PC2, color = Condition)) +
                    #geom_point(size = 3, alpha = 0.8) +
                    #theme_bw() +
                    #theme(legend.position = "bottom") +
                    #labs(title = "PCA - Before Batch Correction")
                    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 10, width = 10)

                    #p1 + p2
                    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 5, width = 10)

                    ######################################################## PRINT PCAs

                    ## Calculate the median profile
                    expr_pb_corr_disease <- subset(pseudo_bulk_corr, Condition == condition)
                    unique(expr_pb_corr_disease$Condition)
                    
                    expr_pb_corr_pn <- subset(pseudo_bulk_corr, Condition == "PN")
                    unique(expr_pb_corr_pn$Condition)
                    
                    disease_expr <- expr_pb_corr_disease@assays$RNA@layers$counts
                    print(dim(disease_expr))
                    rownames(disease_expr) <- rownames(expr_pb_corr_disease)
                    
                    pn_expr <- expr_pb_corr_pn@assays$RNA@layers$counts
                    print(dim(pn_expr))
                    rownames(pn_expr) <- rownames(expr_pb_corr_pn)
                    
                    disease_expr <- as.matrix(disease_expr)
                    pn_expr <- as.matrix(pn_expr)

                    disease_expr <- apply(disease_expr, 1, median)
                    length(disease_expr)
                    
                    pn_expr <- apply(pn_expr, 1, median)
                    length(pn_expr)

                    ## Can we concatonate the median expression for PCA?
                    expr_pb_corr <- pseudo_bulk_corr@assays$RNA@layers$counts
                    expr_pb_corr_df <- data.frame(expr_pb_corr)
                    expr_pb_corr_df$median_disease <- disease_expr
                    expr_pb_corr_df$median_pn <- pn_expr
                    
                    expr_pb_corrected_filtered <- expr_pb_corr_df[apply(expr_pb_corr_df, 1, var) > 0, ]
                    pca_corr <- prcomp(t(expr_pb_corrected_filtered), scale. = TRUE)

                    sample_id <- data.frame(pseudo_bulk_corr$orig.ident)
                    sample_id <- rbind(
                        sample_id,
                        data.frame(pseudo_bulk_corr.orig.ident = "median_disease"),
                        data.frame(pseudo_bulk_corr.orig.ident = "median_control")
                        )
                    tail(sample_id)

                    disease_status <- data.frame(pseudo_bulk_corr$Condition)
                    disease_status <- rbind(
                        disease_status,
                        data.frame(pseudo_bulk_corr.Condition = "median_disease"),
                        data.frame(pseudo_bulk_corr.Condition = "median_control")
                        )
                    tail(disease_status)


                    pca_corr_df <- data.frame(
                    Sample = colnames(expr_pb_corrected_filtered),
                    PC1 = pca_corr$x[,1],
                    PC2 = pca_corr$x[,2],
                    orig.ident = sample_id$pseudo_bulk_corr.orig.ident,
                    Condition = disease_status$pseudo_bulk_corr.Condition
                    )

                    p2 <- ggplot(pca_corr_df, aes(x = PC1, y = PC2, color = Condition)) +
                    geom_point(size = 3) +
                    theme_bw() +
                    theme(legend.position = "none") +
                    scale_colour_manual(values = c("#febdb9", "darkblue", "darkred", 'lightblue')) 
                    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp_TxD_',celltype,'.pdf'), height = 5, width = 6)

    
                    divergence_score <- calculate_divergence_score(disease_expr, pn_expr)


                    ## Dataframe
                    score_results <- data.frame(
                    celltype = celltype,
                    divergence_score = divergence_score
                    )
                    score_results$Condition <- condition

                    fill <- rbind(fill, score_results)
                }
                }


            write.csv(fill, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/Pineda_TxD.csv')
    ##

    ###########################
    ## Compute performance
    ###########################
    ## code HVGs
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/All_HVG_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, run) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
        
        summary_stats <- data.frame(summary_stats)

        ## Set validation and evakuation columns.
        summary_stats$run_type[summary_stats$run == "validation"] <- "Validation"
        summary_stats$run_type[summary_stats$run == "evaluation_unbalanced"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_1X"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_2X"] <- "Evaluation"

        ## Create a validation dataset
        summary_stats_val <- subset(summary_stats, run_type == "Validation")

        ## Create a evaluation dataset
        summary_stats_eval <- subset(summary_stats, run_type == "Evaluation")

        summary_stats_eval_best <- summary_stats_eval %>%
            group_by(celltype) %>%
            slice_max(order_by = mean_accuracy, n = 1) %>%
            ungroup()

        ## Create an overall accuracy dataframe
        overall_accuracy <- rbind(summary_stats_val, summary_stats_eval_best)
        overall_accuracy <- data.frame(overall_accuracy)
        overall_accuracy$run_type[overall_accuracy$run == "validation"] <- "Validation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_unbalanced"] <- "Evaluation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_1X"] <- "Evaluation"
        overall_accuracy$run_type[overall_accuracy$run == "evaluation_2X"] <- "Evaluation"
        overall_accuracy$key <- paste0(overall_accuracy$celltype, '_', overall_accuracy$run )
        overall_accuracy$group <- "Overall"

        table(overall_accuracy$celltype)

        overall_accuracy_best <- overall_accuracy


        #############
        
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/All_HVG_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add region column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(region = case_when(
                str_detect(donor, "_BA4") ~ "BA4",
                str_detect(donor, "_BA9") ~ "BA9",
                TRUE ~ NA_character_
            ))

        ## Calculate sample wise mean
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype, run, region) %>%
            summarise(
                sample_mean_accuracy = mean(test_accuracy_all),        
                .groups = "drop"               
            )

        ## add group column
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "ALS",
                str_ends(donor, "_2") ~ "FTLD",
                str_ends(donor, "_0") ~ "control",
                TRUE ~ NA_character_
            ))

        BA4_bind_all <- data.frame(BA4_bind_all)
        BA4_bind_all$key <- paste0(BA4_bind_all$celltype, '_', BA4_bind_all$run )
        BA4_bind_all_best <- subset(BA4_bind_all, key %in% overall_accuracy_best$key)
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all_best %>%
            group_by(celltype, run, group) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        summary_stats <- data.frame(summary_stats)
        

        summary_stats$run_type[summary_stats$run == "validation"] <- "Validation"
        summary_stats$run_type[summary_stats$run == "evaluation_unbalanced"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_1X"] <- "Evaluation"
        summary_stats$run_type[summary_stats$run == "evaluation_2X"] <- "Evaluation"

        summary_stats$merge_key <- paste0(summary_stats$run_type, '_', summary_stats$group, '_', summary_stats$celltype )

        summary_stats <- subset(summary_stats, run == "validation")
        summary_stats <- subset(summary_stats, group != "control")

        summary_stats$merge_key <- paste0( summary_stats$group, '_', summary_stats$celltype )

        summary_stats_HVG <- summary_stats

    ##

    
    ###########################
    ## Compute correlation
    ###########################
    TxD <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Pineda_TxD.csv')
    TxD <- TxD[-1,]
    TxD$merge_key <- paste0( TxD$Condition, '_', TxD$celltype )
    TxD <- subset(TxD, celltype != "T_Cell")
    nrow(TxD) == nrow(summary_stats_HVG)

    merge_HVG <- merge(summary_stats_HVG, TxD, by = "merge_key")
    nrow(merge_HVG)

    colors <- c(
        "L3_L5" = "#f3c300",
        "L2_L3" = "#f38400",
        "L4_L6" = "#a1caf1",
        "L4_L5" = "#be0032",
        "L5_L6" = "#c2b280",
        "L5" = "#008856",
        "L6" = "#2b3d26",
        "PV" = "#e25822",
        "5HT3aR" = "#e68fac",
        "Rosehip" = "#0067a5",
        "SOM" = "#f99379",
        "Oligo" = "#604e97",
        "Astro" = "#875692",
        "OPC" = "#f6a600",
        "Micro" = "#b3446c",
        "T_Cell" = "#dcd300",
        "Mural" = "#882d17",
        "Endo" = "#8db600",
        "Fibro" = "#654522")
    
    ggplot(merge_HVG, aes(x = divergence_score/1000, y = mean_accuracy)) + 
        theme_bw() + 
        #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
        geom_point(aes(colour = celltype.x)) +
        theme(
            legend.position = "right",
            legend.key.size = unit(0.2, "cm"),      
            legend.text = element_text(size = 6),   
            legend.title = element_text(size = 7),  
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title = element_text(face = "bold", size = 8),
            strip.text = element_text(size = 8, face = "bold")
        ) +
        scale_colour_manual(values = colors) +
        scale_shape_manual(values = c(16,15)) +
        ylab("Accuracy") +
        xlab("N cells") +
        facet_grid(.~group) +
        stat_cor(
            method = "pearson",
            label.x = 2.75,  # X-coordinate in data units
            label.y = 0.01, # Y-coordinate in data units
            size = 2.5) #+
        #scale_x_log10(
        #    breaks = c(1, 1e1, 1e2, 1e3, 1e4, 1e5),
        #    labels = c("1", "1e1", "1e2", "1e3", "1e4", "1e5")
        #) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 5, width = 5)





    ##
    
    
    ## Method 1 
        # --------------------------------------------------
        # 1. Covariate correction at single-cell level
        # --------------------------------------------------

        # Example covariates in meta.data
        # (replace with what you have: batch, sex, nCount_RNA, percent.mt, etc.)
        vars_to_regress <- c("batch", "nCount_RNA", "percent.mt")

        seu <- ScaleData(seu, vars.to.regress = vars_to_regress)

        # --------------------------------------------------
        # 2. Aggregate into pseudo-bulk profiles
        # --------------------------------------------------
        pseudo_bulk <- AggregateExpression(
            seu_sub,
            group.by = c("CellType", "Sample_ID", "Condition"),
            assays = "RNA",
            slot = "scale.data",   # <-- use covariate-corrected values
            return.seurat = TRUE
        )
        # Extract corrected expression matrix
        expr_mat <- GetAssayData(pseudo_bulk, slot = "scale.data", assay = "RNA")

        # Metadata for each pseudo-bulk profile
        pb_meta <- pseudo_bulk@meta.data

        # --------------------------------------------------
        # 3. Compute median profiles and scores
        # --------------------------------------------------

        celltypes <- unique(pb_meta$CellType)
        conditions <- c("ALS", "FTLD")

        results <- list()

        for (ct in celltypes) {
            message("Processing cell type: ", ct)
            
            # Subset this cell type
            idx <- pb_meta$CellType == ct
            expr_ct <- expr_mat[, idx, drop = FALSE]
            meta_ct <- pb_meta[idx, ]
            
            # Split PN vs ALS vs FTLD
            pn_idx   <- meta_ct$Condition == "PN"
            als_idx  <- meta_ct$Condition == "ALS"
            ftld_idx <- meta_ct$Condition == "FTLD"
            
            pn_median   <- apply(expr_ct[, pn_idx, drop = FALSE], 1, median)
            als_median  <- apply(expr_ct[, als_idx, drop = FALSE], 1, median)
            ftld_median <- apply(expr_ct[, ftld_idx, drop = FALSE], 1, median)
            
            # ------------------------------------------------
            # 4. Compute scores
            # ------------------------------------------------
            
            # TxD = mean absolute change vs PN
            txd_als  <- mean(abs(als_median  - pn_median))
            txd_ftld <- mean(abs(ftld_median - pn_median))
            
            # Divergence = Euclidean distance vs PN
            div_als  <- sqrt(sum((als_median  - pn_median)^2))
            div_ftld <- sqrt(sum((ftld_median - pn_median)^2))
            
            results[[ct]] <- data.frame(
                celltype = ct,
                Condition = c("ALS", "FTLD"),
                TxD = c(txd_als, txd_ftld),
                Divergence = c(div_als, div_ftld)
            )
        }

        scores <- do.call(rbind, results)
        rownames(scores) <- NULL

        print(scores)
    ##
    
    

        
    
    
    
    
    
    
    
    
    
    
    
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

 