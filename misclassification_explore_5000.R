salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

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


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ HVG BA4 and BA9 together
## code

    ###########################
    ## HVG diagnosis
    ###########################

    ## code
        ## read in file
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


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "ALS", "FTLD")) +
            scale_y_discrete(labels = c("Control", "ALS", "FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_diagnosis_prop_missclassification_5000.csv')

    ##

    ###########################
    ## HVG genetic
    ###########################

    ## code
        ## read in file
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

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Create group label
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
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
        df_avg_props$true_class[df_avg_props$group == "Sporadic"] <- "1"
        df_avg_props$true_class[df_avg_props$group == "C9orf72"] <- "2"
        df_avg_props$true_class[df_avg_props$group == "control"] <- "0"

        df_avg_props$predicted_class <- df_avg_props$class

        df_avg_props$celltype <- factor(df_avg_props$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "Sporadic", "C9orf72")) +
            scale_y_discrete(labels = c("Control", "Sporadic", "C9orf72"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_genetic_prop_missclassification_5000.csv')

    ##

    ###########################
    ## HVG disease
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

        
        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD")) +
            scale_y_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_disease_prop_missclassification_5000.csv')
    ##

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ NMF BA4 and BA9 together
## code

    ###########################
    ## NMF diagnosis
    ###########################

    ## code
        ## read in file
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


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "ALS", "FTLD")) +
            scale_y_discrete(labels = c("Control", "ALS", "FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_diagnosis_prop_missclassification_5000.csv')

    ##

    ###########################
    ## NMF genetic
    ###########################

    ## code
        ## read in file
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

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Create group label
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
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
        df_avg_props$true_class[df_avg_props$group == "Sporadic"] <- "1"
        df_avg_props$true_class[df_avg_props$group == "C9orf72"] <- "2"
        df_avg_props$true_class[df_avg_props$group == "control"] <- "0"

        df_avg_props$predicted_class <- df_avg_props$class

        df_avg_props$celltype <- factor(df_avg_props$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "Sporadic", "C9orf72")) +
            scale_y_discrete(labels = c("Control", "Sporadic", "C9orf72"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_genetic_prop_missclassification_5000.csv')

    ##

    ###########################
    ## NMF disease
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

        
        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD")) +
            scale_y_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_disease_prop_missclassification_5000.csv')
    ##

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ HVG BA4
## code

    ###########################
    ## HVG diagnosis
    ###########################

    ## code
        ## read in file
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


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "ALS", "FTLD")) +
            scale_y_discrete(labels = c("Control", "ALS", "FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_diagnosis_prop_missclassification_5000.csv')

    ##

    ###########################
    ## HVG genetic
    ###########################

    ## code
        ## read in file
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

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Create group label
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
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
        df_avg_props$true_class[df_avg_props$group == "Sporadic"] <- "1"
        df_avg_props$true_class[df_avg_props$group == "C9orf72"] <- "2"
        df_avg_props$true_class[df_avg_props$group == "control"] <- "0"

        df_avg_props$predicted_class <- df_avg_props$class

        df_avg_props$celltype <- factor(df_avg_props$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "Sporadic", "C9orf72")) +
            scale_y_discrete(labels = c("Control", "Sporadic", "C9orf72"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_genetic_prop_missclassification_5000.csv')

    ##

    ###########################
    ## HVG disease
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

        
        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD")) +
            scale_y_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_disease_prop_missclassification_5000.csv')
    ##

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ NMF BA4
## code

    ###########################
    ## NMF diagnosis
    ###########################

    ## code
        ## read in file
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


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "ALS", "FTLD")) +
            scale_y_discrete(labels = c("Control", "ALS", "FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_diagnosis_prop_missclassification_5000.csv')

    ##

    ###########################
    ## NMF genetic
    ###########################

    ## code
        ## read in file
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

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Create group label
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
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
        df_avg_props$true_class[df_avg_props$group == "Sporadic"] <- "1"
        df_avg_props$true_class[df_avg_props$group == "C9orf72"] <- "2"
        df_avg_props$true_class[df_avg_props$group == "control"] <- "0"

        df_avg_props$predicted_class <- df_avg_props$class

        df_avg_props$celltype <- factor(df_avg_props$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "Sporadic", "C9orf72")) +
            scale_y_discrete(labels = c("Control", "Sporadic", "C9orf72"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_genetic_prop_missclassification_5000.csv')

    ##

    ###########################
    ## NMF disease
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

        
        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD")) +
            scale_y_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_disease_prop_missclassification_5000.csv')
    ##

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ HVG BA9
## code

    ###########################
    ## HVG diagnosis
    ###########################

    ## code
        ## read in file
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


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "ALS", "FTLD")) +
            scale_y_discrete(labels = c("Control", "ALS", "FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_diagnosis_prop_missclassification_5000.csv')

    ##

    ###########################
    ## HVG genetic
    ###########################

    ## code
        ## read in file
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

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Create group label
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
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
        df_avg_props$true_class[df_avg_props$group == "Sporadic"] <- "1"
        df_avg_props$true_class[df_avg_props$group == "C9orf72"] <- "2"
        df_avg_props$true_class[df_avg_props$group == "control"] <- "0"

        df_avg_props$predicted_class <- df_avg_props$class

        df_avg_props$celltype <- factor(df_avg_props$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "Sporadic", "C9orf72")) +
            scale_y_discrete(labels = c("Control", "Sporadic", "C9orf72"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_genetic_prop_missclassification_5000.csv')

    ##

    ###########################
    ## HVG disease
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

        
        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD")) +
            scale_y_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_disease_prop_missclassification_5000.csv')
    ##

##


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ NMF BA9
## code

    ###########################
    ## NMF diagnosis
    ###########################

    ## code
        ## read in file
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


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "ALS", "FTLD")) +
            scale_y_discrete(labels = c("Control", "ALS", "FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_diagnosis_prop_missclassification_5000.csv')

    ##

    ###########################
    ## NMF genetic
    ###########################

    ## code
        ## read in file
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

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

        ## Create group label
        BA4_bind_all <- BA4_bind_all %>%
            mutate(group = case_when(
                str_ends(donor, "_1") ~ "Sporadic",
                str_ends(donor, "_2") ~ "C9orf72",
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
        df_avg_props$true_class[df_avg_props$group == "Sporadic"] <- "1"
        df_avg_props$true_class[df_avg_props$group == "C9orf72"] <- "2"
        df_avg_props$true_class[df_avg_props$group == "control"] <- "0"

        df_avg_props$predicted_class <- df_avg_props$class

        df_avg_props$celltype <- factor(df_avg_props$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))


        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "Sporadic", "C9orf72")) +
            scale_y_discrete(labels = c("Control", "Sporadic", "C9orf72"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_genetic_prop_missclassification_5000.csv')

    ##

    ###########################
    ## NMF disease
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

        
        df_avg_props$colour <- "black"
        df_avg_props$colour[df_avg_props$mean_proportion > 0.7] = "white"

        temp <- ggplot(df_avg_props, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
            theme_bw() + 
            #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
            geom_tile() +
            #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
            theme(
                legend.position = "none",
                legend.key.size = unit(0.2, "cm"),      
                legend.text = element_text(size = 6),   
                legend.title = element_text(size = 7),  
                panel.grid = element_blank(),
                axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title = element_text(face = "bold", size = 8),
                strip.text = element_text(size = 8, face = "bold")
            ) +
            facet_nested(.~ celltype, scales = "free_x", space = "free_x") +
            #scale_fill_viridis_c(direction = -1) +
            scale_colour_manual(values = rev(c("white", "black"))) +
            scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
            #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
            #scale_colour_manual(values = "black") +
            #scale_shape_manual(values = c(16,15)) +
            ylab("Predicted label") +
            xlab("True label") +
            scale_x_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD")) +
            scale_y_discrete(labels = c("Control", "SALS", "C9ALS", "SFTLD",  "C9FTLD"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 1.6, width = 13)

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_disease_prop_missclassification_5000.csv')
    ##

##

########################
## Mega plot
########################
## code 
    ## Read in summary files. 
    HVG_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_diagnosis_prop_missclassification_5000.csv', sep = ",")
    HVG_diagnosis$model = "Diagnosis"
    HVG_diagnosis$region = "M/FCx"
    HVG_diagnosis$input = "HVG"

    HVG_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_genetic_prop_missclassification_5000.csv', sep = ",")
    HVG_genetic$model = "Genetic group"
    HVG_genetic$region = "M/FCx"
    HVG_genetic$input = "HVG"

    HVG_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_disease_prop_missclassification_5000.csv', sep = ",")
    HVG_disease$model = "Disease group"
    HVG_disease$region = "M/FCx"
    HVG_disease$input = "HVG"


    NMF_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_diagnosis_prop_missclassification_5000.csv', sep = ",")
    NMF_diagnosis$model = "Diagnosis"
    NMF_diagnosis$region = "M/FCx"
    NMF_diagnosis$input = "NMF"

    NMF_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_genetic_prop_missclassification_5000.csv', sep = ",")
    NMF_genetic$model = "Genetic group"
    NMF_genetic$region = "M/FCx"
    NMF_genetic$input = "NMF"

    NMF_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_disease_prop_missclassification_5000.csv', sep = ",")
    NMF_disease$model = "Disease group"
    NMF_disease$region = "M/FCx"
    NMF_disease$input = "NMF"


    BA4_HVG_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_diagnosis_prop_missclassification_5000.csv', sep = ",")
    BA4_HVG_diagnosis$model = "Diagnosis"
    BA4_HVG_diagnosis$region = "MCx"
    BA4_HVG_diagnosis$input = "HVG"

    BA4_HVG_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_genetic_prop_missclassification_5000.csv', sep = ",")
    BA4_HVG_genetic$model = "Genetic group"
    BA4_HVG_genetic$region = "MCx"
    BA4_HVG_genetic$input = "HVG"

    BA4_HVG_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_disease_prop_missclassification_5000.csv', sep = ",")
    BA4_HVG_disease$model = "Disease group"
    BA4_HVG_disease$region = "MCx"
    BA4_HVG_disease$input = "HVG"


    BA9_HVG_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_diagnosis_prop_missclassification_5000.csv', sep = ",")
    BA9_HVG_diagnosis$model = "Diagnosis"
    BA9_HVG_diagnosis$region = "FCx"
    BA9_HVG_diagnosis$input = "HVG"

    BA9_HVG_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_genetic_prop_missclassification_5000.csv', sep = ",")
    BA9_HVG_genetic$model = "Genetic group"
    BA9_HVG_genetic$region = "FCx"
    BA9_HVG_genetic$input = "HVG"

    BA9_HVG_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_disease_prop_missclassification_5000.csv', sep = ",")
    BA9_HVG_disease$model = "Disease group"
    BA9_HVG_disease$region = "FCx"
    BA9_HVG_disease$input = "HVG"


    BA4_NMF_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_diagnosis_prop_missclassification_5000.csv', sep = ",")
    BA4_NMF_diagnosis$model = "Diagnosis"
    BA4_NMF_diagnosis$region = "MCx"
    BA4_NMF_diagnosis$input = "NMF"

    BA4_NMF_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_genetic_prop_missclassification_5000.csv', sep = ",")
    BA4_NMF_genetic$model = "Genetic group"
    BA4_NMF_genetic$region = "MCx"
    BA4_NMF_genetic$input = "NMF"

    BA4_NMF_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_disease_prop_missclassification_5000.csv', sep = ",")
    BA4_NMF_disease$model = "Disease group"
    BA4_NMF_disease$region = "MCx"
    BA4_NMF_disease$input = "NMF"


    BA9_NMF_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_diagnosis_prop_missclassification_5000.csv', sep = ",")
    BA9_NMF_diagnosis$model = "Diagnosis"
    BA9_NMF_diagnosis$region = "FCx"
    BA9_NMF_diagnosis$input = "NMF"

    BA9_NMF_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_genetic_prop_missclassification_5000.csv', sep = ",")
    BA9_NMF_genetic$model = "Genetic group"
    BA9_NMF_genetic$region = "FCx"
    BA9_NMF_genetic$input = "NMF"

    BA9_NMF_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_disease_prop_missclassification_5000.csv', sep = ",")
    BA9_NMF_disease$model = "Disease group"
    BA9_NMF_disease$region = "FCx"
    BA9_NMF_disease$input = "NMF"

    bind_total <- rbind(HVG_diagnosis, HVG_genetic, HVG_disease,
                        NMF_diagnosis, NMF_genetic, NMF_disease,
                        BA4_HVG_diagnosis, BA4_HVG_genetic, BA4_HVG_disease,
                        BA9_HVG_diagnosis, BA9_HVG_genetic, BA9_HVG_disease,
                        BA4_NMF_diagnosis, BA4_NMF_genetic, BA4_NMF_disease,
                        BA9_NMF_diagnosis, BA9_NMF_genetic, BA9_NMF_disease)

    ## SET FACTOR LEVELS
    bind_total$region <- factor(bind_total$region, levels = c("M/FCx", "MCx", "FCx"))
    bind_total$celltype <- factor(bind_total$celltype, levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC'))

    ## Plot diagnosis
    bind_total_model <- subset(bind_total, model == "Diagnosis")
    bind_total_model$true_class[bind_total_model$true_class == "0"] <- "Control"
    bind_total_model$true_class[bind_total_model$true_class == "1"] <- "ALS"
    bind_total_model$true_class[bind_total_model$true_class == "2"] <- "FTLD"

    bind_total_model$predicted_class[bind_total_model$predicted_class == "0"] <- "Control"
    bind_total_model$predicted_class[bind_total_model$predicted_class == "1"] <- "ALS"
    bind_total_model$predicted_class[bind_total_model$predicted_class == "2"] <- "FTLD"

    bind_total_model$predicted_class <- factor(bind_total_model$predicted_class, levels = c("Control", "ALS", "FTLD"))
    bind_total_model$true_class <- factor(bind_total_model$true_class, levels = c("Control", "ALS", "FTLD"))

    temp <- ggplot(bind_total_model, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
                theme_bw() + 
                #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
                geom_tile() +
                #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
                theme(
                    legend.position = "none",
                    legend.key.size = unit(0.2, "cm"),      
                    legend.text = element_text(size = 6),   
                    legend.title = element_text(size = 7),  
                    panel.grid = element_blank(),
                    axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                    axis.text.y = element_text(size = 8, colour = "black"),
                    axis.title = element_text(face = "bold", size = 8),
                    strip.text = element_text(size = 8, face = "bold")
                ) +
                facet_nested( input + region ~ celltype) +
                #scale_fill_viridis_c(direction = -1) +
                scale_colour_manual(values = rev(c("white", "black"))) +
                scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
                #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
                #scale_colour_manual(values = "black") +
                #scale_shape_manual(values = c(16,15)) +
                ylab("Predicted label") +
                xlab("True label") 

            ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 4, width = 13)


    ## Plot genetic
    bind_total_model <- subset(bind_total, model == "Genetic group")
    bind_total_model$true_class[bind_total_model$true_class == "0"] <- "Control"
    bind_total_model$true_class[bind_total_model$true_class == "1"] <- "Sporadic"
    bind_total_model$true_class[bind_total_model$true_class == "2"] <- "C9orf72"

    bind_total_model$predicted_class[bind_total_model$predicted_class == "0"] <- "Control"
    bind_total_model$predicted_class[bind_total_model$predicted_class == "1"] <- "Sporadic"
    bind_total_model$predicted_class[bind_total_model$predicted_class == "2"] <- "C9orf72"

    bind_total_model$predicted_class <- factor(bind_total_model$predicted_class, levels = c("Control", "Sporadic", "C9orf72"))
    bind_total_model$true_class <- factor(bind_total_model$true_class, levels = c("Control", "Sporadic", "C9orf72"))

    temp <- ggplot(bind_total_model, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
                theme_bw() + 
                #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
                geom_tile() +
                #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
                theme(
                    legend.position = "none",
                    legend.key.size = unit(0.2, "cm"),      
                    legend.text = element_text(size = 6),   
                    legend.title = element_text(size = 7),  
                    panel.grid = element_blank(),
                    axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                    axis.text.y = element_text(size = 8, colour = "black"),
                    axis.title = element_text(face = "bold", size = 8),
                    strip.text = element_text(size = 8, face = "bold")
                ) +
                facet_nested( input + region ~ celltype) +
                #scale_fill_viridis_c(direction = -1) +
                scale_colour_manual(values = rev(c("white", "black"))) +
                scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
                #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
                #scale_colour_manual(values = "black") +
                #scale_shape_manual(values = c(16,15)) +
                ylab("Predicted label") +
                xlab("True label") 

            ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 4, width = 13)


    ## Plot disease
    bind_total_model <- subset(bind_total, model == "Disease group")
    bind_total_model$true_class[bind_total_model$true_class == "0"] <- "Control"
    bind_total_model$true_class[bind_total_model$true_class == "1"] <- "SALS"
    bind_total_model$true_class[bind_total_model$true_class == "2"] <- "C9ALS"
    bind_total_model$true_class[bind_total_model$true_class == "3"] <- "SFTLD"
    bind_total_model$true_class[bind_total_model$true_class == "4"] <- "C9FTLD"

    bind_total_model$predicted_class[bind_total_model$predicted_class == "0"] <- "Control"
    bind_total_model$predicted_class[bind_total_model$predicted_class == "1"] <- "SALS"
    bind_total_model$predicted_class[bind_total_model$predicted_class == "2"] <- "C9ALS"
    bind_total_model$predicted_class[bind_total_model$predicted_class == "3"] <- "SFTLD"
    bind_total_model$predicted_class[bind_total_model$predicted_class == "4"] <- "C9FTLD"

    bind_total_model$predicted_class <- factor(bind_total_model$predicted_class, levels = c("Control", "SALS", "C9ALS", "SFTLD", "C9FTLD"))
    bind_total_model$true_class <- factor(bind_total_model$true_class, levels = c("Control", "SALS", "C9ALS", "SFTLD", "C9FTLD"))

    temp <- ggplot(bind_total_model, aes(x = true_class, y = predicted_class, fill = mean_proportion)) + 
                theme_bw() + 
                #geom_hline(yintercept = 0.6, colour = "red", linetype = "dashed") +
                geom_tile() +
                #geom_text(aes(label = round(mean_proportion, 2), colour = colour), size = 2.25) +
                theme(
                    legend.position = "none",
                    legend.key.size = unit(0.2, "cm"),      
                    legend.text = element_text(size = 6),   
                    legend.title = element_text(size = 7),  
                    panel.grid = element_blank(),
                    axis.text.x = element_text(size = 8, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
                    axis.text.y = element_text(size = 8, colour = "black"),
                    axis.title = element_text(face = "bold", size = 8),
                    strip.text = element_text(size = 8, face = "bold")
                ) +
                facet_nested( input + region ~ celltype) +
                #scale_fill_viridis_c(direction = -1) +
                scale_colour_manual(values = rev(c("white", "black"))) +
                scale_fill_gradient(low = "white", high = "darkred", limits = c(0, 1)) +   # <-- white = low, darkred = high
                #scale_size_continuous(range = c(0.1, 6)) +   # <- make lower end smaller
                #scale_colour_manual(values = "black") +
                #scale_shape_manual(values = c(16,15)) +
                ylab("Predicted label") +
                xlab("True label") 

            ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 4.4, width = 13)
##