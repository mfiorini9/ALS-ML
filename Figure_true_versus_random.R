salloc -A def-groulea --time=0-8 -c 1 --mem=10g

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


WE ARE USING THE NORMALIZED CELL COUNTS.
TRUE SUMMARY COUNTS COME FROM: Figure_model_accuracy_generalizable_normalize_counts.R
RANDOM SUMMARY COUNTS COME FROM: Figure_model_accuracy_generalizable_normalize_counts_random.R

######################
## HVG: BA4 + BA9 TRUE
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5000_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_HVG <- performance
##

######################
## HVG: BA4 + BA9 RANDOM
######################
## code
    # Read in file
    performance_random <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_5000_cell_random_label.csv')
    
    # Only retain overall performance
    performance_random <- subset(performance_random, group == "Overall")

    # Retain columns of interest
    performance_random <- performance_random %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance_random)

    performance_HVG_random <- performance_random
##

######################
## NMF: BA4 + BA9 TRUE
######################
## code
    # Read in file
    performance_NMF <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_5000_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance_NMF <- subset(performance_NMF, group == "Overall")

    # Retain columns of interest
    performance_NMF <- performance_NMF %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance_NMF)
##

######################
## NMF: BA4 + BA9 RANDOM
######################
## code
    # Read in file
    performance_NMF_random <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_5000_cell_random_label.csv')
    
    # Only retain overall performance
    performance_NMF_random <- subset(performance_NMF_random, group == "Overall")

    # Retain columns of interest
    performance_NMF_random <- performance_NMF_random %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance_NMF_random)
##


######################
## HVG: BA4 TRUE
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_5000_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_HVG_BA4 <- performance
##


######################
## HVG: BA4 RANDOM
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_5000_cell_random_label.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_HVG_BA4_random <- performance
##


######################
## NMF: BA4 TRUE
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_5000_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_NMF_BA4 <- performance
##


######################
## NMF: BA4 RANDOM
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_5000_cell_random_label.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_NMF_BA4_random <- performance
##


######################
## HVG: BA9 TRUE
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_5000_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_HVG_BA9 <- performance
##

######################
## HVG: BA9 RANDOM
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_5000_cell_random_label.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_HVG_BA9_random <- performance
##

######################
## NMF: BA9 TRUE
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_5000_final_performance_cell_level_top.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_NMF_BA9 <- performance
##

######################
## NMF: BA9 RANDOM
######################
## code
    # Read in file
    performance <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_5000_cell_random_label.csv')
    
    # Only retain overall performance
    performance <- subset(performance, group == "Overall")

    # Retain columns of interest
    performance <- performance %>% dplyr::select(input, celltype, mean_accuracy)
    nrow(performance)

    performance_NMF_BA9_random <- performance
##





######################
## Merge 
######################
## code 
    ## Prep and merge
    performance_HVG$type <- "MCx + FCx"
    performance_HVG$input_type <- "HVGs"

    performance_HVG_random$type <- "MCx + FCx"
    performance_HVG_random$input_type <- "HVGs"

    performance_NMF$type <- "MCx + FCx"
    performance_NMF$input_type <- "NMF"

    performance_NMF_random$type <- "MCx + FCx"
    performance_NMF_random$input_type <- "NMF"

    performance_HVG_BA4$type <- "MCx"
    performance_HVG_BA4$input_type <- "HVGs"

    performance_HVG_BA4_random$type <- "MCx"
    performance_HVG_BA4_random$input_type <- "HVGs"

    performance_HVG_BA9$type <- "FCx"
    performance_HVG_BA9$input_type <- "HVGs"

    performance_HVG_BA9_random$type <- "FCx"
    performance_HVG_BA9_random$input_type <- "HVGs"

    performance_NMF_BA4$type <- "MCx"
    performance_NMF_BA4$input_type <- "NMF"

    performance_NMF_BA4_random$type <- "MCx"
    performance_NMF_BA4_random$input_type <- "NMF"

    performance_NMF_BA9$type <- "FCx"
    performance_NMF_BA9$input_type <- "NMF"

    performance_NMF_BA9_random$type <- "FCx"
    performance_NMF_BA9_random$input_type <- "NMF"

    
    
    performance_total <- rbind(performance_HVG, performance_HVG_random,
                                performance_NMF, performance_NMF_random,
                                performance_HVG_BA4, performance_HVG_BA4_random,
                                performance_HVG_BA9, performance_HVG_BA9_random,
                                performance_NMF_BA4, performance_NMF_BA4_random,
                                performance_NMF_BA9, performance_NMF_BA9_random
                                )
    nrow(performance_total)

    performance_total$input[performance_total$input == "3 Class"] <- "3 class random"
    performance_total$input[performance_total$input == "5 Class"] <- "5 class random"

    performance_total$input <- factor(performance_total$input, levels = c("Diagnosis", "Genetic group", "Disease group", "3 class random", "5 class random"))
##

######################
## Plot performance
######################
## code 
    ## Plot    
    ggplot(performance_total, aes(x = input, y = mean_accuracy, fill = input)) + 
        theme_bw() + 
        geom_boxplot(outlier.shape = NA) +
        geom_vline(xintercept = 3.5, colour = "red", linetype = "dashed", linewidth = 0.25) +
        #geom_line(data = means_df, aes(x = major_group, y = mean_val, group = input),
        #    inherit.aes = FALSE, color = "red", linewidth = 0.4) +
        #stat_summary(fun = mean, geom = "point", shape = 23, size = 1, fill = "red", colour = "red") +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", size = 8),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 6.5, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        #facet_grid(~input, scales = "free_x") +
        facet_nested(.~ input_type, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = c("#ec7014", "#f6a437", "#ffd366", "#758a9b", "#b0b0b0")) +
        scale_y_continuous(lim = c(0,1)) +
        ylab("Accuracy") #+
        #scale_x_discrete(labels = c("Excitatory", "Inhibitory", "Non-neuronal"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2.75, width =2)

    ## Calculate accuracy
    #summary_stats_group <- overall_accuracy %>%
    #    group_by(major_group) %>%
    #    summarise(
    #        median_accuracy = median(mean_accuracy), 
    #        mean_accuracy = mean(mean_accuracy),           # Calculate the median of X0
    #        sd_accuracy = sd(mean_accuracy),                  # Calculate the standard deviation of X0
    #        mad_accuracy = median(abs(mean_accuracy - median(mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
    #        .groups = "drop"                 # Drop the grouping after summarising
    #    )
##

######################
## Prepare fold change
######################
## code

    ## HVG_both_diagnosis
        test <- subset(performance_total, input_type == "HVGs")
        test <- subset(test, type == "MCx + FCx")
        test <- subset(test, input %in% c("Diagnosis", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(Diagnosis / `3 class random`)
        )

        HVG_both_diagnosis <- test
    ##

    ## HVG_both_genetic
        test <- subset(performance_total, input_type == "HVGs")
        test <- subset(test, type == "MCx + FCx")
        test <- subset(test, input %in% c("Genetic group", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Genetic group` / `3 class random`)
        )

        HVG_both_genetic <- test
    ##

    ## HVG_both_disease
        test <- subset(performance_total, input_type == "HVGs")
        test <- subset(test, type == "MCx + FCx")
        test <- subset(test, input %in% c("Disease group", "5 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Disease group` / `5 class random`)
        )

        HVG_both_disease <- test
    ##

    ## HVG_BA4_diagnosis
        test <- subset(performance_total, input_type == "HVGs")
        test <- subset(test, type == "MCx")
        test <- subset(test, input %in% c("Diagnosis", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(Diagnosis / `3 class random`)
        )

        HVG_BA4_diagnosis <- test
    ##

    ## HVG_BA4_genetic
        test <- subset(performance_total, input_type == "HVGs")
        test <- subset(test, type == "MCx")
        test <- subset(test, input %in% c("Genetic group", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Genetic group` / `3 class random`)
        )

        HVG_BA4_genetic <- test
    ##

    ## HVG_BA4_disease
        test <- subset(performance_total, input_type == "HVGs")
        test <- subset(test, type == "MCx")
        test <- subset(test, input %in% c("Disease group", "5 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Disease group` / `5 class random`)
        )

        HVG_BA4_disease <- test
    ##

    ## HVG_BA9_diagnosis
        test <- subset(performance_total, input_type == "HVGs")
        test <- subset(test, type == "FCx")
        test <- subset(test, input %in% c("Diagnosis", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(Diagnosis / `3 class random`)
        )

        HVG_BA9_diagnosis <- test
    ##

    ## HVG_BA9_genetic
        test <- subset(performance_total, input_type == "HVGs")
        test <- subset(test, type == "FCx")
        test <- subset(test, input %in% c("Genetic group", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Genetic group` / `3 class random`)
        )

        HVG_BA9_genetic <- test
    ##

    ## HVG_BA9_disease
        test <- subset(performance_total, input_type == "HVGs")
        test <- subset(test, type == "FCx")
        test <- subset(test, input %in% c("Disease group", "5 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Disease group` / `5 class random`)
        )

        HVG_BA9_disease <- test
    ##

    ## NMF_both_diagnosis
        test <- subset(performance_total, input_type == "NMF")
        test <- subset(test, type == "MCx + FCx")
        test <- subset(test, input %in% c("Diagnosis", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(Diagnosis / `3 class random`)
        )

        NMF_both_diagnosis <- test
    ##

    ## NMF_both_genetic
        test <- subset(performance_total, input_type == "NMF")
        test <- subset(test, type == "MCx + FCx")
        test <- subset(test, input %in% c("Genetic group", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Genetic group` / `3 class random`)
        )

        NMF_both_genetic <- test
    ##

    ## NMF_both_disease
        test <- subset(performance_total, input_type == "NMF")
        test <- subset(test, type == "MCx + FCx")
        test <- subset(test, input %in% c("Disease group", "5 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Disease group` / `5 class random`)
        )

        NMF_both_disease <- test
    ##

    ## NMF_BA4_diagnosis
        test <- subset(performance_total, input_type == "NMF")
        test <- subset(test, type == "MCx")
        test <- subset(test, input %in% c("Diagnosis", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(Diagnosis / `3 class random`)
        )

        NMF_BA4_diagnosis <- test
    ##

    ## NMF_BA4_genetic
        test <- subset(performance_total, input_type == "NMF")
        test <- subset(test, type == "MCx")
        test <- subset(test, input %in% c("Genetic group", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Genetic group` / `3 class random`)
        )

        NMF_BA4_genetic <- test
    ##

    ## NMF_BA4_disease
        test <- subset(performance_total, input_type == "NMF")
        test <- subset(test, type == "MCx")
        test <- subset(test, input %in% c("Disease group", "5 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Disease group` / `5 class random`)
        )

        NMF_BA4_disease <- test
    ##

    ## NMF_BA9_diagnosis
        test <- subset(performance_total, input_type == "NMF")
        test <- subset(test, type == "FCx")
        test <- subset(test, input %in% c("Diagnosis", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(Diagnosis / `3 class random`)
        )

        NMF_BA9_diagnosis <- test
    ##

    ## NMF_BA9_genetic
        test <- subset(performance_total, input_type == "NMF")
        test <- subset(test, type == "FCx")
        test <- subset(test, input %in% c("Genetic group", "3 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Genetic group` / `3 class random`)
        )

        NMF_BA9_genetic <- test
    ##

    ## NMF_BA9_disease
        test <- subset(performance_total, input_type == "NMF")
        test <- subset(test, type == "FCx")
        test <- subset(test, input %in% c("Disease group", "5 class random"))

        test <- test %>%
        select(input, celltype, mean_accuracy) %>%
        pivot_wider(names_from = input, values_from = mean_accuracy) %>%
        mutate(
            log2FC = log2(`Disease group` / `5 class random`)
        )

        NMF_BA9_disease <- test
    ##

##


######################
## Plot fold change
######################
## code
    df_names <- c(
    "HVG_both_diagnosis", "HVG_both_genetic", "HVG_both_disease",
    "HVG_BA4_diagnosis",  "HVG_BA4_genetic",  "HVG_BA4_disease",
    "HVG_BA9_diagnosis",  "HVG_BA9_genetic",  "HVG_BA9_disease",
    "NMF_both_diagnosis", "NMF_both_genetic", "NMF_both_disease",
    "NMF_BA4_diagnosis",  "NMF_BA4_genetic",  "NMF_BA4_disease",
    "NMF_BA9_diagnosis",  "NMF_BA9_genetic",  "NMF_BA9_disease"
    )

    # loop through and bind
    df_all <- lapply(df_names, function(nm) {
    df <- get(nm)  # pull dataframe by name
    
    parts <- str_split(nm, "_")[[1]]
    tibble(df,
            input  = parts[1],
            region = parts[2],
            model  = parts[3])
    }) %>%
    bind_rows()

    df_all <- data.frame(df_all)

    df_all <- df_all %>% dplyr::select(celltype, log2FC, input, region, model)

    ## Calculate accuracy
    summary_stats_group <- df_all %>%
        group_by(celltype) %>%
        summarise(
            mean_accuracy = mean(log2FC, na.rm = TRUE),           # Calculate the median of X0
            .groups = "drop"                 # Drop the grouping after summarising
    )

    df_ordered <- summary_stats_group %>%
        arrange(desc(mean_accuracy))

    df_all$celltype <- factor(df_all$celltype, levels = rev(unique(df_ordered$celltype)))

    ## Plot
    ggplot(df_all, aes(x = model, y = celltype, fill = log2FC)) +
        theme_bw() + 
        geom_tile() +
        #geom_text(aes(label = sig), size = 3, colour = "black") +
        theme(
            legend.position = "right",
            panel.grid = element_blank(),
            axis.text.x = element_text(
            colour = "black",
            size = 8,
            angle = 90,
            hjust = 1,   # aligns the end of text with tick
            vjust = 0.5    # pushes text closer to axis
            ),
            #axis.ticks.x = element_blank(),
            axis.text.y = element_text(colour = "black", size = 8),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 6.5, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
            facet_nested(.~ input + region, scales = "free_x", space = "free_x") +
        ylab("Accuracy") +
        #scale_size_continuous(range = c(1, 5))+
        scale_x_discrete(expand = c(0,0)) +
        scale_y_discrete(expand = c(0,0)) +
        scale_fill_gradientn(colors = brewer.pal(9, "YlGn"), na.value = "darkgrey" )

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3.5, width = 4)
##

######################
## Plot correlation against TxD
######################
## code All cells diagnosis
    # Read in file
    TxD_1 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_ALS_FTLD_control_python.csv')
    
    # Retain columns of interest
    TxD_1 <- TxD_1 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_1 <- unique(TxD_1)

    # Add input info
    TxD_1$input <- "Diagnosis"

    # Add type info
    TxD_1$type <- "All cells"


    # change column names
    colnames(TxD_1) <- c("celltype", "Total_divergence", "input", "type" )
##

## code All cells genetic
    # Read in file
    TxD_2 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_C9_sporadic_control_python.csv')
    
    # Retain columns of interest
    TxD_2 <- TxD_2 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_2 <- unique(TxD_2)

    # Add input info
    TxD_2$input <- "Genetic group"

    # Add type info
    TxD_2$type <- "All cells"

    # change column names
    colnames(TxD_2) <- c("celltype", "Total_divergence", "input", "type" )
##

## code All cells disease
    # Read in file
    TxD_3 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_subtypes_control_python.csv')
    
    # Retain columns of interest
    TxD_3 <- TxD_3 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_3 <- unique(TxD_3)

    # Add input info
    TxD_3$input <- "Disease group"

    # Add type info
    TxD_3$type <- "All cells"

    # change column names
    colnames(TxD_3) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA4 diagnosis
    # Read in file
    TxD_1_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_ALS_FTLD_control_python_BA4.csv')
    
    # Retain columns of interest
    TxD_1_BA4 <- TxD_1_BA4 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_1_BA4 <- unique(TxD_1_BA4)

    # Add input info
    TxD_1_BA4$input <- "Diagnosis"

    # Add type info
    TxD_1_BA4$type <- "BA4"


    # change column names
    colnames(TxD_1_BA4) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA4 genetic
    # Read in file
    TxD_2_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_sporadic_C9_control_python_BA4.csv')
    
    # Retain columns of interest
    TxD_2_BA4 <- TxD_2_BA4 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_2_BA4 <- unique(TxD_2_BA4)

    # Add input info
    TxD_2_BA4$input <- "Genetic group"

    # Add type info
    TxD_2_BA4$type <- "BA4"

    # change column names
    colnames(TxD_2_BA4) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA4 disease
    # Read in file
    TxD_3_BA4 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_subtypes_control_python_BA4.csv')
    
    # Retain columns of interest
    TxD_3_BA4 <- TxD_3_BA4 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_3_BA4 <- unique(TxD_3_BA4)

    # Add input info
    TxD_3_BA4$input <- "Disease group"

    # Add type info
    TxD_3_BA4$type <- "BA4"

    # change column names
    colnames(TxD_3_BA4) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA9 diagnosis
    # Read in file
    TxD_1_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_ALS_FTLD_control_python_BA9.csv')
    
    # Retain columns of interest
    TxD_1_BA9 <- TxD_1_BA9 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_1_BA9 <- unique(TxD_1_BA9)

    # Add input info
    TxD_1_BA9$input <- "Diagnosis"

    # Add type info
    TxD_1_BA9$type <- "BA9"


    # change column names
    colnames(TxD_1_BA9) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA9 genetic
    # Read in file
    TxD_2_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_sporadic_C9_control_python_BA9.csv')
    
    # Retain columns of interest
    TxD_2_BA9 <- TxD_2_BA9 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_2_BA9 <- unique(TxD_2_BA9)

    # Add input info
    TxD_2_BA9$input <- "Genetic group"

    # Add type info
    TxD_2_BA9$type <- "BA9"

    # change column names
    colnames(TxD_2_BA9) <- c("celltype", "Total_divergence", "input", "type" )
##

## code BA9 disease
    # Read in file
    TxD_3_BA9 <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/TxD_subtypes_control_python_BA9.csv')
    
    # Retain columns of interest
    TxD_3_BA9 <- TxD_3_BA9 %>% dplyr::select(CellType, Total_divergence)

    # Remove duplicated rows
    TxD_3_BA9 <- unique(TxD_3_BA9)

    # Add input info
    TxD_3_BA9$input <- "Disease group"

    # Add type info
    TxD_3_BA9$type <- "BA9"

    # change column names
    colnames(TxD_3_BA9) <- c("celltype", "Total_divergence", "input", "type" )
##

## Code: bind and merge
    # bind
    TxD_all <- rbind(TxD_1, TxD_2, TxD_3,
                    TxD_1_BA4, TxD_2_BA4, TxD_3_BA4,
                    TxD_1_BA9, TxD_2_BA9, TxD_3_BA9)
    
##

    ## Calculate accuracy
    summary_stats_group <- TxD_all %>%
        group_by(celltype) %>%
        summarise(
            mean_TxD = mean(Total_divergence, na.rm = TRUE),           # Calculate the median of X0
            .groups = "drop"                 # Drop the grouping after summarising
    )

    df_ordered_TxD <- summary_stats_group %>%
        arrange(desc(mean_TxD))



df_ordered_all <- merge(df_ordered, df_ordered_TxD, by = "celltype")

ggplot(df_ordered_all, aes(x = mean_TxD, y = mean_accuracy)) +
    geom_point() +
    stat_cor(
        method = "pearson",
        #label.x = 2.75,  # X-coordinate in data units
        #label.y = 0.01, # Y-coordinate in data units
        size = 2.5) 
    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3.5, width = 4)

