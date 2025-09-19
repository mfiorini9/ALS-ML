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
SUMMARY COUNTS COME FROM: Figure_model_accuracy_generalizable_normalize_counts.R

######################
## HVG: BA4 + BA9
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
## NMF: BA4 + BA9
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
## HVG: BA4
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
## NMF: BA4
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
## HVG: BA9
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
## NMF: BA9
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
## Merge performance and plot
######################
## code 
    ## Prep and merge
    performance_HVG$type <- "MCx + FCx"
    performance_HVG$input_type <- "HVGs"

    performance_NMF$type <- "MCx + FCx"
    performance_NMF$input_type <- "NMF"

    performance_HVG_BA4$type <- "MCx"
    performance_HVG_BA4$input_type <- "HVGs"

    performance_HVG_BA9$type <- "FCx"
    performance_HVG_BA9$input_type <- "HVGs"

    performance_NMF_BA4$type <- "MCx"
    performance_NMF_BA4$input_type <- "NMF"

    performance_NMF_BA9$type <- "FCx"
    performance_NMF_BA9$input_type <- "NMF"

    performance_total <- rbind(performance_HVG, performance_NMF, performance_HVG_BA4, performance_HVG_BA9, performance_NMF_BA4, performance_NMF_BA9)
    nrow(performance_total)
    performance_total$celltype[performance_total$celltype == "Rose"] <- "Rosehip"

    ## assign broad levels
    performance_total$major_group[performance_total$celltype == "L2_L3"] <- "Excitatory"
    performance_total$major_group[performance_total$celltype == "L3_L5"] <- "Excitatory"
    performance_total$major_group[performance_total$celltype == "L4_L5"] <- "Excitatory"
    performance_total$major_group[performance_total$celltype == "L4_L6"] <- "Excitatory"
    performance_total$major_group[performance_total$celltype == "L5"] <- "Excitatory"
    performance_total$major_group[performance_total$celltype == "L5_L6"] <- "Excitatory"
    performance_total$major_group[performance_total$celltype == "L6"] <- "Excitatory"

    performance_total$major_group[performance_total$celltype == "5HT3aR"] <- "Inhibitory"
    performance_total$major_group[performance_total$celltype == "PV"] <- "Inhibitory"
    performance_total$major_group[performance_total$celltype == "Rosehip"] <- "Inhibitory"
    performance_total$major_group[performance_total$celltype == "SOM"] <- "Inhibitory"

    performance_total$major_group[performance_total$celltype == "Astro"] <- "non-neuronal"
    performance_total$major_group[performance_total$celltype == "Endo"] <- "non-neuronal"
    performance_total$major_group[performance_total$celltype == "Fibro"] <- "non-neuronal"
    performance_total$major_group[performance_total$celltype == "Micro"] <- "non-neuronal"
    performance_total$major_group[performance_total$celltype == "Mural"] <- "non-neuronal"
    performance_total$major_group[performance_total$celltype == "Oligo"] <- "non-neuronal"
    performance_total$major_group[performance_total$celltype == "OPC"] <- "non-neuronal"

    ## Plot    
    ggplot(performance_total, aes(x = major_group, y = mean_accuracy, fill = major_group)) + 
        theme_bw() + 
        geom_boxplot(outlier.shape = NA) +
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
        facet_nested(.~ input_type + type, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = c("Excitatory" = "#00A8A8",
                                "Inhibitory"  = "#9AD3DA",
                                "non-neuronal" = "#FF924D")) +
        scale_y_continuous(lim = c(0,1)) +
        ylab("Accuracy") #+
        #scale_x_discrete(labels = c("Excitatory", "Inhibitory", "Non-neuronal"))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2.75, width =3.5)



    ### Calculate accuracy
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