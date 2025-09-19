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

###########################
## Dataset mapping
###########################
## code 
    donors_pineda <- c( "101_1", "132_1","108_1","131_1","114_1","140_1","209_2", "240_2",
                        "127_1", "216_2", "106_1", "135_1","311_0", "234_2", "307_0", "228_2", 
                        "301_0", "317_0", "319_0","129_1", "217_2", "138_1" , "130_1", "222_2", 
                        "206_2", "323_0", "238_2", "102_1", "322_0", "303_0", "134_1", "236_2",
                        "318_0", "122_1", "210_2", "112_1", "137_1", "110_1", "111_1", "128_1",
                        "214_2", "120_1", "118_1", "232_2", "328_0", "226_2", "229_2", "324_0",
                        "235_2", "103_1", "309_0", "204_2", "113_1", "306_0", "207_2", "304_0",
                        "324_0", "235_2", "103_1", "309_0", "204_2", "113_1", "306_0", "207_2",
                        "304_0", "239_2", "225_2", "115_1", "202_2", "218_2", "239_2", "225_2",
                        "115_1", "202_2", "218_2", "126_1", "116_1", "136_1", "117_1", "139_1",
                        "124_1", "302_0", "325_0", "133_1", "211_2", "104_1", "237_2")
    
    donors_li <- c("C3_0", "A2_1", "A1_1", "F4_2","F5_2","C4_0", "F3_2", "C5_0", "C2_0","C1_0", "A5_1",
                    "A6_1", "A4_1", "C6_0", "F1_2", "A3_1", "F2_2")

    donors_limone <- c("ALS_MotorCortex.012318_FC17.RDS_1", "Control_MotorCortex.121417_FC11.RDS_0",
                    "Control_MotorCortex.012218_FC13.RDS_0", "ALS_MotorCortex.012218_FC21.RDS_1", 
                    "ALS_MotorCortex.121417_FC19.RDS_1", "Control_MotorCortex.012218_FC12.RDS_0")
    
    
    donors_pineda <- sub("_[^_]+$", "", donors_pineda)
    donors_li <- sub("_[^_]+$", "", donors_li)
    donors_limone <- sub("_[^_]+$", "", donors_limone)
##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ BA4 and BA9 together
## code

                                                          
    ###########################
    ## HVG: Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Diagnosis"
        summary_stats_diagnosis <- BA4_bind_all
    ##

    ###########################
    ## HVG: Genetic group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Genetic"
        summary_stats_genetic <- BA4_bind_all
    ##

    ###########################
    ## HVG: Disease group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Genetic"
        summary_stats_genetic <- BA4_bind_all
    ##

    ###########################
    ## HVG: Merge 
    ###########################
    ## code
        total_bind_HVG <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_disease)
        total_bind_HVG$input_type <- "HVG"
    ##

    ###########################
    ## NMF: Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Diagnosis"
        summary_stats_diagnosis <- BA4_bind_all
    ##

    ###########################
    ## NMF: Genetic group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Genetic"
        summary_stats_genetic <- BA4_bind_all
    ##

    ###########################
    ## NMF: Disease group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        BA4_bind_all$input <- "Disease"
        summary_stats_disease <- BA4_bind_all
    ##

    ###########################
    ## NMF: Merge 
    ###########################
    ## code
        total_bind_NMF <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_disease)
        total_bind_NMF$input_type <- "NMF"
    ##
   
    ###########################
    ## Merge and plot
    ###########################
    ## code
        total_bind <- rbind(total_bind_HVG, total_bind_NMF)

        total_bind_all <- total_bind

        total_bind <- data.frame(total_bind)
        
        total_bind$dataset <- factor(total_bind$dataset, levels = c("Pineda et al.", "Li et al.", "Limone et al."))
        total_bind$input <- factor(total_bind$input, levels = c("Diagnosis", "Genetic", "Disease"))

        ggplot(total_bind, aes(x = dataset, y = test_accuracy_all, fill = dataset)) + 
            theme_bw() + 
            geom_boxplot(outlier.shape = NA) +
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
            facet_nested(.~ input_type + input, scales = "free_x", space = "free_x") +
            scale_fill_manual(values = c("Pineda et al." = "#FF7F00",
                                    "Li et al."  = "#FDBF6F",
                                    "Limone et al." = "#B15928")) +
            ylab("Accuracy") +
            scale_x_discrete(labels = c("Pineda et al.", "Li et al.", "Limone et al."))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width =4)

        ## Calculate dataset-specific overall accuracy
        #overall_summary_stats <- total_bind %>%
        #    group_by(dataset, input) %>%
        #    summarise(
        #        median_accuracy = median(sample_mean_accuracy), 
        #        mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
        #        sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
        #        mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
        #        .groups = "drop"                 # Drop the grouping after summarising
        #    )
    ##
##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ BA4 
## code
                                                     
    ###########################
    ## HVG: Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_diagnosis.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Diagnosis"
        summary_stats_diagnosis <- BA4_bind_all
    ##

    ###########################
    ## HVG: Genetic group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_genetic.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Genetic"
        summary_stats_genetic <- BA4_bind_all
    ##

    ###########################
    ## HVG: Disease group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_disease.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        BA4_bind_all$input <- "Disease"
        summary_stats_disease <- BA4_bind_all
    ##

    ###########################
    ## HVG: Merge 
    ###########################
    ## code
        total_bind_HVG <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_disease)
        total_bind_HVG$input_type <- "HVG"
    ##

    ###########################
    ## NMF: Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_diagnosis.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Diagnosis"
        summary_stats_diagnosis <- BA4_bind_all
    ##

    ###########################
    ## NMF: Genetic group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_genetic.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Genetic"
        summary_stats_genetic <- BA4_bind_all
    ##

    ###########################
    ## NMF: Disease group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_disease.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        BA4_bind_all$input <- "Disease"
        summary_stats_disease <- BA4_bind_all
    ##

    ###########################
    ## NMF: Merge 
    ###########################
    ## code
        total_bind_NMF <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_disease)
        total_bind_NMF$input_type <- "NMF"
    ##
   
    ###########################
    ## Merge and plot
    ###########################
    ## code
        total_bind <- rbind(total_bind_HVG, total_bind_NMF)

        total_bind_BA4 <- total_bind
        
        total_bind <- data.frame(total_bind)
        
        total_bind$dataset <- factor(total_bind$dataset, levels = c("Pineda et al.", "Li et al.", "Limone et al."))
        total_bind$input <- factor(total_bind$input, levels = c("Diagnosis", "Genetic", "Disease"))

        ggplot(total_bind, aes(x = dataset, y = test_accuracy_all, fill = dataset)) + 
            theme_bw() + 
            geom_boxplot(outlier.shape = NA) +
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
            facet_nested(.~ input_type + input, scales = "free_x", space = "free_x") +
            scale_fill_manual(values = c("Pineda et al." = "#FF7F00",
                                    "Li et al."  = "#FDBF6F",
                                    "Limone et al." = "#B15928")) +
            ylab("Accuracy") +
            scale_x_discrete(labels = c("Pineda et al.", "Li et al.", "Limone et al."))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width =4)

        ## Calculate dataset-specific overall accuracy
        #overall_summary_stats <- total_bind %>%
        #    group_by(dataset, input) %>%
        #    summarise(
        #        median_accuracy = median(sample_mean_accuracy), 
        #        mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
        #        sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
        #        mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
        #        .groups = "drop"                 # Drop the grouping after summarising
        #    )
    ##
##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ BA9
## code
                                                     
    ###########################
    ## HVG: Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_diagnosis.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Diagnosis"
        summary_stats_diagnosis <- BA4_bind_all
    ##

    ###########################
    ## HVG: Genetic group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_genetic.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Genetic"
        summary_stats_genetic <- BA4_bind_all
    ##

    ###########################
    ## HVG: Disease group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_disease.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        BA4_bind_all$input <- "Disease"
        summary_stats_disease <- BA4_bind_all
    ##

    ###########################
    ## HVG: Merge 
    ###########################
    ## code
        total_bind_HVG <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_disease)
        total_bind_HVG$input_type <- "HVG"
    ##

    ###########################
    ## NMF: Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_diagnosis.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Diagnosis"
        summary_stats_diagnosis <- BA4_bind_all
    ##

    ###########################
    ## NMF: Genetic group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_genetic.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        BA4_bind_all$input <- "Genetic"
        summary_stats_genetic <- BA4_bind_all
    ##

    ###########################
    ## NMF: Disease group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_disease.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_BA.*", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset, n_cells) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Only retain samples with at least 25 cells
        #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

        BA4_bind_all$input <- "Disease"
        summary_stats_disease <- BA4_bind_all
    ##

    ###########################
    ## NMF: Merge 
    ###########################
    ## code
        total_bind_NMF <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_disease)
        total_bind_NMF$input_type <- "NMF"
    ##
   
    ###########################
    ## Merge and plot
    ###########################
    ## code
        total_bind <- rbind(total_bind_HVG, total_bind_NMF)

        total_bind_BA9 <- total_bind
        
        total_bind <- data.frame(total_bind)
        
        total_bind$dataset <- factor(total_bind$dataset, levels = c("Pineda et al.", "Li et al.", "Limone et al."))
        total_bind$input <- factor(total_bind$input, levels = c("Diagnosis", "Genetic", "Disease"))

        ggplot(total_bind, aes(x = dataset, y = test_accuracy_all, fill = dataset)) + 
            theme_bw() + 
            geom_boxplot(outlier.shape = NA) +
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
            facet_nested(.~ input_type + input, scales = "free_x", space = "free_x") +
            scale_fill_manual(values = c("Pineda et al." = "#FF7F00",
                                    "Li et al."  = "#FDBF6F",
                                    "Limone et al." = "#B15928")) +
            ylab("Accuracy") +
            scale_x_discrete(labels = c("Pineda et al.", "Li et al.", "Limone et al."))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 3, width =4)

        ## Calculate dataset-specific overall accuracy
        #overall_summary_stats <- total_bind %>%
        #    group_by(dataset, input) %>%
        #    summarise(
        #        median_accuracy = median(sample_mean_accuracy), 
        #        mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
        #        sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
        #        mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
        #        .groups = "drop"                 # Drop the grouping after summarising
        #    )
    ##
##


###########################
## Plot all together
###########################
## code
    total_bind <- rbind(total_bind_all, total_bind_BA4, total_bind_BA9 )

    overall_summary_stats <- total_bind %>%
        group_by(donor, dataset, input_type) %>%
        summarise(
            median_accuracy = median(test_accuracy_all), 
            mean_accuracy = mean(test_accuracy_all),           # Calculate the median of X0
            sd_accuracy = sd(test_accuracy_all),                  # Calculate the standard deviation of X0
            mad_accuracy = median(abs(test_accuracy_all - median(test_accuracy_all))),  # Calculate MAD: median of absolute deviations from the median
            .groups = "drop"                 # Drop the grouping after summarising
        )

    overall_summary_stats$dataset <- factor(overall_summary_stats$dataset, levels = c("Pineda et al.", "Li et al.", "Limone et al."))

    ggplot(overall_summary_stats, aes(x = dataset, y = mean_accuracy, fill = dataset)) + 
        theme_bw() + 
        geom_boxplot(outlier.shape = NA) +
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
        facet_nested(.~ input_type, scales = "free_x", space = "free_x") +
        scale_fill_manual(values = c("Pineda et al." = "#FF7F00",
                                "Li et al."  = "#FDBF6F",
                                "Limone et al." = "#B15928")) +
        ylab("Accuracy") +
        scale_x_discrete(labels = c("Pineda et al.", "Li et al.", "Limone et al.")) +
        scale_y_continuous(lim = c(0,1)) #+
        #stat_compare_means(
        #    method = "wilcox.test", 
        #    comparisons = list(c("Pineda et al.", "Li et al."), 
        #                    c("Pineda et al.", "Limone et al."), 
        #                    c("Li et al.", "Limone et al.")),
        #    label = "p.signif",
        #    size = 1.5
        #)

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 2.75, width =2)


##

####################################################################################################################################

NOT USING BEYOND THIS POINT.

total_bind_all
total_bind_BA9 <- total_bind







Need to calculate cell wise.
Need to calculate sample wise.
Need to calculate sample wise cell aggregate. 

Consider a bar plot here.


Not using beyond this point. 

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ BA4 and BA9 together -- 5,000 cell 
## code
    ###########################
    ## Dataset mapping
    ###########################
    ## code 
        donors_pineda <- c( "101_1", "132_1","108_1","131_1","114_1","140_1","209_2", "240_2",
                            "127_1", "216_2", "106_1", "135_1","311_0", "234_2", "307_0", "228_2", 
                            "301_0", "317_0", "319_0","129_1", "217_2", "138_1" , "130_1", "222_2", 
                            "206_2", "323_0", "238_2", "102_1", "322_0", "303_0", "134_1", "236_2",
                            "318_0", "122_1", "210_2", "112_1", "137_1", "110_1", "111_1", "128_1",
                            "214_2", "120_1", "118_1", "232_2", "328_0", "226_2", "229_2", "324_0",
                            "235_2", "103_1", "309_0", "204_2", "113_1", "306_0", "207_2", "304_0",
                            "324_0", "235_2", "103_1", "309_0", "204_2", "113_1", "306_0", "207_2",
                            "304_0", "239_2", "225_2", "115_1", "202_2", "218_2", "239_2", "225_2",
                            "115_1", "202_2", "218_2", "126_1", "116_1", "136_1", "117_1", "139_1",
                            "124_1", "302_0", "325_0", "133_1", "211_2", "104_1", "237_2")
        
        donors_li <- c("C3_0", "A2_1", "A1_1", "F4_2","F5_2","C4_0", "F3_2", "C5_0", "C2_0","C1_0", "A5_1",
                        "A6_1", "A4_1", "C6_0", "F1_2", "A3_1", "F2_2")

        donors_limone <- c("ALS_MotorCortex.012318_FC17.RDS_1", "Control_MotorCortex.121417_FC11.RDS_0",
                        "Control_MotorCortex.012218_FC13.RDS_0", "ALS_MotorCortex.012218_FC21.RDS_1", 
                        "ALS_MotorCortex.121417_FC19.RDS_1", "Control_MotorCortex.012218_FC12.RDS_0")
        
        
        donors_pineda <- sub("_[^_]+$", "", donors_pineda)
        donors_li <- sub("_[^_]+$", "", donors_li)
        donors_limone <- sub("_[^_]+$", "", donors_limone)
    ##
                                                          
    ###########################
    ## HVG: Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_diagnosis_5000_cell_true_label.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(sample_mean_accuracy),        
        #        .groups = "drop"               
        #    )

        BA4_bind_all$input <- "Diagnosis"
        summary_stats_diagnosis <- BA4_bind_all
    ##

    ###########################
    ## HVG: Genetic group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_genetic_5000_cell_true_label.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        ## Calculate sample wise mean
        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(sample_mean_accuracy),        
        #        .groups = "drop"               
        #    )
        
        BA4_bind_all$input <- "Genetic"
        summary_stats_genetic <- BA4_bind_all
    ##

    ###########################
    ## HVG: Disease group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_disease_5000_cell_true_label.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        ## Calculate sample wise mean
        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(sample_mean_accuracy),        
        #        .groups = "drop"               
        #    )

        BA4_bind_all$input <- "Disease"
        summary_stats_disease <- BA4_bind_all
    ##

    ###########################
    ## HVG: Merge 
    ###########################
    ## code
        total_bind_HVG <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_disease)
        total_bind_HVG$input <- "HVG"
    ##

    ###########################
    ## NMF: Diagnosis
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        ## Calculate sample wise mean
        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Calculate sample wise mean
        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(sample_mean_accuracy),        
        #        .groups = "drop"               
        #    )
        
        BA4_bind_all$input <- "Diagnosis"
        summary_stats_diagnosis <- BA4_bind_all
    ##

    ###########################
    ## NMF: Genetic group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        ## Calculate sample wise mean
        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Calculate sample wise mean
        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(sample_mean_accuracy),        
        #        .groups = "drop"               
        #    )
        
        BA4_bind_all$input <- "Genetic"
        summary_stats_genetic <- BA4_bind_all
    ##

    ###########################
    ## NMF: Disease group
    ###########################
    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Add dataset info
        BA4_bind_all$donor <- sub("_[^_]+$", "", BA4_bind_all$donor)

        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_pineda] <- "Pineda et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_li] <- "Li et al."
        BA4_bind_all$dataset[BA4_bind_all$donor %in% donors_limone] <- "Limone et al."

        unique(BA4_bind_all$dataset)

        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 3, with_ties = FALSE) %>%
            ungroup()

        ## Only retain samples with at least 25 cells
        BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)
        
        ## Calculate sample wise mean
        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(donor, celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(test_accuracy_all),        
        #        .groups = "drop"               
        #    )

        ## Calculate sample wise mean
        ## Calculate sample wise mean
        #BA4_bind_all <- BA4_bind_all %>%
        #    group_by(celltype, dataset) %>%
        #    summarise(
        #        sample_mean_accuracy = mean(sample_mean_accuracy),        
        #        .groups = "drop"               
        #    )

        BA4_bind_all$input <- "Disease"
        summary_stats_disease <- BA4_bind_all
    ##

    ###########################
    ## NMF: Merge 
    ###########################
    ## code
        total_bind_NMF <- rbind(summary_stats_diagnosis, summary_stats_genetic, summary_stats_disease)
        total_bind_NMF$input <- "NMF"
    ##
   
    ###########################
    ## Merge and plot
    ###########################
    ## code
        total_bind <- rbind(total_bind_HVG, total_bind_NMF)
        
        total_bind$dataset <- factor(total_bind$dataset, levels = c("Pineda et al.", "Li et al.", "Limone et al."))

        ggplot(total_bind, aes(x = dataset, y = test_accuracy_all, fill = dataset)) + 
            theme_bw() + 
            geom_boxplot(outlier.shape = NA) +
            theme(
                legend.position = "none",
                panel.grid = element_blank(),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, colour = "black"),
                axis.text.y = element_text(size = 8, colour = "black"),
                axis.title.x = element_blank(),
                axis.title.y = element_text(face = "bold")
            ) +
            facet_grid(~input, scales = "free_x") +
            scale_fill_manual(values = c("Pineda et al." = "#FF7F00",
                                    "Li et al."  = "#FDBF6F",
                                    "Limone et al." = "#B15928")) +
            ylab("Accuracy") +
            scale_x_discrete(labels = c("Pineda et al.", "Li et al.", "Limone et al."))

        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3, width =2.25)

        ## Calculate dataset-specific overall accuracy
        overall_summary_stats <- total_bind %>%
            group_by(dataset, input) %>%
            summarise(
                median_accuracy = median(sample_mean_accuracy), 
                mean_accuracy = mean(sample_mean_accuracy),           # Calculate the median of X0
                sd_accuracy = sd(sample_mean_accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(sample_mean_accuracy - median(sample_mean_accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )
    ##
##
