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
## HVGs BA4 + BA9 together. 
###############################################################################################


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ ALS v FTLD v Control

## Code: true labels
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
        group_by(donor, celltype, run, region) %>%
        summarise(
            sample_mean_accuracy = mean(test_accuracy_all),        
            .groups = "drop"               
        )

    ## Only retain samples with at least 25 cells
    #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

    BA4_bind_all$input <- "True"
    BA4_bind_all_true <- BA4_bind_all
##

## Code: random labels
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
        slice_min(test_accuracy_all, n = 3, with_ties = FALSE) %>%
        ungroup()

    ## Calculate sample wise mean
    BA4_bind_all <- BA4_bind_all %>%
        group_by(donor, celltype, run, region) %>%
        summarise(
            sample_mean_accuracy = mean(test_accuracy_all),        
            .groups = "drop"               
        )

    BA4_bind_all$input <- "Random"
    BA4_bind_all_random <- BA4_bind_all

##
    
## Code: Wilcoxon test
    celltype_list <- c("L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6",
                   "PV", "5HT3aR", "Rosehip", "SOM",
                   "Oligo", "Astro", "OPC", "Micro", "Mural", "Endo", "Fibro", "L5")

    results_list <- list()

    for(par_celltype in celltype_list){
        BA4_bind_all_random_lim <- subset(BA4_bind_all_random, celltype == par_celltype)
        BA4_bind_all_true_lim <- subset(BA4_bind_all_true, celltype == par_celltype)

        BA4_bind_all_random_lim <- subset(BA4_bind_all_random_lim, donor %in% BA4_bind_all_true_lim$donor)

        acc_true <- BA4_bind_all_true_lim$sample_mean_accuracy
        acc_rand <- BA4_bind_all_random_lim$sample_mean_accuracy

        # Skip if there are not enough paired samples
        if(length(acc_true) < 2 || length(acc_rand) < 2){
            results_list[[par_celltype]] <- data.frame(
            celltype = par_celltype,
            pval = NA,
            effect_size = NA,
            mean_true_lab = mean(acc_true, na.rm = TRUE),
            mean_rand_lab = mean(acc_rand, na.rm = TRUE)
            )
            next
        }

        # Wilcoxon signed-rank test
        wilcox_res <- wilcox.test(acc_true, acc_rand, paired = TRUE, exact = FALSE)
        pval <- wilcox_res$p.value

        # Effect size: Cohen's d on differences
        differences <- acc_true - acc_rand
        cohens_d <- mean(differences) / sd(differences)

        results_list[[par_celltype]] <- data.frame(
            celltype = par_celltype,
            pval = pval,
            effect_size = cohens_d,  # or rb for rank-biserial
            mean_true_lab = mean(acc_true),
            mean_rand_lab = mean(acc_rand)
        )
    }

    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL
    #results_df$diff <- results_df$mean_true_lab - results_df$mean_rand_lab

    results_df$pval_adj <- p.adjust(results_df$pval, method = "BH")
    results_df$pval_bonf <- p.adjust(results_df$pval, method = "bonferroni")

##

## Code: permutation testing
    unique(BA4_bind_all_random$donor) == unique(BA4_bind_all_true$donor)
    
    #par_celltype = "L2_L3"

    celltype_list <- c("L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6",
    "PV", "5HT3aR", "Rosehip", "SOM",
    "Oligo", "Astro", "OPC", "Micro", "Mural", "Endo", "Fibro", "L5")

    results_list <- list()
    
    for(par_celltype in celltype_list){
        BA4_bind_all_random_lim <- subset(BA4_bind_all_random, celltype == par_celltype)
        BA4_bind_all_true_lim <- subset(BA4_bind_all_true, celltype == par_celltype)

        BA4_bind_all_random_lim <- subset(BA4_bind_all_random_lim, donor %in% BA4_bind_all_true_lim$donor )
        unique(BA4_bind_all_random_lim$donor) == unique(BA4_bind_all_true_lim$donor)
    
        nrow(BA4_bind_all_random_lim) == nrow(BA4_bind_all_true_lim)
        
        acc_true <- BA4_bind_all_true_lim$sample_mean_accuracy
        acc_rand <- BA4_bind_all_random_lim$sample_mean_accuracy

        mean_true <- mean(acc_true)
        mean_rand <- mean(acc_rand)

        ## Compute observed difference
        differences <- acc_true - acc_rand
        obs_diff <- mean(differences)
        obs_diff

        ## Permutation tet
        set.seed(123)          # reproducibility
        n_perm <- 10000        # number of permutations
        perm_diffs <- numeric(n_perm)

        for (i in 1:n_perm) {
        flips <- sample(c(1, -1), size = length(differences), replace = TRUE)
        perm_diffs[i] <- mean(differences * flips)
        }

    
        pval <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)
        pval

        ## Effect size
        cohens_d <- mean(differences) / sd(differences)
        cohens_d

        ## Store results
        results_list[[par_celltype]] <- data.frame(
            celltype = par_celltype,
            pval = pval,
            cohens_d = cohens_d,
            mean_true_lab = mean_true,
            mean_rand_lab = mean_rand
        )
    }

    # Combine into one dataframe
    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL

    results_df

    results_df$diff <- results_df$mean_true_lab - results_df$mean_rand_lab


    ## Viualization
    hist(perm_diffs, breaks = 50, main = "Null distribution (paired permutation test)",
    xlab = "Mean difference under null")
    abline(v = obs_diff, col = "red", lwd = 2)
##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Disease group v control
## Code: true labels
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
        slice_max(test_accuracy_all, n = 2, with_ties = FALSE) %>%
        ungroup()

    ## Calculate sample wise mean
    BA4_bind_all <- BA4_bind_all %>%
        group_by(donor, celltype, run, region) %>%
        summarise(
            sample_mean_accuracy = mean(test_accuracy_all),        
            .groups = "drop"               
        )

    ## Only retain samples with at least 25 cells
    #BA4_bind_all <- subset(BA4_bind_all, n_cells > 10)

    BA4_bind_all$input <- "True"
    BA4_bind_all_true <- BA4_bind_all
##

## Code: random labels
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
        slice_min(test_accuracy_all, n = 2, with_ties = FALSE) %>%
        ungroup()

    ## Calculate sample wise mean
    BA4_bind_all <- BA4_bind_all %>%
        group_by(donor, celltype, run, region) %>%
        summarise(
            sample_mean_accuracy = mean(test_accuracy_all),        
            .groups = "drop"               
        )

    BA4_bind_all$input <- "Random"
    BA4_bind_all_random <- BA4_bind_all

##
    
## Code: Wilcoxon test
    celltype_list <- c("L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6",
                   "PV", "5HT3aR", "Rosehip", "SOM",
                   "Oligo", "Astro", "OPC", "Micro", "Mural", "Endo", "Fibro", "L5")

    results_list <- list()

    for(par_celltype in celltype_list){
        BA4_bind_all_random_lim <- subset(BA4_bind_all_random, celltype == par_celltype)
        BA4_bind_all_true_lim <- subset(BA4_bind_all_true, celltype == par_celltype)

        BA4_bind_all_random_lim <- subset(BA4_bind_all_random_lim, donor %in% BA4_bind_all_true_lim$donor)

        acc_true <- BA4_bind_all_true_lim$sample_mean_accuracy
        acc_rand <- BA4_bind_all_random_lim$sample_mean_accuracy

        # Skip if there are not enough paired samples
        if(length(acc_true) < 2 || length(acc_rand) < 2){
            results_list[[par_celltype]] <- data.frame(
            celltype = par_celltype,
            pval = NA,
            effect_size = NA,
            mean_true_lab = mean(acc_true, na.rm = TRUE),
            mean_rand_lab = mean(acc_rand, na.rm = TRUE)
            )
            next
        }

        # Wilcoxon signed-rank test
        wilcox_res <- wilcox.test(acc_true, acc_rand, paired = TRUE, exact = FALSE)
        pval <- wilcox_res$p.value

        # Effect size: Cohen's d on differences
        differences <- acc_true - acc_rand
        cohens_d <- mean(differences) / sd(differences)

        results_list[[par_celltype]] <- data.frame(
            celltype = par_celltype,
            pval = pval,
            effect_size = cohens_d,  # or rb for rank-biserial
            mean_true_lab = mean(acc_true),
            mean_rand_lab = mean(acc_rand)
        )
    }

    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL
    #results_df$diff <- results_df$mean_true_lab - results_df$mean_rand_lab

    results_df$pval_adj <- p.adjust(results_df$pval, method = "BH")
    results_df$pval_bonf <- p.adjust(results_df$pval, method = "bonferroni")

##

## Code: permutation testing
    unique(BA4_bind_all_random$donor) == unique(BA4_bind_all_true$donor)
    
    #par_celltype = "L2_L3"

    celltype_list <- c("L3_L5", "L2_L3", "L4_L6", "L4_L5", "L5_L6", "L6",
    "PV", "5HT3aR", "Rosehip", "SOM",
    "Oligo", "Astro", "OPC", "Micro", "Mural", "Endo", "Fibro", "L5")

    results_list <- list()
    
    for(par_celltype in celltype_list){
        BA4_bind_all_random_lim <- subset(BA4_bind_all_random, celltype == par_celltype)
        BA4_bind_all_true_lim <- subset(BA4_bind_all_true, celltype == par_celltype)

        BA4_bind_all_random_lim <- subset(BA4_bind_all_random_lim, donor %in% BA4_bind_all_true_lim$donor )
        unique(BA4_bind_all_random_lim$donor) == unique(BA4_bind_all_true_lim$donor)
    
        nrow(BA4_bind_all_random_lim) == nrow(BA4_bind_all_true_lim)
        
        acc_true <- BA4_bind_all_true_lim$sample_mean_accuracy
        acc_rand <- BA4_bind_all_random_lim$sample_mean_accuracy

        mean_true <- mean(acc_true)
        mean_rand <- mean(acc_rand)

        ## Compute observed difference
        differences <- acc_true - acc_rand
        obs_diff <- mean(differences)
        obs_diff

        ## Permutation tet
        set.seed(123)          # reproducibility
        n_perm <- 10000        # number of permutations
        perm_diffs <- numeric(n_perm)

        for (i in 1:n_perm) {
        flips <- sample(c(1, -1), size = length(differences), replace = TRUE)
        perm_diffs[i] <- mean(differences * flips)
        }

    
        pval <- (sum(abs(perm_diffs) >= abs(obs_diff)) + 1) / (n_perm + 1)
        pval

        ## Effect size
        cohens_d <- mean(differences) / sd(differences)
        cohens_d

        ## Store results
        results_list[[par_celltype]] <- data.frame(
            celltype = par_celltype,
            pval = pval,
            cohens_d = cohens_d,
            mean_true_lab = mean_true,
            mean_rand_lab = mean_rand
        )
    }

    # Combine into one dataframe
    results_df <- do.call(rbind, results_list)
    rownames(results_df) <- NULL

    results_df

    results_df$diff <- results_df$mean_true_lab - results_df$mean_rand_lab


    ## Viualization
    hist(perm_diffs, breaks = 50, main = "Null distribution (paired permutation test)",
    xlab = "Mean difference under null")
    abline(v = obs_diff, col = "red", lwd = 2)
##

