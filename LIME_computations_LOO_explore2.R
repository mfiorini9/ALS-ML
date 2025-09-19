## Options
1. Do we do absolute value of LIME?
2. Do we weigh by percent expression or not?
3. Do we take case and control or only case?
4. Make sure we are only using the correctly classified cells -- I think we did this at the LIME level in the original script
5. Do we want to remove subject with insufficient number of cells?
6. Shold compare list to original list



screen -S percent_expressed_LOO
salloc -A def-tdurcan --time=0-8 -c 1 --mem=10g

module load StdEnv/2023
module load r/4.4.0
R


# Load required libraries
#install.packages("lubridate", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("tidyverse", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("cluster", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("factoextra", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("reshape2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("umap", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("purrr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("ashr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### SFTLD BA4

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_SFTLD_BA4_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_SFTLD_BA4_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_SFTLD_BA4_2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_SFTLD_BA4_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA4"
    par_status = "SFTLD"
    par_remove_group = c("SALS", "C9ALS", "C9FTLD")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, .9995), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, .9995)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_.9995_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control; no delta threshold.  -- USE THIS ONE
    
    ## code 
    
        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        nano LIME_LOO_explore_optimal_test_SFTLD_BA4_2_1.sh

        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #!/bin/bash  
        #SBATCH --account=def-sfarhan
        #SBATCH --time=01-00:00           # time (DD-HH:MM)
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=40g          # memory per cor
        #SBATCH --job-name=LIME_LOO_explore_optimal_test_SFTLD_BA4_2_1
        #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
        #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

        module load StdEnv/2023
        module load r/4.4.0

        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_SFTLD_BA4_2_1.R


        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        nano LIME_LOO_explore_optimal_test_SFTLD_BA4_2_1.R


        library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
        library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
        library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
        library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

        par_brain_region = "BA4"
        par_status = "SFTLD"
        par_remove_group = c("SALS", "C9ALS", "C9FTLD")
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        

        ## code 
            for (celltype2 in celltype_list ){    
            # Load your LIME data (replace with your actual CSV file)
            lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
            lime_df$test <- str_count(lime_df$feature, ' ')

            lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
            lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
            
            express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

            ## Normalize the percent expressed value
            express_fill_df <- express_fill_df %>%
            group_by(subject) %>%
            mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
            ungroup()

            express_fill_df <- data.frame(express_fill_df)

            ## We need to take the aboslute value of LIME
            lime_df$importance <- abs(lime_df$importance)
            nrow(lime_df)

            ## merge lime with percent expressed
            lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

            ## compute importance weighted by expression
            lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

            ## Compute mean importance across cells from each subject
            donor_gene_matrix <- lime_df2 %>%
            group_by(test_donor, feature) %>%
            summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
            pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

            dim(donor_gene_matrix)

            # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
            lime_matrix <- as.data.frame(donor_gene_matrix)
            n_subjects <- nrow(lime_matrix)
            gene_names <- colnames(lime_matrix)

            # Store significant genes in each LOSO fold
            loso_sig_genes <- list()
            loso_mean_delta_lime <- list()

            for (i in 1:n_subjects) {
            # Exclude subject i
            test_data <- lime_matrix[-i, ]  # Remaining subjects
            n_remaining <- nrow(test_data)

            # Global gene mean for each subject (i.e., across genes)
            test_data <- test_data[, sapply(test_data, is.numeric)]
            global_gene_means <- rowMeans(test_data)

            # Build a matrix of (LIME value - global mean) for each gene
            delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

            # Now test if the deltas are significantly > 0
                pvals <- sapply(delta_matrix, function(delta_vals) {
                wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
                })

            # Adjust for multiple testing
            adj_pvals <- p.adjust(pvals, method = "BH")

            # Store significant genes
            sig_genes <- names(adj_pvals[adj_pvals < 0.05])
            loso_sig_genes[[i]] <- sig_genes

            loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
            }

            # Intersect across LOSO folds
            robust_genes <- reduce(loso_sig_genes, intersect)
            length(robust_genes)

            ##############################################################################

            ## Calculate mean Delta across LOO
            delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

            # Now calculate mean for each gene (i.e., column)
            mean_deltas <- colMeans(delta_df)

            ## Non sig deltas
            nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
            nonsig_deltas_df <- data.frame(nonsig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
            hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
            dev.off()

            threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


            ## Sig deltas
            sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
            sig_deltas_df <- data.frame(sig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
            hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = threshold, col = "red", lty = 2)
            dev.off()

            sig_deltas_df$feature <- rownames(sig_deltas_df)
            refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
            nrow(refined_genes)


            ############################################################################## 
            mean_deltas_df <- data.frame(
            gene = names(mean_deltas),
            mean_delta = as.numeric(mean_deltas)
            )

            mean_deltas_sorted <- mean_deltas_df %>%
            arrange(desc(mean_delta))

            mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
            mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
            print(nrow(mean_deltas_sorted))

            write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")

            }
    ##
    
    
    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Not taking absolute value; weighted by percent expression; case and control
    
    
    ## code 
        for (celltype2 in celltype_list ){    
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
        
        express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        #lime_df$importance <- abs(lime_df$importance)
        #nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
        print(nrow(mean_deltas_sorted))

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_no_abs_case_control_narval_2.csv"), sep = ",")

        }
    ##

##


####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### SFTLD BA9

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_SFTLD_BA9_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_SFTLD_BA9_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_SFTLD_BA9_2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_SFTLD_BA9_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA9"
    par_status = "SFTLD"
    par_remove_group = c("SALS", "C9ALS", "C9FTLD")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control 
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, .9995), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, .9995)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_.9995_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control; no delta threshold. -- USE THIS ONE
    
    ## code 
    
        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        nano LIME_LOO_explore_optimal_test_SFTLD_BA9_2_1.sh

        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #!/bin/bash  
        #SBATCH --account=def-sfarhan
        #SBATCH --time=01-00:00           # time (DD-HH:MM)
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=40g          # memory per cor
        #SBATCH --job-name=LIME_LOO_explore_optimal_test_SFTLD_BA9_2_1
        #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
        #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

        module load StdEnv/2023
        module load r/4.4.0

        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_SFTLD_BA9_2_1.R


        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        nano LIME_LOO_explore_optimal_test_SFTLD_BA9_2_1.R


        library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
        library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
        library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
        library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

        par_brain_region = "BA9"
        par_status = "SFTLD"
        par_remove_group = c("SALS", "C9ALS", "C9FTLD")
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        

        ## code 
            for (celltype2 in celltype_list ){    
            # Load your LIME data (replace with your actual CSV file)
            lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
            lime_df$test <- str_count(lime_df$feature, ' ')

            lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
            lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
            
            express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

            ## Normalize the percent expressed value
            express_fill_df <- express_fill_df %>%
            group_by(subject) %>%
            mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
            ungroup()

            express_fill_df <- data.frame(express_fill_df)

            ## We need to take the aboslute value of LIME
            lime_df$importance <- abs(lime_df$importance)
            nrow(lime_df)

            ## merge lime with percent expressed
            lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

            ## compute importance weighted by expression
            lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

            ## Compute mean importance across cells from each subject
            donor_gene_matrix <- lime_df2 %>%
            group_by(test_donor, feature) %>%
            summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
            pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

            dim(donor_gene_matrix)

            # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
            lime_matrix <- as.data.frame(donor_gene_matrix)
            n_subjects <- nrow(lime_matrix)
            gene_names <- colnames(lime_matrix)

            # Store significant genes in each LOSO fold
            loso_sig_genes <- list()
            loso_mean_delta_lime <- list()

            for (i in 1:n_subjects) {
            # Exclude subject i
            test_data <- lime_matrix[-i, ]  # Remaining subjects
            n_remaining <- nrow(test_data)

            # Global gene mean for each subject (i.e., across genes)
            test_data <- test_data[, sapply(test_data, is.numeric)]
            global_gene_means <- rowMeans(test_data)

            # Build a matrix of (LIME value - global mean) for each gene
            delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

            # Now test if the deltas are significantly > 0
                pvals <- sapply(delta_matrix, function(delta_vals) {
                wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
                })

            # Adjust for multiple testing
            adj_pvals <- p.adjust(pvals, method = "BH")

            # Store significant genes
            sig_genes <- names(adj_pvals[adj_pvals < 0.05])
            loso_sig_genes[[i]] <- sig_genes

            loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
            }

            # Intersect across LOSO folds
            robust_genes <- reduce(loso_sig_genes, intersect)
            length(robust_genes)

            ##############################################################################

            ## Calculate mean Delta across LOO
            delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

            # Now calculate mean for each gene (i.e., column)
            mean_deltas <- colMeans(delta_df)

            ## Non sig deltas
            nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
            nonsig_deltas_df <- data.frame(nonsig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
            hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
            dev.off()

            threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


            ## Sig deltas
            sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
            sig_deltas_df <- data.frame(sig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
            hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = threshold, col = "red", lty = 2)
            dev.off()

            sig_deltas_df$feature <- rownames(sig_deltas_df)
            refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
            nrow(refined_genes)


            ############################################################################## 
            mean_deltas_df <- data.frame(
            gene = names(mean_deltas),
            mean_delta = as.numeric(mean_deltas)
            )

            mean_deltas_sorted <- mean_deltas_df %>%
            arrange(desc(mean_delta))

            mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
            mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
            print(nrow(mean_deltas_sorted))

            write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")

            }
    ##
    
    
    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Not taking absolute value; weighted by percent expression; case and control
    
    
    ## code 
        for (celltype2 in celltype_list ){    
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
        
        express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        #lime_df$importance <- abs(lime_df$importance)
        #nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
        print(nrow(mean_deltas_sorted))

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_no_abs_case_control_narval_2.csv"), sep = ",")

        }
    ##

##


####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### C9FTLD BA4

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_C9FTLD_BA4_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_C9FTLD_BA4_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_C9FTLD_BA4_2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_C9FTLD_BA4_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA4"
    par_status = "C9FTLD"
    par_remove_group = c("SALS", "C9ALS", "SFTLD")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, .9995), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, .9995)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_.9995_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control; no delta threshold.  -- USE THIS ONE
    
    ## code 
    
        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        nano LIME_LOO_explore_optimal_test_C9FTLD_BA4_2_1.sh

        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #!/bin/bash  
        #SBATCH --account=def-sfarhan
        #SBATCH --time=01-00:00           # time (DD-HH:MM)
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=40g          # memory per cor
        #SBATCH --job-name=LIME_LOO_explore_optimal_test_C9FTLD_BA4_2_1
        #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
        #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

        module load StdEnv/2023
        module load r/4.4.0

        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_C9FTLD_BA4_2_1.R


        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        nano LIME_LOO_explore_optimal_test_C9FTLD_BA4_2_1.R


        library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
        library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
        library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
        library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

        par_brain_region = "BA4"
        par_status = "C9FTLD"
        par_remove_group = c("SALS", "C9ALS", "SFTLD")
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        

        ## code 
            for (celltype2 in celltype_list ){    
            # Load your LIME data (replace with your actual CSV file)
            lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
            lime_df$test <- str_count(lime_df$feature, ' ')

            lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
            lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
            
            express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

            ## Normalize the percent expressed value
            express_fill_df <- express_fill_df %>%
            group_by(subject) %>%
            mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
            ungroup()

            express_fill_df <- data.frame(express_fill_df)

            ## We need to take the aboslute value of LIME
            lime_df$importance <- abs(lime_df$importance)
            nrow(lime_df)

            ## merge lime with percent expressed
            lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

            ## compute importance weighted by expression
            lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

            ## Compute mean importance across cells from each subject
            donor_gene_matrix <- lime_df2 %>%
            group_by(test_donor, feature) %>%
            summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
            pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

            dim(donor_gene_matrix)

            # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
            lime_matrix <- as.data.frame(donor_gene_matrix)
            n_subjects <- nrow(lime_matrix)
            gene_names <- colnames(lime_matrix)

            # Store significant genes in each LOSO fold
            loso_sig_genes <- list()
            loso_mean_delta_lime <- list()

            for (i in 1:n_subjects) {
            # Exclude subject i
            test_data <- lime_matrix[-i, ]  # Remaining subjects
            n_remaining <- nrow(test_data)

            # Global gene mean for each subject (i.e., across genes)
            test_data <- test_data[, sapply(test_data, is.numeric)]
            global_gene_means <- rowMeans(test_data)

            # Build a matrix of (LIME value - global mean) for each gene
            delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

            # Now test if the deltas are significantly > 0
                pvals <- sapply(delta_matrix, function(delta_vals) {
                wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
                })

            # Adjust for multiple testing
            adj_pvals <- p.adjust(pvals, method = "BH")

            # Store significant genes
            sig_genes <- names(adj_pvals[adj_pvals < 0.05])
            loso_sig_genes[[i]] <- sig_genes

            loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
            }

            # Intersect across LOSO folds
            robust_genes <- reduce(loso_sig_genes, intersect)
            length(robust_genes)

            ##############################################################################

            ## Calculate mean Delta across LOO
            delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

            # Now calculate mean for each gene (i.e., column)
            mean_deltas <- colMeans(delta_df)

            ## Non sig deltas
            nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
            nonsig_deltas_df <- data.frame(nonsig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
            hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
            dev.off()

            threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


            ## Sig deltas
            sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
            sig_deltas_df <- data.frame(sig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
            hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = threshold, col = "red", lty = 2)
            dev.off()

            sig_deltas_df$feature <- rownames(sig_deltas_df)
            refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
            nrow(refined_genes)


            ############################################################################## 
            mean_deltas_df <- data.frame(
            gene = names(mean_deltas),
            mean_delta = as.numeric(mean_deltas)
            )

            mean_deltas_sorted <- mean_deltas_df %>%
            arrange(desc(mean_delta))

            mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
            mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
            print(nrow(mean_deltas_sorted))

            write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")

            }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Not taking absolute value; weighted by percent expression; case and control
    
    ## code 
        for (celltype2 in celltype_list ){    
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
        
        express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        #lime_df$importance <- abs(lime_df$importance)
        #nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
        print(nrow(mean_deltas_sorted))

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_no_abs_case_control_narval_2.csv"), sep = ",")

        }
    ##

##


####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### C9FTLD BA9

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_C9FTLD_BA9_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_C9FTLD_BA9_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_C9FTLD_BA9_2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_C9FTLD_BA9_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA9"
    par_status = "C9FTLD"
    par_remove_group = c("SALS", "C9ALS", "SFTLD")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control 
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, .9995), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, .9995)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_.9995_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control; no delta threshold. -- USE THIS ONE
    
    ## code 
    
        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        nano LIME_LOO_explore_optimal_test_C9FTLD_BA9_2_1.sh

        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #!/bin/bash  
        #SBATCH --account=def-sfarhan
        #SBATCH --time=01-00:00           # time (DD-HH:MM)
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=40g          # memory per cor
        #SBATCH --job-name=LIME_LOO_explore_optimal_test_C9FTLD_BA9_2_1
        #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
        #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

        module load StdEnv/2023
        module load r/4.4.0

        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_C9FTLD_BA9_2_1.R


        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        nano LIME_LOO_explore_optimal_test_C9FTLD_BA9_2_1.R


        library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
        library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
        library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
        library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

        par_brain_region = "BA9"
        par_status = "C9FTLD"
        par_remove_group = c("SALS", "C9ALS", "SFTLD")
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        

        ## code 
            for (celltype2 in celltype_list ){    
            # Load your LIME data (replace with your actual CSV file)
            lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
            lime_df$test <- str_count(lime_df$feature, ' ')

            lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
            lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
            
            express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

            ## Normalize the percent expressed value
            express_fill_df <- express_fill_df %>%
            group_by(subject) %>%
            mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
            ungroup()

            express_fill_df <- data.frame(express_fill_df)

            ## We need to take the aboslute value of LIME
            lime_df$importance <- abs(lime_df$importance)
            nrow(lime_df)

            ## merge lime with percent expressed
            lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

            ## compute importance weighted by expression
            lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

            ## Compute mean importance across cells from each subject
            donor_gene_matrix <- lime_df2 %>%
            group_by(test_donor, feature) %>%
            summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
            pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

            dim(donor_gene_matrix)

            # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
            lime_matrix <- as.data.frame(donor_gene_matrix)
            n_subjects <- nrow(lime_matrix)
            gene_names <- colnames(lime_matrix)

            # Store significant genes in each LOSO fold
            loso_sig_genes <- list()
            loso_mean_delta_lime <- list()

            for (i in 1:n_subjects) {
            # Exclude subject i
            test_data <- lime_matrix[-i, ]  # Remaining subjects
            n_remaining <- nrow(test_data)

            # Global gene mean for each subject (i.e., across genes)
            test_data <- test_data[, sapply(test_data, is.numeric)]
            global_gene_means <- rowMeans(test_data)

            # Build a matrix of (LIME value - global mean) for each gene
            delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

            # Now test if the deltas are significantly > 0
                pvals <- sapply(delta_matrix, function(delta_vals) {
                wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
                })

            # Adjust for multiple testing
            adj_pvals <- p.adjust(pvals, method = "BH")

            # Store significant genes
            sig_genes <- names(adj_pvals[adj_pvals < 0.05])
            loso_sig_genes[[i]] <- sig_genes

            loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
            }

            # Intersect across LOSO folds
            robust_genes <- reduce(loso_sig_genes, intersect)
            length(robust_genes)

            ##############################################################################

            ## Calculate mean Delta across LOO
            delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

            # Now calculate mean for each gene (i.e., column)
            mean_deltas <- colMeans(delta_df)

            ## Non sig deltas
            nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
            nonsig_deltas_df <- data.frame(nonsig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
            hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
            dev.off()

            threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


            ## Sig deltas
            sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
            sig_deltas_df <- data.frame(sig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
            hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = threshold, col = "red", lty = 2)
            dev.off()

            sig_deltas_df$feature <- rownames(sig_deltas_df)
            refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
            nrow(refined_genes)


            ############################################################################## 
            mean_deltas_df <- data.frame(
            gene = names(mean_deltas),
            mean_delta = as.numeric(mean_deltas)
            )

            mean_deltas_sorted <- mean_deltas_df %>%
            arrange(desc(mean_delta))

            mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
            mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
            print(nrow(mean_deltas_sorted))

            write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")

            }
    ##


    ####################################################
    ####################################################
    ####################################################
    #################################################### Not taking absolute value; weighted by percent expression; case and control
    
    
    ## code 
        for (celltype2 in celltype_list ){    
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
        
        express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        #lime_df$importance <- abs(lime_df$importance)
        #nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
        print(nrow(mean_deltas_sorted))

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_no_abs_case_control_narval_2.csv"), sep = ",")

        }
    ##

##


####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### C9ALS BA4

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_C9ALS_BA4_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_C9ALS_BA4_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_C9ALS_BA4_2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_C9ALS_BA4_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA4"
    par_status = "C9ALS"
    par_remove_group = c("SALS", "C9FTLD", "SFTLD")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control 
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, .9995), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, .9995)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_.9995_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control; no delta threshold. -- USE THIS ONE
    
    ## code 
    
        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        nano LIME_LOO_explore_optimal_test_C9ALS_BA4_2_1.sh

        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #!/bin/bash  
        #SBATCH --account=def-sfarhan
        #SBATCH --time=01-00:00           # time (DD-HH:MM)
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=40g          # memory per cor
        #SBATCH --job-name=LIME_LOO_explore_optimal_test_C9ALS_BA4_2_1
        #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
        #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

        module load StdEnv/2023
        module load r/4.4.0

        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_C9ALS_BA4_2_1.R


        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        nano LIME_LOO_explore_optimal_test_C9ALS_BA4_2_1.R


        library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
        library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
        library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
        library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

        par_brain_region = "BA4"
        par_status = "C9ALS"
        par_remove_group = c("SALS", "C9FTLD", "SFTLD")
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        

        ## code 
            for (celltype2 in celltype_list ){    
            # Load your LIME data (replace with your actual CSV file)
            lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
            lime_df$test <- str_count(lime_df$feature, ' ')

            lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
            lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
            
            express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

            ## Normalize the percent expressed value
            express_fill_df <- express_fill_df %>%
            group_by(subject) %>%
            mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
            ungroup()

            express_fill_df <- data.frame(express_fill_df)

            ## We need to take the aboslute value of LIME
            lime_df$importance <- abs(lime_df$importance)
            nrow(lime_df)

            ## merge lime with percent expressed
            lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

            ## compute importance weighted by expression
            lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

            ## Compute mean importance across cells from each subject
            donor_gene_matrix <- lime_df2 %>%
            group_by(test_donor, feature) %>%
            summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
            pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

            dim(donor_gene_matrix)

            # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
            lime_matrix <- as.data.frame(donor_gene_matrix)
            n_subjects <- nrow(lime_matrix)
            gene_names <- colnames(lime_matrix)

            # Store significant genes in each LOSO fold
            loso_sig_genes <- list()
            loso_mean_delta_lime <- list()

            for (i in 1:n_subjects) {
            # Exclude subject i
            test_data <- lime_matrix[-i, ]  # Remaining subjects
            n_remaining <- nrow(test_data)

            # Global gene mean for each subject (i.e., across genes)
            test_data <- test_data[, sapply(test_data, is.numeric)]
            global_gene_means <- rowMeans(test_data)

            # Build a matrix of (LIME value - global mean) for each gene
            delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

            # Now test if the deltas are significantly > 0
                pvals <- sapply(delta_matrix, function(delta_vals) {
                wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
                })

            # Adjust for multiple testing
            adj_pvals <- p.adjust(pvals, method = "BH")

            # Store significant genes
            sig_genes <- names(adj_pvals[adj_pvals < 0.05])
            loso_sig_genes[[i]] <- sig_genes

            loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
            }

            # Intersect across LOSO folds
            robust_genes <- reduce(loso_sig_genes, intersect)
            length(robust_genes)

            ##############################################################################

            ## Calculate mean Delta across LOO
            delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

            # Now calculate mean for each gene (i.e., column)
            mean_deltas <- colMeans(delta_df)

            ## Non sig deltas
            nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
            nonsig_deltas_df <- data.frame(nonsig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
            hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
            dev.off()

            threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


            ## Sig deltas
            sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
            sig_deltas_df <- data.frame(sig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
            hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = threshold, col = "red", lty = 2)
            dev.off()

            sig_deltas_df$feature <- rownames(sig_deltas_df)
            refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
            nrow(refined_genes)


            ############################################################################## 
            mean_deltas_df <- data.frame(
            gene = names(mean_deltas),
            mean_delta = as.numeric(mean_deltas)
            )

            mean_deltas_sorted <- mean_deltas_df %>%
            arrange(desc(mean_delta))

            mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
            mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
            print(nrow(mean_deltas_sorted))

            write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")

            }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Not taking absolute value; weighted by percent expression; case and control
    
    
    ## code 
        for (celltype2 in celltype_list ){    
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
        
        express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        #lime_df$importance <- abs(lime_df$importance)
        #nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
        print(nrow(mean_deltas_sorted))

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_no_abs_case_control_narval_2.csv"), sep = ",")

        }
    ##

##


####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### C9ALS BA9

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_C9ALS_BA9_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_C9ALS_BA9_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_C9ALS_BA9_2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_C9ALS_BA9_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA9"
    par_status = "C9ALS"
    par_remove_group = c("SALS", "C9FTLD", "SFTLD")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control 
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, .9995), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, .9995)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_.9995_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control; no delta threshold. -- USE THIS ONE
    
    ## code 
    
        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        nano LIME_LOO_explore_optimal_test_C9ALS_BA9_2_1.sh

        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #!/bin/bash  
        #SBATCH --account=def-sfarhan
        #SBATCH --time=01-00:00           # time (DD-HH:MM)
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=40g          # memory per cor
        #SBATCH --job-name=LIME_LOO_explore_optimal_test_C9ALS_BA9_2_1
        #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
        #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

        module load StdEnv/2023
        module load r/4.4.0

        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_C9ALS_BA9_2_1.R


        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        nano LIME_LOO_explore_optimal_test_C9ALS_BA9_2_1.R


        library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
        library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
        library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
        library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

        par_brain_region = "BA9"
        par_status = "C9ALS"
        par_remove_group = c("SALS", "C9FTLD", "SFTLD")
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        

        ## code 
            for (celltype2 in celltype_list ){    
            # Load your LIME data (replace with your actual CSV file)
            lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
            lime_df$test <- str_count(lime_df$feature, ' ')

            lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
            lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
            
            express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

            ## Normalize the percent expressed value
            express_fill_df <- express_fill_df %>%
            group_by(subject) %>%
            mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
            ungroup()

            express_fill_df <- data.frame(express_fill_df)

            ## We need to take the aboslute value of LIME
            lime_df$importance <- abs(lime_df$importance)
            nrow(lime_df)

            ## merge lime with percent expressed
            lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

            ## compute importance weighted by expression
            lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

            ## Compute mean importance across cells from each subject
            donor_gene_matrix <- lime_df2 %>%
            group_by(test_donor, feature) %>%
            summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
            pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

            dim(donor_gene_matrix)

            # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
            lime_matrix <- as.data.frame(donor_gene_matrix)
            n_subjects <- nrow(lime_matrix)
            gene_names <- colnames(lime_matrix)

            # Store significant genes in each LOSO fold
            loso_sig_genes <- list()
            loso_mean_delta_lime <- list()

            for (i in 1:n_subjects) {
            # Exclude subject i
            test_data <- lime_matrix[-i, ]  # Remaining subjects
            n_remaining <- nrow(test_data)

            # Global gene mean for each subject (i.e., across genes)
            test_data <- test_data[, sapply(test_data, is.numeric)]
            global_gene_means <- rowMeans(test_data)

            # Build a matrix of (LIME value - global mean) for each gene
            delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

            # Now test if the deltas are significantly > 0
                pvals <- sapply(delta_matrix, function(delta_vals) {
                wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
                })

            # Adjust for multiple testing
            adj_pvals <- p.adjust(pvals, method = "BH")

            # Store significant genes
            sig_genes <- names(adj_pvals[adj_pvals < 0.05])
            loso_sig_genes[[i]] <- sig_genes

            loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
            }

            # Intersect across LOSO folds
            robust_genes <- reduce(loso_sig_genes, intersect)
            length(robust_genes)

            ##############################################################################

            ## Calculate mean Delta across LOO
            delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

            # Now calculate mean for each gene (i.e., column)
            mean_deltas <- colMeans(delta_df)

            ## Non sig deltas
            nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
            nonsig_deltas_df <- data.frame(nonsig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
            hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
            dev.off()

            threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


            ## Sig deltas
            sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
            sig_deltas_df <- data.frame(sig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
            hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = threshold, col = "red", lty = 2)
            dev.off()

            sig_deltas_df$feature <- rownames(sig_deltas_df)
            refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
            nrow(refined_genes)


            ############################################################################## 
            mean_deltas_df <- data.frame(
            gene = names(mean_deltas),
            mean_delta = as.numeric(mean_deltas)
            )

            mean_deltas_sorted <- mean_deltas_df %>%
            arrange(desc(mean_delta))

            mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
            mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
            print(nrow(mean_deltas_sorted))

            write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")

            }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Not taking absolute value; weighted by percent expression; case and control
    
    
    ## code 
        for (celltype2 in celltype_list ){    
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
        
        express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        #lime_df$importance <- abs(lime_df$importance)
        #nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
        print(nrow(mean_deltas_sorted))

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_no_abs_case_control_narval_2.csv"), sep = ",")

        }
    ##

##


####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### SALS BA4

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_SALS_BA4_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_SALS_BA4_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_SALS_BA4_2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_SALS_BA4_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA4"
    par_status = "SALS"
    par_remove_group = c("C9ALS", "C9FTLD", "SFTLD")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control 
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, .9995), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, .9995)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_.9995_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control; no delta threshold. -- USE THIS ONE
    
    ## code 
    
        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        nano LIME_LOO_explore_optimal_test_SALS_BA4_2_1.sh

        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #!/bin/bash  
        #SBATCH --account=def-sfarhan
        #SBATCH --time=01-00:00           # time (DD-HH:MM)
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=40g          # memory per cor
        #SBATCH --job-name=LIME_LOO_explore_optimal_test_SALS_BA4_2_1
        #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
        #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

        module load StdEnv/2023
        module load r/4.4.0

        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_SALS_BA4_2_1.R


        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        nano LIME_LOO_explore_optimal_test_SALS_BA4_2_1.R

        library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
        library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
        library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
        library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

        par_brain_region = "BA4"
        par_status = "SALS"
        par_remove_group = c("SALS", "C9FTLD", "SFTLD")
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        

        ## code 
            for (celltype2 in celltype_list ){    
            # Load your LIME data (replace with your actual CSV file)
            lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
            lime_df$test <- str_count(lime_df$feature, ' ')

            lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
            lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
            
            express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

            ## Normalize the percent expressed value
            express_fill_df <- express_fill_df %>%
            group_by(subject) %>%
            mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
            ungroup()

            express_fill_df <- data.frame(express_fill_df)

            ## We need to take the aboslute value of LIME
            lime_df$importance <- abs(lime_df$importance)
            nrow(lime_df)

            ## merge lime with percent expressed
            lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

            ## compute importance weighted by expression
            lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

            ## Compute mean importance across cells from each subject
            donor_gene_matrix <- lime_df2 %>%
            group_by(test_donor, feature) %>%
            summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
            pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

            dim(donor_gene_matrix)

            # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
            lime_matrix <- as.data.frame(donor_gene_matrix)
            n_subjects <- nrow(lime_matrix)
            gene_names <- colnames(lime_matrix)

            # Store significant genes in each LOSO fold
            loso_sig_genes <- list()
            loso_mean_delta_lime <- list()

            for (i in 1:n_subjects) {
            # Exclude subject i
            test_data <- lime_matrix[-i, ]  # Remaining subjects
            n_remaining <- nrow(test_data)

            # Global gene mean for each subject (i.e., across genes)
            test_data <- test_data[, sapply(test_data, is.numeric)]
            global_gene_means <- rowMeans(test_data)

            # Build a matrix of (LIME value - global mean) for each gene
            delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

            # Now test if the deltas are significantly > 0
                pvals <- sapply(delta_matrix, function(delta_vals) {
                wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
                })

            # Adjust for multiple testing
            adj_pvals <- p.adjust(pvals, method = "BH")

            # Store significant genes
            sig_genes <- names(adj_pvals[adj_pvals < 0.05])
            loso_sig_genes[[i]] <- sig_genes

            loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
            }

            # Intersect across LOSO folds
            robust_genes <- reduce(loso_sig_genes, intersect)
            length(robust_genes)

            ##############################################################################

            ## Calculate mean Delta across LOO
            delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

            # Now calculate mean for each gene (i.e., column)
            mean_deltas <- colMeans(delta_df)

            ## Non sig deltas
            nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
            nonsig_deltas_df <- data.frame(nonsig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
            hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
            dev.off()

            threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


            ## Sig deltas
            sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
            sig_deltas_df <- data.frame(sig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
            hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = threshold, col = "red", lty = 2)
            dev.off()

            sig_deltas_df$feature <- rownames(sig_deltas_df)
            refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
            nrow(refined_genes)


            ############################################################################## 
            mean_deltas_df <- data.frame(
            gene = names(mean_deltas),
            mean_delta = as.numeric(mean_deltas)
            )

            mean_deltas_sorted <- mean_deltas_df %>%
            arrange(desc(mean_delta))

            mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
            mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
            print(nrow(mean_deltas_sorted))

            write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")

            }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Not taking absolute value; weighted by percent expression; case and control
    
    
    ## code 
        for (celltype2 in celltype_list ){    
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
        
        express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        #lime_df$importance <- abs(lime_df$importance)
        #nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
        print(nrow(mean_deltas_sorted))

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_no_abs_case_control_narval_2.csv"), sep = ",")

        }
    ##

##



####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### SALS BA9

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_SALS_BA9_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_SALS_BA9_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_SALS_BA9_2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_SALS_BA9_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA9"
    par_status = "SALS"
    par_remove_group = c("C9ALS", "C9FTLD", "SFTLD")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control 
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, .9995), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, .9995)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_.9995_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control; no delta threshold. -- USE THIS ONE
    
    ## code 
    
        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        nano LIME_LOO_explore_optimal_test_SALS_BA9_2_1.sh

        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #!/bin/bash  
        #SBATCH --account=def-sfarhan
        #SBATCH --time=01-00:00           # time (DD-HH:MM)
        #SBATCH --cpus-per-task=1
        #SBATCH --mem=40g          # memory per cor
        #SBATCH --job-name=LIME_LOO_explore_optimal_test_SALS_BA9_2_1
        #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
        #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

        module load StdEnv/2023
        module load r/4.4.0

        cd /home/fiorini9/scratch/machine_learning_ALS/scripts

        Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_SALS_BA9_2_1.R


        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        nano LIME_LOO_explore_optimal_test_SALS_BA9_2_1.R

        library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
        library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
        library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
        library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

        par_brain_region = "BA9"
        par_status = "SALS"
        par_remove_group = c("SALS", "C9FTLD", "SFTLD")
        par_prep = "CombatSeq"
            
        celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
                "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
                "Micro",  "Mural",   "Endo",    "Fibro", "L5")
        

        ## code 
            for (celltype2 in celltype_list ){    
            # Load your LIME data (replace with your actual CSV file)
            lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
            lime_df$test <- str_count(lime_df$feature, ' ')

            lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
            lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
            
            express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

            ## Normalize the percent expressed value
            express_fill_df <- express_fill_df %>%
            group_by(subject) %>%
            mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
            ungroup()

            express_fill_df <- data.frame(express_fill_df)

            ## We need to take the aboslute value of LIME
            lime_df$importance <- abs(lime_df$importance)
            nrow(lime_df)

            ## merge lime with percent expressed
            lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

            ## compute importance weighted by expression
            lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

            ## Compute mean importance across cells from each subject
            donor_gene_matrix <- lime_df2 %>%
            group_by(test_donor, feature) %>%
            summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
            pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

            dim(donor_gene_matrix)

            # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
            lime_matrix <- as.data.frame(donor_gene_matrix)
            n_subjects <- nrow(lime_matrix)
            gene_names <- colnames(lime_matrix)

            # Store significant genes in each LOSO fold
            loso_sig_genes <- list()
            loso_mean_delta_lime <- list()

            for (i in 1:n_subjects) {
            # Exclude subject i
            test_data <- lime_matrix[-i, ]  # Remaining subjects
            n_remaining <- nrow(test_data)

            # Global gene mean for each subject (i.e., across genes)
            test_data <- test_data[, sapply(test_data, is.numeric)]
            global_gene_means <- rowMeans(test_data)

            # Build a matrix of (LIME value - global mean) for each gene
            delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

            # Now test if the deltas are significantly > 0
                pvals <- sapply(delta_matrix, function(delta_vals) {
                wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
                })

            # Adjust for multiple testing
            adj_pvals <- p.adjust(pvals, method = "BH")

            # Store significant genes
            sig_genes <- names(adj_pvals[adj_pvals < 0.05])
            loso_sig_genes[[i]] <- sig_genes

            loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
            }

            # Intersect across LOSO folds
            robust_genes <- reduce(loso_sig_genes, intersect)
            length(robust_genes)

            ##############################################################################

            ## Calculate mean Delta across LOO
            delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

            # Now calculate mean for each gene (i.e., column)
            mean_deltas <- colMeans(delta_df)

            ## Non sig deltas
            nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
            nonsig_deltas_df <- data.frame(nonsig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
            hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
            dev.off()

            threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


            ## Sig deltas
            sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
            sig_deltas_df <- data.frame(sig_deltas)
            pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
            hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
                xlab = "Mean Delta-LIME", col = "gray")
            abline(v = threshold, col = "red", lty = 2)
            dev.off()

            sig_deltas_df$feature <- rownames(sig_deltas_df)
            refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
            nrow(refined_genes)


            ############################################################################## 
            mean_deltas_df <- data.frame(
            gene = names(mean_deltas),
            mean_delta = as.numeric(mean_deltas)
            )

            mean_deltas_sorted <- mean_deltas_df %>%
            arrange(desc(mean_delta))

            mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
            mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
            print(nrow(mean_deltas_sorted))

            write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")

            }
    ##

    ####################################################
    ####################################################
    ####################################################
    #################################################### Not taking absolute value; weighted by percent expression; case and control
    
    
    ## code 
        for (celltype2 in celltype_list ){    
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed
        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
        
        express_fill_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"))  # <-- Update path if needed

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        #lime_df$importance <- abs(lime_df$importance)
        #nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)


        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )
        print(nrow(mean_deltas_sorted))

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_no_abs_case_control_narval_2.csv"), sep = ",")

        }
    ##

##

####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### All ALS BA4

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_all_ALS_BA4_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_all_ALS_BA4_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_all_ALS_BA4_2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_all_ALS_BA4_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA4"
    par_status = "All ALS"
    
    par_remove_group = c("C9FTLD", "SFTLD")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control -- USE THIS ONE
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed

        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)

        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

##

####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### All ALS BA9

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_all_ALS_BA9_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_all_ALS_BA9_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_all_ALS_BA9_2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_all_ALS_BA9_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA9"
    par_status = "All ALS"
    
    par_remove_group = c("C9FTLD", "SFTLD")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control -- USE THIS ONE
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed

        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)

        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

##

####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### All FTLD BA9

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_all_FTLD_BA9_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_all_FTLD_BA9_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_all_FTLD_BA9_2.R

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_all_FTLD_BA9_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA9"
    par_status = "All FTLD"
    
    par_remove_group = c("C9ALS", "SALS")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control -- USE THIS ONE
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed

        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)

        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

##


####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
####################################################
#################################################### All FTLD BA4

## code 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano LIME_LOO_explore_optimal_test_all_FTLD_BA4_2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=03-00:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=100g          # memory per cor
    #SBATCH --job-name=LIME_LOO_explore_optimal_test_all_FTLD_BA4_2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2023
    module load r/4.4.0

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/LIME_LOO_explore_optimal_test_all_FTLD_BA4_2.R

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano LIME_LOO_explore_optimal_test_all_FTLD_BA4_2.R


    library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
    library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
    library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
    library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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

    #celltype2 = "SOM"

    par_brain_region = "BA4"
    par_status = "All FTLD"
    
    par_remove_group = c("C9ALS", "SALS")
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")

    
    ####################################################
    ####################################################
    ####################################################
    #################################################### Taking absolute value; weighted by percent expression; case and control -- USE THIS ONE
    
    ## code 
        for (celltype2 in celltype_list ){
        # Load your LIME data (replace with your actual CSV file)
        lime_df <- read_csv(paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2_LOO_LIME.csv"))  # <-- Update path if needed

        lime_df$test <- str_count(lime_df$feature, ' ')

        lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
        lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])

        seu <-readRDS(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_",par_brain_region,"_",celltype2,"_int.rds")) 
        ncol(seu) 
        DefaultAssay(seu)
        print(unique(seu@meta.data$CellType)) ##
        unique(seu@meta.data$Group) ##

        ## Only retain group interest
        ## Subset to remove group that we do not want
        unique(seu@meta.data$Group)
            
            xx <- unique(seu@meta.data$Group)
            xx <- xx[!(xx %in% c(par_remove_group))]

            Idents(seu) <- "Group"
            seu=subset(seu,idents=xx)

        unique(seu@meta.data$Group) 
        nrow(seu@meta.data)

        unique_features <- unique(lime_df$feature)

        ## Need to separate by subject
        unique(seu@meta.data$orig.ident)
        unique(lime_df$test_donor)

        unique(seu@meta.data$orig.ident) %in% unique(lime_df$test_donor)

        # Initialize an empty list to store results
        express_fill_list <- list()

        # Get unique subjects from Seurat metadata
        subjects <- unique(seu@meta.data$orig.ident)
        genes <- unique_features  # your list of genes

        # Loop over each subject
        for (subject in subjects) {
        
        # Set idents
        Idents(seu) <- "orig.ident"
        
        # Get cells for this subject
        subject_cells <- WhichCells(seu, idents = subject)  # or filter with meta.data
        
        # Loop over each gene
        for (gene in genes) {
            
            # Extract expression for this gene in subject's cells
            expression <- seu@assays$RNA@data[gene, subject_cells]
            
            # Count cells where expression > 0
            cells_expressing <- expression[expression > 0]
            
            # Calculate percent of cells expressing the gene
            percent_expressed <- length(cells_expressing) / length(subject_cells)
            
            # Store the result
            express_fill_list[[paste(subject, gene, sep = "_")]] <- data.frame(
            subject = subject,
            gene = gene,
            percent_expressed = percent_expressed
            )
        }
        }

        # Combine all results into a single data frame
        express_fill_df <- do.call(rbind, express_fill_list)
        rownames(express_fill_df) <- NULL  # Optional: clean row names

        write.csv(express_fill_df, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LOO_percent_expressed_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_narval_2.csv"), sep = ",")

        ## Normalize the percent expressed value
        express_fill_df <- express_fill_df %>%
        group_by(subject) %>%
        mutate(normalized_percent_expressed = percent_expressed / sum(percent_expressed)) %>%
        ungroup()

        express_fill_df <- data.frame(express_fill_df)

        ## We need to take the aboslute value of LIME
        lime_df$importance <- abs(lime_df$importance)
        nrow(lime_df)

        ## merge lime with percent expressed
        lime_df2 <- merge(lime_df, express_fill_df, by.x = c("test_donor", "feature"), by.y = c("subject", "gene"))

        ## compute importance weighted by expression
        lime_df2$weighted_importance <- lime_df2$importance*lime_df2$normalized_percent_expressed

        ## Compute mean importance across cells from each subject
        donor_gene_matrix <- lime_df2 %>%
        group_by(test_donor, feature) %>%
        summarize(mean_weighted_importance = mean(weighted_importance), .groups = "drop") %>%
        pivot_wider(names_from = feature, values_from = mean_weighted_importance, values_fill = 0)

        dim(donor_gene_matrix)

        # Assume `lime_matrix` is a data frame: rows = subjects, columns = genes
        lime_matrix <- as.data.frame(donor_gene_matrix)
        n_subjects <- nrow(lime_matrix)
        gene_names <- colnames(lime_matrix)

        # Store significant genes in each LOSO fold
        loso_sig_genes <- list()
        loso_mean_delta_lime <- list()

        for (i in 1:n_subjects) {
        # Exclude subject i
        test_data <- lime_matrix[-i, ]  # Remaining subjects
        n_remaining <- nrow(test_data)

        # Global gene mean for each subject (i.e., across genes)
        test_data <- test_data[, sapply(test_data, is.numeric)]
        global_gene_means <- rowMeans(test_data)

        # Build a matrix of (LIME value - global mean) for each gene
        delta_matrix <- sweep(test_data, 1, global_gene_means, "-")

        # Now test if the deltas are significantly > 0
            pvals <- sapply(delta_matrix, function(delta_vals) {
            wilcox.test(delta_vals, mu = 0, alternative = "greater")$p.value
            })

        # Adjust for multiple testing
        adj_pvals <- p.adjust(pvals, method = "BH")

        # Store significant genes
        sig_genes <- names(adj_pvals[adj_pvals < 0.05])
        loso_sig_genes[[i]] <- sig_genes

        loso_mean_delta_lime[[i]] <- colMeans(delta_matrix)
        }

        # Intersect across LOSO folds
        robust_genes <- reduce(loso_sig_genes, intersect)
        length(robust_genes)

        ##############################################################################

        ## Calculate mean Delta across LOO
        delta_df <- as.data.frame(do.call(rbind, loso_mean_delta_lime))

        # Now calculate mean for each gene (i.e., column)
        mean_deltas <- colMeans(delta_df)

        ## Non sig deltas
        nonsig_deltas <- mean_deltas[!names(mean_deltas) %in% robust_genes]
        nonsig_deltas_df <- data.frame(nonsig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf", width = 7, height = 5)
        hist(nonsig_deltas_df$nonsig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = quantile(nonsig_deltas_df$nonsig_deltas, 0), col = "red", lty = 2)
        dev.off()

        threshold <- quantile(nonsig_deltas_df$nonsig_deltas, 0)

        ## Sig deltas
        sig_deltas <- mean_deltas[names(mean_deltas) %in% robust_genes]
        sig_deltas_df <- data.frame(sig_deltas)
        pdf("/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf", width = 7, height = 5)
        hist(sig_deltas_df$sig_deltas, breaks = 50, main = "Effect Size Distribution of Non-Significant Genes",
            xlab = "Mean Delta-LIME", col = "gray")
        abline(v = threshold, col = "red", lty = 2)
        dev.off()

        sig_deltas_df$feature <- rownames(sig_deltas_df)
        refined_genes <- sig_deltas_df[sig_deltas_df$sig_deltas > threshold, ]
        nrow(refined_genes)


        ############################################################################## 
        mean_deltas_df <- data.frame(
        gene = names(mean_deltas),
        mean_delta = as.numeric(mean_deltas)
        )

        mean_deltas_sorted <- mean_deltas_df %>%
        arrange(desc(mean_delta))

        mean_deltas_sorted$row <- 1:nrow(mean_deltas_sorted)
        mean_deltas_sorted <- subset(mean_deltas_sorted,gene %in%  refined_genes$feature )

        write.csv(mean_deltas_sorted, file = paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",par_prep,"_",par_status,"_",par_brain_region,"_",celltype2,"_0_abs_case_control_narval_2.csv"), sep = ",")


        ##############################################################################
        }
    ##

##









################################################################################################ Not using beyond this point. 
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################



module load StdEnv/2023
module load r/4.4.0
R

# Load required libraries
#install.packages("lubridate", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("tidyverse", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("cluster", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("factoextra", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("reshape2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("umap", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("purrr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#install.packages("ashr", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

library(ashr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(lubridate, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(purrr, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization
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



## code to see how many genes we have for each cell type

    par_brain_region = "BA4"
    par_status = "SALS"
    par_prep = "CombatSeq"
        
    celltype_list = c("L3_L5", "L2_L3",   "L4_L6",   "L4_L5",   "L5_L6",   "L6",           
            "PV",      "5HT3aR",  "Rosehip", "SOM",     "Oligo",   "Astro",   "OPC",    
            "Micro",  "Mural",   "Endo",    "Fibro", "L5")


    df_list <- list()
    
    for (celltype2 in celltype_list) {
    file_path <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_optimal_delta_",
                        par_prep, "_", par_status, "_", par_brain_region, "_", 
                        celltype2, "_0.9995_abs_case_control_narval_2.csv")
    
    # Read in as character to avoid type mismatch
    temp_df <- read_csv(file_path, col_types = cols(.default = col_character())) %>%
        mutate(celltype = celltype2)
    
    df_list[[celltype2]] <- temp_df
    }

    # Combine all dataframes
    merged_df <- bind_rows(df_list)

    table(merged_df$celltype)

    length(unique(merged_df$celltype))
##