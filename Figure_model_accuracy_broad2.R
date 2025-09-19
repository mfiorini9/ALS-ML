## We are running this in Narval
salloc -A def-grouleau --time=0-8 -c 1 --mem=50g

module load StdEnv/2023
module load r/4.4.0
R


## Load libraries
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
#install.packages("RColorBrewer", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Number of highly variable genes
## code
    
    #BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_SALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9ALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SFTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all$condition[BA4_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c("Ex", 'In','Glia', 'Vasc'))

        ggplot(BA4_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_nHVGs_broad.pdf', height = 2, width = 5)
    ##

    #BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_SALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9ALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SFTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')
        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            HVGs <- as.numeric(temp$n_HVGs)
            result_df <- data.frame(HVGs = HVGs)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all$condition[BA9_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c("Ex", 'In','Glia', 'Vasc'))

        ggplot(BA9_bind_all, aes(x = condition, y = HVGs, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n HVGs") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_nHVGs_broad.pdf', height = 2, width = 5)
    ##

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Number of cells
## code
    
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
            result_df <- data.frame(n_cells = n_cells)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_SALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9ALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SFTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c('Ex', 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################
    ## code
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all$condition[BA4_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA4_bind_all$condition <- factor(BA4_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA4_bind_all$celltype <- factor(BA4_bind_all$celltype, levels = c('Ex','In','Glia','Vasc'))


        ggplot(BA4_bind_all, aes(x = condition, y = n_cells, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n cells") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_nCells_broad.pdf', height = 2, width = 5)
    ##

    ## BA9

    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
            result_df <- data.frame(n_cells = n_cells)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df
        }

        final_df_SALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9ALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SFTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c('Ex', 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file_1 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/test_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_1 <- read.delim(file_1, header = T, sep = ",")
            length_1 <- nrow(temp_1)

            file_2 <- paste0("/home/fiorini9/scratch/machine_learning_ALS/model_outs/train_meta_",condition,"_",region,"_",i,"_pineda_narval.csv")
            temp_2 <- read.delim(file_2, header = T, sep = ",")
            length_2 <- nrow(temp_2)

            n_cells = length_1 + length_2
                        result_df <- data.frame(n_cells = n_cells)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################
    ## code
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all$condition[BA9_bind_all$condition == "all_conditions"] <- "All conditions"

        ## set factor levels
        BA9_bind_all$condition <- factor(BA9_bind_all$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        BA9_bind_all$celltype <- factor(BA9_bind_all$celltype, levels = c('Ex','In','Glia','Vasc'))


        ggplot(BA9_bind_all, aes(x = condition, y = n_cells, fill = condition, colour = condition)) + 
        theme_bw() + 
        geom_bar(stat = 'identity') +
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("n cells") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black"))
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_nCells_broad.pdf', height = 2, width = 5)
    ##

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ DNN 5CV classification balanced accuracies (Supplemental)
## code
    
    ## BA4
    
    ###########################
    ## SALS BA4
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA4
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9ALS_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SFTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA4
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA4
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA4"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA4 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA4
    ###########################

    ## code for barplot
        BA4_bind_all <- rbind(final_df_SALS_BA4, final_df_C9ALS_BA4, final_df_SFTLD_BA4, final_df_C9FTLD_BA4, final_df_all_conditions_BA4)
        BA4_bind_all$condition[BA4_bind_all$condition == "all_conditions"] <- "All conditions"

        ## Only retain the 5CV values, remove test (we will plot test in the main figure)
        BA4_bind_all <- subset(BA4_bind_all, test == "validation" )
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA4_bind_all %>%
            group_by(celltype, condition) %>%
            summarise(
                median_accuracy = median(accuracy), 
                mean_accuracy = mean(accuracy),           # Calculate the median of X0
                sd_accuracy = sd(accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(accuracy - median(accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$condition <- factor(summary_stats$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c("Ex", 'In', 'Vasc', 'Glia' ))

        ## Plot
        ggplot(summary_stats, aes(x = condition, y = mean_accuracy, fill = condition, colour = condition)) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Balanced\naccuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA4_5CV_balanced_acc_barplot_broad.pdf', height = 2, width = 5)
    ##

    ## BA9
    
    ###########################
    ## SALS BA9
    ###########################
    ## code
        ## define variables
        condition = "SALS" #one ofSALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9ALS BA9
    ###########################
    ## code
        ## define variables
        condition = "C9ALS" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9ALS_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## SFTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "SFTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_SFTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## C9FTLD BA9
    ###########################
    ## code
        ## define variables
        condition = "C9FTLD" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_C9FTLD_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## all_conditions BA9
    ###########################
    ## code
        ## define variables
        condition = "all_conditions" #one of SALS, C9ALS, SFTLD, C9FTLD, or all_conditions
        region = "BA9"
        cell_type_list <- c("Ex", 'In', 'Vasc', 'Glia')

        result_list <- list()

        for(i in unique(cell_type_list)){
            file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/model_outs/model_report_CombatSeq_',condition,'_',region,'_',i,'_narval_2.csv')
            temp <- read.delim(file, header = T, sep = ",")
            accuracies <- as.numeric(unlist(strsplit(temp$five_fold_accuracies, ";")))
            all_accuracies <- c(accuracies, temp$test_accuracy)
            test_column <- c(rep("validation", length(accuracies)), "test")
            result_df <- data.frame(accuracy = all_accuracies, test = test_column)
            result_df$celltype <- i
            result_df$condition <- condition
            result_df$region <- region

            result_list[[i]] <- result_df

        }

        final_df_all_conditions_BA9 <- do.call(rbind, result_list)
    ##

    ###########################
    ## merge and plot BA9
    ###########################

    ## code for barplot
        BA9_bind_all <- rbind(final_df_SALS_BA9, final_df_C9ALS_BA9, final_df_SFTLD_BA9, final_df_C9FTLD_BA9, final_df_all_conditions_BA9)
        BA9_bind_all$condition[BA9_bind_all$condition == "all_conditions"] <- "All conditions"

        ## Only retain the 5CV values, remove test (we will plot test in the main figure)
        BA9_bind_all <- subset(BA9_bind_all, test == "validation" )
        
        ## Calculate the mean and standard deviation
        summary_stats <- BA9_bind_all %>%
            group_by(celltype, condition) %>%
            summarise(
                median_accuracy = median(accuracy), 
                mean_accuracy = mean(accuracy),           # Calculate the median of X0
                sd_accuracy = sd(accuracy),                  # Calculate the standard deviation of X0
                mad_accuracy = median(abs(accuracy - median(accuracy))),  # Calculate MAD: median of absolute deviations from the median
                .groups = "drop"                 # Drop the grouping after summarising
            )

        ## set factor levels
        summary_stats$condition <- factor(summary_stats$condition, levels = c("SALS", "C9ALS", "SFTLD", "C9FTLD", "All conditions"))
        summary_stats$celltype <- factor(summary_stats$celltype, levels = c("Ex", 'In', 'Vasc', 'Glia' ))

        ## Plot
        ggplot(summary_stats, aes(x = condition, y = mean_accuracy, fill = condition, colour = condition)) + 
        geom_errorbar(aes(ymin=mean_accuracy - sd_accuracy, ymax=mean_accuracy + sd_accuracy), width=.2,
            position=position_dodge(.9), colour = "black") +
        geom_bar(stat="identity", position=position_dodge()) +
        geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey") +
        theme_bw() + 
        theme(
            legend.position = "none",
            panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size =8, colour = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold")
        ) +
        facet_grid(~celltype, scales = "free_x") +
        ylab("Balanced\naccuracy") +
        scale_fill_manual(values = c("orange", "red", "blue", "purple", "black")) +
        scale_colour_manual(values = c("orange", "red", "blue", "purple", "black")) #+
        #coord_cartesian(ylim=c(0.5,1.0)) 
        ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/BA9_5CV_balanced_acc_barplot_broad.pdf', height = 2, width = 5)
    ##

##