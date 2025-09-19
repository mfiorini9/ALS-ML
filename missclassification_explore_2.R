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
    ## HVG diagnosis -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "ALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_diagnosis_prop_missclassification_2.csv')

    ##

    ###########################
    ## HVG genetic -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "Sporadic"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9orf72"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "Sporadic"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9orf72"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_genetic_prop_missclassification_2.csv')

    ##

    ###########################
    ## HVG disease -- done
    ###########################

    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_All_HVGs_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "SALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9ALS"
        df_avg_props$true_class[df_avg_props$true_class == "3"] <- "SFTLD"
        df_avg_props$true_class[df_avg_props$true_class == "4"] <- "C9FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "SALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "3"] <- "SFTLD"
        df_avg_props$pred_class[df_avg_props$pred_class == "4"] <- "C9FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_disease_prop_missclassification_2.csv')
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
    ## NMF diagnosis -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_ALS_FTLD_control_combined.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "ALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_diagnosis_prop_missclassification_2.csv')

    ##

    ###########################
    ## NMF genetic -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_LOSO_combat_generalizable_sporadic_C9_control_combined.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "Sporadic"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9orf72"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "Sporadic"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9orf72"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_genetic_prop_missclassification_2.csv')

    ##

    ###########################
    ## NMF disease -- done
    ###########################

    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/total_donor_NMF_fix_HVG_generalizable_subtype_combined.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "SALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9ALS"
        df_avg_props$true_class[df_avg_props$true_class == "3"] <- "SFTLD"
        df_avg_props$true_class[df_avg_props$true_class == "4"] <- "C9FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "SALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "3"] <- "SFTLD"
        df_avg_props$pred_class[df_avg_props$pred_class == "4"] <- "C9FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_disease_prop_missclassification_2.csv')
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
    ## HVG diagnosis -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_diagnosis.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "ALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_diagnosis_prop_missclassification_2.csv')

    ##

    ###########################
    ## HVG genetic -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_genetic.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "Sporadic"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9orf72"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "Sporadic"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9orf72"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_genetic_prop_missclassification_2.csv')

    ##

    ###########################
    ## HVG disease -- done
    ###########################

    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA4_disease.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "SALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9ALS"
        df_avg_props$true_class[df_avg_props$true_class == "3"] <- "SFTLD"
        df_avg_props$true_class[df_avg_props$true_class == "4"] <- "C9FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "SALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "3"] <- "SFTLD"
        df_avg_props$pred_class[df_avg_props$pred_class == "4"] <- "C9FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_disease_prop_missclassification_2.csv')
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
    ## NMF diagnosis -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_diagnosis.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "ALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_diagnosis_prop_missclassification_2.csv')

    ##

    ###########################
    ## NMF genetic -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_genetic.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "Sporadic"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9orf72"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "Sporadic"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9orf72"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_genetic_prop_missclassification_2.csv')

    ##

    ###########################
    ## NMF disease -- done
    ###########################

    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA4_disease.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "SALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9ALS"
        df_avg_props$true_class[df_avg_props$true_class == "3"] <- "SFTLD"
        df_avg_props$true_class[df_avg_props$true_class == "4"] <- "C9FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "SALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "3"] <- "SFTLD"
        df_avg_props$pred_class[df_avg_props$pred_class == "4"] <- "C9FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_disease_prop_missclassification_2.csv')
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
    ## HVG diagnosis -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_diagnosis.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "ALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_diagnosis_prop_missclassification_2.csv')

    ##

    ###########################
    ## HVG genetic -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_genetic.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "Sporadic"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9orf72"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "Sporadic"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9orf72"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_genetic_prop_missclassification_2.csv')

    ##

    ###########################
    ## HVG disease -- done
    ###########################

    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_BA9_disease.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "SALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9ALS"
        df_avg_props$true_class[df_avg_props$true_class == "3"] <- "SFTLD"
        df_avg_props$true_class[df_avg_props$true_class == "4"] <- "C9FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "SALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "3"] <- "SFTLD"
        df_avg_props$pred_class[df_avg_props$pred_class == "4"] <- "C9FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_disease_prop_missclassification_2.csv')
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
    ## NMF diagnosis -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_diagnosis.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "ALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_diagnosis_prop_missclassification_2.csv')

    ##

    ###########################
    ## NMF genetic -- done
    ###########################

    ## code
        ## read in file
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_genetic.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "Sporadic"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9orf72"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "Sporadic"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9orf72"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_genetic_prop_missclassification_2.csv')

    ##

    ###########################
    ## NMF disease -- done
    ###########################

    ## code
        BA4_bind_all <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_BA9_disease.csv')
        
        ## Retain best
        BA4_bind_all <- BA4_bind_all %>%
            group_by(donor, celltype) %>%
            slice_max(test_accuracy_all, n = 1, with_ties = FALSE) %>%
            ungroup()

        BA4_bind_all <- data.frame(BA4_bind_all)
        nrow(BA4_bind_all)

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
        BA4_bind_all <- BA4_bind_all %>% dplyr::select(donor, celltype, group, test_accuracy_all, true_class, counts_pred)


        #BA4_bind_all <- subset(BA4_bind_all, celltype == "L2_L3")

        df_avg_props <- BA4_bind_all %>%
            # split counts_pred into separate rows
            separate_rows(counts_pred, sep = ";") %>%
            # remove quotes and spaces
            mutate(counts_pred = str_trim(counts_pred),
                    pred_class = str_extract(counts_pred, "\\d+"),
                    count = as.numeric(str_extract(counts_pred, "\\d+$"))) %>%
            # group and compute proportions
            group_by(celltype, true_class, pred_class) %>%
            summarise(total = sum(count), .groups = "drop") %>%
            group_by(celltype, true_class) %>%
            mutate(proportion = total / sum(total)) %>%
            ungroup()

        ## prepare plot columns
        df_avg_props$true_class[df_avg_props$true_class == "0"] <- "Control"
        df_avg_props$true_class[df_avg_props$true_class == "1"] <- "SALS"
        df_avg_props$true_class[df_avg_props$true_class == "2"] <- "C9ALS"
        df_avg_props$true_class[df_avg_props$true_class == "3"] <- "SFTLD"
        df_avg_props$true_class[df_avg_props$true_class == "4"] <- "C9FTLD"

        df_avg_props$pred_class[df_avg_props$pred_class == "0"] <- "Control"
        df_avg_props$pred_class[df_avg_props$pred_class == "1"] <- "SALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "2"] <- "C9ALS"
        df_avg_props$pred_class[df_avg_props$pred_class == "3"] <- "SFTLD"
        df_avg_props$pred_class[df_avg_props$pred_class == "4"] <- "C9FTLD"

        write.csv(df_avg_props, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_disease_prop_missclassification_2.csv')
    ##

##

########################
## Mega heatmap plot
########################
## code 
    ## Read in summary files. 
    HVG_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_diagnosis_prop_missclassification_2.csv', sep = ",")
    HVG_diagnosis$model = "Diagnosis"
    HVG_diagnosis$region = "M/FCx"
    HVG_diagnosis$input = "HVG"

    HVG_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_genetic_prop_missclassification_2.csv', sep = ",")
    HVG_genetic$model = "Genetic group"
    HVG_genetic$region = "M/FCx"
    HVG_genetic$input = "HVG"

    HVG_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_disease_prop_missclassification_2.csv', sep = ",")
    HVG_disease$model = "Disease group"
    HVG_disease$region = "M/FCx"
    HVG_disease$input = "HVG"


    NMF_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_diagnosis_prop_missclassification_2.csv', sep = ",")
    NMF_diagnosis$model = "Diagnosis"
    NMF_diagnosis$region = "M/FCx"
    NMF_diagnosis$input = "NMF"

    NMF_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_genetic_prop_missclassification_2.csv', sep = ",")
    NMF_genetic$model = "Genetic group"
    NMF_genetic$region = "M/FCx"
    NMF_genetic$input = "NMF"

    NMF_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_disease_prop_missclassification_2.csv', sep = ",")
    NMF_disease$model = "Disease group"
    NMF_disease$region = "M/FCx"
    NMF_disease$input = "NMF"


    BA4_HVG_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_diagnosis_prop_missclassification_2.csv', sep = ",")
    BA4_HVG_diagnosis$model = "Diagnosis"
    BA4_HVG_diagnosis$region = "MCx"
    BA4_HVG_diagnosis$input = "HVG"

    BA4_HVG_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_genetic_prop_missclassification_2.csv', sep = ",")
    BA4_HVG_genetic$model = "Genetic group"
    BA4_HVG_genetic$region = "MCx"
    BA4_HVG_genetic$input = "HVG"

    BA4_HVG_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_disease_prop_missclassification_2.csv', sep = ",")
    BA4_HVG_disease$model = "Disease group"
    BA4_HVG_disease$region = "MCx"
    BA4_HVG_disease$input = "HVG"


    BA9_HVG_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_diagnosis_prop_missclassification_2.csv', sep = ",")
    BA9_HVG_diagnosis$model = "Diagnosis"
    BA9_HVG_diagnosis$region = "FCx"
    BA9_HVG_diagnosis$input = "HVG"

    BA9_HVG_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_genetic_prop_missclassification_2.csv', sep = ",")
    BA9_HVG_genetic$model = "Genetic group"
    BA9_HVG_genetic$region = "FCx"
    BA9_HVG_genetic$input = "HVG"

    BA9_HVG_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_disease_prop_missclassification_2.csv', sep = ",")
    BA9_HVG_disease$model = "Disease group"
    BA9_HVG_disease$region = "FCx"
    BA9_HVG_disease$input = "HVG"


    BA4_NMF_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_diagnosis_prop_missclassification_2.csv', sep = ",")
    BA4_NMF_diagnosis$model = "Diagnosis"
    BA4_NMF_diagnosis$region = "MCx"
    BA4_NMF_diagnosis$input = "NMF"

    BA4_NMF_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_genetic_prop_missclassification_2.csv', sep = ",")
    BA4_NMF_genetic$model = "Genetic group"
    BA4_NMF_genetic$region = "MCx"
    BA4_NMF_genetic$input = "NMF"

    BA4_NMF_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_disease_prop_missclassification_2.csv', sep = ",")
    BA4_NMF_disease$model = "Disease group"
    BA4_NMF_disease$region = "MCx"
    BA4_NMF_disease$input = "NMF"


    BA9_NMF_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_diagnosis_prop_missclassification_2.csv', sep = ",")
    BA9_NMF_diagnosis$model = "Diagnosis"
    BA9_NMF_diagnosis$region = "FCx"
    BA9_NMF_diagnosis$input = "NMF"

    BA9_NMF_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_genetic_prop_missclassification_2.csv', sep = ",")
    BA9_NMF_genetic$model = "Genetic group"
    BA9_NMF_genetic$region = "FCx"
    BA9_NMF_genetic$input = "NMF"

    BA9_NMF_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_disease_prop_missclassification_2.csv', sep = ",")
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

    bind_total_model$pred_class <- factor(bind_total_model$pred_class, levels = c("Control", "ALS", "FTLD"))
    bind_total_model$true_class <- factor(bind_total_model$true_class, levels = c("Control", "ALS", "FTLD"))

    temp <- ggplot(bind_total_model, aes(y = true_class, x = proportion,  fill = pred_class)) + 
        theme_bw() + 
        geom_col() +
        theme(
            legend.position = "none",
            legend.key.size = unit(0.2, "cm"),      
            legend.text = element_text(size = 6),   
            legend.title = element_text(size = 7),  
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title = element_text(face = "bold", size = 8),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested( input + region ~ celltype) +
        scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
        ylab("True label") +
        xlab("Proportion of DNN classifications") +
        scale_fill_manual(values = c("ALS" = "darkred",
                                    "FTLD"  = "darkblue",
                                    "Control" = "#339966")) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 4, width = 13)


    ## Plot genetic
    bind_total_model <- subset(bind_total, model == "Genetic group")

    bind_total_model$pred_class <- factor(bind_total_model$pred_class, levels = c("Control", "Sporadic", "C9orf72"))
    bind_total_model$true_class <- factor(bind_total_model$true_class, levels = c("Control", "Sporadic", "C9orf72"))

    temp <- ggplot(bind_total_model, aes(y = true_class, x = proportion,  fill = pred_class)) + 
        theme_bw() + 
        geom_col() +
        theme(
            legend.position = "none",
            legend.key.size = unit(0.2, "cm"),      
            legend.text = element_text(size = 6),   
            legend.title = element_text(size = 7),  
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = "black"),
            axis.title = element_text(face = "bold", size = 8),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested( input + region ~ celltype) +
        scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
        ylab("True label") +
        xlab("Proportion of DNN classifications") +
        scale_fill_manual(values = c("Sporadic" = "#ec7014",
                                    "C9orf72"  = "#b0b0b0",
                                    "Control" = "#339966")) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 4, width = 13)


    ## Plot disease
    bind_total_model <- subset(bind_total, model == "Disease group")

    bind_total_model$pred_class <- factor(bind_total_model$pred_class, levels = c("Control", "SALS", "C9ALS", "SFTLD", "C9FTLD"))
    bind_total_model$true_class <- factor(bind_total_model$true_class, levels = c("Control", "SALS", "C9ALS", "SFTLD", "C9FTLD"))

    temp <- ggplot(bind_total_model, aes(y = true_class, x = proportion,  fill = pred_class)) + 
        theme_bw() + 
        geom_col() +
        theme(
            legend.position = "none",
            legend.key.size = unit(0.2, "cm"),      
            legend.text = element_text(size = 6),   
            legend.title = element_text(size = 7),  
            panel.grid = element_blank(),
            axis.text.x = element_text(size = 8, colour = "black"),
            axis.text.y = element_text(size = 8, colour = c("#339966","orange", "red", "blue", "purple" )),
            axis.title = element_text(face = "bold", size = 8),
            strip.background = element_rect(fill="lightgrey", colour = "white"),
            strip.text = element_text(size= 8, face="bold", colour = "black"),
            strip.placement = "outside"
        ) +
        facet_nested( input + region ~ celltype) +
        scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
        ylab("True label") +
        xlab("Proportion of DNN classifications") +
        scale_fill_manual(values = c("SALS" = "orange",
                                    "C9ALS"  = "red",
                                    "SFTLD" = "blue",
                                    "C9FTLD"  = "purple",
                                    "Control" = "#339966")) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 4.4, width = 13)

##

########################
## Quantify L5 explore. 
########################

## Read in summary files. 
    HVG_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_diagnosis_prop_missclassification_2.csv', sep = ",")
    HVG_diagnosis$model = "Diagnosis"
    HVG_diagnosis$region = "M/FCx"
    HVG_diagnosis$input = "HVG"

    HVG_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_genetic_prop_missclassification_2.csv', sep = ",")
    HVG_genetic$model = "Genetic group"
    HVG_genetic$region = "M/FCx"
    HVG_genetic$input = "HVG"

    HVG_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/HVG_disease_prop_missclassification_2.csv', sep = ",")
    HVG_disease$model = "Disease group"
    HVG_disease$region = "M/FCx"
    HVG_disease$input = "HVG"


    NMF_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_diagnosis_prop_missclassification_2.csv', sep = ",")
    NMF_diagnosis$model = "Diagnosis"
    NMF_diagnosis$region = "M/FCx"
    NMF_diagnosis$input = "NMF"

    NMF_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_genetic_prop_missclassification_2.csv', sep = ",")
    NMF_genetic$model = "Genetic group"
    NMF_genetic$region = "M/FCx"
    NMF_genetic$input = "NMF"

    NMF_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/NMF_disease_prop_missclassification_2.csv', sep = ",")
    NMF_disease$model = "Disease group"
    NMF_disease$region = "M/FCx"
    NMF_disease$input = "NMF"


    BA4_HVG_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_diagnosis_prop_missclassification_2.csv', sep = ",")
    BA4_HVG_diagnosis$model = "Diagnosis"
    BA4_HVG_diagnosis$region = "MCx"
    BA4_HVG_diagnosis$input = "HVG"

    BA4_HVG_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_genetic_prop_missclassification_2.csv', sep = ",")
    BA4_HVG_genetic$model = "Genetic group"
    BA4_HVG_genetic$region = "MCx"
    BA4_HVG_genetic$input = "HVG"

    BA4_HVG_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_HVG_disease_prop_missclassification_2.csv', sep = ",")
    BA4_HVG_disease$model = "Disease group"
    BA4_HVG_disease$region = "MCx"
    BA4_HVG_disease$input = "HVG"


    BA9_HVG_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_diagnosis_prop_missclassification_2.csv', sep = ",")
    BA9_HVG_diagnosis$model = "Diagnosis"
    BA9_HVG_diagnosis$region = "FCx"
    BA9_HVG_diagnosis$input = "HVG"

    BA9_HVG_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_genetic_prop_missclassification_2.csv', sep = ",")
    BA9_HVG_genetic$model = "Genetic group"
    BA9_HVG_genetic$region = "FCx"
    BA9_HVG_genetic$input = "HVG"

    BA9_HVG_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_HVG_disease_prop_missclassification_2.csv', sep = ",")
    BA9_HVG_disease$model = "Disease group"
    BA9_HVG_disease$region = "FCx"
    BA9_HVG_disease$input = "HVG"


    BA4_NMF_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_diagnosis_prop_missclassification_2.csv', sep = ",")
    BA4_NMF_diagnosis$model = "Diagnosis"
    BA4_NMF_diagnosis$region = "MCx"
    BA4_NMF_diagnosis$input = "NMF"

    BA4_NMF_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_genetic_prop_missclassification_2.csv', sep = ",")
    BA4_NMF_genetic$model = "Genetic group"
    BA4_NMF_genetic$region = "MCx"
    BA4_NMF_genetic$input = "NMF"

    BA4_NMF_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA4_NMF_disease_prop_missclassification_2.csv', sep = ",")
    BA4_NMF_disease$model = "Disease group"
    BA4_NMF_disease$region = "MCx"
    BA4_NMF_disease$input = "NMF"


    BA9_NMF_diagnosis <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_diagnosis_prop_missclassification_2.csv', sep = ",")
    BA9_NMF_diagnosis$model = "Diagnosis"
    BA9_NMF_diagnosis$region = "FCx"
    BA9_NMF_diagnosis$input = "NMF"

    BA9_NMF_genetic <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_genetic_prop_missclassification_2.csv', sep = ",")
    BA9_NMF_genetic$model = "Genetic group"
    BA9_NMF_genetic$region = "FCx"
    BA9_NMF_genetic$input = "NMF"

    BA9_NMF_disease <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/model_outs/BA9_NMF_disease_prop_missclassification_2.csv', sep = ",")
    BA9_NMF_disease$model = "Disease group"
    BA9_NMF_disease$region = "FCx"
    BA9_NMF_disease$input = "NMF"

    bind_total <- rbind(HVG_diagnosis, HVG_genetic, HVG_disease,
                        NMF_diagnosis, NMF_genetic, NMF_disease,
                        BA4_HVG_diagnosis, BA4_HVG_genetic, BA4_HVG_disease,
                        BA9_HVG_diagnosis, BA9_HVG_genetic, BA9_HVG_disease,
                        BA4_NMF_diagnosis, BA4_NMF_genetic, BA4_NMF_disease,
                        BA9_NMF_diagnosis, BA9_NMF_genetic, BA9_NMF_disease)

    ## Extract misclassifiation performance

    bind_total_model <- subset(bind_total, model == "Disease group")

    # Define excitatory subtypes
    exc_neurons <- c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6')

    # Filter to C9FTLD true class and C9ALS predicted
    df_ftld_to_als <- bind_total_model %>%
        filter(true_class == "C9FTLD", pred_class == "C9ALS", celltype %in% exc_neurons)

    # Summarize proportions per excitatory subtype
    misclass_summary <- df_ftld_to_als %>%
        select(celltype, total, proportion)

    print(misclass_summary)

    ## Build contingency table for L5 vs others
    # L5 counts
    l5_mis <- sum(df_ftld_to_als %>% filter(celltype == "L5") %>% pull(total))
    l5_total <- sum(bind_total_model %>% filter(celltype == "L5", true_class == "C9FTLD") %>% pull(total))
    l5_correct <- l5_total - l5_mis

    # Other excitatory counts
    others_mis <- sum(df_ftld_to_als %>% filter(celltype != "L5") %>% pull(total))
    others_total <- sum(bind_total_model %>% filter(celltype %in% exc_neurons & celltype != "L5", 
                                    true_class == "C9FTLD") %>% pull(total))
    others_correct <- others_total - others_mis

    # Contingency table
    tab <- matrix(c(l5_mis, l5_correct,
                    others_mis, others_correct),
                nrow = 2, byrow = TRUE,
                dimnames = list(
                    Group = c("L5", "Other_excitatory"),
                    Classification = c("Misclassified_as_C9ALS", "Not_C9ALS")
                ))

    print(tab)

    # Statistical test 
    chisq_res <- chisq.test(tab)
    fisher_res <- fisher.test(tab)  # safer if counts are small

    print(chisq_res)
    print(fisher_res)

    # --- 4. Visualization ---
    df_ftld_to_als %>%
    ggplot(aes(x = celltype, y = proportion)) +
    geom_col(fill = "darkred") +
    theme_bw() +
    ylab("Proportion misclassified as C9ALS") +
    xlab("Excitatory subtype (C9FTLD true class)") +
    ggtitle("C9FTLD --> C9ALS misclassification by cell type")

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 4, width = 13)

##








    

