salloc -A def-grouleau --time=0-1 -c 1 --mem=10g

## The outputs from LIME_outputs_explore3.R will go into this script.

## The results will be placed here
#/home/fiorini9/scratch/machine_learning_ALS/ALSKP_rare_burden


module load StdEnv/2023
module load r/4.4.0
R

#setwd("C:/Users/Alli_/Documents/PostDoc/Functional Domain_ALS_Collaboration")
#install.packages('data.table', lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
#install.packages('R.utils', lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
#install.packages('dplyr', lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
#install.packages('splitstackshape', lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
#install.packages('purrr', lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
#install.packages('logistf', lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
#install.packages('stringr', lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")

library(data.table, lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
library(R.utils, lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
library(dplyr, lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
library(splitstackshape, lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
library(purrr, lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
library(logistf, lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
library(stringr, lib ="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
library(ggplot2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
library(ggrepel, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R") #good
library(labeling, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

ALSKP <- fread('/home/fiorini9/projects/def-sfarhan/AALS/ALS_cohort_genetic_data/2024_10_03_ALSKP_BayesDel_AlphaMiss_SpliceAI_REVEL.csv', sep=',', header=TRUE)

#ALSKP <- subset(ALSKP, gene == "MICA")
#ALSKP <- subset(ALSKP, consequence == "missense" | consequence == "lof")
#write.csv(ALSKP, '/home/fiorini9/projects/def-sfarhan/AALS/ALS_cohort_genetic_data/MICA_2024_10_03_ALSKP_BayesDel_AlphaMiss_SpliceAI_REVEL.csv')

ALSKP$REVEL

ALSKP$isAffected <- ifelse(ALSKP$case_control == "Case", 1, 0)
ALSKP <- ALSKP %>% mutate(SpliceAI_Overall = pmax(SpliceAI_pred_DS_AG, SpliceAI_pred_DS_AL, SpliceAI_pred_DS_DG, SpliceAI_pred_DS_DL))

table(ALSKP$case_control)
table(ALSKP$isAffected)

############################
#### check carrier counts
############################
ALSKP %>%
  group_by(case_control) %>%
  summarise(unique_samples = n_distinct(sample))

############################
#### Allele frequency -- already filtered to AF < 0.01
############################
head(ALSKP$allele_frequency)
max(ALSKP$allele_frequency)


############################
####Gene Burden
############################

Gene_Burden_Counts <- function(Data){
  syn <- Data %>% 
    filter(consequence == "synonymous") %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("syn" = "n")
  miss <- Data %>% 
    filter(consequence == "missense") %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("miss" = "n")
  miss_path <- Data %>% 
    filter(consequence == "missense") %>%
    filter(am_class == "pathogenic" | REVEL >= 0.75) %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("miss_path" = "n")
  
  any_path <- Data %>% 
    filter(am_class == "pathogenic" | REVEL >= 0.75) %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("any_path" = "n")
  
  splicing <- Data %>% 
    filter(SpliceAI_Overall >= 0.5) %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("splicing" = "n")
  
  lof <- Data %>% 
    filter(consequence == "lof") %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("lof" = "n")
  
  lof_splice_anypath <- Data %>% 
    filter(am_class == "pathogenic" | REVEL >= 0.75 | SpliceAI_Overall >= 0.5 | consequence == "lof") %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("lof_splice_anypath" = "n")
  
  output <- merge(syn, miss, all = TRUE) %>% merge(miss_path, all = TRUE) %>% 
    merge(any_path, all = TRUE) %>% merge(splicing, all = TRUE) %>% merge(lof, all = TRUE) %>% merge(lof_splice_anypath, all = TRUE)

  #output <- merge(syn, miss, all = TRUE) %>% merge(lof, all = TRUE)
  return(output) 
}

Gene_Burden_Counts_alphaM <- function(Data){
  syn <- Data %>% 
    filter(consequence == "synonymous") %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("syn" = "n")
  miss <- Data %>% 
    filter(consequence == "missense") %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("miss" = "n")
  miss_path <- Data %>% 
    filter(consequence == "missense") %>%
    filter(am_class == "pathogenic") %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("miss_path" = "n")
  
  any_path <- Data %>% 
    filter(am_class == "pathogenic" ) %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("any_path" = "n")
  
  splicing <- Data %>% 
    filter(SpliceAI_Overall >= 0.5) %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("splicing" = "n")
  
  lof <- Data %>% 
    filter(consequence == "lof") %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("lof" = "n")
  
  lof_splice_anypath <- Data %>% 
    filter(am_class == "pathogenic" | SpliceAI_Overall >= 0.5 | consequence == "lof") %>%
    group_by(gene, case_control) %>%
    count() %>%
    rename("lof_splice_anypath" = "n")
  
  output <- merge(syn, miss, all = TRUE) %>% merge(miss_path, all = TRUE) %>% 
    merge(any_path, all = TRUE) %>% merge(splicing, all = TRUE) %>% merge(lof, all = TRUE) %>% merge(lof_splice_anypath, all = TRUE)

  #output <- merge(syn, miss, all = TRUE) %>% merge(lof, all = TRUE)
  return(output) 
}

ALSKP_Gene_Counts <- Gene_Burden_Counts(ALSKP)

all_combinations <- expand.grid(gene = unique(ALSKP$gene),
                                case_control = unique(ALSKP$case_control))

ALSKP_Gene_Counts <- merge(ALSKP_Gene_Counts, all_combinations, by = c("gene", "case_control"), all.y = TRUE)

ALSKP_Gene_Counts[is.na(ALSKP_Gene_Counts)] <- 0

ALSKP_Gene_Counts <- ALSKP_Gene_Counts %>% filter(gene != "")

casevscontrol_counts <- merge(ALSKP_Gene_Counts %>% filter( case_control == "Case") %>% 
                                                            rename(syn_case = syn, miss_case = miss, miss_path_case = miss_path,
                                                                   any_path_case = any_path, splicing_case = splicing, lof_case = lof,
                                                                   lof_splice_anypath_case = lof_splice_anypath) %>%
                                                            select(!c(case_control)),
                                     ALSKP_Gene_Counts %>% filter(case_control == "Control") %>% 
                                       rename(syn_control = syn, miss_control = miss, miss_path_control = miss_path,
                                              any_path_control = any_path, splicing_control = splicing, lof_control = lof,
                                              lof_splice_anypath_control = lof_splice_anypath) %>%
                                       select(!c(case_control)))

## fisher syn
fisher_syn <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$syn_case[i]
    control_carrier <- df$syn_control[i]
    case_noncarrier <- 3864 - df$syn_case[i]
    control_noncarrier <- 7839 - df$syn_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$syn_pval <- pval
  df$syn_OR <- OR
  df$syn_lowci <- lowci
  df$syn_highci <- highci
  df[order(df$syn_pval),]
}

## miss fisher test 
fisher_miss <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$miss_case[i]
    control_carrier <- df$miss_control[i]
    case_noncarrier <- 3864 - df$miss_case[i]
    control_noncarrier <- 7839 - df$miss_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$miss_pval <- pval
  df$miss_OR <- OR
  df$miss_lowci <- lowci
  df$miss_highci <- highci
  df[order(df$miss_pval),]
}

## miss_path fisher test 
fisher_miss_path <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$miss_path_case[i]
    control_carrier <- df$miss_path_control[i]
    case_noncarrier <- 3864 - df$miss_path_case[i]
    control_noncarrier <- 7839 - df$miss_path_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$miss_path_pval <- pval
  df$miss_path_OR <- OR
  df$miss_path_lowci <- lowci
  df$miss_path_highci <- highci
  df[order(df$miss_path_pval),]
}

## fisher any path
fisher_any_path <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$any_path_case[i]
    control_carrier <- df$any_path_control[i]
    case_noncarrier <- 3864 - df$any_path_case[i]
    control_noncarrier <- 7839 - df$any_path_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$any_path_pval <- pval
  df$any_path_OR <- OR
  df$any_path_lowci <- lowci
  df$any_pathp_highci <- highci
  df[order(df$any_path_pval),]
}

## fisher splicing
fisher_splicing <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$splicing_case[i]
    control_carrier <- df$splicing_control[i]
    case_noncarrier <- 3864 - df$splicing_case[i]
    control_noncarrier <- 7839 - df$splicing_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$splicing_pval <- pval
  df$splicing_OR <- OR
  df$splicing_lowci <- lowci
  df$splicing_highci <- highci
  df[order(df$splicing_pval),]
}

## lof fisher test 
fisher_lof <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$lof_case[i]
    control_carrier <- df$lof_control[i]
    case_noncarrier <- 3864 - df$lof_case[i]
    control_noncarrier <- 7839 - df$lof_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$lof_pval <- pval
  df$lof_OR <- OR
  df$lof_lowci <- lowci
  df$lof_highci <- highci
  df[order(df$lof_pval),]
}

## lof_splice_anypath fisher test 
fisher_lof_splice_anypath <- function(mydata){
  df <- data.frame(mydata)
  pval <- NA
  OR <- NA
  lowci <- NA
  highci <- NA
  for (i in 1:nrow(df)) {
    case_carrier <- df$lof_splice_anypath_case[i]
    control_carrier <- df$lof_splice_anypath_control[i]
    case_noncarrier <- 3864 - df$lof_splice_anypath_case[i]
    control_noncarrier <- 7839 - df$lof_splice_anypath_control[i]
    tbl <- matrix(c(case_carrier,case_noncarrier,control_carrier,control_noncarrier),2,2)
    mod <- fisher.test(tbl)
    pval[i] <- mod$p.value
    OR[i] <- mod$estimate
    lowci[i] <- mod$conf.int[1]
    highci[i] <- mod$conf.int[2]
    print(i)
  }
  df$lof_splice_anypath_pval <- pval
  df$lof_splice_anypath_OR <- OR
  df$lof_splice_anypath_lowci <- lowci
  df$lof_splice_anypath_highci <- highci
  df[order(df$lof_splice_anypath_pval),]
}

## the first row is gene = "0", for some reason, we will remove this. 
nrow(casevscontrol_counts)

casevscontrol_counts <- subset(casevscontrol_counts, gene != "TTN")

#merge fishers per gene set together 
#Enrichment_fisher_ALSKP <- merge(fisher_syn(casevscontrol_counts), fisher_miss(casevscontrol_counts), all = TRUE) %>%
#  merge(fisher_lof(casevscontrol_counts), all = TRUE)

Enrichment_fisher_ALSKP <- merge(fisher_syn(casevscontrol_counts), fisher_miss(casevscontrol_counts), all = TRUE) %>%
  merge(fisher_miss_path(casevscontrol_counts), all = TRUE) %>%
  merge(fisher_any_path(casevscontrol_counts), all = TRUE) %>%
  merge(fisher_splicing(casevscontrol_counts), all = TRUE) %>%
  merge(fisher_lof(casevscontrol_counts), all = TRUE) %>%
  merge(fisher_lof_splice_anypath(casevscontrol_counts), all = TRUE)


##########################
## Save csv file
##########################
write.csv(Enrichment_fisher_ALSKP,'/home/fiorini9/scratch/machine_learning_ALS/ALSKP_rare_burden/ALSKP_gene_burden_all_genes.csv')


##########################
## READ IN GENES OF INTEREST and subset the file to only include those genes
##########################
Enrichment_fisher_ALSKP <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/ALSKP_rare_burden/ALSKP_gene_burden_all_genes.csv', sep = ',') 

LIME_optimal_sets <- read.delim('/home/fiorini9/scratch/machine_learning_ALS/LIME_gene_set_outputs/temp_pineda_LIME_lambda_0.00025_optimal_gene_sets.csv', sep = ',')
unique_features <- unique(LIME_optimal_sets$feature)

Enrichment_fisher_ALSKP_LIME <- subset(Enrichment_fisher_ALSKP, gene %in% unique_features)

length(unique(Enrichment_fisher_ALSKP_LIME$gene))

write.csv(Enrichment_fisher_ALSKP_LIME,'/home/fiorini9/scratch/machine_learning_ALS/ALSKP_rare_burden/temp_ALSKP_gene_burden_LIME_lambda_0.00025_optimal_gene_sets.csv')
