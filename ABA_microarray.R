
module load StdEnv/2023
module load r/4.4.0
R

library(Seurat)
library(MAST)
library(ggplot2)
#library(ggpubr)
library(tidyr)
library(dplyr)
#library(SeuratData, lib="/home/fiorini9/R/x86_64-pc-linux-gnu-library/4.4")
library(SeuratDisk)
library(DESeq2)

/home/fiorini9/scratch/machine_learning_NPJ_review/machine_learning/ABA_microarray/XU180V1OH012_0zscore

columns <- read.csv('/home/fiorini9/scratch/machine_learning_NPJ_review/machine_learning/ABA_microarray/XU180V1OH012_0zscore/Columns.csv')
expression <- read.csv('/home/fiorini9/scratch/machine_learning_NPJ_review/machine_learning/ABA_microarray/XU180V1OH012_0zscore/Expression.csv')

dim(columns)
dim(expression)

unique(columns$structure_name)
unique(columns$donor_id)
unique(columns$top_level_structure_name)
#substantia nigra

keep <- which(columns$structure_name == "substantia nigra")
keep_mod <- keep + 1

sn <- expression[,c(111,280,449,618,787,956)]
non_sn <- expression[,-c(111,280,449,618,787,956)]
dim(sn)
dim(non_sn)

sn_t <- t(sn)
sn_t <- data.frame(sn_t)
sn_t <- sn_t[-1,]
sn_t <- data.frame(sn_t)
sn_t$region <- "sn"
colnames(sn_t) <- c("expression", "region")

non_sn_t <- t(non_sn)
non_sn_t <- data.frame(non_sn_t)
non_sn_t <- non_sn_t[-1,]
non_sn_t <- data.frame(non_sn_t)
non_sn_t$region <- "non_sn"
colnames(non_sn_t) <- c("expression", "region")

bind <- rbind(non_sn_t, sn_t)

head(bind)

ggplot(bind, aes(y = expression, x = region)) +
geom_violin() +
theme_classic()
ggsave('/home/fiorini9/scratch/machine_learning_NPJ_review/temp.pdf', height = 4, width = 3)



