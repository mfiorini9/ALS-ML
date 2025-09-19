LIME_LOO_values_CombatSeq_SFTLD_BA9_Rosehip_narval_2_LOO_LIME.csv


salloc -A def-sfarhan --time=0-6 -c 1 --mem=40g

module load StdEnv/2023
module load r/4.4.0
R


## load libraries
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
library(labeling, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")


# Load required libraries
install.packages("lubridate", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
install.packages("tidyverse", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
install.packages("cluster", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
install.packages("factoextra", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
install.packages("reshape2", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
install.packages("umap", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(cluster, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")      # silhouette
library(factoextra, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")   # clustering helpers
library(reshape2, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")     # for reshaping
library(umap, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")         # optional: visualization

# Load your LIME data (replace with your actual CSV file)
lime_df <- read_csv("/home/fiorini9/scratch/machine_learning_ALS/model_outs/LIME_LOO_values_CombatSeq_SFTLD_BA9_Rosehip_narval_2_LOO_LIME.csv")  # <-- Update path if needed

# Clean feature names (optional)
lime_df$test <- str_count(lime_df$feature, ' ')

lime_df$feature[lime_df$test == 4] <- sub(".* .* (.*) .* .*", "\\1", lime_df$feature[lime_df$test == 4])
lime_df$feature[lime_df$test == 2] <- sub("(.*) .* .*", "\\1", lime_df$feature[lime_df$test == 2])
        

# Aggregate: Average LIME importance per donor per gene
donor_gene_matrix <- lime_df %>%
  group_by(test_donor, feature) %>%
  summarize(mean_importance = mean(importance), .groups = "drop") %>%
  pivot_wider(names_from = feature, values_from = mean_importance, values_fill = 0)

# Save donor IDs and make matrix
donor_ids <- donor_gene_matrix$test_donor
lime_mat <- donor_gene_matrix %>% select(-test_donor)
rownames(lime_mat) <- donor_ids

# Rank genes by average importance (absolute)
ranked_genes <- lime_df %>%
  group_by(feature) %>%
  summarize(mean_abs_importance = mean(abs(importance))) %>%
  arrange(desc(mean_abs_importance)) %>%
  pull(feature)

# Function to compute average silhouette score
evaluate_clustering <- function(mat, k = 2) {
  dist_mat <- dist(scale(mat))
  clustering <- kmeans(mat, centers = k, nstart = 10)
  sil <- silhouette(clustering$cluster, dist_mat)
  mean(sil[, 3])
}

# Try different top-N gene sets
gene_counts <- seq(10, 3000, by = 10)
silhouette_scores <- numeric(length(gene_counts))

for (i in seq_along(gene_counts)) {
  top_genes <- ranked_genes[1:gene_counts[i]]
  sub_mat <- lime_mat[, top_genes, drop = FALSE]
  silhouette_scores[i] <- evaluate_clustering(sub_mat)
}

# Plot: Number of genes vs silhouette score
plot(gene_counts, silhouette_scores, type = "b", pch = 19,
     xlab = "Number of Top LIME Genes",
     ylab = "Average Silhouette Score (k = 2)",
     main = "Donor Similarity in LIME Space")

# Optional: Visualize donor similarity with UMAP
set.seed(123)
top_genes_final <- ranked_genes[1:50]  # ← adjust if you pick a different cutoff
umap_out <- umap(scale(lime_mat[, top_genes_final]))

umap_df <- as.data.frame(umap_out$layout)
umap_df$donor_id <- donor_ids

ggplot(umap_df, aes(V1, V2, label = donor_id)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  ggtitle("UMAP of Donors Based on LIME Profiles")