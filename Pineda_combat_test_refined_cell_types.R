salloc -A def-sfarhan --time=0-8 -c 1 --mem=200g
module load StdEnv/2020 
module load r/4.2.2 
R

## code ALS BA4
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano combat_ALS_BA4_refined.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#!/bin/bash  
#SBATCH --account=def-sfarhan
#SBATCH --time=02-00:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=200g          # memory per cor
#SBATCH --job-name=combat_ALS_BA4_refined
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load r/4.2.2 

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/combat_ALS_BA4_refined.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano combat_ALS_BA4_refined.R

    ## load libraries
    library(sva)
    library(Seurat)
    library(SeuratDisk, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2)
    library(Matrix)
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    # Load the Seurat object
    seu <- readRDS("/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA4_merged.rds")
    unique(seu@meta.data$full_label)

    ## See how many cells we have for each cell type
    table(seu@meta.data$full_label)
    celltype <- unique(seu@meta.data$full_label)

    for(seu_celltype in unique(celltype)){
        print(seu_celltype)
        ## select only cell type of interest
        xx <- unique(seu@meta.data$full_label)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu) <- "full_label"
        seu_lim=subset(seu,idents=xx)
        unique(seu_lim@meta.data$full_label)
        dim(seu_lim)

        ## we need to remove samples that only have one cell.
        num_cells <- data.frame(table(seu_lim@meta.data$Donor))
        num_cells <- subset(num_cells, Freq > 1)
        keep_donor <- unique(num_cells$Var1)

        Idents(seu_lim) <- "Donor"
        seu_lim=subset(seu_lim,idents=keep_donor)
        unique(seu_lim@meta.data$Donor)
        table(seu_lim@meta.data$Donor)
        dim(seu_lim)
        
        m = seu_lim@assays$RNA$counts
        adjusted <- ComBat_seq(as.matrix(m), batch=as.numeric(as.factor(seu_lim@meta.data$Donor)))

        ## Add feature and cell names to combat-corrected matrix
        #test <- data.frame(seu_lim@assays$RNA@counts)
        rownames(adjusted) <- rownames(seu_lim)
        colnames(adjusted) <- rownames(seu_lim@meta.data)

        ## Create Seurat assay with combat-generated batch corrected matrix
        assay.v5 <- CreateAssay5Object(counts = as.matrix(adjusted))
        pbmc3k_slim <- CreateSeuratObject(assay.v5)
        pbmc3k_slim <- NormalizeData(pbmc3k_slim)
        pbmc3k_slim <- ScaleData(pbmc3k_slim, verbose = FALSE)
        pbmc3k_slim <- FindVariableFeatures(pbmc3k_slim)
        pbmc3k_slim <- RunPCA(pbmc3k_slim, npcs = 25, verbose = FALSE)
        pbmc3k_slim <- RunUMAP(pbmc3k_slim, dims = 1:25, n.neighbors =45)
        pbmc3k_slim <- RunTSNE(pbmc3k_slim, dims = 1:25)
        dim(pbmc3k_slim)

        ## Add meta data
        df <-data.frame(seu_lim@meta.data)
        pbmc3k_slim <- AddMetaData(pbmc3k_slim, df)

        ## Print UMAP
        p1 <- DimPlot(pbmc3k_slim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "none",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        
        p2 <- DimPlot(pbmc3k_slim, reduction = "umap", group.by = "Group", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "bottom",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

        ## print plots for un-integrated objects
        xx <- unique(seu@meta.data$full_label)
        xx <- xx[xx %in% c(seu_celltype)]

        Idents(seu) <- "full_label"
        seu_lim=subset(seu,idents=xx)
        unique(seu_lim@meta.data$full_label)

        ## Create Seurat assay with out combat correction
        seu_lim <- NormalizeData(seu_lim)
        seu_lim <- ScaleData(seu_lim, verbose = FALSE)
        seu_lim <- FindVariableFeatures(seu_lim)
        seu_lim <- RunPCA(seu_lim, npcs = 25, verbose = FALSE)
        seu_lim <- RunUMAP(seu_lim, dims = 1:25, n.neighbors =45)
        seu_lim <- RunTSNE(seu_lim, dims = 1:25)

        ## Print UMAP
        p3 <- DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "none",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

        p4 <- DimPlot(seu_lim, reduction = "umap", group.by = "Group", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "bottom",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

        ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Combat_UMAP_ALS_BA4_',seu_celltype,'.pdf'), height = 10, width = 10)

        ## Save combat Seurat object
        assay.v5 <- CreateAssayObject(counts = as.matrix(adjusted))
        pbmc3k_slim <- CreateSeuratObject(assay.v5)
        
        #pbmc3k_slim[["RNA"]] <- as(object = pbmc3k_slim[["RNA"]], Class = "Assay")

        df <-data.frame(seu_lim@meta.data)
        pbmc3k_slim <- AddMetaData(pbmc3k_slim, df)
        
        SaveH5Seurat(pbmc3k_slim, 
                filename = paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_ALS_BA4_",seu_celltype,"_int.h5Seurat"),
                overwrite = T,
                verbose = T
                )

        Convert(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_ALS_BA4_",seu_celltype,"_int.h5Seurat"), dest = "h5ad", assay = "RNA")
    }
##

## code ALS BA9
    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano combat_ALS_BA9_refined.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#!/bin/bash  
#SBATCH --account=def-sfarhan
#SBATCH --time=02-00:00           # time (DD-HH:MM)
#SBATCH --cpus-per-task=1
#SBATCH --mem=200g          # memory per cor
#SBATCH --job-name=combat_ALS_BA9_refined
#SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
#SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

module load StdEnv/2020 
module load r/4.2.2 
R

cd /home/fiorini9/scratch/machine_learning_ALS/scripts

Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/combat_ALS_BA9_refined.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano combat_ALS_BA9_refined.R

    ## load libraries
    library(sva)
    library(Seurat)
    library(SeuratDisk, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2)
    library(Matrix)
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    # Load the Seurat object
    seu <- readRDS("/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA9_merged.rds")
    unique(seu@meta.data$full_label)

    ## See how many cells we have for each cell type
    table(seu@meta.data$full_label)
    celltype <- unique(seu@meta.data$full_label)

    seu_celltype <- "Astro"
    for(seu_celltype in unique(celltype)){
        print(seu_celltype)
        ## select only cell type of interest
        xx <- unique(seu@meta.data$full_label)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu) <- "full_label"
        seu_lim=subset(seu,idents=xx)
        unique(seu_lim@meta.data$full_label)
        dim(seu_lim)

        ## we need to remove samples that only have one cell.
        num_cells <- data.frame(table(seu_lim@meta.data$Donor))
        num_cells <- subset(num_cells, Freq > 1)
        keep_donor <- unique(num_cells$Var1)

        Idents(seu_lim) <- "Donor"
        seu_lim=subset(seu_lim,idents=keep_donor)
        unique(seu_lim@meta.data$Donor)
        table(seu_lim@meta.data$Donor)
        dim(seu_lim)
        
        m = seu_lim@assays$RNA$counts
        adjusted <- ComBat_seq(as.matrix(m), batch=as.numeric(as.factor(seu_lim@meta.data$Donor)))

        ## Add feature and cell names to combat-corrected matrix
        #test <- data.frame(seu_lim@assays$RNA@counts)
        rownames(adjusted) <- rownames(seu_lim)
        colnames(adjusted) <- rownames(seu_lim@meta.data)

        ## Create Seurat assay with combat-generated batch corrected matrix
        assay.v5 <- CreateAssay5Object(counts = as.matrix(adjusted))
        pbmc3k_slim <- CreateSeuratObject(assay.v5)
        pbmc3k_slim <- NormalizeData(pbmc3k_slim)
        pbmc3k_slim <- ScaleData(pbmc3k_slim, verbose = FALSE)
        pbmc3k_slim <- FindVariableFeatures(pbmc3k_slim)
        pbmc3k_slim <- RunPCA(pbmc3k_slim, npcs = 25, verbose = FALSE)
        pbmc3k_slim <- RunUMAP(pbmc3k_slim, dims = 1:25, n.neighbors =45)
        pbmc3k_slim <- RunTSNE(pbmc3k_slim, dims = 1:25)
        dim(pbmc3k_slim)

        ## Add meta data
        df <-data.frame(seu_lim@meta.data)
        pbmc3k_slim <- AddMetaData(pbmc3k_slim, df)

        ## Print UMAP
        p1 <- DimPlot(pbmc3k_slim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "none",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        
        p2 <- DimPlot(pbmc3k_slim, reduction = "umap", group.by = "Group", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "bottom",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

        ## print plots for un-integrated objects
        xx <- unique(seu@meta.data$full_label)
        xx <- xx[xx %in% c(seu_celltype)]

        Idents(seu) <- "full_label"
        seu_lim=subset(seu,idents=xx)
        unique(seu_lim@meta.data$full_label)

        ## Create Seurat assay with out combat correction
        seu_lim <- NormalizeData(seu_lim)
        seu_lim <- ScaleData(seu_lim, verbose = FALSE)
        seu_lim <- FindVariableFeatures(seu_lim)
        seu_lim <- RunPCA(seu_lim, npcs = 25, verbose = FALSE)
        seu_lim <- RunUMAP(seu_lim, dims = 1:25, n.neighbors =45)
        seu_lim <- RunTSNE(seu_lim, dims = 1:25)

        ## Print UMAP
        p3 <- DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "none",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

        p4 <- DimPlot(seu_lim, reduction = "umap", group.by = "Group", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "bottom",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

        ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Combat_UMAP_ALS_BA9_',seu_celltype,'.pdf'), height = 10, width = 10)

        ## Save combat Seurat object
        assay.v5 <- CreateAssayObject(counts = as.matrix(adjusted))
        pbmc3k_slim <- CreateSeuratObject(assay.v5)
        
        #pbmc3k_slim[["RNA"]] <- as(object = pbmc3k_slim[["RNA"]], Class = "Assay")

        df <-data.frame(seu_lim@meta.data)
        pbmc3k_slim <- AddMetaData(pbmc3k_slim, df)
        
        SaveH5Seurat(pbmc3k_slim, 
                filename = paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_ALS_BA9_",seu_celltype,"_int.h5Seurat"),
                overwrite = T,
                verbose = T
                )

        Convert(paste0("/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_ALS_BA9_",seu_celltype,"_int.h5Seurat"), dest = "h5ad", assay = "RNA")
    }
##
