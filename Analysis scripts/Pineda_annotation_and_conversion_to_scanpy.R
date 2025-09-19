salloc -A def-sfarhan --time=0-8 -c 1 --mem=200g

module load StdEnv/2020 
module load r/4.2.2 
R

/home/fiorini9/scratch/machine_learning_ALS/base_objects

################################################################################################################
################################################################################################################
################################################################################################################ Pineda ALSB4 annotation and conversion to AD
## code
    ## libraries
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(tidyr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(stringr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(ggrepel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(SeuratData, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(SeuratDisk, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    ## new Seurat
    seu_pineda <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4/step6/objs6/seu_step6.rds')
    DefaultAssay(seu_pineda)

    str(seu_pineda@meta.data)
    unique(seu_pineda@meta.data$CellType)

    ## base UMAP plot with simple cell types
    p1 <- DimPlot(seu_pineda, reduction = "umap", group.by = "CellType", label = T, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf',sep=""))

    ## base UMAP plot with disease status
    p1 <- DimPlot(seu_pineda, reduction = "umap", group.by = "Donor", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf',sep=""))

    ## base UMAP plot with disease status
    p1 <- DimPlot(seu_pineda, reduction = "umap", group.by = "Group", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf',sep=""))
    
    ## base UMAP plot with subtypes
    unique(seu_pineda@meta.data$full_label)

    p1 <- DimPlot(seu_pineda, reduction = "umap", group.by = "full_label", label = F, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") +
    ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf',sep=""))

    ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/base_UMAP_scrnabox_ALS_BA4.pdf',sep=""))
    
    ## convert to Ann data
    DefaultAssay(seu_pineda)
    DefaultAssay(seu_pineda) <- "RNA"
    DefaultAssay(seu_pineda)

    ## fix the meta data columns
    df <- data.frame(seu_pineda@meta.data)
    colnames(df)
    df$Cell_Type <- df$predicted.id
    df$Sample_ID <- df$Sample_ID
    # Find rows where 'Disease' column contains exact match "ALS"
    rows_with_ALS <- grep(disease, df$Sample_ID)
    df$Disease_Status <- "ctrl"
    df$Disease_Status[grep(disease, df$Sample_ID) ] <- disease

    table(df$Sample_ID, df$Disease_Status)


    seu_pineda <- AddMetaData(seu_pineda,df)
    ncol(seu_pineda) #107479

    p1 <- DimPlot(seu_pineda, reduction = "umap", group.by = "Cell_Type", label = T, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/Final_UMAP_scrnabox_ALS_BA4.pdf',sep=""))

    assays(seu_pineda)

    SaveH5Seurat(seu_pineda, filename = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B4.h5Seurat")
    Convert("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B4.h5Seurat", dest = "h5ad")

    max(seu_kam_lim@assays$RNA@data)

    #Fix raw
    seu_pineda <-LoadH5Seurat("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B4.h5Seurat")  

    new_seurat_object <- CreateSeuratObject(counts = seu_pineda@assays$RNA@counts)  
    new_seurat_object@meta.data <- seu_pineda@meta.data

    rm Pineda_ALS_B4.h5Seurat
    rm Pineda_ALS_B4.h5ad
    
    SaveH5Seurat(new_seurat_object, filename = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B4.h5Seurat")
    Convert("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B4.h5Seurat", dest = "h5ad")

##


################################################################################################################
################################################################################################################
################################################################################################################ Pineda scrnabox_ALS_BA9 annotation and conversion to AD
## code
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    scrnabox_ALS_BA9


    cd /home/fiorini9/scratch/machine_learning/scripts

    nano scrnabox_ALS_BA9_annot.sh

    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=00-05:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=200g          # memory per cor
    #SBATCH --job-name=scrnabox_ALS_BA9_annot
    #SBATCH --error=/home/fiorini9/scratch/machine_learning/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning/scripts

    Rscript /home/fiorini9/scratch/machine_learning/scripts/scrnabox_ALS_BA9_annot.R

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
    nano scrnabox_ALS_BA9_annot.R

    disease = "ALS"


    ## libraries
    library(Seurat)
    library(ggplot2)
    library(Seurat)
    library(ggplot2)
    library(ggpubr, lib="/home/fiorini9/scratch/mjff/MF/practice_ensemblux/R")
    library(ComplexUpset, lib="/home/fiorini9/scratch/mjff/MF/practice_ensemblux/R")
    library(tidyr)
    library(stringr)
    library(dplyr)
    library(ggrepel)
    library(phateR)
    library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Seurat)
    library(SeuratData, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(SeuratDisk, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    ## new Seurat
    seu_pineda <- readRDS('/home/fiorini9/scratch/machine_learning/scRNAbox_pineda/scrnabox_ALS_BA4/step6/objs6/seu_step6.rds')
    DefaultAssay(seu_pineda)

    str(seu_pineda@meta.data)
    unique(seu_pineda@meta.data$CellType)

    ## base UMAP plot
    p1 <- DimPlot(seu_pineda, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/base_UMAP_scrnabox_ALS_BA9.pdf',sep=""))

    ## reference based annotation with KAM old object

    ## set user defined clustering resolution
    seu_pineda <- SetIdent(seu_pineda, value = 'RNA_snn_res.0.8')


    ###################################
    ## new Kamath (with peri) as reference
    ###################################
    reference0 <- LoadH5Seurat("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Kamath.h5Seurat")  
    DefaultAssay(reference0) <- "RNA" ## new code
    unique(reference0@meta.data$Disease_Status)
    unique(reference0@meta.data$Cell_Type)
    
    ## downsample the Kamath reference
    ncol(reference0)
    # Step 1: Get the cell names or indices
    cell_names <- colnames(reference0)

    # Step 2: Generate random indices
    num_random_cells <- 150000
    random_indices <- sample(length(cell_names), num_random_cells)

    # Step 3: Subset the Seurat object based on random indices
    reference0_subset <- reference0[, random_indices]

    ## remove whole object
    rm(reference0)

    # perform standard preprocessing on reference object
    reference0_subset<- NormalizeData(reference0_subset)
    reference0_subset <- FindVariableFeatures(reference0_subset)
    reference0_subset<- ScaleData(reference0_subset)
    reference0_subset <- RunPCA(object = reference0_subset, assay = "RNA", npcs = 20)

    ## find transfer anchors between reference and query Seurat objects
    transfer.anchors <- FindTransferAnchors(reference = reference0_subset, query = seu_pineda, dims = 1:20, reference.reduction = "pca")

    ## add reference-based annotations to the qeury object
    eval(parse(text = paste('predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0_subset$',"Cell_Type" ,',dims = 1:',20,')', sep='')))
    seu_pineda <- AddMetaData(object = seu_pineda, metadata = predictions)

    # Add metadata column for reference object
    seu_pineda$temp_temp_2 <- seu_pineda@meta.data$predicted.id
    name_meta <- names(seu_pineda@meta.data) 
    length <- length(name_meta)
    name_meta[length] <- paste("Kamath", "_predictions", sep = "")
    names(seu_pineda@meta.data) <- name_meta

    ## Print a umap projection showing the predicted cell types on the query object 
    reference0_subset <- RunUMAP(reference0_subset, dims = 1:20, reduction = "pca", return.model = TRUE)
    seu_pineda <- MapQuery(anchorset = transfer.anchors, reference = reference0_subset, query = seu_pineda,
        refdata = list(celltype = "Cell_Type"), reference.reduction = "pca", reduction.model = "umap")
    p1 <- DimPlot(reference0_subset, reduction = "umap", group.by = "Cell_Type", label = TRUE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() + ggtitle("Reference annotations")
    p2 <- DimPlot(seu_pineda, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE, raster=FALSE) + NoLegend() + ggtitle("Query transferred labels")
    p1 + p2
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/base_UMAP_ref_scrnabox_ALS_BA9.pdf',sep=""), height = 10, width =10)  
    
    ## convert to Ann data
    DefaultAssay(seu_pineda)
    DefaultAssay(seu_pineda) <- "RNA"
    DefaultAssay(seu_pineda)

    ## fix the meta data columns
    df <- data.frame(seu_pineda@meta.data)
    colnames(df)
    df$Cell_Type <- df$predicted.id
    df$Sample_ID <- df$Sample_ID
    # Find rows where 'Disease' column contains exact match "ALS"
    rows_with_ALS <- grep(disease, df$Sample_ID)
    df$Disease_Status <- "ctrl"
    df$Disease_Status[grep(disease, df$Sample_ID) ] <- disease

    table(df$Sample_ID, df$Disease_Status)


    seu_pineda <- AddMetaData(seu_pineda,df)
    ncol(seu_pineda) #107479

    p1 <- DimPlot(seu_pineda, reduction = "umap", group.by = "Cell_Type", label = T, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.text = element_text( size =10),
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste('/home/fiorini9/scratch/machine_learning/temp_figures/Final_UMAP_scrnabox_ALS_BA9.pdf',sep=""))

    SaveH5Seurat(seu_pineda, filename = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B9.h5Seurat")
    Convert("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B9.h5Seurat", dest = "h5ad")

    max(seu_kam_lim@assays$RNA@data)

    #Fix raw
    seu_pineda <-LoadH5Seurat("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B9.h5Seurat")  

    new_seurat_object <- CreateSeuratObject(counts = seu_pineda@assays$RNA@counts)  
    new_seurat_object@meta.data <- seu_pineda@meta.data

    rm Pineda_ALS_B4.h5Seurat
    rm Pineda_ALS_B4.h5ad
    
    SaveH5Seurat(new_seurat_object, filename = "/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B9.h5Seurat")
    Convert("/home/fiorini9/scratch/machine_learning/base_seurat_objects/Pineda_ALS_B9.h5Seurat", dest = "h5ad")
##