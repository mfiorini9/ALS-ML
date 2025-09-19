salloc -A def-sfarhan --time=0-8 -c 1 --mem=40g

module load StdEnv/2023
module load r/4.4.0
R

## code

    ###########################
    ## load libraries
    ###########################
    ##
    library(harmony)
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
    library(scCustomize, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(sva)
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
        
    ##

    ###########################
    ## Test
    ###########################
    ## read in Pineda
    seu_pineda <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_all_celltypes_lim_narval.rds')
    dim(seu_pineda)
    n_cells <- 100000 
    set.seed(123)  # for reproducibility
    cells_to_keep <- sample(colnames(seu_pineda), n_cells)
    seu_pineda <- subset(seu_pineda, cells = cells_to_keep)
    dim(seu_pineda)
    seu_pineda@meta.data$Dataset <- "Pineda"
    seu_pineda@meta.data$Sex[seu_pineda@meta.data$Sex == "FALSE"] <- "F"
    seu_pineda@meta.data$Sex[seu_pineda@meta.data$Sex == "M"] <- "M"

    unique(seu_pineda@meta.data$CellType)
    unique(seu_pineda@meta.data$Region)
    unique(seu_pineda@meta.data$Dataset)
    unique(seu_pineda@meta.data$Group)
    unique(seu_pineda@meta.data$Sex)
    
    ## Read in Li
    seu_li <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_BA4_BA9_all_celltypes_lim_narval.rds')
    dim(seu_li)
    n_cells <- 100000 #100000
    set.seed(123)  # for reproducibility
    cells_to_keep <- sample(colnames(seu_li), n_cells)
    seu_li <- subset(seu_li, cells = cells_to_keep)
    dim(seu_li)
    seu_li@meta.data$Dataset <- "Li"

    seu_li@meta.data$Region[seu_li@meta.data$Region == "motor cortex"] <- "BA4"
    seu_li@meta.data$Region[seu_li@meta.data$Region == "medial frontal cortex"] <- "BA9"
    seu_li@meta.data$Group[seu_li@meta.data$Group == "C9-ALS"] <- "C9ALS"
    seu_li@meta.data$Group[seu_li@meta.data$Group == "C9-FTD"] <- "C9FTD"
    seu_li@meta.data$Group[seu_li@meta.data$Group == "Control"] <- "PN"
    seu_li@meta.data$Sex <- seu_li@meta.data$sex
    
    unique(seu_li@meta.data$CellType)
    unique(seu_li@meta.data$Region)
    unique(seu_li@meta.data$Dataset)
    unique(seu_li@meta.data$Group)
    unique(seu_li@meta.data$Sex)

    ## Read in Limone
    seu_limone <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_BA4_all_celltypes_lim_narval.rds')
    dim(seu_limone)
    seu_limone@meta.data$Dataset <- "Limone"

    seu_limone@meta.data$Group[seu_limone@meta.data$Group == "ALS"] <- "SALS"
    seu_limone@meta.data$Group[seu_limone@meta.data$Group == "Control"] <- "PN"
    seu_limone@meta.data$Sex[seu_limone@meta.data$Sample_ID == "ALS_MotorCortex.012218_FC21"] <- "F"
    seu_limone@meta.data$Sex[seu_limone@meta.data$Sample_ID == "ALS_MotorCortex.121417_FC19"] <- "M"
    seu_limone@meta.data$Sex[seu_limone@meta.data$Sample_ID == "Control_MotorCortex.012218_FC13"] <- "M"
    seu_limone@meta.data$Sex[seu_limone@meta.data$Sample_ID == "ALS_MotorCortex.012318_FC17"] <- "F"
    seu_limone@meta.data$Sex[seu_limone@meta.data$Sample_ID == "Control_MotorCortex.012218_FC12"] <- "M"
    seu_limone@meta.data$Sex[seu_limone@meta.data$Sample_ID == "Control_MotorCortex.121417_FC11"] <- "F"

    unique(seu_limone@meta.data$CellType)
    unique(seu_limone@meta.data$Region)
    unique(seu_limone@meta.data$Dataset)
    unique(seu_limone@meta.data$Group)
    unique(seu_limone@meta.data$Sex)

    ## get common genes
    keep_genes <- Reduce(intersect, list(
    rownames(seu_pineda),
    rownames(seu_li),
    rownames(seu_limone)
    ))

    ## Subset to only include common genes
    seu_pineda <- subset(seu_pineda, features = keep_genes)
    seu_li <- subset(seu_li, features = keep_genes)
    seu_limone <- subset(seu_limone, features = keep_genes)
    dim(seu_pineda)
    dim(seu_li)
    dim(seu_limone)

    
    
    ############ LOOP STARTS HERE
    start by subsetting to cell type on interest. 
    
    
    
    
    
    
    
    ## Subset relevant metadata columns
    meta_cols <- c("CellType", "Sample_ID", "Group", "Region", "Donor", 
               "Dataset", "nCount_RNA", "nFeature_RNA", "Sex")

    seu_pineda@meta.data <- seu_pineda@meta.data[, meta_cols]
    seu_li@meta.data <- seu_li@meta.data[, meta_cols]
    seu_limone@meta.data <- seu_limone@meta.data[, meta_cols]

    ## create joint metadata
    meta_data_total <- rbind(data.frame(seu_pineda@meta.data), data.frame(seu_li@meta.data), data.frame(seu_limone@meta.data))
    nrow(meta_data_total)
    unique(meta_data_total$Sex)
    write.csv(meta_data_total, '/home/fiorini9/scratch/machine_learning_ALS/model_outs/integrated_metadata.csv')
    
    # 1. Put Seurat objects in a list
    seu_list <- list(Pineda = seu_pineda, Li = seu_li, Limone = seu_limone)

    # 2. Normalize + find variable features for each dataset
    seu_list <- lapply(seu_list, function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    return(x)
    })

    # 3. Find integration anchors
    anchors <- FindIntegrationAnchors(object.list = seu_list, dims = 1:30)

    # 4. Integrate data
    seu_combined <- IntegrateData(anchorset = anchors, dims = 1:30)

    saveRDS(seu_combined, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/integrated_pineda_li_limone_BA4_BA9_lim.RDS')

    #########
    seu_combined <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/integrated_pineda_li_limone_BA4_BA9_lim.RDS')

    meta_data_total <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/integrated_metadata.csv')
    rownames(meta_data_total) <- meta_data_total$X
    
    ## Add metadata so that we get sample Sex.
    seu_combined <- AddMetaData(seu_combined, meta_data_total) 
    str(seu_combined@meta.data)

    
    # 5. Standard workflow on integrated data
    DefaultAssay(seu_combined) <- "integrated"
    seu_combined <- ScaleData(seu_combined, verbose = FALSE)
    seu_combined <- RunPCA(seu_combined, npcs = 50, verbose = FALSE)
    seu_combined <- RunUMAP(seu_combined, dims = 1:40)
    #DimPlot(seu_combined, reduction = "umap", group.by = "Dataset")

    #################### 
    # UMAP plots
    ####################
    colors <- c(
            "L3_L5" = "#f3c300",
            "L2_L3" = "#f38400",
            "L4_L6" = "#a1caf1",
            "L4_L5" = "#be0032",
            "L5_L6" = "#c2b280",
            "L5" = "#008856",
            "L6" = "#2b3d26",
            "PV" = "#e25822",
            "5HT3aR" = "#e68fac",
            "Rosehip" = "#0067a5",
            "SOM" = "#f99379",
            "Oligo" = "#604e97",
            "Astro" = "#875692",
            "OPC" = "#f6a600",
            "Micro" = "#b3446c",
            "T_Cell" = "#dcd300",
            "Mural" = "#882d17",
            "Endo" = "#8db600",
            "Fibro" = "#654522")
    
    str(seu_combined@meta.data)

    ## By Celltype
    DimPlot(seu_combined, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE, pt.size = 0.1) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values =  colors)

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 15, width = 15)

    ## By Region
    spectral_base <- brewer.pal(12, "Paired")

    DimPlot(seu_combined, reduction = "umap", group.by = "Region", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE, pt.size = 0.1) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values =  c( "#E31A1C", "#1B9E77"))

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 15, width = 15)

    ## By Sex
    spectral_base <- brewer.pal(12, "Paired")

    DimPlot(seu_combined, reduction = "umap", group.by = "Sex", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE, pt.size = 0.1) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values =  c( "#FB9A99", "#1F78B4"))

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp1.pdf', height = 15, width = 15)

    ## By Dataset
    DimPlot(seu_combined, reduction = "umap", group.by = "Dataset", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE, pt.size = 0.1) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c("Pineda" = "#FF7F00",
                                "Li"  = "#FDBF6F" ,  
                                "Limone" = "#B15928")) 


    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 15, width = 15)
    
    ## By Group
    DimPlot(seu_combined, reduction = "umap", group.by = "Group", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE, pt.size = 0.1) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c("SALS" = "orange",
                                "C9ALS"  = "red",
                                "SFTLD" = "blue",
                                "C9FTLD"  = "purple",
                                "PN" = "#339966")) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 15, width = 15)
    
    ## By Donor
    
    spectral_base <- brewer.pal(12, "Paired")

    # Interpolate to 96 colors in Lab space for smooth perceptual differences
    spectral_96 <- colorRampPalette(spectral_base, space = "Lab")(96)

    DimPlot(seu_combined, reduction = "umap", group.by = "Donor", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE, pt.size = 0.1) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = rev(spectral_96)) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 15, width = 15)
    
##