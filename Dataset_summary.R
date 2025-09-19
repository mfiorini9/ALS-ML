/home/fiorini9/nearline/ctb-grouleau/fiorini9/segpy_AALS



salloc -A def-sfarhan --time=0-8 -c 1 --mem=40g

module load StdEnv/2023
module load r/4.4.0
R

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Print dataframe from Python objects
## code
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/home/fiorini9/projects/def-sfarhan/fiorini9/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

    python  

    import pandas as pd
    import scanpy as sc


    ###################################
    # Pineda
    ###################################
    ## code
        cell_type_list = ['L3_L5', 'L2_L3', 'L4_L6', 'L4_L5', 'L5_L6', 'L6', 'PV', 
                        '5HT3aR', 'Rosehip', 'SOM', 'Oligo', 'Astro', 'OPC', 
                        'Micro', 'Mural', 'Endo', 'Fibro', 'L5']

        len(cell_type_list)
        out_file = "/home/fiorini9/scratch/machine_learning_ALS/model_outs/Pineda_combined_metadata.csv"
        first = True  # track if it's the first write

        for cell_type in cell_type_list:
            for region in ["BA4", "BA9"]:
                print(cell_type)
                path = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/combat_{region}_{cell_type}_int.h5ad"
                adata = sc.read_h5ad(path, backed="r")   # backed mode → doesn’t fully load into memory
                
                # convert obs to DataFrame
                df = adata.obs.copy()
                df["cell_type"] = cell_type
                df["region"] = region
                
                # write to CSV incrementally
                if first:
                    df.to_csv(out_file, mode="w", index=True, header=True)
                    first = False
                else:
                    df.to_csv(out_file, mode="a", index=True, header=False)
                
                # cleanup
                del adata, df

    ##

    ###################################
    # Li
    ###################################
    ## code 
        cell_type_list = ['L3_L5', 'L2_L3', 'L4_L6', 'L4_L5', 'L5_L6', 'L6', 'PV', 
                        '5HT3aR', 'Rosehip', 'SOM', 'Oligo', 'Astro', 'OPC', 
                        'Micro', 'Mural', 'Endo', 'Fibro', 'L5']

        len(cell_type_list)
        out_file = "/home/fiorini9/scratch/machine_learning_ALS/model_outs/Li_combined_metadata.csv"
        first = True  # track if it's the first write

        for cell_type in cell_type_list:
            for region in ["BA4", "BA9"]:
                print(cell_type)
                path = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_combat_{region}_{cell_type}_int.h5ad"
                adata = sc.read_h5ad(path, backed="r")   # backed mode → doesn’t fully load into memory
                
                # convert obs to DataFrame
                df = adata.obs.copy()
                df["cell_type"] = cell_type
                df["region"] = region
                
                # write to CSV incrementally
                if first:
                    df.to_csv(out_file, mode="w", index=True, header=True)
                    first = False
                else:
                    df.to_csv(out_file, mode="a", index=True, header=False)
                
                # cleanup
                del adata, df
    ##

    ###################################
    # Limone
    ###################################
    ## Limone
        cell_type_list = ['L3_L5', 'L2_L3', 'L4_L6', 'L4_L5', 'L5_L6', 'L6', 'PV', 
                        '5HT3aR', 'Rosehip', 'SOM', 'Oligo', 'Astro', 'OPC', 
                        'Micro', 'Mural', 'Endo', 'Fibro', 'L5']

        len(cell_type_list)
        out_file = "/home/fiorini9/scratch/machine_learning_ALS/model_outs/Limone_combined_metadata.csv"
        first = True  # track if it's the first write

        for cell_type in cell_type_list:
            for region in ["BA4"]:
                print(cell_type)
                path = f"/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_combat_{region}_{cell_type}_int.h5ad"
                adata = sc.read_h5ad(path, backed="r")   # backed mode → doesn’t fully load into memory
                
                # convert obs to DataFrame
                df = adata.obs.copy()
                df["cell_type"] = cell_type
                df["region"] = region
                
                # write to CSV incrementally
                if first:
                    df.to_csv(out_file, mode="w", index=True, header=True)
                    first = False
                else:
                    df.to_csv(out_file, mode="a", index=True, header=False)
                
                # cleanup
                del adata, df
    ##
##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ UMAPS

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
    
    # 1. Put your Seurat objects in a list
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


################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Combined datasets: Figure 1
## code 
    module load StdEnv/2023
    module load r/4.4.0
    R

    ## Load libraries
    library(MAST)
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
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)
    library(colorspace)

    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)

    ############################################
    ## Read in metadata files
    ############################################
    ## Pineda
    pineda_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Pineda_combined_metadata.csv')
    pineda_meta <- pineda_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, CellType, Region)
    pineda_meta$dataset <- "Pineda et al."
    colnames(pineda_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    
    ## Li
    li_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Li_combined_metadata.csv')
    li_meta$Condition[li_meta$Group == "C9-ALS"] <- "ALS"
    li_meta$Condition[li_meta$Group == "C9-FTD"] <- "FTLD"
    li_meta$Condition[li_meta$Group == "Control"] <- "Control"
    li_meta <- li_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, sex, Condition, Group, cell_type, region)
    li_meta$dataset <- "Li et al."
    colnames(li_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")
    
    ## Limone
    limone_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Limone_combined_metadata.csv')
    colnames(limone_meta) 
    unique(limone_meta$Sample_ID)
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012218_FC21"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.121417_FC19"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC13"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012318_FC17"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC12"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.121417_FC11"] <- "F"
    limone_meta$Condition <- limone_meta$Group
    limone_meta$Group[limone_meta$Group == "ALS"] <- "SALS"
    limone_meta$Group[limone_meta$Group == "Control"] <- "Control"
    limone_meta <- limone_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, cell_type, region)
    limone_meta$dataset <- "Limone et al."
    colnames(limone_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    ## Bind and standardize
    total_meta <- rbind(pineda_meta, li_meta, limone_meta)
    total_meta$Sex[total_meta$Sex == "FALSE"] <- "F"
    total_meta$Condition[total_meta$Condition == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "C9-ALS"] <- "C9ALS"
    total_meta$Group[total_meta$Group == "C9-FTD"] <- "C9FTLD"
    
    unique(total_meta$Donor)
    unique(total_meta$Sex)
    unique(total_meta$Condition)
    unique(total_meta$Group)
    unique(total_meta$Cell_Type)
    unique(total_meta$Region)
    unique(total_meta$Dataset)

    ############################################
    ## Set cell types levels
    ############################################
    cell_type_levels <- c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC')
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

    ############################################
    ## Plot number of cells per cell type
    ############################################
    n_cell <- data.frame(table(total_meta$Cell_Type))
    colnames(n_cell) <- c("cell_type", "Freq")
    n_cell$cell_type <- factor(n_cell$cell_type, levels =  rev(cell_type_levels))

    ## Plot
    n_cell_plot <-  ggplot(n_cell, aes(x = Freq, y = cell_type, fill = cell_type)) + 
    theme_classic() + 
    geom_bar(stat = 'identity') +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size =8, colour = "black"),
        axis.text.y = element_text(size =8, colour = "black"),
        axis.title = element_blank(),
    ) +
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = "white") +
    scale_x_continuous(breaks = c(0, 50000, 100000, 150000), labels = c("0", "5e4", "10e4", "15e4") )
    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Dataset prop
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    unique(total_meta$Dataset)

    spectral_base <- brewer.pal(12, "Paired")

    dataset_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Dataset )) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("Pineda et al." = "#FF7F00",
                                "Li et al."  = "#FDBF6F" ,  
                                "Limone et al." = "#B15928")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Brain region prop
    ############################################
    spectral_base <- brewer.pal(11, "Paired")

    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Region <- factor(total_meta$Region , levels =  rev(c("BA4", "BA9")))
    unique(total_meta$Region)

    region_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Region )) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("BA9" = "#1B9E77",
                                "BA4"  = "#E31A1C")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Diagnosis prop
    ############################################
    #total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    #total_meta$Condition <- factor(total_meta$Condition , levels =  rev(c("ALS", "FTLD", "Control")))
    #unique(total_meta$Condition)

    #condition_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Condition )) + 
    #geom_bar(position = "fill") +   # normalize counts to proportions
    #coord_flip() +                  # flip so proportions are on x-axis
    #theme_classic() + 
    #theme(
    #    legend.position = "none",
    #    panel.grid = element_blank(),
    #    axis.text.x = element_text(size = 8, colour = "black"),
    #    axis.text.y = element_blank(),
    #    axis.ticks.y = element_blank(),
    #    axis.title = element_blank()
    #) +
    #scale_fill_manual(values = c("ALS" = "darkred",
    #                            "FTLD"  = "darkblue",
    #                            "Control"  = "#339966")) +
    #scale_colour_manual(values = "white") +
    #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Disease group prop
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Group <- factor(total_meta$Group , levels =  rev(c("SALS", "C9ALS", "SFTLD", "C9FTLD", "Control")))
    unique(total_meta$Group)

    group_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Group)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("SALS" = "orange",
                                "C9ALS"  = "red",
                                "SFTLD" = "blue",
                                "C9FTLD"  = "purple",
                                "Control" = "#339966")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Sex prop
    ############################################
    spectral_base <- brewer.pal(12, "Paired")
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Sex <- factor(total_meta$Sex , levels =  rev(c("M", "F")))
    unique(total_meta$Sex)

    sex_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Sex)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("M" = "#1F78B4",
                                "F"  = "#FB9A99")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)


    ############################################
    ## Donor
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    unique(total_meta$Donor)
    
    spectral_base <- brewer.pal(12, "Paired")

    # Interpolate to 96 colors in Lab space for smooth perceptual differences
    spectral_96 <- colorRampPalette(spectral_base, space = "Lab")(96)

    # Print all 96 colors
    print(spectral_96)

    donor_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Donor)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = spectral_96) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## UMI
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))

    UMI_plot <- ggplot(total_meta, aes(y = Cell_Type, x = log(nCount_RNA), fill = Cell_Type, colour = NULL)) + 
    geom_violin() +   # normalize counts to proportions
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = colors) #+
    #scale_colour_manual(values = colors) #+
    #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

    ############################################
    ## Gene
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))

    feature_plot <- ggplot(total_meta, aes(y = Cell_Type, x = log(nFeature_RNA), fill = Cell_Type, colour = NULL)) + 
    geom_violin() +   # normalize counts to proportions
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = colors) #+
    #scale_colour_manual(values = colors) #+
    #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)


    ############################################
    ## Arrange
    ############################################
    ggarrange(n_cell_plot, 
                dataset_prop_plot, 
                region_prop_plot, 
                sex_prop_plot, 
                group_prop_plot, 
                donor_prop_plot, 
                UMI_plot,
                feature_plot,
                nrow = 1, ncol =8,
                widths = c(1.2,1,1,1,1,1,1,1),
                align = "h")
    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 3.5, width = 12)


##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Pineda supplemental figure
## code 
    module load StdEnv/2023
    module load r/4.4.0
    R

    ## Load libraries
    library(MAST)
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
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)
    library(colorspace)

    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)

    ############################################
    ## Read in metadata files
    ############################################
    ## Pineda
    pineda_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Pineda_combined_metadata.csv')
    pineda_meta <- pineda_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, CellType, Region)
    pineda_meta$dataset <- "Pineda et al."
    colnames(pineda_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    
    ## Li
    li_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Li_combined_metadata.csv')
    li_meta$Condition[li_meta$Group == "C9-ALS"] <- "ALS"
    li_meta$Condition[li_meta$Group == "C9-FTD"] <- "FTLD"
    li_meta$Condition[li_meta$Group == "Control"] <- "Control"
    li_meta <- li_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, sex, Condition, Group, cell_type, region)
    li_meta$dataset <- "Li et al."
    colnames(li_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")
    
    ## Limone
    limone_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Limone_combined_metadata.csv')
    colnames(limone_meta) 
    unique(limone_meta$Sample_ID)
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012218_FC21"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.121417_FC19"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC13"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012318_FC17"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC12"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.121417_FC11"] <- "F"
    limone_meta$Condition <- limone_meta$Group
    limone_meta$Group[limone_meta$Group == "ALS"] <- "SALS"
    limone_meta$Group[limone_meta$Group == "Control"] <- "Control"
    limone_meta <- limone_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, cell_type, region)
    limone_meta$dataset <- "Limone et al."
    colnames(limone_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    ## Bind and standardize
    total_meta <- rbind(pineda_meta, li_meta, limone_meta)
    total_meta$Sex[total_meta$Sex == "FALSE"] <- "F"
    total_meta$Condition[total_meta$Condition == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "C9-ALS"] <- "C9ALS"
    total_meta$Group[total_meta$Group == "C9-FTD"] <- "C9FTLD"
    
    unique(total_meta$Donor)
    unique(total_meta$Sex)
    unique(total_meta$Condition)
    unique(total_meta$Group)
    unique(total_meta$Cell_Type)
    unique(total_meta$Region)
    unique(total_meta$Dataset)

    ############################################
    ## Subset dataset of interest
    ############################################
    total_meta <- subset(total_meta, Dataset == "Pineda et al.")
    unique(total_meta$Dataset)

    ############################################
    ## Set cell types levels
    ############################################
    cell_type_levels <- c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC')
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

    ############################################
    ## Plot number of cells per cell type
    ############################################
    n_cell <- data.frame(table(total_meta$Cell_Type))
    colnames(n_cell) <- c("cell_type", "Freq")
    n_cell$cell_type <- factor(n_cell$cell_type, levels =  rev(cell_type_levels))

    ## Plot
    n_cell_plot <-  ggplot(n_cell, aes(x = Freq, y = cell_type, fill = cell_type)) + 
    theme_classic() + 
    geom_bar(stat = 'identity') +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size =8, colour = "black"),
        axis.text.y = element_text(size =8, colour = "black"),
        axis.title = element_blank(),
    ) +
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = "white") +
    scale_x_continuous(breaks = c(0, 40000, 80000), labels = c("0", "4e4", "8e4") )
    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

    ############################################
    ## Brain region prop
    ############################################
    spectral_base <- brewer.pal(11, "Paired")

    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Region <- factor(total_meta$Region , levels =  rev(c("BA4", "BA9")))
    unique(total_meta$Region)

    region_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Region )) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("BA9" = "#1B9E77",
                                "BA4"  = "#E31A1C")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Disease group prop
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Group <- factor(total_meta$Group , levels =  rev(c("SALS", "C9ALS", "SFTLD", "C9FTLD", "Control")))
    unique(total_meta$Group)

    group_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Group)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("SALS" = "orange",
                                "C9ALS"  = "red",
                                "SFTLD" = "blue",
                                "C9FTLD"  = "purple",
                                "Control" = "#339966")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Sex prop
    ############################################
    spectral_base <- brewer.pal(12, "Paired")
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Sex <- factor(total_meta$Sex , levels =  rev(c("M", "F")))
    unique(total_meta$Sex)

    sex_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Sex)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("M" = "#1F78B4",
                                "F"  = "#FB9A99")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Donor
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    unique(total_meta$Donor)
    
    spectral_base <- brewer.pal(12, "Paired")

    number <- length(unique(total_meta$Donor))

    # Interpolate to 96 colors in Lab space for smooth perceptual differences
    spectral_96 <- colorRampPalette(spectral_base, space = "Lab")(number)

    # Print all 96 colors
    print(spectral_96)

    donor_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Donor)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = spectral_96) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## UMI
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))

    UMI_plot <- ggplot(total_meta, aes(y = Cell_Type, x = log(nCount_RNA), fill = Cell_Type, colour = NULL)) + 
    geom_violin() +   # normalize counts to proportions
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = colors) #+
    #scale_colour_manual(values = colors) #+
    #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

    ############################################
    ## Gene
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))

    feature_plot <- ggplot(total_meta, aes(y = Cell_Type, x = log(nFeature_RNA), fill = Cell_Type, colour = NULL)) + 
    geom_violin() +   # normalize counts to proportions
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = colors) #+
    #scale_colour_manual(values = colors) #+
    #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)


    ############################################
    ## Arrange
    ############################################
    ggarrange(n_cell_plot, 
                region_prop_plot, 
                sex_prop_plot, 
                group_prop_plot, 
                donor_prop_plot, 
                UMI_plot,
                feature_plot,
                nrow = 1, ncol =7,
                widths = c(1.2,1,1,1,1,1,1),
                align = "h")
    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 10)


##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Li supplemental figure
## code 
    module load StdEnv/2023
    module load r/4.4.0
    R

    ## Load libraries
    library(MAST)
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
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)
    library(colorspace)

    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)

    ############################################
    ## Read in metadata files
    ############################################
    ## Pineda
    pineda_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Pineda_combined_metadata.csv')
    pineda_meta <- pineda_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, CellType, Region)
    pineda_meta$dataset <- "Pineda et al."
    colnames(pineda_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    
    ## Li
    li_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Li_combined_metadata.csv')
    li_meta$Condition[li_meta$Group == "C9-ALS"] <- "ALS"
    li_meta$Condition[li_meta$Group == "C9-FTD"] <- "FTLD"
    li_meta$Condition[li_meta$Group == "Control"] <- "Control"
    li_meta <- li_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, sex, Condition, Group, cell_type, region)
    li_meta$dataset <- "Li et al."
    colnames(li_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")
    
    ## Limone
    limone_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Limone_combined_metadata.csv')
    colnames(limone_meta) 
    unique(limone_meta$Sample_ID)
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012218_FC21"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.121417_FC19"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC13"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012318_FC17"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC12"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.121417_FC11"] <- "F"
    limone_meta$Condition <- limone_meta$Group
    limone_meta$Group[limone_meta$Group == "ALS"] <- "SALS"
    limone_meta$Group[limone_meta$Group == "Control"] <- "Control"
    limone_meta <- limone_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, cell_type, region)
    limone_meta$dataset <- "Limone et al."
    colnames(limone_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    ## Bind and standardize
    total_meta <- rbind(pineda_meta, li_meta, limone_meta)
    total_meta$Sex[total_meta$Sex == "FALSE"] <- "F"
    total_meta$Condition[total_meta$Condition == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "C9-ALS"] <- "C9ALS"
    total_meta$Group[total_meta$Group == "C9-FTD"] <- "C9FTLD"
    
    unique(total_meta$Donor)
    unique(total_meta$Sex)
    unique(total_meta$Condition)
    unique(total_meta$Group)
    unique(total_meta$Cell_Type)
    unique(total_meta$Region)
    unique(total_meta$Dataset)

    ############################################
    ## Subset dataset of interest
    ############################################
    total_meta <- subset(total_meta, Dataset == "Li et al.")
    unique(total_meta$Dataset)

    ############################################
    ## Set cell types levels
    ############################################
    cell_type_levels <- c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC')
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

    ############################################
    ## Plot number of cells per cell type
    ############################################
    n_cell <- data.frame(table(total_meta$Cell_Type))
    colnames(n_cell) <- c("cell_type", "Freq")
    n_cell$cell_type <- factor(n_cell$cell_type, levels =  rev(cell_type_levels))

    ## Plot
    n_cell_plot <-  ggplot(n_cell, aes(x = Freq, y = cell_type, fill = cell_type)) + 
    theme_classic() + 
    geom_bar(stat = 'identity') +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size =8, colour = "black"),
        axis.text.y = element_text(size =8, colour = "black"),
        axis.title = element_blank(),
    ) +
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = "white") +
    scale_x_continuous(breaks = c(0, 20000, 40000), labels = c("0", "2e4", "4e4") )
    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

    ############################################
    ## Brain region prop
    ############################################
    spectral_base <- brewer.pal(11, "Paired")

    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Region <- factor(total_meta$Region , levels =  rev(c("BA4", "BA9")))
    unique(total_meta$Region)

    region_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Region )) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("BA9" = "#1B9E77",
                                "BA4"  = "#E31A1C")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Disease group prop
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Group <- factor(total_meta$Group , levels =  rev(c("SALS", "C9ALS", "SFTLD", "C9FTLD", "Control")))
    unique(total_meta$Group)

    group_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Group)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("SALS" = "orange",
                                "C9ALS"  = "red",
                                "SFTLD" = "blue",
                                "C9FTLD"  = "purple",
                                "Control" = "#339966")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Sex prop
    ############################################
    spectral_base <- brewer.pal(12, "Paired")
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Sex <- factor(total_meta$Sex , levels =  rev(c("M", "F")))
    unique(total_meta$Sex)

    sex_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Sex)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("M" = "#1F78B4",
                                "F"  = "#FB9A99")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Donor
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    unique(total_meta$Donor)
    
    spectral_base <- brewer.pal(12, "Paired")

    number <- length(unique(total_meta$Donor))

    # Interpolate to 96 colors in Lab space for smooth perceptual differences
    spectral_96 <- colorRampPalette(spectral_base, space = "Lab")(number)

    # Print all 96 colors
    print(spectral_96)

    donor_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Donor)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = spectral_96) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## UMI
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))

    UMI_plot <- ggplot(total_meta, aes(y = Cell_Type, x = log(nCount_RNA), fill = Cell_Type, colour = NULL)) + 
    geom_violin() +   # normalize counts to proportions
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = colors) #+
    #scale_colour_manual(values = colors) #+
    #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

    ############################################
    ## Gene
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))

    feature_plot <- ggplot(total_meta, aes(y = Cell_Type, x = log(nFeature_RNA), fill = Cell_Type, colour = NULL)) + 
    geom_violin() +   # normalize counts to proportions
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = colors) #+
    #scale_colour_manual(values = colors) #+
    #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)


    ############################################
    ## Arrange
    ############################################
    ggarrange(n_cell_plot, 
                region_prop_plot, 
                sex_prop_plot, 
                group_prop_plot, 
                donor_prop_plot, 
                UMI_plot,
                feature_plot,
                nrow = 1, ncol =7,
                widths = c(1.2,1,1,1,1,1,1),
                align = "h")
    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 10)


##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Limone supplemental figure
## code 
    module load StdEnv/2023
    module load r/4.4.0
    R

    ## Load libraries
    library(MAST)
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
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)
    library(colorspace)

    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)

    ############################################
    ## Read in metadata files
    ############################################
    ## Pineda
    pineda_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Pineda_combined_metadata.csv')
    pineda_meta <- pineda_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, CellType, Region)
    pineda_meta$dataset <- "Pineda et al."
    colnames(pineda_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    
    ## Li
    li_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Li_combined_metadata.csv')
    li_meta$Condition[li_meta$Group == "C9-ALS"] <- "ALS"
    li_meta$Condition[li_meta$Group == "C9-FTD"] <- "FTLD"
    li_meta$Condition[li_meta$Group == "Control"] <- "Control"
    li_meta <- li_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, sex, Condition, Group, cell_type, region)
    li_meta$dataset <- "Li et al."
    colnames(li_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")
    
    ## Limone
    limone_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Limone_combined_metadata.csv')
    colnames(limone_meta) 
    unique(limone_meta$Sample_ID)
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012218_FC21"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.121417_FC19"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC13"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012318_FC17"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC12"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.121417_FC11"] <- "F"
    limone_meta$Condition <- limone_meta$Group
    limone_meta$Group[limone_meta$Group == "ALS"] <- "SALS"
    limone_meta$Group[limone_meta$Group == "Control"] <- "Control"
    limone_meta <- limone_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, cell_type, region)
    limone_meta$dataset <- "Limone et al."
    colnames(limone_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    ## Bind and standardize
    total_meta <- rbind(pineda_meta, li_meta, limone_meta)
    total_meta$Sex[total_meta$Sex == "FALSE"] <- "F"
    total_meta$Condition[total_meta$Condition == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "C9-ALS"] <- "C9ALS"
    total_meta$Group[total_meta$Group == "C9-FTD"] <- "C9FTLD"
    
    unique(total_meta$Donor)
    unique(total_meta$Sex)
    unique(total_meta$Condition)
    unique(total_meta$Group)
    unique(total_meta$Cell_Type)
    unique(total_meta$Region)
    unique(total_meta$Dataset)

    ############################################
    ## Subset dataset of interest
    ############################################
    total_meta <- subset(total_meta, Dataset == "Limone et al.")
    unique(total_meta$Dataset)

    ############################################
    ## Set cell types levels
    ############################################
    cell_type_levels <- c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC')
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

    ############################################
    ## Plot number of cells per cell type
    ############################################
    n_cell <- data.frame(table(total_meta$Cell_Type))
    colnames(n_cell) <- c("cell_type", "Freq")
    n_cell$cell_type <- factor(n_cell$cell_type, levels =  rev(cell_type_levels))

    ## Plot
    n_cell_plot <-  ggplot(n_cell, aes(x = Freq, y = cell_type, fill = cell_type)) + 
    theme_classic() + 
    geom_bar(stat = 'identity') +
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size =8, colour = "black"),
        axis.text.y = element_text(size =8, colour = "black"),
        axis.title = element_blank(),
    ) +
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = "white") +
    scale_x_continuous(breaks = c(0, 5000, 10000), labels = c("0", "5e3", "10e3") )
    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

    ############################################
    ## Brain region prop
    ############################################
    spectral_base <- brewer.pal(11, "Paired")

    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Region <- factor(total_meta$Region , levels =  rev(c("BA4", "BA9")))
    unique(total_meta$Region)

    region_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Region )) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("BA9" = "#1B9E77",
                                "BA4"  = "#E31A1C")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Disease group prop
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Group <- factor(total_meta$Group , levels =  rev(c("SALS", "C9ALS", "SFTLD", "C9FTLD", "Control")))
    unique(total_meta$Group)

    group_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Group)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("SALS" = "orange",
                                "C9ALS"  = "red",
                                "SFTLD" = "blue",
                                "C9FTLD"  = "purple",
                                "Control" = "#339966")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Sex prop
    ############################################
    spectral_base <- brewer.pal(12, "Paired")
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    total_meta$Sex <- factor(total_meta$Sex , levels =  rev(c("M", "F")))
    unique(total_meta$Sex)

    sex_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Sex)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = c("M" = "#1F78B4",
                                "F"  = "#FB9A99")) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## Donor
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))
    unique(total_meta$Donor)
    
    spectral_base <- brewer.pal(12, "Paired")

    number <- length(unique(total_meta$Donor))

    # Interpolate to 96 colors in Lab space for smooth perceptual differences
    spectral_96 <- colorRampPalette(spectral_base, space = "Lab")(number)

    # Print all 96 colors
    print(spectral_96)

    donor_prop_plot <- ggplot(total_meta, aes(x = Cell_Type, fill = Donor)) + 
    geom_bar(position = "fill") +   # normalize counts to proportions
    coord_flip() +                  # flip so proportions are on x-axis
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = spectral_96) +
    scale_colour_manual(values = "white") +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)

    ############################################
    ## UMI
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))

    UMI_plot <- ggplot(total_meta, aes(y = Cell_Type, x = log(nCount_RNA), fill = Cell_Type, colour = NULL)) + 
    geom_violin() +   # normalize counts to proportions
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = colors) #+
    #scale_colour_manual(values = colors) #+
    #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 6)

    ############################################
    ## Gene
    ############################################
    total_meta$Cell_Type <- factor(total_meta$Cell_Type , levels =  rev(cell_type_levels))

    feature_plot <- ggplot(total_meta, aes(y = Cell_Type, x = log(nFeature_RNA), fill = Cell_Type, colour = NULL)) + 
    geom_violin() +   # normalize counts to proportions
    theme_classic() + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank()
    ) +
    scale_fill_manual(values = colors) #+
    #scale_colour_manual(values = colors) #+
    #scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1") ) 

    #ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 4, width = 1.75)


    ############################################
    ## Arrange
    ############################################
    ggarrange(n_cell_plot, 
                region_prop_plot, 
                sex_prop_plot, 
                group_prop_plot, 
                donor_prop_plot, 
                UMI_plot,
                feature_plot,
                nrow = 1, ncol =7,
                widths = c(1.2,1,1,1,1,1,1),
                align = "h")
    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 10)


##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Supplemental cell counts by disease group heatmap
## code 
    module load StdEnv/2023
    module load r/4.4.0
    R

    ## Load libraries
    library(MAST)
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
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)
    library(colorspace)

    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)

    ############################################
    ## Read in metadata files
    ############################################
    ## Pineda
    pineda_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Pineda_combined_metadata.csv')
    pineda_meta <- pineda_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, CellType, Region)
    pineda_meta$dataset <- "Pineda et al."
    colnames(pineda_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    
    ## Li
    li_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Li_combined_metadata.csv')
    li_meta$Condition[li_meta$Group == "C9-ALS"] <- "ALS"
    li_meta$Condition[li_meta$Group == "C9-FTD"] <- "FTLD"
    li_meta$Condition[li_meta$Group == "Control"] <- "Control"
    li_meta <- li_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, sex, Condition, Group, cell_type, region)
    li_meta$dataset <- "Li et al."
    colnames(li_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")
    
    ## Limone
    limone_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Limone_combined_metadata.csv')
    colnames(limone_meta) 
    unique(limone_meta$Sample_ID)
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012218_FC21"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.121417_FC19"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC13"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012318_FC17"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC12"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.121417_FC11"] <- "F"
    limone_meta$Condition <- limone_meta$Group
    limone_meta$Group[limone_meta$Group == "ALS"] <- "SALS"
    limone_meta$Group[limone_meta$Group == "Control"] <- "Control"
    limone_meta <- limone_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, cell_type, region)
    limone_meta$dataset <- "Limone et al."
    colnames(limone_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    ## Bind and standardize
    total_meta <- rbind(pineda_meta, li_meta, limone_meta)
    total_meta$Sex[total_meta$Sex == "FALSE"] <- "F"
    total_meta$Condition[total_meta$Condition == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "C9-ALS"] <- "C9ALS"
    total_meta$Group[total_meta$Group == "C9-FTD"] <- "C9FTLD"
    
    unique(total_meta$Donor)
    unique(total_meta$Sex)
    unique(total_meta$Condition)
    unique(total_meta$Group)
    unique(total_meta$Cell_Type)
    unique(total_meta$Region)
    unique(total_meta$Dataset)

    ############################################
    ## Only keep li and limone
    ############################################
    total_meta_eval <- subset(total_meta, Dataset %in% c("Li et al.", "Limone et al."))
    total_meta_eval$Dataset <- "Li et al. and Limone et al."

    total_meta <- rbind(total_meta, total_meta_eval)

    ############################################
    ## Set cell types levels
    ############################################
    cell_type_levels <- c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC')


    ############################################
    ## Plot number of cells per cell type
    ############################################
    n_cell <- as.data.frame(table(
        total_meta$Dataset,
        total_meta$Group,
        total_meta$Cell_Type
        ))
    colnames(n_cell) <- c("Dataset", "Group", "Cell_Type", "Freq")
    n_cell$Cell_Type <- factor(n_cell$Cell_Type, levels =  rev(cell_type_levels))
    n_cell$Dataset <- factor(n_cell$Dataset, levels =  c("Pineda et al.", "Li et al.", "Limone et al.", "Li et al. and Limone et al."))
    n_cell$Group <- factor(n_cell$Group, levels =  c("SALS", "C9ALS", "SFTLD", "C9FTLD", "Control"))
    n_cell <- subset(n_cell, Freq != 0)

    n_cell$Cat[n_cell$Dataset == "Pineda et al."] <- "Pineda et al.\n(training and LOSO validation)"
    n_cell$Cat[n_cell$Dataset == "Li et al. and Limone et al."] <- "Li et al. and Limone et al.\n(LOSO evaluation)"
    n_cell$Cat[n_cell$Dataset == "Li et al."] <- "Li et al.\n(LOSO evaluation)"
    n_cell$Cat[n_cell$Dataset == "Limone et al."] <- "Limone et al.\n(LOSO evaluation)"

    n_cell$Cat <- factor(n_cell$Cat, levels =  c("Pineda et al.\n(training and LOSO validation)", 
                                                "Li et al.\n(LOSO evaluation)", 
                                                "Limone et al.\n(LOSO evaluation)", 
                                                "Li et al. and Limone et al.\n(LOSO evaluation)"))

    n_cell$fill <- "white"
    n_cell$fill[n_cell$Freq < 75 & n_cell$Cat == "Li et al. and Limone et al.\n(LOSO evaluation)"] <- "red"

    ## Plot
    n_cell_plot <-  ggplot(n_cell, aes(x = Group, y = Cell_Type, label = Freq, fill = fill, colour = "black")) + 
    theme_bw() + 
    geom_tile(stat = 'identity') +
    geom_text(size = 2.5) + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size =8, colour = c("black")),
        axis.text.y = element_text(size =8, colour = "black"),
        axis.title = element_blank(),
    ) +
    scale_fill_manual(values = c("red","white")) +
    scale_colour_manual(values = "black") +
    facet_grid(. ~ Cat, scales = "free_x", space = "free_x")

    #scale_fill_manual(values = colors) +
    #scale_colour_manual(values = "white") +
    #scale_x_continuous(breaks = c(0, 5000, 10000), labels = c("0", "5e3", "10e3") )
    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 10)

##

################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################
################################################################################################ Supplemental cell counts by broad disease group heatmap
## code 
    module load StdEnv/2023
    module load r/4.4.0
    R

    ## Load libraries
    library(MAST)
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
    library(tidyverse, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)
    library(colorspace)

    library(tidyverse, lib="/home/fiorini9/scracth/R")
    library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
    library(readr)
    library(reshape2)

    ############################################
    ## Read in metadata files
    ############################################
    ## Pineda
    pineda_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Pineda_combined_metadata.csv')
    pineda_meta <- pineda_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, CellType, Region)
    pineda_meta$dataset <- "Pineda et al."
    colnames(pineda_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    
    ## Li
    li_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Li_combined_metadata.csv')
    li_meta$Condition[li_meta$Group == "C9-ALS"] <- "ALS"
    li_meta$Condition[li_meta$Group == "C9-FTD"] <- "FTLD"
    li_meta$Condition[li_meta$Group == "Control"] <- "Control"
    li_meta <- li_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, sex, Condition, Group, cell_type, region)
    li_meta$dataset <- "Li et al."
    colnames(li_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")
    
    ## Limone
    limone_meta <- read.csv('/home/fiorini9/scratch/machine_learning_ALS/model_outs/Limone_combined_metadata.csv')
    colnames(limone_meta) 
    unique(limone_meta$Sample_ID)
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012218_FC21"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.121417_FC19"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC13"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "ALS_MotorCortex.012318_FC17"] <- "F"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.012218_FC12"] <- "M"
    limone_meta$Sex[limone_meta$Sample_ID == "Control_MotorCortex.121417_FC11"] <- "F"
    limone_meta$Condition <- limone_meta$Group
    limone_meta$Group[limone_meta$Group == "ALS"] <- "SALS"
    limone_meta$Group[limone_meta$Group == "Control"] <- "Control"
    limone_meta <- limone_meta %>% dplyr::select(nCount_RNA, nFeature_RNA, Sample_ID, Donor, Sex, Condition, Group, cell_type, region)
    limone_meta$dataset <- "Limone et al."
    colnames(limone_meta) <- c("nCount_RNA", "nFeature_RNA", "Sample_ID", "Donor", "Sex", "Condition", "Group", "Cell_Type", "Region", "Dataset")

    ## Bind and standardize
    total_meta <- rbind(pineda_meta, li_meta, limone_meta)
    total_meta$Sex[total_meta$Sex == "FALSE"] <- "F"
    total_meta$Condition[total_meta$Condition == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "PN"] <- "Control"
    total_meta$Group[total_meta$Group == "C9-ALS"] <- "C9ALS"
    total_meta$Group[total_meta$Group == "C9-FTD"] <- "C9FTLD"
    
    unique(total_meta$Donor)
    unique(total_meta$Sex)
    unique(total_meta$Condition)
    unique(total_meta$Group)
    unique(total_meta$Cell_Type)
    unique(total_meta$Region)
    unique(total_meta$Dataset)

    total_meta$Group[total_meta$Group == "SALS"] <- "ALS"
    total_meta$Group[total_meta$Group == "C9ALS"] <- "ALS"
    total_meta$Group[total_meta$Group == "SFTLD"] <- "FTLD"
    total_meta$Group[total_meta$Group == "C9FTLD"] <- "FTLD"
    
    ############################################
    ## Only keep li and limone
    ############################################
    total_meta_eval <- subset(total_meta, Dataset %in% c("Li et al.", "Limone et al."))
    total_meta_eval$Dataset <- "Li et al. and Limone et al."

    total_meta <- rbind(total_meta, total_meta_eval)

    ############################################
    ## Set cell types levels
    ############################################
    cell_type_levels <- c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC')


    ############################################
    ## Plot number of cells per cell type
    ############################################
    n_cell <- as.data.frame(table(
        total_meta$Dataset,
        total_meta$Group,
        total_meta$Cell_Type
        ))
    colnames(n_cell) <- c("Dataset", "Group", "Cell_Type", "Freq")
    n_cell$Cell_Type <- factor(n_cell$Cell_Type, levels =  rev(cell_type_levels))
    n_cell$Dataset <- factor(n_cell$Dataset, levels =  c("Pineda et al.", "Li et al.", "Limone et al.", "Li et al. and Limone et al."))
    n_cell$Group <- factor(n_cell$Group, levels =  c("ALS", "FTLD",  "Control"))
    n_cell <- subset(n_cell, Freq != 0)

    n_cell$Cat[n_cell$Dataset == "Pineda et al."] <- "Pineda et al.\n(training and LOSO validation)"
    n_cell$Cat[n_cell$Dataset == "Li et al. and Limone et al."] <- "Li et al. and Limone et al.\n(LOSO evaluation)"
    n_cell$Cat[n_cell$Dataset == "Li et al."] <- "Li et al.\n(LOSO evaluation)"
    n_cell$Cat[n_cell$Dataset == "Limone et al."] <- "Limone et al.\n(LOSO evaluation)"

    n_cell$Cat <- factor(n_cell$Cat, levels =  c("Pineda et al.\n(training and LOSO validation)", 
                                                "Li et al.\n(LOSO evaluation)", 
                                                "Limone et al.\n(LOSO evaluation)", 
                                                "Li et al. and Limone et al.\n(LOSO evaluation)"))

    n_cell$fill <- "white"
    n_cell$fill[n_cell$Freq < 50 & n_cell$Cat == "Li et al. and Limone et al.\n(LOSO evaluation)"] <- "red"

    ## Plot
    n_cell_plot <-  ggplot(n_cell, aes(x = Group, y = Cell_Type, label = Freq, fill = fill, colour = "black")) + 
    theme_bw() + 
    geom_tile(stat = 'identity') +
    geom_text(size = 2.5) + 
    theme(
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(size =8, colour = c("black")),
        axis.text.y = element_text(size =8, colour = "black"),
        axis.title = element_blank(),
    ) +
    scale_fill_manual(values = c("red","white")) +
    scale_colour_manual(values = "black") +
    facet_grid(. ~ Cat, scales = "free_x", space = "free_x")

    #scale_fill_manual(values = colors) +
    #scale_colour_manual(values = "white") +
    #scale_x_continuous(breaks = c(0, 5000, 10000), labels = c("0", "5e3", "10e3") )
    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf', height = 2.5, width = 10)

    write.csv(n_cell, '/home/fiorini9/scratch/machine_learning_ALS/outs/cell_count_broad_disease_group.csv')

##