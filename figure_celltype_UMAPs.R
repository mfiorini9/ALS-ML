## We will run this in Beluga then transfer base models and figures to Narval

###################################################################################################################################################################################################### Beluga code
######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################

################################################################## ## Merge seurat objects and print main figure UMAPS
##################################################################
################################################################## Performed this in Beluga. 
##################################################################
##################################################################
################################################################## Pineda
## code
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=300g

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano UMAPs_script2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-grouleau
    #SBATCH --time=00-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=300g          # memory per cor
    #SBATCH --job-name=UMAPs_script2
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/UMAPs_script2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano UMAPs_script2.R

    ## code
    library(harmony)
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2)
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(scCustomize, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggrepel)
    library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")


    #############################
    # Excitatory neurons
    #############################
    ## code
        ## read in BA4
        seu_BA4_ex <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA4.rds')

        ## only retain excitatory class
        unique(seu_BA4_ex@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu_BA4_ex@meta.data$Region) #"BA4"
        seu_celltype <- c("Ex")
        xx <- unique(seu_BA4_ex@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA4_ex) <- "CellClass"
        seu_BA4_ex=subset(seu_BA4_ex,idents=xx)
        unique(seu_BA4_ex@meta.data$CellClass)
        dim(seu_BA4_ex)

        ## read in BA9
        seu_BA9_ex <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA9.rds')

        ## only retain excitatory class
        unique(seu_BA9_ex@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu_BA9_ex@meta.data$Region) #"BA9"
        seu_celltype <- c("Ex")
        xx <- unique(seu_BA9_ex@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA9_ex) <- "CellClass"
        seu_BA9_ex=subset(seu_BA9_ex,idents=xx)
        unique(seu_BA9_ex@meta.data$CellClass)
        dim(seu_BA9_ex)

        ## merge
        DefaultAssay(seu_BA4_ex) <- "RNA"
        DefaultAssay(seu_BA9_ex) <- "RNA"

        seu_list <- list(seu_BA4_ex = seu_BA4_ex, seu_BA9_ex = seu_BA9_ex)

        seu_merge_ex <- Merge_Seurat_List(
        list_seurat = seu_list,
        add.cell.ids = NULL,
        merge.data = TRUE,
        project = "MergeSeurat"
        )

        unique(seu_merge_ex@meta.data$Region)
        unique(seu_merge_ex@meta.data$Condition)

        rm(seu_BA4_ex)
        rm(seu_BA9_ex)

        ## Run Harmony
        seu_merge_ex <- ScaleData(seu_merge_ex)    
        seu_merge_ex<- FindVariableFeatures(seu_merge_ex, selection.method = "vst", nfeatures = 2500)
        seu_merge_ex <- RunPCA(seu_merge_ex, npcs = 30)
        seu_merge_ex <- RunHarmony(seu_merge_ex, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_merge_ex <- RunUMAP(seu_merge_ex, reduction = "harmony", dims = 1:15)
        seu_merge_ex <- FindNeighbors(seu_merge_ex, reduction = "harmony", dims = 1:15)
        seu_merge_ex <- FindClusters(seu_merge_ex)

        ## By Celltype
        DimPlot(seu_merge_ex, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "L3_L5" = "#f3c300",
            "L2_L3" = "#f38400",
            "L4_L6" = "#a1caf1",
            "L4_L5" = "#be0032",
            "L5_L6" = "#c2b280",
            "L5" = "#008856",
            "L6" = "#2b3d26"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        ncol(seu_merge_ex)

        ## fix that one cluster
        #DimPlot(seu_merge_ex, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3,repel = TRUE, raster=FALSE) +
        #theme(axis.text.x = element_text(size =10),
        #axis.text.y = element_text(size =10),
        #axis.title = element_text(face="bold", size =10),
        #legend.position = "right",
        #plot.title = element_text( size =10)) + 
        #xlab("UMAP1") + ylab("UMAP2")
        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf'))


        #seu_merge_ex@meta.data$mod_subtype_clusters <- seu_merge_ex@meta.data$full_label
        #seu_merge_ex@meta.data$mod_subtype_clusters[seu_merge_ex@meta.data$seurat_clusters == 4] <- "Ex.L4_L5.RORB_POU3F2"

        ## By Subtype
        #DimPlot(seu_merge_ex, reduction = "umap", group.by = "mod_subtype_clusters", label = FALSE, label.size = 2, repel = TRUE, raster=FALSE) +
        #theme_void() +
        #theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        #axis.title = element_blank(),
        #legend.position = "none",
        #plot.title = element_blank()) + 
        #xlab("UMAP1") + ylab("UMAP2") +
        #scale_color_manual(values = c(
        #    "Ex.L3_L5.CUX2_RORB" = "#f3c300",
        #    "Ex.L2_L3.CUX2_RASGRF2"="#875692",
        #    "Ex.L4_L6.RORB_LRRK1"="#f38400",
        #    "Ex.L4_L5.RORB_POU3F2"="#a1caf1",
        #    "Ex.L4_L5.RORB_FOXO1"="#be0032",
        #    "Ex.L5_L6.THEMIS_NR4A2"="#c2b280",
        #    "Ex.L5_L6.THEMIS_TMEM233"="#848482",
        #    "Ex.L4_L6.RORB_ADGRL4"="#008856",
        #    "Ex.L6.TLE4_MEGF11"="#e68fac",
        #    "Ex.L3_L5.SCN4B_NEFH"="#0067a5",
        #    "Ex.L6.TLE4_SEMA3D"="#f99379",
        #    "Ex.L5.VAT1L_EYA4"="#604e97",
        #    "Ex.L5.PCP4_NXPH2"="#f6a600",
        #   "Ex.L5.VAT1L_THSD4"="#b3446c",
        #    "Ex.L6.TLE4_CCBE1"="#dcd300"))

        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'), height = 5, width = 5)
        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Pineda_excitatory_UMAP.pdf'), height = 5, width = 5)

        saveRDS(seu_merge_ex, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_Excitatory_neurons.rds')
        rm(seu_BA9_ex)
    ##

    #############################
    # Inhibitory neurons
    #############################
    ## code
        ## read in BA4
        seu_BA4_in <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA4.rds')

        ## only retain excitatory class
        unique(seu_BA4_in@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu_BA4_in@meta.data$Region) #"BA4"
        seu_celltype <- c("In")
        xx <- unique(seu_BA4_in@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA4_in) <- "CellClass"
        seu_BA4_in=subset(seu_BA4_in,idents=xx)
        unique(seu_BA4_in@meta.data$CellClass)
        dim(seu_BA4_in)

        ## read in BA9
        seu_BA9_in <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA9.rds')

        ## only retain excitatory class
        unique(seu_BA9_in@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu_BA9_in@meta.data$Region) #"BA9"
        seu_celltype <- c("In")
        xx <- unique(seu_BA9_in@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA9_in) <- "CellClass"
        seu_BA9_in=subset(seu_BA9_in,idents=xx)
        unique(seu_BA9_in@meta.data$CellClass)
        dim(seu_BA9_in)

        ## merge
        DefaultAssay(seu_BA4_in) <- "RNA"
        DefaultAssay(seu_BA9_in) <- "RNA"

        seu_list <- list(seu_BA4_in = seu_BA4_in, seu_BA9_in = seu_BA9_in)

        seu_merge_in <- Merge_Seurat_List(
        list_seurat = seu_list,
        add.cell.ids = NULL,
        merge.data = TRUE,
        project = "MergeSeurat"
        )

        unique(seu_merge_in@meta.data$Region)
        unique(seu_merge_in@meta.data$Condition)

        rm(seu_BA4_in)
        rm(seu_BA9_in)

        ## Run Harmony
        seu_merge_in <- ScaleData(seu_merge_in)    
        seu_merge_in<- FindVariableFeatures(seu_merge_in, selection.method = "vst", nfeatures = 2500)
        seu_merge_in <- RunPCA(seu_merge_in, npcs = 15)
        seu_merge_in <- RunHarmony(seu_merge_in, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_merge_in <- RunUMAP(seu_merge_in, reduction = "harmony", dims = 1:15)
        seu_merge_in <- FindNeighbors(seu_merge_in, reduction = "harmony", dims = 1:15)
        seu_merge_in <- FindClusters(seu_merge_in)

        ## By Celltype
        DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "PV" = "#e25822",
            "5HT3aR" = "#e68fac",
            "Rosehip" = "#0067a5",
            "SOM" = "#f99379"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        ncol(seu_merge_in)

        ## By subtype
        #DimPlot(seu_merge_in, reduction = "umap", group.by = "full_label", label = TRUE, label.size = 3,repel = TRUE, raster=FALSE) +
        #theme(axis.text.x = element_text(size =10),
        #axis.text.y = element_text(size =10),
        #axis.title = element_text(face="bold", size =10),
        #legend.position = "right",
        #plot.title = element_text( size =10)) + 
        #xlab("UMAP1") + ylab("UMAP2")
        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'), height = 5, width = 7)

        ## fix that one cluster
        #DimPlot(seu_merge_in, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3,repel = TRUE, raster=FALSE) +
        #theme(axis.text.x = element_text(size =10),
        #axis.text.y = element_text(size =10),
        #axis.title = element_text(face="bold", size =10),
        #legend.position = "right",
        #plot.title = element_text( size =10)) + 
        #xlab("UMAP1") + ylab("UMAP2")
        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf'))

        #seu_merge_in@meta.data$mod_subtype_clusters <- seu_merge_in@meta.data$full_label
        #seu_merge_in@meta.data$mod_subtype_clusters[seu_merge_in@meta.data$seurat_clusters == 4] <- "Ex.L4_L5.RORB_POU3F2"

        #unique(seu_merge_in@meta.data$full_label)
        ## By Subtype
        #DimPlot(seu_merge_in, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 2, repel = TRUE, raster=FALSE) +
        #theme_void() +
        #theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        #axis.title = element_blank(),
        #legend.position = "none",
        #plot.title = element_blank()) + 
        #xlab("UMAP1") + ylab("UMAP2") +
        #scale_color_manual(values = c(
        #    "In.PV.PVALB_CEMIP" = "#f3c300",
        #    "In.5HT3aR.VIP_HTR2C"="#875692",
        #    "In.5HT3aR.VIP_CLSTN2"="#f38400",
        #    "In.5HT3aR.VIP_LAMA3"="#a1caf1",
        #    "In.Rosehip.LAMP5_CA3"="#be0032",
        #    "In.5HT3aR.DISC1_CCK"="#c2b280",
        #    "In.PV.PVALB_PTHLH"="#848482",
        #    "In.5HT3aR.CDH4_CCK"="#008856",
        #    "In.SOM.SST_GALNT14"="#e68fac",
        #    "In.SOM.SST_ADAMTS19"="#0067a5",
        #    "In.Rosehip.LAMP5_PMEPA1"="#f99379",
        #    "In.SOM.SST_BRINP3"="#604e97",
        #    "In.5HT3aR.DISC1_RELN"="#f6a600",
        #    "In.SOM.SST_NPY"="#b3446c",
        #    "In.PV.PVALB_MYBPC1"="#dcd300",
        #   "In.5HT3aR.CDH4_SCGN"="#882d17"))

        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'), height = 5, width = 5)
        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Pineda_inhibitory_UMAP.pdf'), height = 5, width = 5)


      ## Cell type marker genes

        saveRDS(seu_merge_in, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_inhibitory_neurons.rds')
        rm(seu_merge_in)
    ##

    #############################
    # Non-neuronal cells
    #############################
    ## code
        ## read in BA4
        seu_BA4_in <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA4.rds')

        ## only retain excitatory class
        unique(seu_BA4_in@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu_BA4_in@meta.data$Region) #"BA4"
        seu_celltype <- c("Glia", "Vasc")
        xx <- unique(seu_BA4_in@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA4_in) <- "CellClass"
        seu_BA4_in=subset(seu_BA4_in,idents=xx)
        unique(seu_BA4_in@meta.data$CellClass)
        dim(seu_BA4_in)

        ## read in BA9
        seu_BA9_in <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA9.rds')

        ## only retain excitatory class
        unique(seu_BA9_in@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu_BA9_in@meta.data$Region) #"BA9"
        seu_celltype <- c("Glia", "Vasc")
        xx <- unique(seu_BA9_in@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA9_in) <- "CellClass"
        seu_BA9_in=subset(seu_BA9_in,idents=xx)
        unique(seu_BA9_in@meta.data$CellClass)
        dim(seu_BA9_in)

        ## merge
        DefaultAssay(seu_BA4_in) <- "RNA"
        DefaultAssay(seu_BA9_in) <- "RNA"

        seu_list <- list(seu_BA4_in = seu_BA4_in, seu_BA9_in = seu_BA9_in)

        seu_merge_in <- Merge_Seurat_List(
        list_seurat = seu_list,
        add.cell.ids = NULL,
        merge.data = TRUE,
        project = "MergeSeurat"
        )

        unique(seu_merge_in@meta.data$Region)
        unique(seu_merge_in@meta.data$Condition)

        rm(seu_BA4_in)
        rm(seu_BA9_in)

        ## Run Harmony
        seu_merge_in <- ScaleData(seu_merge_in)    
        seu_merge_in<- FindVariableFeatures(seu_merge_in, selection.method = "vst", nfeatures = 2500)
        seu_merge_in <- RunPCA(seu_merge_in, npcs = 20)
        seu_merge_in <- RunHarmony(seu_merge_in, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_merge_in <- RunUMAP(seu_merge_in, reduction = "harmony", dims = 1:20)
        seu_merge_in <- FindNeighbors(seu_merge_in, reduction = "harmony", dims = 1:20)
        seu_merge_in <- FindClusters(seu_merge_in)

        ## By Celltype
        DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "Oligo" = "#604e97",
            "Astro" = "#875692",
            "OPC" = "#f6a600",
            "Micro" = "#b3446c",
            "T_Cell" = "#dcd300",
            "Mural" = "#882d17",
            "Endo" = "#8db600",
            "Fibro" = "#654522"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        ncol(seu_merge_in)

        ## By subtype
        #DimPlot(seu_merge_in, reduction = "umap", group.by = "full_label", label = TRUE, label.size = 3,repel = TRUE, raster=FALSE) +
        #theme(axis.text.x = element_text(size =10),
        #axis.text.y = element_text(size =10),
        #axis.title = element_text(face="bold", size =10),
        #legend.position = "right",
        #plot.title = element_text( size =10)) + 
        #xlab("UMAP1") + ylab("UMAP2")
        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'), height = 5, width = 7)

        ## fix that one cluster
        #DimPlot(seu_merge_in, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3,repel = TRUE, raster=FALSE) +
        #theme(axis.text.x = element_text(size =10),
        #axis.text.y = element_text(size =10),
        #axis.title = element_text(face="bold", size =10),
        #legend.position = "right",
        #plot.title = element_text( size =10)) + 
        #xlab("UMAP1") + ylab("UMAP2")
        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp2.pdf'))

        #seu_merge_in@meta.data$mod_subtype_clusters <- seu_merge_in@meta.data$full_label
        #seu_merge_in@meta.data$mod_subtype_clusters[seu_merge_in@meta.data$seurat_clusters == 4] <- "Ex.L4_L5.RORB_POU3F2"

        #unique(seu_merge_in@meta.data$full_label)
        ## By Subtype
        #DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 2, repel = TRUE, raster=FALSE) +
        #theme_void() +
        #theme(axis.text.x = element_blank(),
        #axis.text.y = element_blank(),
        #axis.title = element_blank(),
        #legend.position = "none",
        #plot.title = element_blank()) + 
        #xlab("UMAP1") + ylab("UMAP2") +
        #scale_color_manual(values = c(
        #    "Glia.Oligo" = "#f3c300",
        #    "Glia.Astro.GFAP-neg"="#875692",
        #    "Glia.Astro.GFAP-pos"="#f38400",
        #    "Glia.OPC"="#a1caf1",
        #    "Glia.Micro"="#be0032",
        #    "Vasc.T_Cell"="#c2b280",
        #   "Vasc.Mural.Pericyte"="#848482",
        #    "Vasc.Endo.Capillary"="#008856",
        #    "Vasc.Fibro.CLMP_PDGFRA"="#e68fac",
        #    "Vasc.Endo.Venous"="#0067a5",
        #    "Vasc.Mural.SMC"="#f99379",
        #    "Vasc.Endo.Arterial"="#604e97",
        #    "Vasc.Fibro.CLMP_KCNMA1"="#f6a600"))

        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'), height = 5, width = 5)
        #ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Pineda_non_neuron_UMAP.pdf'), height = 5, width = 5)


      ## Cell type marker genes

        saveRDS(seu_merge_in, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_non_neurons.rds')
    ##
##


################################################################## ## Create a single merged object to perform joint analyses -- this has BA4 and BA9 all cell types, subsetted to only include 300 000 cells
##################################################################
################################################################## Performed this in Beluga. 
##################################################################
##################################################################
################################################################## Pineda

## code

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano merged_all_object.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-grouleau
    #SBATCH --time=00-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=4500g          # memory per cor
    #SBATCH --job-name=merged_all_object
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/merged_all_object.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano merged_all_object.R

    ## load libraries
    library(harmony)
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2)
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(scCustomize, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggrepel)
    library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")


    ##############
    ## Excitatory neurons
    ##############
    
    seu_ex <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_Excitatory_neurons.rds') ## This object was created in the above chunk of this document
    ncol(seu_ex)
    # subset to 132000
    # Step 1: Get the total number of cells in your Seurat object
    total_cells <- ncol(seu_ex)

    # Step 2: Randomly sample 132,000 cell names
    set.seed(42)  # Set seed for reproducibility (optional)
    cells_to_keep <- sample(colnames(seu_ex), 132000)

    # Step 3: Subset the Seurat object
    seu_ex_lim <- subset(seu_ex, cells = cells_to_keep)
    ncol(seu_ex_lim)
    rm(seu_ex)

    ##############
    ## Inhibitory neurons
    ##############

    seu_in <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_inhibitory_neurons.rds') ## This object was created in the above chunk of this document
    ncol(seu_in)
    # subset to 54000
    # Step 1: Get the total number of cells in your Seurat object
    total_cells <- ncol(seu_in)

    # Step 2: Randomly sample 54000 cell names
    set.seed(42)  # Set seed for reproducibility (optional)
    cells_to_keep <- sample(colnames(seu_in), 54000)

    # Step 3: Subset the Seurat object
    seu_in_lim <- subset(seu_in, cells = cells_to_keep)
    ncol(seu_in_lim)
    rm(seu_in)

    ##############
    ## non neurons
    ##############
    
    seu_non <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_non_neurons.rds') ## This object was created in the above chunk of this document
    ncol(seu_non)
    # subset to 99000
    # Step 1: Get the total number of cells in your Seurat object
    total_cells <- ncol(seu_non)

    # Step 2: Randomly sample 99000 cell names
    set.seed(42)  # Set seed for reproducibility (optional)
    cells_to_keep <- sample(colnames(seu_non), 99000)

    # Step 3: Subset the Seurat object
    seu_non_lim <- subset(seu_non, cells = cells_to_keep)
    ncol(seu_non_lim)
    rm(seu_non)

    ##############
    ## merge
    ##############
    DefaultAssay(seu_ex_lim) <- "RNA"
    DefaultAssay(seu_in_lim) <- "RNA"
    DefaultAssay(seu_non_lim) <- "RNA"

    seu_list <- list(seu_ex_lim = seu_ex_lim, seu_in_lim = seu_in_lim, seu_non_lim = seu_non_lim )

    seu_merge_in <- Merge_Seurat_List(
    list_seurat = seu_list,
    add.cell.ids = NULL,
    merge.data = TRUE,
    project = "MergeSeurat"
    )

    unique(seu_merge_in@meta.data$Region)
    unique(seu_merge_in@meta.data$Condition)
    unique(seu_merge_in@meta.data$CellType)

    rm(seu_ex_lim)
    rm(seu_in_lim)
    rm(seu_non_lim)

    ## Run Harmony
    seu_merge_in <- ScaleData(seu_merge_in)    
    seu_merge_in<- FindVariableFeatures(seu_merge_in, selection.method = "vst", nfeatures = 2500)
    seu_merge_in <- RunPCA(seu_merge_in, npcs = 30)
    seu_merge_in <- RunHarmony(seu_merge_in, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_merge_in <- RunUMAP(seu_merge_in, reduction = "harmony", dims = 1:30)
    seu_merge_in <- FindNeighbors(seu_merge_in, reduction = "harmony", dims = 1:30)
    seu_merge_in <- FindClusters(seu_merge_in)

    ## By Celltype
    DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme_void()+
    theme(axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_blank()) + 
    xlab("UMAP1") + ylab("UMAP2") 

    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/UMAP_all_celltypes_temp.pdf'))
    ncol(seu_merge_in)

    ## Save the merged all object
    saveRDS(seu_merge_in, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_all_celltypes_lim_narval.rds')


##


################################################################## ## Merge seurat objects and print main figure UMAPS
##################################################################
################################################################## Performed this in Beluga. 
##################################################################
##################################################################
################################################################## Li
## code
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=300g

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano UMAPs_script_Li.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-grouleau
    #SBATCH --time=00-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=300g          # memory per cor
    #SBATCH --job-name=UMAPs_script_Li
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/UMAPs_script_Li.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano UMAPs_script_Li.R

    ## code
    library(harmony)
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2)
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(scCustomize, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggrepel)
    library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")


    #############################
    # Excitatory neurons
    #############################
    ## code
        ## read in BA4
        seu_BA4_ex <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_All_conditions_BA4.rds')
        ncol(seu_BA4_ex)

        ## only retain excitatory class
        unique(seu_BA4_ex@meta.data$Region) #"BA4"
        unique(seu_BA4_ex@meta.data$CellType) #"BA4
        unique(seu_BA4_ex@meta.data$annotation_major_cell_type) #"BA4
        
        ## Add major CellClass metadata
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L2_L3"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L3_L5"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L4_L5"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L4_L6"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L5"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L5_L6"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L6"] <- "Ex"

        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "5HT3aR"] <- "In"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "PV"] <- "In"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Rosehip"] <- "In"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "SOM"] <- "In"

        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Astro"] <- "Glia"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Endo"] <- "Vasc"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Fibro"] <- "Vasc"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Micro"] <- "Glia"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Mural"] <- "Vasc"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Oligo"] <- "Glia"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "OPC"] <- "Glia"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "T_Cell"] <- "Vasc"

        unique(seu_BA4_ex@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

        seu_celltype <- c("Ex")
        xx <- unique(seu_BA4_ex@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA4_ex) <- "CellClass"
        seu_BA4_ex=subset(seu_BA4_ex,idents=xx)
        unique(seu_BA4_ex@meta.data$CellClass)
        dim(seu_BA4_ex)

        ## read in BA9
        seu_BA9_ex <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_All_conditions_BA9.rds')
        ncol(seu_BA9_ex)

        ## only retain excitatory class
        unique(seu_BA9_ex@meta.data$Region) #"BA9"
        unique(seu_BA9_ex@meta.data$CellType) #"BA9
        unique(seu_BA9_ex@meta.data$annotation_major_cell_type) #"BA9
        
        ## Add major CellClass metadata
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "L2_L3"] <- "Ex"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "L3_L5"] <- "Ex"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "L4_L5"] <- "Ex"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "L4_L6"] <- "Ex"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "L5"] <- "Ex"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "L5_L6"] <- "Ex"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "L6"] <- "Ex"

        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "5HT3aR"] <- "In"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "PV"] <- "In"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "Rosehip"] <- "In"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "SOM"] <- "In"

        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "Astro"] <- "Glia"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "Endo"] <- "Vasc"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "Fibro"] <- "Vasc"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "Micro"] <- "Glia"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "Mural"] <- "Vasc"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "Oligo"] <- "Glia"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "OPC"] <- "Glia"
        seu_BA9_ex@meta.data$CellClass[seu_BA9_ex@meta.data$CellType == "T_Cell"] <- "Vasc"

        unique(seu_BA9_ex@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

        seu_celltype <- c("Ex")
        xx <- unique(seu_BA9_ex@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA9_ex) <- "CellClass"
        seu_BA9_ex=subset(seu_BA9_ex,idents=xx)
        unique(seu_BA9_ex@meta.data$CellClass)
        dim(seu_BA9_ex)

        ## merge
        DefaultAssay(seu_BA4_ex) <- "RNA"
        DefaultAssay(seu_BA9_ex) <- "RNA"

        seu_list <- list(seu_BA4_ex = seu_BA4_ex, seu_BA9_ex = seu_BA9_ex)

        seu_merge_ex <- Merge_Seurat_List(
        list_seurat = seu_list,
        add.cell.ids = NULL,
        merge.data = TRUE,
        project = "MergeSeurat"
        )

        unique(seu_merge_ex@meta.data$Region)
        unique(seu_merge_ex@meta.data$Group)

        rm(seu_BA4_ex)
        rm(seu_BA9_ex)

        ## Run Harmony
        seu_merge_ex <- ScaleData(seu_merge_ex)    
        seu_merge_ex<- FindVariableFeatures(seu_merge_ex, selection.method = "vst", nfeatures = 2500)
        seu_merge_ex <- RunPCA(seu_merge_ex, npcs = 15)
        seu_merge_ex <- RunHarmony(seu_merge_ex, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_merge_ex <- RunUMAP(seu_merge_ex, reduction = "harmony", dims = 1:15)
        seu_merge_ex <- FindNeighbors(seu_merge_ex, reduction = "harmony", dims = 1:15)
        seu_merge_ex <- FindClusters(seu_merge_ex)

        ## By Celltype
        DimPlot(seu_merge_ex, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "L3_L5" = "#f3c300",
            "L2_L3" = "#f38400",
            "L4_L6" = "#a1caf1",
            "L4_L5" = "#be0032",
            "L5_L6" = "#c2b280",
            "L5" = "#008856",
            "L6" = "#2b3d26"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        ncol(seu_merge_ex)

        ## By Seurat cluster
        DimPlot(seu_merge_ex, reduction = "umap", group.by = "seurat_clusters", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") 

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        
        
        ## Cluster 15 did not match to any cell types from Pineda et al. we will remove it.
        seu_merge_ex <- subset(seu_merge_ex, seurat_clusters != 15)

        ## By Celltype
        DimPlot(seu_merge_ex, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "L3_L5" = "#f3c300",
            "L2_L3" = "#f38400",
            "L4_L6" = "#a1caf1",
            "L4_L5" = "#be0032",
            "L5_L6" = "#c2b280",
            "L5" = "#008856",
            "L6" = "#2b3d26"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        ncol(seu_merge_ex)

        saveRDS(seu_merge_ex, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_BA4_BA9_Excitatory_neurons.rds')
        rm(seu_BA9_ex)
    ##

    #############################
    # Inhibitory neurons
    #############################
    ## code
        ## read in BA4
        seu_BA4_in <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_All_conditions_BA4.rds')
        ncol(seu_BA4_in)

        ## only retain Inhibitory class
        unique(seu_BA4_in@meta.data$Region) #"BA4"
        unique(seu_BA4_in@meta.data$CellType) #"BA4
        unique(seu_BA4_in@meta.data$annotation_major_cell_type) #"BA4
        
        ## Add major CellClass metadata
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L2_L3"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L3_L5"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L4_L5"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L4_L6"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L5"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L5_L6"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L6"] <- "Ex"

        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "5HT3aR"] <- "In"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "PV"] <- "In"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Rosehip"] <- "In"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "SOM"] <- "In"

        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Astro"] <- "Glia"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Endo"] <- "Vasc"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Fibro"] <- "Vasc"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Micro"] <- "Glia"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Mural"] <- "Vasc"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Oligo"] <- "Glia"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "OPC"] <- "Glia"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "T_Cell"] <- "Vasc"

        unique(seu_BA4_in@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

        seu_celltype <- c("In")
        xx <- unique(seu_BA4_in@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA4_in) <- "CellClass"
        seu_BA4_in=subset(seu_BA4_in,idents=xx)
        unique(seu_BA4_in@meta.data$CellClass)
        dim(seu_BA4_in)

        ## read in BA9
        seu_BA9_in <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_All_conditions_BA9.rds')
        ncol(seu_BA9_in)

        ## only retain Inhibitory class
        unique(seu_BA9_in@meta.data$Region) #"BA9"
        unique(seu_BA9_in@meta.data$CellType) #"BA9
        unique(seu_BA9_in@meta.data$annotation_major_cell_type) #"BA9
        
        ## Add major CellClass metadata
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "L2_L3"] <- "Ex"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "L3_L5"] <- "Ex"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "L4_L5"] <- "Ex"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "L4_L6"] <- "Ex"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "L5"] <- "Ex"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "L5_L6"] <- "Ex"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "L6"] <- "Ex"

        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "5HT3aR"] <- "In"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "PV"] <- "In"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "Rosehip"] <- "In"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "SOM"] <- "In"

        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "Astro"] <- "Glia"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "Endo"] <- "Vasc"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "Fibro"] <- "Vasc"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "Micro"] <- "Glia"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "Mural"] <- "Vasc"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "Oligo"] <- "Glia"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "OPC"] <- "Glia"
        seu_BA9_in@meta.data$CellClass[seu_BA9_in@meta.data$CellType == "T_Cell"] <- "Vasc"

        unique(seu_BA9_in@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

        seu_celltype <- c("In")
        xx <- unique(seu_BA9_in@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA9_in) <- "CellClass"
        seu_BA9_in=subset(seu_BA9_in,idents=xx)
        unique(seu_BA9_in@meta.data$CellClass)
        dim(seu_BA9_in)

        ## merge
        DefaultAssay(seu_BA4_in) <- "RNA"
        DefaultAssay(seu_BA9_in) <- "RNA"

        seu_list <- list(seu_BA4_in = seu_BA4_in, seu_BA9_in = seu_BA9_in)

        seu_merge_in <- Merge_Seurat_List(
        list_seurat = seu_list,
        add.cell.ids = NULL,
        merge.data = TRUE,
        project = "MergeSeurat"
        )

        unique(seu_merge_in@meta.data$Region)
        unique(seu_merge_in@meta.data$Group)

        rm(seu_BA4_in)
        rm(seu_BA9_in)

        ## Run Harmony
        seu_merge_in <- ScaleData(seu_merge_in)    
        seu_merge_in<- FindVariableFeatures(seu_merge_in, selection.method = "vst", nfeatures = 2500)
        seu_merge_in <- RunPCA(seu_merge_in, npcs = 15)
        seu_merge_in <- RunHarmony(seu_merge_in, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_merge_in <- RunUMAP(seu_merge_in, reduction = "harmony", dims = 1:15)
        seu_merge_in <- FindNeighbors(seu_merge_in, reduction = "harmony", dims = 1:15)
        seu_merge_in <- FindClusters(seu_merge_in)

        ## By Celltype
        DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "PV" = "#e25822",
            "5HT3aR" = "#e68fac",
            "Rosehip" = "#0067a5",
            "SOM" = "#f99379"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'), height = 4, width = 4)
        ncol(seu_merge_in)

        saveRDS(seu_merge_in, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_BA4_BA9_Inhibitory_neurons.rds')
        rm(seu_BA9_in)
    ##

    #############################
    # Non-neuronal cells
    #############################
    ## code
        ## read in BA4
        seu_BA4_non <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_All_conditions_BA4.rds')
        ncol(seu_BA4_non)

        ## only retain Inhibitory class
        unique(seu_BA4_non@meta.data$Region) #"BA4"
        unique(seu_BA4_non@meta.data$CellType) #"BA4
        unique(seu_BA4_non@meta.data$annotation_major_cell_type) #"BA4
        
        ## Add major CellClass metadata
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L2_L3"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L3_L5"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L4_L5"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L4_L6"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L5"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L5_L6"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L6"] <- "Ex"

        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "5HT3aR"] <- "In"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "PV"] <- "In"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Rosehip"] <- "In"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "SOM"] <- "In"

        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Astro"] <- "Glia"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Endo"] <- "Vasc"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Fibro"] <- "Vasc"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Micro"] <- "Glia"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Mural"] <- "Vasc"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Oligo"] <- "Glia"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "OPC"] <- "Glia"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "T_Cell"] <- "Vasc"

        unique(seu_BA4_non@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

        seu_celltype <- c("Glia", "Vasc")
        xx <- unique(seu_BA4_non@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA4_non) <- "CellClass"
        seu_BA4_non=subset(seu_BA4_non,idents=xx)
        unique(seu_BA4_non@meta.data$CellClass)
        dim(seu_BA4_non)

        ## read in BA9
        seu_BA9_non <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_All_conditions_BA9.rds')
        ncol(seu_BA9_non)

        ## only retain Inhibitory class
        unique(seu_BA9_non@meta.data$Region) #"BA9"
        unique(seu_BA9_non@meta.data$CellType) #"BA9
        unique(seu_BA9_non@meta.data$annotation_major_cell_type) #"BA9
        
        ## Add major CellClass metadata
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "L2_L3"] <- "Ex"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "L3_L5"] <- "Ex"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "L4_L5"] <- "Ex"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "L4_L6"] <- "Ex"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "L5"] <- "Ex"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "L5_L6"] <- "Ex"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "L6"] <- "Ex"

        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "5HT3aR"] <- "In"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "PV"] <- "In"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "Rosehip"] <- "In"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "SOM"] <- "In"

        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "Astro"] <- "Glia"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "Endo"] <- "Vasc"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "Fibro"] <- "Vasc"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "Micro"] <- "Glia"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "Mural"] <- "Vasc"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "Oligo"] <- "Glia"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "OPC"] <- "Glia"
        seu_BA9_non@meta.data$CellClass[seu_BA9_non@meta.data$CellType == "T_Cell"] <- "Vasc"

        unique(seu_BA9_non@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

        seu_celltype <- c("Glia", "Vasc")
        xx <- unique(seu_BA9_non@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA9_non) <- "CellClass"
        seu_BA9_non=subset(seu_BA9_non,idents=xx)
        unique(seu_BA9_non@meta.data$CellClass)
        dim(seu_BA9_non)

        ## merge
        DefaultAssay(seu_BA4_non) <- "RNA"
        DefaultAssay(seu_BA9_non) <- "RNA"

        seu_list <- list(seu_BA4_non = seu_BA4_non, seu_BA9_non = seu_BA9_non)

        seu_merge_in <- Merge_Seurat_List(
        list_seurat = seu_list,
        add.cell.ids = NULL,
        merge.data = TRUE,
        project = "MergeSeurat"
        )

        unique(seu_merge_in@meta.data$Region)
        unique(seu_merge_in@meta.data$Group)

        rm(seu_BA4_non)
        rm(seu_BA9_non)

        ## Run Harmony
        seu_merge_in <- ScaleData(seu_merge_in)    
        seu_merge_in<- FindVariableFeatures(seu_merge_in, selection.method = "vst", nfeatures = 2500)
        seu_merge_in <- RunPCA(seu_merge_in, npcs = 15)
        seu_merge_in <- RunHarmony(seu_merge_in, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_merge_in <- RunUMAP(seu_merge_in, reduction = "harmony", dims = 1:15)
        seu_merge_in <- FindNeighbors(seu_merge_in, reduction = "harmony", dims = 1:15)
        seu_merge_in <- FindClusters(seu_merge_in)

        ## By Celltype
        DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "Oligo" = "#604e97",
            "Astro" = "#875692",
            "OPC" = "#f6a600",
            "Micro" = "#b3446c",
            "T_Cell" = "#dcd300",
            "Mural" = "#882d17",
            "Endo" = "#8db600",
            "Fibro" = "#654522"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        ncol(seu_merge_in)

        ## By Seurat cluster
        DimPlot(seu_merge_in, reduction = "umap", group.by = "seurat_clusters", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") 

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

        ## Cluster 15 did not match to any cell types from Pineda et al. we will remove it.
        seu_merge_in <- subset(seu_merge_in, seurat_clusters != 15)
        seu_merge_in <- subset(seu_merge_in, seurat_clusters != 16)

        seu_merge_in@meta.data$CellType[seu_merge_in@meta.data$seurat_clusters == 13] <- "Endo"
        seu_merge_in@meta.data$CellType[seu_merge_in@meta.data$seurat_clusters == 14] <- "Mural"

        ## By Celltype
        DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "Oligo" = "#604e97",
            "Astro" = "#875692",
            "OPC" = "#f6a600",
            "Micro" = "#b3446c",
            "T_Cell" = "#dcd300",
            "Mural" = "#882d17",
            "Endo" = "#8db600",
            "Fibro" = "#654522"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        ncol(seu_merge_in)

        saveRDS(seu_merge_in, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_BA4_BA9_non_neurons.rds')
        rm(seu_BA9_non)
    ##

##

################################################################## ## Create a single merged object to perform joint analyses -- this has BA4 and BA9 all cell types
##################################################################
################################################################## Performed this in Beluga. 
##################################################################
##################################################################
################################################################## Li

## code

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano merged_all_object_Li.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-grouleau
    #SBATCH --time=00-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=4500g          # memory per cor
    #SBATCH --job-name=merged_all_object_Li
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/merged_all_object_Li.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano merged_all_object_Li.R

    ## load libraries
    library(harmony)
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2)
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(scCustomize, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggrepel)
    library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")


    ##############
    ## Excitatory neurons
    ##############
    
    seu_ex_lim <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_BA4_BA9_Excitatory_neurons.rds') ## This object was created in the above chunk of this document
    ncol(seu_ex_lim)

    ##############
    ## Inhibitory neurons
    ##############

    seu_in_lim <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_BA4_BA9_Inhibitory_neurons.rds') ## This object was created in the above chunk of this document
    ncol(seu_in_lim)

    ##############
    ## non neurons
    ##############
    
    seu_non_lim <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_BA4_BA9_non_neurons.rds') ## This object was created in the above chunk of this document
    ncol(seu_non_lim)

    ##############
    ## merge
    ##############
    DefaultAssay(seu_ex_lim) <- "RNA"
    DefaultAssay(seu_in_lim) <- "RNA"
    DefaultAssay(seu_non_lim) <- "RNA"

    seu_list <- list(seu_ex_lim = seu_ex_lim, seu_in_lim = seu_in_lim, seu_non_lim = seu_non_lim )

    seu_merge_in <- Merge_Seurat_List(
    list_seurat = seu_list,
    add.cell.ids = NULL,
    merge.data = TRUE,
    project = "MergeSeurat"
    )

    unique(seu_merge_in@meta.data$Region)
    unique(seu_merge_in@meta.data$Condition)
    unique(seu_merge_in@meta.data$CellType)

    rm(seu_ex_lim)
    rm(seu_in_lim)
    rm(seu_non_lim)

    ## Run Harmony
    seu_merge_in <- ScaleData(seu_merge_in)    
    seu_merge_in<- FindVariableFeatures(seu_merge_in, selection.method = "vst", nfeatures = 2500)
    seu_merge_in <- RunPCA(seu_merge_in, npcs = 30)
    seu_merge_in <- RunHarmony(seu_merge_in, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_merge_in <- RunUMAP(seu_merge_in, reduction = "harmony", dims = 1:30)
    seu_merge_in <- FindNeighbors(seu_merge_in, reduction = "harmony", dims = 1:30)
    seu_merge_in <- FindClusters(seu_merge_in)

    ## By Celltype
    DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme_void()+
    theme(axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_blank()) + 
    xlab("UMAP1") + ylab("UMAP2") 

    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/UMAP_all_celltypes_temp.pdf'))
    ncol(seu_merge_in)

    ## Save the merged all object
    saveRDS(seu_merge_in, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_BA4_BA9_all_celltypes_lim_narval.rds')


##


################################################################## ## Merge seurat objects and print main figure UMAPS
##################################################################
################################################################## Performed this in Beluga. 
##################################################################
##################################################################
################################################################## Limone
## code
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano UMAPs_script_Li.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-grouleau
    #SBATCH --time=00-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=300g          # memory per cor
    #SBATCH --job-name=UMAPs_script_Li
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/UMAPs_script_Li.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano UMAPs_script_Li.R

    ## code
    library(harmony)
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2)
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(scCustomize, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggrepel)
    library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")


    #############################
    # Excitatory neurons
    #############################
    ## code
        ## read in BA4
        seu_BA4_ex <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_All_conditions_BA4.rds')
        ncol(seu_BA4_ex)

        ## only retain excitatory class
        unique(seu_BA4_ex@meta.data$Region) #"BA4"
        unique(seu_BA4_ex@meta.data$CellType) #"BA4
        
        ## Add major CellClass metadata
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L2_L3"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L3_L5"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L4_L5"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L4_L6"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L5"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L5_L6"] <- "Ex"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "L6"] <- "Ex"

        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "5HT3aR"] <- "In"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "PV"] <- "In"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Rosehip"] <- "In"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "SOM"] <- "In"

        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Astro"] <- "Glia"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Endo"] <- "Vasc"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Fibro"] <- "Vasc"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Micro"] <- "Glia"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Mural"] <- "Vasc"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "Oligo"] <- "Glia"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "OPC"] <- "Glia"
        seu_BA4_ex@meta.data$CellClass[seu_BA4_ex@meta.data$CellType == "T_Cell"] <- "Vasc"

        unique(seu_BA4_ex@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

        seu_celltype <- c("Ex")
        xx <- unique(seu_BA4_ex@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA4_ex) <- "CellClass"
        seu_BA4_ex=subset(seu_BA4_ex,idents=xx)
        unique(seu_BA4_ex@meta.data$CellClass)
        dim(seu_BA4_ex)

        seu_merge_ex <- seu_BA4_ex

        unique(seu_merge_ex@meta.data$Region)
        unique(seu_merge_ex@meta.data$Group)

        rm(seu_BA4_ex)

        ## Run Harmony
        seu_merge_ex <- ScaleData(seu_merge_ex)    
        seu_merge_ex<- FindVariableFeatures(seu_merge_ex, selection.method = "vst", nfeatures = 2500)
        seu_merge_ex <- RunPCA(seu_merge_ex, npcs = 35)
        seu_merge_ex <- RunHarmony(seu_merge_ex, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_merge_ex <- RunUMAP(seu_merge_ex, reduction = "harmony", dims = 1:35)
        seu_merge_ex <- FindNeighbors(seu_merge_ex, reduction = "harmony", dims = 1:35)
        seu_merge_ex <- FindClusters(seu_merge_ex)

        ## By Celltype
        DimPlot(seu_merge_ex, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "L3_L5" = "#f3c300",
            "L2_L3" = "#f38400",
            "L4_L6" = "#a1caf1",
            "L4_L5" = "#be0032",
            "L5_L6" = "#c2b280",
            "L5" = "#008856",
            "L6" = "#2b3d26"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        ncol(seu_merge_ex)

        saveRDS(seu_merge_ex, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_BA4_Excitatory_neurons.rds')
        rm(seu_merge_ex)
    ##

    #############################
    # Inhibitory neurons
    #############################
    ## code
        ## read in BA4
        seu_BA4_in <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_All_conditions_BA4.rds')
        ncol(seu_BA4_in)

        ## only retain Inhibitory class
        unique(seu_BA4_in@meta.data$Region) #"BA4"
        unique(seu_BA4_in@meta.data$CellType) #"BA4
        
        ## Add major CellClass metadata
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L2_L3"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L3_L5"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L4_L5"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L4_L6"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L5"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L5_L6"] <- "Ex"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "L6"] <- "Ex"

        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "5HT3aR"] <- "In"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "PV"] <- "In"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Rosehip"] <- "In"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "SOM"] <- "In"

        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Astro"] <- "Glia"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Endo"] <- "Vasc"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Fibro"] <- "Vasc"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Micro"] <- "Glia"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Mural"] <- "Vasc"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "Oligo"] <- "Glia"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "OPC"] <- "Glia"
        seu_BA4_in@meta.data$CellClass[seu_BA4_in@meta.data$CellType == "T_Cell"] <- "Vasc"

        unique(seu_BA4_in@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

        seu_celltype <- c("In")
        xx <- unique(seu_BA4_in@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA4_in) <- "CellClass"
        seu_BA4_in=subset(seu_BA4_in,idents=xx)
        unique(seu_BA4_in@meta.data$CellClass)
        dim(seu_BA4_in)

        seu_merge_in <- seu_BA4_in

        unique(seu_merge_in@meta.data$Region)
        unique(seu_merge_in@meta.data$Group)

        rm(seu_BA4_in)

        ## Run Harmony
        seu_merge_in <- ScaleData(seu_merge_in)    
        seu_merge_in<- FindVariableFeatures(seu_merge_in, selection.method = "vst", nfeatures = 2500)
        seu_merge_in <- RunPCA(seu_merge_in, npcs = 35)
        seu_merge_in <- RunHarmony(seu_merge_in, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_merge_in <- RunUMAP(seu_merge_in, reduction = "harmony", dims = 1:35)
        seu_merge_in <- FindNeighbors(seu_merge_in, reduction = "harmony", dims = 1:35)
        seu_merge_in <- FindClusters(seu_merge_in)

        ## By Celltype
        DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "PV" = "#e25822",
            "5HT3aR" = "#e68fac",
            "Rosehip" = "#0067a5",
            "SOM" = "#f99379"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        ncol(seu_merge_in)

        saveRDS(seu_merge_in, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_BA4_Inhibitory_neurons.rds')
        rm(seu_merge_in)
    ##

    #############################
    # Non-neuronal cells
    #############################
    ## code
        ## read in BA4
        seu_BA4_non <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_All_conditions_BA4.rds')
        ncol(seu_BA4_non)

        ## only retain Inhibitory class
        unique(seu_BA4_non@meta.data$Region) #"BA4"
        unique(seu_BA4_non@meta.data$CellType) #"BA4
        
        ## Add major CellClass metadata
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L2_L3"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L3_L5"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L4_L5"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L4_L6"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L5"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L5_L6"] <- "Ex"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "L6"] <- "Ex"

        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "5HT3aR"] <- "In"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "PV"] <- "In"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Rosehip"] <- "In"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "SOM"] <- "In"

        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Astro"] <- "Glia"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Endo"] <- "Vasc"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Fibro"] <- "Vasc"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Micro"] <- "Glia"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Mural"] <- "Vasc"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "Oligo"] <- "Glia"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "OPC"] <- "Glia"
        seu_BA4_non@meta.data$CellClass[seu_BA4_non@meta.data$CellType == "T_Cell"] <- "Vasc"

        unique(seu_BA4_non@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

        seu_celltype <- c("Glia", "Vasc")
        xx <- unique(seu_BA4_non@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu_BA4_non) <- "CellClass"
        seu_BA4_non=subset(seu_BA4_non,idents=xx)
        unique(seu_BA4_non@meta.data$CellClass)
        dim(seu_BA4_non)

        seu_merge_in <- seu_BA4_non

        unique(seu_merge_in@meta.data$Region)
        unique(seu_merge_in@meta.data$Group)
        unique(seu_merge_in@meta.data$CellClass)

        rm(seu_BA4_non)

        ## Run Harmony
        seu_merge_in <- ScaleData(seu_merge_in)    
        seu_merge_in<- FindVariableFeatures(seu_merge_in, selection.method = "vst", nfeatures = 2500)
        seu_merge_in <- RunPCA(seu_merge_in, npcs = 35)
        seu_merge_in <- RunHarmony(seu_merge_in, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_merge_in <- RunUMAP(seu_merge_in, reduction = "harmony", dims = 1:35)
        seu_merge_in <- FindNeighbors(seu_merge_in, reduction = "harmony", dims = 1:35)
        seu_merge_in <- FindClusters(seu_merge_in)

        ## By Celltype
        DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme_void()+
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "Oligo" = "#604e97",
            "Astro" = "#875692",
            "OPC" = "#f6a600",
            "Micro" = "#b3446c",
            "T_Cell" = "#dcd300",
            "Mural" = "#882d17",
            "Endo" = "#8db600",
            "Fibro" = "#654522"
            ))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
        ncol(seu_merge_in)

        saveRDS(seu_merge_in, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_BA4_non_neurons.rds')
        rm(seu_BA9_non)
    ##

##

################################################################## ## Create a single merged object to perform joint analyses -- this has BA4 and BA9 all cell types
##################################################################
################################################################## Performed this in Beluga. 
##################################################################
##################################################################
################################################################## Limone

## code

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    nano merged_all_object_Li.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #!/bin/bash  
    #SBATCH --account=def-grouleau
    #SBATCH --time=00-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=4500g          # memory per cor
    #SBATCH --job-name=merged_all_object_Li
    #SBATCH --error=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.err
    #SBATCH --output=/home/fiorini9/scratch/machine_learning_ALS/scripts/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning_ALS/scripts

    Rscript /home/fiorini9/scratch/machine_learning_ALS/scripts/merged_all_object_Li.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano merged_all_object_Li.R

    ## load libraries
    library(harmony)
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2)
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(scCustomize, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggrepel)
    library(RColorBrewer, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")


    ##############
    ## Excitatory neurons
    ##############
    
    seu_ex_lim <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_BA4_Excitatory_neurons.rds') ## This object was created in the above chunk of this document
    ncol(seu_ex_lim)
    unique(seu_ex_lim@meta.data$CellType)

    ##############
    ## Inhibitory neurons
    ##############

    seu_in_lim <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_BA4_Inhibitory_neurons.rds') ## This object was created in the above chunk of this document
    ncol(seu_in_lim)
    unique(seu_in_lim@meta.data$CellType)

    ##############
    ## non neurons
    ##############
    
    seu_non_lim <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_BA4_non_neurons.rds') ## This object was created in the above chunk of this document
    ncol(seu_non_lim)
    unique(seu_non_lim@meta.data$CellType)

    ##############
    ## merge
    ##############
    DefaultAssay(seu_ex_lim) <- "RNA"
    DefaultAssay(seu_in_lim) <- "RNA"
    DefaultAssay(seu_non_lim) <- "RNA"

    seu_list <- list(seu_ex_lim = seu_ex_lim, seu_in_lim = seu_in_lim, seu_non_lim = seu_non_lim )

    seu_merge_in <- Merge_Seurat_List(
    list_seurat = seu_list,
    add.cell.ids = NULL,
    merge.data = TRUE,
    project = "MergeSeurat"
    )

    unique(seu_merge_in@meta.data$Region)
    unique(seu_merge_in@meta.data$Group)
    unique(seu_merge_in@meta.data$CellType)

    rm(seu_ex_lim)
    rm(seu_in_lim)
    rm(seu_non_lim)

    ## Run Harmony
    seu_merge_in <- ScaleData(seu_merge_in)    
    seu_merge_in<- FindVariableFeatures(seu_merge_in, selection.method = "vst", nfeatures = 2500)
    seu_merge_in <- RunPCA(seu_merge_in, npcs = 35)
    seu_merge_in <- RunHarmony(seu_merge_in, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_merge_in <- RunUMAP(seu_merge_in, reduction = "harmony", dims = 1:35)
    seu_merge_in <- FindNeighbors(seu_merge_in, reduction = "harmony", dims = 1:35)
    seu_merge_in <- FindClusters(seu_merge_in)

    ## By Celltype
    DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme_void()+
    theme(axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    legend.position = "none",
    plot.title = element_blank()) + 
    xlab("UMAP1") + ylab("UMAP2") 

    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/UMAP_all_celltypes_temp.pdf'))
    ncol(seu_merge_in)

    ## Save the merged all object
    saveRDS(seu_merge_in, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_BA4_all_celltypes_lim_narval.rds')


##



###################################################################################################################################################################################################### Narval code
######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################


## transfer objects from beluga to Narval


## run this in Narval

salloc -A def-grouleau --time=0-4 -c 1 --mem=300g

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
library(scCustomize, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(RColorBrewer, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(viridis, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(SingleCellExperiment, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(MetaNeighbor, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
library(data.table, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")

#install.packages("data.table", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#library(BiocManager, lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#BiocManager::install("SingleCellExperiment", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")
#BiocManager::install("MetaNeighbor", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R")




################################################################## ## Marker gene dotplot 
##################################################################
##################################################################
##################################################################
##################################################################
################################################################## Pineda

## Code 
    ## Notes from pineda manuscript: # Taxonomy of annotated subtypes. Values represent average log-transformed marker expression normalized to a maximum per column.
    
    seu_ex <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_all_celltypes_lim_narval.rds') ## This object was created in the above chunk of this document

    ## features to plot
    par_select_features_list= c("SLC17A7", "GAD1", "CUX2", "THEMIS", "RORB", "FEZF2", "TLE4", "CCN2", "SST", "NPY", "PVALB", "LAMP5", "CNR1", "RELN", 'VIP', 'CALB2', 'AQP4', 'GFAP', 'PLP1', 'VCAN', 'CTSS', 'CD8A', 'CLDN5', 'PDGFRB', 'BICC1' )

    ## Kam 
    DefaultAssay(seu_ex) <- "RNA"

    ## dotplot to get scaled expression
    dot_plt <- DotPlot(seu_ex, features = par_select_features_list, group.by = 'CellType')
    
    dot_plt2 <- data.frame(dot_plt$data)
    dot_plt2$id <- factor(dot_plt2$id, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
    dot_plt2$features.plot <- factor(dot_plt2$features.plot, levels = c("SLC17A7", "GAD1", "CUX2", "THEMIS", "RORB", "FEZF2", "TLE4", "CCN2", "SST", "NPY", "PVALB", "LAMP5", "CNR1", "RELN", 'VIP', 'CALB2', 'AQP4', 'GFAP', 'PLP1', 'VCAN', 'CTSS', 'CD8A', 'CLDN5', 'PDGFRB', 'BICC1'))
    
    ## Add major group labels so that we can subset. 
    dot_plt2$major_group[dot_plt2$id == "L2_L3"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L3_L5"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L4_L5"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L4_L6"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L5"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L5_L6"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L6"] <- "Excitatory\nneurons"

    dot_plt2$major_group[dot_plt2$id == "5HT3aR"] <- "Inhibitory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "PV"] <- "Inhibitory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "Rosehip"] <- "Inhibitory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "SOM"] <- "Inhibitory\nneurons"

    dot_plt2$major_group[dot_plt2$id == "Astro"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Endo"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Fibro"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Micro"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Mural"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Oligo"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "OPC"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "T_Cell"] <- "Non-neuronal\ncells"
    
    
    ## heatmap
    kam_plot <- ggplot(dot_plt2, aes(x = features.plot , y = id, fill = avg.exp.scaled)) + 
    theme_bw() + 
    geom_tile() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, colour = "black", face = "italic"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = "black"),
        axis.title = element_text(face="bold", size =12),
        legend.text = element_text( size =10),
        legend.title = element_text( size= 10),
        legend.position = "bottom"
    ) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_viridis() +
    ylab("") + xlab("") 
    ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf' ,sep=""),width = 4, height = 4)
##


################################################################## ## Marker gene dotplot 
##################################################################
##################################################################
##################################################################
##################################################################
################################################################## Li

## Code 
    ## Notes from pineda manuscript: # Taxonomy of annotated subtypes. Values represent average log-transformed marker expression normalized to a maximum per column.
    
    seu_ex <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_BA4_BA9_all_celltypes_lim_narval.rds') ## This object was created in the above chunk of this document

    ## features to plot
    par_select_features_list= c("SLC17A7", "GAD1", "CUX2", "THEMIS", "RORB", "FEZF2", "TLE4", "CCN2", "SST", "NPY", "PVALB", "LAMP5", "CNR1", "RELN", 'VIP', 'CALB2', 'AQP4', 'GFAP', 'PLP1', 'VCAN', 'CTSS', 'CD8A', 'CLDN5', 'PDGFRB', 'BICC1' )

    ## Kam 
    DefaultAssay(seu_ex) <- "RNA"

    ## dotplot to get scaled expression
    dot_plt <- DotPlot(seu_ex, features = par_select_features_list, group.by = 'CellType')
    
    dot_plt2 <- data.frame(dot_plt$data)
    dot_plt2$id <- factor(dot_plt2$id, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
    dot_plt2$features.plot <- factor(dot_plt2$features.plot, levels = c("SLC17A7", "GAD1", "CUX2", "THEMIS", "RORB", "FEZF2", "TLE4", "CCN2", "SST", "NPY", "PVALB", "LAMP5", "CNR1", "RELN", 'VIP', 'CALB2', 'AQP4', 'GFAP', 'PLP1', 'VCAN', 'CTSS', 'CD8A', 'CLDN5', 'PDGFRB', 'BICC1'))
    
    ## Add major group labels so that we can subset. 
    dot_plt2$major_group[dot_plt2$id == "L2_L3"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L3_L5"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L4_L5"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L4_L6"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L5"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L5_L6"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L6"] <- "Excitatory\nneurons"

    dot_plt2$major_group[dot_plt2$id == "5HT3aR"] <- "Inhibitory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "PV"] <- "Inhibitory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "Rosehip"] <- "Inhibitory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "SOM"] <- "Inhibitory\nneurons"

    dot_plt2$major_group[dot_plt2$id == "Astro"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Endo"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Fibro"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Micro"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Mural"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Oligo"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "OPC"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "T_Cell"] <- "Non-neuronal\ncells"
    
    
    ## heatmap
    kam_plot <- ggplot(dot_plt2, aes(x = features.plot , y = id, fill = avg.exp.scaled)) + 
    theme_bw() + 
    geom_tile() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, colour = "black", face = "italic"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = "black"),
        axis.title = element_text(face="bold", size =12),
        legend.text = element_text( size =10),
        legend.title = element_text( size= 10),
        legend.position = "bottom"
    ) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_viridis() +
    ylab("") + xlab("") 
    #facet_grid(major_group~., scales = "free_y")
    ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf' ,sep=""),width = 4, height = 4)
##

################################################################## ## Marker gene dotplot 
##################################################################
##################################################################
##################################################################
##################################################################
################################################################## Limone

## Code 
    ## Notes from pineda manuscript: # Taxonomy of annotated subtypes. Values represent average log-transformed marker expression normalized to a maximum per column.
    
    seu_ex <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_BA4_all_celltypes_lim_narval.rds') ## This object was created in the above chunk of this document

    ## features to plot
    par_select_features_list= c("SLC17A7", "GAD1", "CUX2", "THEMIS", "RORB", "FEZF2", "TLE4", "CCN2", "SST", "NPY", "PVALB", "LAMP5", "CNR1", "RELN", 'VIP', 'CALB2', 'AQP4', 'GFAP', 'PLP1', 'VCAN', 'CTSS', 'CD8A', 'CLDN5', 'PDGFRB', 'BICC1' )

    ## Kam 
    DefaultAssay(seu_ex) <- "RNA"

    ## dotplot to get scaled expression
    dot_plt <- DotPlot(seu_ex, features = par_select_features_list, group.by = 'CellType')
    
    dot_plt2 <- data.frame(dot_plt$data)
    dot_plt2$id <- factor(dot_plt2$id, levels = rev(c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell'  )))
    dot_plt2$features.plot <- factor(dot_plt2$features.plot, levels = c("SLC17A7", "GAD1", "CUX2", "THEMIS", "RORB", "FEZF2", "TLE4", "CCN2", "SST", "NPY", "PVALB", "LAMP5", "CNR1", "RELN", 'VIP', 'CALB2', 'AQP4', 'GFAP', 'PLP1', 'VCAN', 'CTSS', 'CD8A', 'CLDN5', 'PDGFRB', 'BICC1'))
    
    ## Add major group labels so that we can subset. 
    dot_plt2$major_group[dot_plt2$id == "L2_L3"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L3_L5"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L4_L5"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L4_L6"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L5"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L5_L6"] <- "Excitatory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "L6"] <- "Excitatory\nneurons"

    dot_plt2$major_group[dot_plt2$id == "5HT3aR"] <- "Inhibitory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "PV"] <- "Inhibitory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "Rosehip"] <- "Inhibitory\nneurons"
    dot_plt2$major_group[dot_plt2$id == "SOM"] <- "Inhibitory\nneurons"

    dot_plt2$major_group[dot_plt2$id == "Astro"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Endo"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Fibro"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Micro"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Mural"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "Oligo"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "OPC"] <- "Non-neuronal\ncells"
    dot_plt2$major_group[dot_plt2$id == "T_Cell"] <- "Non-neuronal\ncells"
    
    
    ## heatmap
    kam_plot <- ggplot(dot_plt2, aes(x = features.plot , y = id, fill = avg.exp.scaled)) + 
    theme_bw() + 
    geom_tile() +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5, colour = "black", face = "italic"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = "black"),
        axis.title = element_text(face="bold", size =12),
        legend.text = element_text( size =10),
        legend.title = element_text( size= 10),
        legend.position = "bottom"
    ) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    scale_fill_viridis() +
    ylab("") + xlab("") 
    #facet_grid(major_group~., scales = "free_y")
    ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf' ,sep=""),width = 4, height = 4)
##


################################################################## ## Metaneighbour
##################################################################
##################################################################
##################################################################
##################################################################
################################################################## Pineda/Li/Limone (enventually)

## Code 
    ## Notes from pineda manuscript: # Taxonomy of annotated subtypes. Values represent average log-transformed marker expression normalized to a maximum per column.
    
    ## Import Pineda and subset down to 100,000 cells
    seu_pineda <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_all_celltypes_lim_narval.rds') ## This object was created in the above chunk of this document
    set.seed(42)  # Set seed for reproducibility (optional)
    cells_to_keep <- sample(colnames(seu_pineda), 100000)
    seu_pineda <- subset(seu_pineda, cells = cells_to_keep)
    ncol(seu_pineda)
    DefaultAssay(seu_pineda) <- "RNA"
    
    ## Import Li
    seu_li <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_BA4_BA9_all_celltypes_lim_narval.rds') ## This object was created in the above chunk of this document
    DefaultAssay(seu_li) <- "RNA"

    ## Import Limone
    seu_limone <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_BA4_all_celltypes_lim_narval.rds') ## This object was created in the above chunk of this document
    DefaultAssay(seu_limone) <- "RNA"

    ##########################################
    ## Pineda and Li
    ##########################################
    ## code
        ## train model with kam
        seu_pineda <- FindVariableFeatures(seu_pineda)

        var_features <- VariableFeatures(seu_pineda)
        Idents(seu_pineda) <- seu_pineda$CellType
        seu_pineda_SCE <- as.SingleCellExperiment(seu_pineda)
        seu_pineda_SCE <- seu_pineda_SCE[var_features,]

        pretrained_model_major_cluster = MetaNeighbor::trainModel(
        var_genes = rownames(seu_pineda_SCE),
        dat = seu_pineda_SCE,
        study_id = rep("Lister", dim(seu_pineda_SCE)[2]),
        cell_type = seu_pineda_SCE$CellType)

        aurocs = MetaNeighborUS(
        trained_model = pretrained_model_major_cluster, dat = as.SingleCellExperiment(seu_li),
        study_id = rep("integrated", dim(seu_li)[2]), 
        cell_type = seu_li@meta.data$CellType,
        fast_version = TRUE
        )

        ## aurocs
        aurocs_lim <- aurocs

        # Get lower triangle of the correlation matrix
        get_lower_tri<-function(cormat){
            cormat[upper.tri(cormat)] <- NA
            return(cormat)
        }
        # Get upper triangle of the correlation matrix
        get_upper_tri <- function(cormat){
            cormat[lower.tri(cormat)]<- NA
            return(cormat)
        }

        order <- rev(c('Lister|L2_L3','Lister|L3_L5','Lister|L4_L5','Lister|L4_L6', 'Lister|L5', 'Lister|L5_L6', 'Lister|L6', "Lister|5HT3aR", 'Lister|PV', 'Lister|Rosehip',  'Lister|SOM', 'Lister|Astro', 'Lister|Endo', 'Lister|Fibro', 'Lister|Micro', 'Lister|Mural', 'Lister|Oligo', 'Lister|OPC','Lister|T_Cell' ))
        aurocs_lim_ordered <- aurocs_lim[,rev(order)]

        order <- rev(c( 'integrated|L2_L3','integrated|L3_L5','integrated|L4_L5','integrated|L4_L6', 'integrated|L5', 'integrated|L5_L6', 'integrated|L6', "integrated|5HT3aR", 'integrated|PV', 'integrated|Rosehip',  'integrated|SOM', 'integrated|Astro', 'integrated|Endo', 'integrated|Fibro', 'integrated|Micro', 'integrated|Mural', 'integrated|Oligo', 'integrated|OPC','integrated|T_Cell'))
        aurocs_lim_ordered <- aurocs_lim_ordered[rev(order),]


        upper_tri <- get_upper_tri(aurocs_lim_ordered)

        melted_cormat <- melt(upper_tri, na.rm = TRUE)

        blups <- brewer.pal(9, "Reds")

        melted_cormat$Var1_new <- sub("^[^|]*\\|", "", melted_cormat$Var1)
        melted_cormat$Var2_new  <- sub("^[^|]*\\|", "", melted_cormat$Var2)
        
        ## reset factor levels
        melted_cormat$Var1_new  <- factor(melted_cormat$Var1_new , levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell' ))
        melted_cormat$Var2_new  <- factor(melted_cormat$Var2_new , levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell' ))

        ## add labels
        melted_cormat$new_value <- ifelse(melted_cormat$Var1_new == melted_cormat$Var2_new, melted_cormat$value, NA)


        # Heatmap
        ggplot(data = melted_cormat, aes(Var2_new, Var1_new, fill = value, labels = new_value))+
        theme_bw() +
        geom_tile(color = "white")+
        geom_text(aes(label = round(new_value, 2)), size = 2, colour = "white") +
        theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
                axis.text.y = element_text( colour = "black"),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                legend.position = "bottom"
            )+
        #coord_fixed() +
        scale_fill_gradientn(colors = blups) + #, values = scales::rescale(c(-1,  0, 3)), na.value = "white")+
        scale_y_discrete(position = "right") +
        labs(fill = "MetaNeighbour AUC") +  # Title for the color scale
        guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, title.theme = element_text(face = "bold")))
        
        ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf',sep=""),width = 4, height = 4.5)
    ##

    ##########################################
    ## Pineda and Limone
    ##########################################
    ## code
        ## train model with kam
        seu_pineda <- FindVariableFeatures(seu_pineda)

        var_features <- VariableFeatures(seu_pineda)
        Idents(seu_pineda) <- seu_pineda$CellType
        seu_pineda_SCE <- as.SingleCellExperiment(seu_pineda)
        seu_pineda_SCE <- seu_pineda_SCE[var_features,]

        pretrained_model_major_cluster = MetaNeighbor::trainModel(
        var_genes = rownames(seu_pineda_SCE),
        dat = seu_pineda_SCE,
        study_id = rep("Lister", dim(seu_pineda_SCE)[2]),
        cell_type = seu_pineda_SCE$CellType)

        aurocs = MetaNeighborUS(
        trained_model = pretrained_model_major_cluster, dat = as.SingleCellExperiment(seu_limone),
        study_id = rep("integrated", dim(seu_limone)[2]), 
        cell_type = seu_limone@meta.data$CellType,
        fast_version = TRUE
        )

        ## aurocs
        aurocs_lim <- aurocs

        # Get lower triangle of the correlation matrix
        get_lower_tri<-function(cormat){
            cormat[upper.tri(cormat)] <- NA
            return(cormat)
        }
        # Get upper triangle of the correlation matrix
        get_upper_tri <- function(cormat){
            cormat[lower.tri(cormat)]<- NA
            return(cormat)
        }

        order <- rev(c('Lister|L2_L3','Lister|L3_L5','Lister|L4_L5','Lister|L4_L6', 'Lister|L5', 'Lister|L5_L6', 'Lister|L6', "Lister|5HT3aR", 'Lister|PV', 'Lister|Rosehip',  'Lister|SOM', 'Lister|Astro', 'Lister|Endo', 'Lister|Fibro', 'Lister|Micro', 'Lister|Mural', 'Lister|Oligo', 'Lister|OPC','Lister|T_Cell' ))
        aurocs_lim_ordered <- aurocs_lim[,rev(order)]

        order <- rev(c( 'integrated|L2_L3','integrated|L3_L5','integrated|L4_L5','integrated|L4_L6', 'integrated|L5', 'integrated|L5_L6', 'integrated|L6', "integrated|5HT3aR", 'integrated|PV', 'integrated|Rosehip',  'integrated|SOM', 'integrated|Astro', 'integrated|Endo', 'integrated|Fibro', 'integrated|Micro', 'integrated|Mural', 'integrated|Oligo', 'integrated|OPC','integrated|T_Cell'))
        aurocs_lim_ordered <- aurocs_lim_ordered[rev(order),]


        upper_tri <- get_upper_tri(aurocs_lim_ordered)

        melted_cormat <- melt(upper_tri, na.rm = TRUE)

        blups <- brewer.pal(9, "Reds")

        melted_cormat$Var1_new <- sub("^[^|]*\\|", "", melted_cormat$Var1)
        melted_cormat$Var2_new  <- sub("^[^|]*\\|", "", melted_cormat$Var2)
        
        ## reset factor levels
        melted_cormat$Var1_new  <- factor(melted_cormat$Var1_new , levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell' ))
        melted_cormat$Var2_new  <- factor(melted_cormat$Var2_new , levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell' ))

        ## add labels
        melted_cormat$new_value <- ifelse(melted_cormat$Var1_new == melted_cormat$Var2_new, melted_cormat$value, NA)


        # Heatmap
        ggplot(data = melted_cormat, aes(Var2_new, Var1_new, fill = value, labels = new_value))+
        theme_bw() +
        geom_tile(color = "white")+
        geom_text(aes(label = round(new_value, 2)), size = 2, colour = "white") +
        theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
                axis.text.y = element_text( colour = "black"),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                legend.position = "bottom"
            )+
        #coord_fixed() +
        scale_fill_gradientn(colors = blups) + #, values = scales::rescale(c(-1,  0, 3)), na.value = "white")+
        scale_y_discrete(position = "right") +
        labs(fill = "MetaNeighbour AUC") +  # Title for the color scale
        guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, title.theme = element_text(face = "bold")))
        
        ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf',sep=""),width = 4, height = 4.5)
    ##

    ##########################################
    ## Li and Limone
    ##########################################
    ## code
        ## train model with kam
        seu_limone <- FindVariableFeatures(seu_limone)

        var_features <- VariableFeatures(seu_limone)
        Idents(seu_limone) <- seu_limone$CellType
        seu_limone_SCE <- as.SingleCellExperiment(seu_limone)
        seu_limone_SCE <- seu_limone_SCE[var_features,]

        pretrained_model_major_cluster = MetaNeighbor::trainModel(
        var_genes = rownames(seu_limone_SCE),
        dat = seu_limone_SCE,
        study_id = rep("Lister", dim(seu_limone_SCE)[2]),
        cell_type = seu_limone_SCE$CellType)

        aurocs = MetaNeighborUS(
        trained_model = pretrained_model_major_cluster, dat = as.SingleCellExperiment(seu_li),
        study_id = rep("integrated", dim(seu_li)[2]), 
        cell_type = seu_li@meta.data$CellType,
        fast_version = TRUE
        )

        ## aurocs
        aurocs_lim <- aurocs

        # Get lower triangle of the correlation matrix
        get_lower_tri<-function(cormat){
            cormat[upper.tri(cormat)] <- NA
            return(cormat)
        }
        # Get upper triangle of the correlation matrix
        get_upper_tri <- function(cormat){
            cormat[lower.tri(cormat)]<- NA
            return(cormat)
        }

        order <- rev(c('Lister|L2_L3','Lister|L3_L5','Lister|L4_L5','Lister|L4_L6', 'Lister|L5', 'Lister|L5_L6', 'Lister|L6', "Lister|5HT3aR", 'Lister|PV', 'Lister|Rosehip',  'Lister|SOM', 'Lister|Astro', 'Lister|Endo', 'Lister|Fibro', 'Lister|Micro', 'Lister|Mural', 'Lister|Oligo', 'Lister|OPC','Lister|T_Cell' ))
        aurocs_lim_ordered <- aurocs_lim[,rev(order)]

        order <- rev(c( 'integrated|L2_L3','integrated|L3_L5','integrated|L4_L5','integrated|L4_L6', 'integrated|L5', 'integrated|L5_L6', 'integrated|L6', "integrated|5HT3aR", 'integrated|PV', 'integrated|Rosehip',  'integrated|SOM', 'integrated|Astro', 'integrated|Endo', 'integrated|Fibro', 'integrated|Micro', 'integrated|Mural', 'integrated|Oligo', 'integrated|OPC','integrated|T_Cell'))
        aurocs_lim_ordered <- aurocs_lim_ordered[rev(order),]
        
                
        upper_tri <- get_upper_tri(aurocs_lim_ordered)

        melted_cormat <- melt(upper_tri, na.rm = TRUE)

        blups <- brewer.pal(9, "Reds")

        melted_cormat$Var1_new <- sub("^[^|]*\\|", "", melted_cormat$Var1)
        melted_cormat$Var2_new  <- sub("^[^|]*\\|", "", melted_cormat$Var2)
        
        ## reset factor levels
        melted_cormat$Var1_new  <- factor(melted_cormat$Var1_new , levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell' ))
        melted_cormat$Var2_new  <- factor(melted_cormat$Var2_new , levels = c('L2_L3','L3_L5','L4_L5','L4_L6', 'L5', 'L5_L6', 'L6', "5HT3aR", 'PV', 'Rosehip',  'SOM', 'Astro', 'Endo', 'Fibro', 'Micro', 'Mural', 'Oligo', 'OPC','T_Cell' ))

        ## add labels
        melted_cormat$new_value <- ifelse(melted_cormat$Var1_new == melted_cormat$Var2_new, melted_cormat$value, NA)


        # Heatmap
        ggplot(data = melted_cormat, aes(Var2_new, Var1_new, fill = value, labels = new_value))+
        theme_bw() +
        geom_tile(color = "white")+
        geom_text(aes(label = round(new_value, 2)), size = 2, colour = "white") +
        theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
                axis.text.y = element_text( colour = "black"),
                axis.title = element_blank(),
                panel.grid = element_blank(),
                legend.position = "bottom"
            )+
        #coord_fixed() +
        scale_fill_gradientn(colors = blups) + #, values = scales::rescale(c(-1,  0, 3)), na.value = "white")+
        scale_y_discrete(position = "right") +
        labs(fill = "MetaNeighbour AUC") +  # Title for the color scale
        guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5, title.theme = element_text(face = "bold")))
        
        ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf',sep=""),width = 4, height = 4.5)
    ##
##












########### OLD CODE


################################################################## ## Print UMAPs in Narval
##################################################################
##################################################################
##################################################################
##################################################################
################################################################## Pineda

## code
    salloc -A def-sfarhan --time=0-4 -c 1 --mem=300g

    module load StdEnv/2023
    module load r/4.4.0
    R

    ## code
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

        #############################
        # Excitatory neurons
        #############################
        ## print UMAP
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Pineda_BA4_BA9_Excitatory_neurons.rds')
        unique(seu@meta.data$CellType)
        
        ## By Cell type
        DimPlot(seu_merge_in, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 2, repel = TRUE, raster=FALSE) +
        theme_void() +
        theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.title = element_blank()) + 
        xlab("UMAP1") + ylab("UMAP2") +
        scale_color_manual(values = c(
            "Glia.Oligo" = "#f3c300",
            "Glia.Astro.GFAP-neg"="#875692",
            "Glia.Astro.GFAP-pos"="#f38400",
            "Glia.OPC"="#a1caf1",
            "Glia.Micro"="#be0032",
            "Vasc.T_Cell"="#c2b280",
            "Vasc.Mural.Pericyte"="#848482",
            "Vasc.Endo.Capillary"="#008856",
            "Vasc.Fibro.CLMP_PDGFRA"="#e68fac",
            "Vasc.Endo.Venous"="#0067a5",
            "Vasc.Mural.SMC"="#f99379",
            "Vasc.Endo.Arterial"="#604e97",
            "Vasc.Fibro.CLMP_KCNMA1"="#f6a600"))

        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'), height = 5, width = 5)
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Pineda_Ex_neuron_UMAP.pdf'), height = 5, width = 5)

        ## Cell type marker genes
##






################################################################## ## Merge with Harmony in beluga
##################################################################
##################################################################
##################################################################
##################################################################
################################################################## Li










######################################################################################################## OLD code. 

#install.packages("scCustomize", lib="/home/fiorini9/projects/def-sfarhan/fiorini9/software/CBIG_R" )

################################################################## ## TXd test
##################################################################
##################################################################
##################################################################
##################################################################
################################################################## Pineda


## TXd temp
    #"The TxD score is a quantification of transcriptional dysregulation. 
    #It represents the change in transcriptome-wide gene expression of each 
    #subtype in disease from its respective PN expression profile. 
    #The divergence score is the Euclidean distance between the median disease 
    #and corresponding PN covariate-corrected, pseudo-bulk expression profiles for each cell type."   
    
    
    ## we can maybe use these (covariat carrected) objects for TXd. 
    #'/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_BA9_',seu_celltype,'_int.rds'

    ##############
    ## Excitatory neurons BA4
    ##############
    
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_BA4_Ex_int.rds')
    ncol(seu)

    
    ## log normalize
    #seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
    #seu <- ScaleData(seu)   

    # Assuming your Seurat object is called `seurat_obj`
    # Subset the Seurat object for disease and PN groups
    seurat_obj_disease <- subset(seu, Group == "C9ALS")
    seurat_obj_pn <- subset(seu, Group == "PN")

    # List of cell types (assuming 'celltype' is in the metadata)
    celltypes <- unique(seu$CellType)

    # Create a function to calculate the pseudo-bulk expression profile (median of gene expression per cell type)
    calculate_pseudo_bulk <- function(seurat_obj, celltype) {
        # Subset the Seurat object by celltype
        #celltype_data <- subset(seurat_obj, subset = celltype == celltype)
        celltype_data <- subset(seurat_obj, CellType == celltype)
        
        # Calculate median expression for each gene in the cell type
        
        #pseudo_bulk_profile <- AverageExpression(celltype_data, return.seurat = TRUE, verbose = FALSE, group.by = "Sample_ID")
        #return(pseudo_bulk_profile)

        pseudo_bulk_profile <- AggregateExpression(celltype_data, return.seurat = TRUE, verbose = FALSE, group.by = "Sample_ID")
        return(pseudo_bulk_profile)
    }

   
                ##################### test
                    #celltype = "Astro"
                    #celltype_data <- subset(seurat_obj_disease, CellType == celltype)
                    
                    # Calculate median expression for each gene in the cell type
                    #pseudo_bulk_profile <- AverageExpression(celltype_data, return.seurat = TRUE, verbose = FALSE, group.by = "Sample_ID")
                    #return(pseudo_bulk_profile)
                ##################### test
   

    # Create empty lists to store pseudo-bulk profiles
    pseudo_bulk_disease <- list()
    pseudo_bulk_pn <- list()

    # Loop through cell types and calculate pseudo-bulk profiles for both disease and PN
    for (celltype in celltypes) {
        pseudo_bulk_disease[[celltype]] <- calculate_pseudo_bulk(seurat_obj_disease, celltype)
        pseudo_bulk_pn[[celltype]] <- calculate_pseudo_bulk(seurat_obj_pn, celltype)
    }

                # Calculate TxD score (change in gene expression)
                #txd_scores <- list()

                #for (celltype in celltypes) {
                #    disease_expr <- pseudo_bulk_disease[[celltype]]@assays$RNA@layers$data ## maybe we use scaled data????????????????????
                #    pn_expr <- pseudo_bulk_pn[[celltype]]@assays$RNA@layers$data
                    
                #    # TxD score is the absolute difference in gene expression between disease and PN
                #    txd_scores[[celltype]] <- abs(disease_expr - pn_expr)
                #}

    # Function to calculate the Euclidean distance
    calculate_divergence_score <- function(disease_expr, pn_expr) {
        # Calculate the Euclidean distance between disease and PN profiles
        distance <- sqrt(sum((disease_expr - pn_expr)^2))
        return(distance)
    }

    # Calculate the divergence score for each cell type
    divergence_scores <- list()

    for (celltype in celltypes) {
        disease_expr <- pseudo_bulk_disease[[celltype]]@assays$RNA@layers$scale.data
        rownames(disease_expr) <- rownames(seurat_obj_disease)
        
        pn_expr <- pseudo_bulk_pn[[celltype]]@assays$RNA@layers$scale.data
        rownames(pn_expr) <- rownames(seurat_obj_pn)
        
        disease_expr <- as.matrix(disease_expr)
        pn_expr <- as.matrix(pn_expr)

        ## calculate the median for each gene
        disease_expr <- apply(disease_expr, 1, median)
        pn_expr <- apply(pn_expr, 1, median)

        # Calculate divergence score
        divergence_scores[[celltype]] <- calculate_divergence_score(disease_expr, pn_expr)
    }


    # Create a data frame with the scores
    score_results <- data.frame(
    celltype = celltypes,
    divergence_score = unlist(divergence_scores)
    )

    # View the results
    print(score_results)
##

## TXd clean
    
    ## NOTES: ########################################################################################
    
    #"The TxD score is a quantification of transcriptional dysregulation. 
    #It represents the change in transcriptome-wide gene expression of each 
    #subtype in disease from its respective PN expression profile. 
    #The divergence score is the Euclidean distance between the median disease 
    #and corresponding PN covariate-corrected, pseudo-bulk expression profiles for each cell type."
    
    ## Functions: ######################################################################################## 

    # Create a function to calculate the pseudo-bulk expression profile (median of gene expression per cell type)
    calculate_pseudo_bulk <- function(seurat_obj, celltype) {
        # Subset the Seurat object by celltype
        celltype_data <- subset(seurat_obj, CellType == celltype)
        
        # Calculate pseudo-bulk expression
        pseudo_bulk_profile <- AggregateExpression(celltype_data, return.seurat = TRUE, verbose = FALSE, group.by = "Sample_ID")
        return(pseudo_bulk_profile)
    }

    # Function to calculate the Euclidean distance
    calculate_divergence_score <- function(disease_expr, pn_expr) {
        distance <- sqrt(sum((disease_expr - pn_expr)^2))
        return(distance)
    }

    complete_workflow_TXd <- function(brain_region, cell_class){
        seu <- readRDS(paste0('/home/fiorini9/scratch/machine_learning_ALS/base_objects/pineda_combat_',brain_region,'_',cell_class,'_int.rds'))
        ncol(seu)

        ## List of cell types
        Groups <- unique(seu$Group)

        for (group in Groups){
        
            ## Subset the Seurat object for disease and PN groups
            seurat_obj_disease <- subset(seu, Group == group)
            seurat_obj_pn <- subset(seu, Group == "PN")

            ## List of cell types
            celltypes <- unique(seu$CellType)

            ## Create empty lists to store pseudo-bulk profiles
            pseudo_bulk_disease <- list()
            pseudo_bulk_pn <- list()

            ## Loop through cell types and calculate pseudo-bulk profiles for both disease and PN
            for (celltype in celltypes) {
                pseudo_bulk_disease[[celltype]] <- calculate_pseudo_bulk(seurat_obj_disease, celltype)
                pseudo_bulk_pn[[celltype]] <- calculate_pseudo_bulk(seurat_obj_pn, celltype)
            }

            ## Calculate the divergence score for each cell type
            divergence_scores <- list()


            for (celltype in celltypes) {
                disease_expr <- pseudo_bulk_disease[[celltype]]@assays$RNA@layers$scale.data
                rownames(disease_expr) <- rownames(seurat_obj_disease)
                
                pn_expr <- pseudo_bulk_pn[[celltype]]@assays$RNA@layers$scale.data
                rownames(pn_expr) <- rownames(seurat_obj_pn)
                
                disease_expr <- as.matrix(disease_expr)
                pn_expr <- as.matrix(pn_expr)

                disease_expr <- apply(disease_expr, 1, median)
                pn_expr <- apply(pn_expr, 1, median)

                divergence_scores[[celltype]] <- calculate_divergence_score(disease_expr, pn_expr)
            }

            ## Dataframe
            score_results <- data.frame(
            celltype = celltypes,
            divergence_score = unlist(divergence_scores)
            )
            score_results$Region <- brain_region
            score_results$Group <- group

            fill <- rbind(fill, score_results)


        }
        return(fill)
    }

    
    
    ## Reset the filler frame: ########################################################################################
    
    fill <- data.frame(
    celltype = "fill",
    divergence_score = 0,
    Region = "fill",
    Group = "fill")
    
    
    ###########################
    ## Excitatory neurons BA4
    ###########################

    ## parameters
    brain_region = "BA4"

    cell_class = "Ex"

    test <- complete_workflow_TXd(brain_region, cell_class)


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Create a data frame with the scores
    score_results <- data.frame(
    celltype = celltypes,
    divergence_score = unlist(divergence_scores)
    )

    # View the results
    print(score_results)



    
    
    
    
    
    
    
    
    
    
    ## Parameters: ########################################################################################   

##
