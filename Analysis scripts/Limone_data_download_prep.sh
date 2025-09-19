#################################################################################################################################### Data download from GEO
################################################################################################################################ 
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/limone_data
    wget -O GSE226753_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE226753&format=file"
    tar -xvf GSE226753_RAW.tar
##

#################################################################################################################################### create barcodes.tsv, features, tsv, and matrix.mtx file from joint .txt file
################################################################################################################################ 
## code
    ## This file is a Matrixmarket in .txt. format, which includes the gene names, cell barcodes, and sparse matrix. 
    ## We need to parse this out into individual files. 

    ## start off by creating separate data for genes, barcodes, and matrix.
    ## read in the sparse matrix. 

    salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    ## load library
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(data.table)

    # Adjust the file path accordingly
    file_path <- "/home/fiorini9/scratch/machine_learning_ALS/limone_data/GSM7083181_MotorCortex_02202018.ICA_res0.1_mingenes400.raw.dge.txt"

    # Read the entire file as lines
    lines <- readLines(file_path)
    
    ########################
    # Genes
    ########################
    genes_index <- grep("^%%GENES", lines)
    genes_index <- setdiff(genes_index, 31)

    line_31 <- lines[31]

    # Extract gene names (columns)
    genes <- unlist(lapply(genes_index, function(i) strsplit(lines[i], "\t")[[1]]))
    line_31_split <- strsplit(line_31, "\t")[[1]]

    length(genes)
    length(line_31_split)
    
    genes_total <- c(genes, line_31_split)
    genes_total <- setdiff(genes_total, "%%GENES")
    length(genes_total) ## 29375 -- > matching dimentions of sparse matrix

    ########################
    ### GET BRACODES
    ########################
    cells_index <- grep("^%%CELL_BARCODES", lines)
    cells <- unlist(lapply(cells_index, function(i) strsplit(lines[i], "\t")[[1]]))
    length(cells)
    cells <- setdiff(cells, "%%CELL_BARCODES")


    ########################
    ### Matrix
    ########################
    matrix_start_index <- cells_index[length(cells_index)] + 1  # Next line after the last '%%CELL_BARCODES'

    matrix_data <- read.table(file_path, skip = matrix_start_index, col.names = c("row", "col", "value"))

    rows <- matrix_data$row
    cols <- matrix_data$col
    values <- matrix_data$value

    sparse_matrix <- sparseMatrix(i = rows, j = cols, x = values, dims = c(length(genes_total), length(cells)))


    ################################################################################################
    ### Now we have to save all of the files independently so that we can read them into Seurat
    ################################################################################################
    # /home/fiorini9/scratch/machine_learning_ALS/limone_data/matrix_files

    ###################
    ## Barcodes
    ###################
    #cells
    writeLines(cells, con = gzfile(paste0('/home/fiorini9/scratch/machine_learning_ALS/limone_data/matrix_files/barcodes.tsv.gz')))

    # check how many we have outside of the pool
    #cells_no_pool <- cells[!grepl("pool", cells)]
    base_names <- sub("_[^_]+$", "", cells)
    unique_base_names <- unique(base_names)

    ###################
    ## Features
    ###################
    #genes_total
    df <- data.frame(
    Column1 = genes_total,
    Column2 = genes_total
    )
    df$type <- "Gene Expression"
    write.table(df, file = gzfile(paste0('/home/fiorini9/scratch/machine_learning_ALS/limone_data/matrix_files/features.tsv.gz')), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )

    ###################
    ## Matrix
    ###################
    writeMM(sparse_matrix, file = "/home/fiorini9/scratch/machine_learning_ALS/limone_data/matrix_files/matrix.mtx")

    ## GO TO FILE AND DIRECTLY GZIP IT.
    cd /home/fiorini9/scratch/machine_learning_ALS/limone_data/matrix_files
    gzip matrix.mtx
##   
    


################################################################################################################################ scrnabox for Limone data
################################################################################################################################ 
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################


####################################
####################################
#################################### Step 0

mkdir /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone
mkdir /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 0 --method SCRNA --container TRUE


####################################
####################################
#################################### Create initial Seurat objects -- we will create distinct objects for each sample so that we can process them independently

## code
    ## load libraries
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")


    ##########################################
    ## Read in total gene expression matrix. 
    ##########################################
    ## test Seurat readin
    sparse_matrix <- Seurat::Read10X(data.dir = paste0('/home/fiorini9/scratch/machine_learning_ALS/limone_data/matrix_files/'))
    seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=1,min.features= 1, project = "all_samples") # only keep cells expressing atleast one gene
    dim(seurat_object)

    ## Modify the metadata to add the information that we know
    head(seurat_object@meta.data)
    meta_df <- data.frame(seurat_object@meta.data)

    ## Sample_ID (make this the orig.ident)
    meta_df$Sample_ID <- sub("_[^_]+$", "", rownames(meta_df))
    unique(meta_df$Sample_ID)

    meta_df$orig.ident <- meta_df$Sample_ID 

    ## Group
    meta_df$Group <- sub("_.*", "", meta_df$Sample_ID)
    table(meta_df$Sample_ID, meta_df$Group)

    ## Add metadata to Seurat object
    seurat_object <- AddMetaData(seurat_object, meta_df)

    ## Subset the Seurat objects according orig.ident
    table(meta_df$orig.ident, meta_df$Sample_ID)
    
    for(i in unique(seurat_object@meta.data$Sample_ID)){
        seu_subset <- subset(seurat_object, subset = Sample_ID == i)
        saveRDS(seu_subset, paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/processed_data/',i,'.RDS'))
    }  
##


####################################
####################################
#################################### Step 2

mkdir /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox/step2/
mkdir /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox/step2/objs/


## code -- here we  process the split objects with Step 2 of scRNAbox

    ## load parameters
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    
    
    ## get sample names
    list <-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/processed_data",sep=""),full.names = TRUE)
    sample_name <-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/processed_data",sep=""))

    for(i in 1:length(sample_name)){
        
        seu <- readRDS(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/processed_data/', sample_name[i] ) )
        
        ## calculate percent mitochondrial
        seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
        seu <- subset(seu, subset = percent.mt < 100)
        print(i)
        
        ## calculate percent ribosomal 
        seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]") 
        seu <- subset(seu, subset = percent.ribo < 100)
        print(i)

        ## Normalize and scale individual Seurat object prior to cell-cycle scoring
        seu <- Seurat::NormalizeData(seu,normalization.method = "LogNormalize" , scale.factor = 10000)
        ## Find variable features and print figure
        seu<- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2500)
        ## scale data
        seu<- ScaleData(seu)   #### we have to see if this works?
        ## perform linear dimensional reduction on individual Seurat objects
        seu <- RunPCA(seu, verbose =FALSE)
                    
      #save seurat object       
      saveRDS(seu,paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox/step2/objs/',sample_name[i], sep=''),compress=TRUE)
      
      ## print violin plot for QC metrics
      Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend()
      ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox/step2/',sample_name[i],".pdf", sep=""))
      
       ## print summary information
        sink(paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox/step2/summary_',sample_name[i],".txt", sep=""))
        cat("Summary of nCount_RNA: \n")
        print(summary(seu$nCount_RNA))
        cat("Summary of nFeature_RNA: \n")
        print(summary(seu$nFeature_RNA))
        cat("Summary of pt_mito: \n")
        print(summary(seu$percent.mt))
        cat("Summary of pt_ribo: \n")
        print(summary(seu$percent.ribo))
        cat("The number of features/genes and number of GEM/barcodes: \n")
        print(dim(seu))
        sink()
        write.csv(colnames(seu[[]]), file= paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox/step2/meta_info_',sample_name[i],".txt", sep=""))
    }

###

####################################
####################################
#################################### Step 3

/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox/step2/objs


export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 3 --method SCRNA --container TRUE


### interactive code
    #!/usr/bin/env Rscript

    ##########################################
    # step3: Quality control and filtering
    ##########################################

    ## load parameters
    args = commandArgs(trailingOnly=TRUE)
    output_dir="/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox"
    r_lib_path="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2"

    ## load library
    .libPaths(r_lib_path)
    packages<-c('Seurat','ggplot2', 'dplyr','Matrix', 'foreach', 'doParallel')
    invisible(lapply(packages, library, character.only = TRUE))

    ## load parameters text file
    source(paste(output_dir,'/job_info/parameters/step3_par.txt',sep=""))

    ## if user has existing Seurat object -- process Seurat objects and create list of distinct objects
    if (exists("par_seurat_object")) {                                                  
        sample_name<-list.files(path = par_seurat_object)
        sample_nameb<-gsub(".rds","",sample_name)
        if(length(sample_name)<1) {
        print("You do not have any existing Seurat object")
        }
    } else {
        sample_name<-list.files(path = paste(output_dir, "/step2/objs2",sep=""))
        sample_nameb<-gsub(".rds","",sample_name)
        if(length(sample_name)<1) {
        print("You do not have any object from step 2 ")
        }
    }

    ## create a list of available seurat objects in 
    for (i in 1:length(sample_name)) {
    if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
        print(c(sample_name[i],"is not R rds"))
    }
    } 

    ## set seed for replicability
    set.seed(1234)

    ## detect number of cores
    numCores <- detectCores()
    cl <- makeCluster(numCores-1)
    registerDoParallel(cl) 

    ###### QC and filtering if users have a preporcessed Seurat object
    i = 1
    foreach (i=1:length(sample_name)) %do% {    
        if (exists("par_seurat_object")) { 
            seu<-readRDS(paste(par_seurat_object, "/", sample_name[i], sep=""))
        } else  {
            seu<-readRDS(paste(output_dir,'/step2/objs2/',sample_name[i], sep=""))
        }
        print(sample_name[i])
        
        ## calculate percent MT and Ribo if users are starting from Step 3 as they may not have this claculated in their Seurat object.
        if (exists("par_seurat_object")) { 
        ## calculate percent MT 
        seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
        ## calculate percent ribo  
        seu[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RP[SL]")    #NEW CODE
        }

        ## fiter according to user-defined thresholds.
        if (exists("par_nFeature_RNA_L")) {
        cat("nFeature_RNA Lower: \n")
        print(par_nFeature_RNA_L)
        seu <- subset(seu, subset = nFeature_RNA > par_nFeature_RNA_L) 
        }
        if (exists("par_nFeature_RNA_U")){
            cat("nFeature_RNA Upper: \n")
            print(par_nFeature_RNA_U)
            seu <- subset(seu, subset = nFeature_RNA < par_nFeature_RNA_U) 
        }
        if (exists("par_nCount_RNA_L")){
            cat("nCount_RNA Lower: \n")
            print(par_nCount_RNA_L)
            seu <- subset(seu, subset = nCount_RNA > par_nCount_RNA_L) 
        }
        if (exists("par_nCount_RNA_U")){
            cat("nCount_RNA Upper: \n")
            print(par_nCount_RNA_U)        
            seu <- subset(seu, subset = nCount_RNA < par_nCount_RNA_U) 
        }
        if (exists("par_mitochondria_percent_L")) {
            cat("Mitochondria_percent Lower: \n")
            print(par_mitochondria_percent_L)        
            seu <- subset(seu, subset = percent.mt > par_mitochondria_percent_L) 
        }
        if (exists("par_mitochondria_percent_U")){
            cat("Mitochondria_percent Upper: \n")
            print(par_mitochondria_percent_U)                
            seu <- subset(seu, subset = percent.mt < par_mitochondria_percent_U) 
        }
        if (exists("par_ribosomal_percent_L")) {
            cat("Ribosomal_percent Lower: \n")
            print(par_ribosomal_percent_L)                        
            seu <- subset(seu, subset = percent.ribo > par_ribosomal_percent_L) 
        }
        if (exists("par_ribosomal_percent_U")) {
            cat("Ribosomal_percent Upper: \n")
            print(par_ribosomal_percent_U)                                
            seu <- subset(seu, subset = percent.ribo < par_ribosomal_percent_U) 
        }

        ## optional: filter out mitochondrial genes
        if (tolower(par_remove_mitochondrial_genes)=='yes') {
        MT_genes <- grep( "^MT-", rownames(seu), value = T)
        counts <- GetAssayData(seu, assay = "RNA")
        counts <- counts[-(which(rownames(counts) %in% MT_genes)),]
        seu[["RNA"]] <- subset(seu[["RNA"]], features = rownames(counts))     
        }

        ## optional: filter out ribosomal genes
        if (tolower(par_remove_ribosomal_genes)=='yes') {
        Ribo_genes <- grep( "^RP[SL]", rownames(seu), value = T)
        counts <- GetAssayData(seu, assay = "RNA")
        counts <- counts[-(which(rownames(counts) %in% Ribo_genes)),]
        seu[["RNA"]] <- subset(seu[["RNA"]], features = rownames(counts))       
        }

        ## remove custom list of genes
        if (exists("par_remove_genes")) {
        counts <- GetAssayData(seu, assay = "RNA")
        counts <- counts[-(which(rownames(counts) %in% par_remove_genes)),]
        seu[["RNA"]] <- subset(seu[["RNA"]], features = rownames(counts))    
        }

        ## normalize after filtering
        seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)
        
        ## find variable features after filtering
        seu<- FindVariableFeatures(seu, selection.method = par_selection.method, nfeatures = par_nfeatures)
        topsel <- head(Seurat::VariableFeatures(seu), par_top)
        write.csv(topsel, file = paste(output_dir,'/step3/info3/most_variable_genes_',sample_nameb[i],'.txt', sep=""), quote = TRUE)
        
        ## print variable features plot
        vf_plot <- Seurat::VariableFeaturePlot(seu)
        Seurat::LabelPoints(plot = vf_plot,points = topsel, repel = TRUE)
        ggsave(paste(output_dir,'/step3/figs3/VariableFeaturePlot_',sample_nameb[i],".pdf",sep=""))

        ##Do not regress out cc genes or custom list
        if (tolower(par_regress_cell_cycle_genes)=='no' & tolower(par_regress_custom_genes)=='no') {
        seu<- ScaleData(seu, verbose = FALSE)
        }

        ## Regress out cc genes only
        if (tolower(par_regress_cell_cycle_genes)=='yes' & tolower(par_regress_custom_genes)=='no') {
        seu<- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
        }

        ## Regress out custom gene list only
        if (tolower(par_regress_custom_genes)=='yes' & tolower(par_regress_cell_cycle_genes)=='no') {
            par_regress_genes <- list(par_regress_genes)
            seu <- AddModuleScore(
            object = seu,
            features = par_regress_genes,
            nbin = 24,
            ctrl = 4,
            k = FALSE,
            name = 'regress_features')    
        seu<- ScaleData(seu, vars.to.regress = "regress_features", verbose = FALSE)
        }

        ## Regress out custom gene list and cc genes
        if (tolower(par_regress_custom_genes)=='yes' & tolower(par_regress_cell_cycle_genes)=='yes') {
            par_regress_genes <- list(par_regress_genes)
            seu <- AddModuleScore(
            object = seu,
            features = par_regress_genes,
            nbin = 24,
            ctrl = 4,
            k = FALSE,
            name = 'regress_features')    
            seu<- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score", "regress_features"), verbose = FALSE)
        }

        ## perform linear dimensional reduction
        seu <- RunPCA(seu, npcs = par_npcs_pca, verbose = FALSE)
                ## print PCA 
                DimPlot(seu, reduction = "pca", raster = FALSE)
                ggsave(paste(output_dir,'/step3/figs3/',"dimplot_pca_",sample_nameb[i],".pdf",sep=""))
                ## print elbow plot
                ElbowPlot(seu, ndims = par_npcs_pca)
                ggsave(paste(output_dir,'/step3/figs3/',"elbowplot_",sample_nameb[i],".pdf",sep=""))
        
        ## save each individual Seurat object as RDS
        saveRDS(seu, paste(output_dir,'/step3/objs3/',sample_nameb[i],".rds", sep=""))
        
        ## print QC violin plot
        Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo","S.Score", "G2M.Score"), pt.size = 0.001,ncol = 3, raster = FALSE) + NoLegend() #new code
        ggsave(paste(output_dir,'/step3/figs3/filtered_QC_vioplot_',sample_nameb[i],".pdf", sep=""))
        
        ## write meta infor available in the Seurat metdata
        write.csv(colnames(seu[[]]), file= paste(output_dir,'/step3/info3/meta_info_',sample_nameb[i],".txt", sep=""))
        
        ## save RNA expression matrix for each individual Seurat object
        if (tolower(par_save_RNA)=='yes') {
        mat <- GetAssayData(object = seu, assay = "RNA", slot = "data")
        writeMM(mat,file= paste(output_dir,'/step3/info3/',sample_nameb[i],"_RNA.txt", sep=""))
        }
        ## save metadata dataframe for each individual Seurat object
        if (tolower(par_save_metadata)=='yes') {
        write.csv(seu[[]], file = paste(output_dir,'/step3/info3/MetaData_',sample_nameb[i],'.txt', sep=""), quote = TRUE)
        }
        ## save summary information for each individual Seurat object
        sink(paste(output_dir,'/step3/info3/summary_',sample_nameb[i],".txt", sep=""))
        cat("Summary of nCount_RNA: \n")
        print(summary(seu$nCount_RNA))
        cat("Summary of nFeature_RNA: \n")
        print(summary(seu$nFeature_RNA))
        cat("Summary of pt_mito: \n")
        print(summary(seu$percent.mt))
        cat("Summary of pt_ribo: \n") 
        print(summary(seu$percent.ribo)) 
        cat("The number of features/genes and number of GEM/barcodes: \n")
        print(dim(seu))
        sink()

    }

    cat(stepp,"has been achieved. Total time:",as.numeric (Sys.time() - start_time, units = "mins"),"minutes\n")
    cat("#####################################\n")

    ## save session information
    writeLines(capture.output(sessionInfo()), paste(output_dir,'/step3/info3/sessionInfo.txt', sep=""))
    if(file.exists("Rplots.pdf")){
        file.remove("Rplots.pdf")
    }

    cat("##########################################################################\n")
    cat(stepp0,"successfully completed. Total time:",as.numeric (Sys.time() - start_time0, units = "mins"),"minutes\n")
    cat("##########################################################################\n")

###

####################################
####################################
#################################### Step 4

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 4 --method SCRNA --container TRUE

####################################
####################################
#################################### Step 5 --> merge

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 5 --method SCRNA --container TRUE

## COPY OVER THE MERGED OBJECT TO SOME OTHER LOCATION OUTSIDE OF STEP 5 RNABOX
cp /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox/step5/objs5/seu_step5.rds /home/fiorini9/scratch/machine_learning_ALS/base_objects/limone_ALS_BA4_merged.rds



####################################
####################################
#################################### We need to transfer the cell type annotations from Pineda to Li. 
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox
    mkdir step7
    cd step7
    mkdir objs7
    mkdir figs7
    mkdir info7

    salloc -A def-sfarhan --time=0-8 -c 1 --mem=150g
    module load StdEnv/2020 
    module load r/4.2.2 
    R

    output_dir="/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox"
    r_lib_path="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2"
    .libPaths(r_lib_path)

    ## load library
    packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix','scCustomize')
    invisible(lapply(packages, library, character.only = TRUE))

    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/limone_ALS_BA4_merged.rds')
    
    colnames(seu@meta.data)

    unique(seu@meta.data$Sample_ID)
    unique(seu@meta.data$orig.ident)
    unique(seu@meta.data$Group)
    dim(seu) #29375 73305

    ##############################################################################
    ## use Pineda data to transfer labels. 
    ## BA4
    ##############################################################################
    ## load reference Seurat object
    par_reference_name = "Pineda_BA4"
    reference0 <-readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA4.rds')
    dim(reference0)
    length(unique(reference0@meta.data$CellType)) #19
    unique(reference0@meta.data$Region)

    set.seed(42)  # For reproducibility (optional)
    sample_cells <- sample(Cells(reference0), 75000)
    reference0 <- subset(reference0, cells = sample_cells)

    dim(reference0)

    length(unique(reference0@meta.data$CellType))

    DefaultAssay(reference0) <- "RNA" ## new code


    # perform standard preprocessing on reference object
    reference0<- NormalizeData(reference0)
    reference0 <- FindVariableFeatures(reference0)
    reference0<- ScaleData(reference0)
    reference0 <- RunPCA(object = reference0, assay = "RNA", npcs = 50)

    ## find transfer anchors between reference and query Seurat objects
    transfer.anchors <- FindTransferAnchors(reference = reference0, query = seu, dims = 1:50, reference.reduction = "pca")

    ## add reference-based annotations to the qeury object
    eval(parse(text = paste('predictions <- TransferData(anchorset = transfer.anchors, refdata = reference0$',"CellType" ,',dims = 1:',50,')', sep='')))
    seu <- AddMetaData(object = seu, metadata = predictions)

    # Add metadata column for reference object
    seu$temp_temp_2 <- seu@meta.data$predicted.id
    name_meta <- names(seu@meta.data) 
    length <- length(name_meta)
    name_meta[length] <- paste(par_reference_name, "_predictions", sep = "")
    names(seu@meta.data) <- name_meta

    ## save query Seurat object with reference annotation predicitions
    saveRDS(seu,paste(output_dir,'/step7/objs7','/seu_step7.rds', sep=''))

    ## save metadata information
    write.csv(colnames(seu[[]]), file= paste(output_dir,'/step7/info7/meta_info_seu_step7',".txt", sep=""))

    ## Print a umap projection showing the predicted cell types on the query object 
    reference0 <- RunUMAP(reference0, dims = 1:50, reduction = "pca", return.model = TRUE)
    seu <- MapQuery(anchorset = transfer.anchors, reference = reference0, query = seu,
        refdata = list(celltype = "CellType"), reference.reduction = "pca", reduction.model = "umap")
    p1 <- DimPlot(reference0, reduction = "umap", group.by = "CellType", label = TRUE, label.size = 3,repel = TRUE, raster=FALSE) + NoLegend() + ggtitle("Reference annotations")
    p2 <- DimPlot(seu, reduction = "umap", group.by = "predicted.id", label = TRUE, label.size = 3, repel = TRUE, raster=FALSE) + NoLegend() + ggtitle("Query transferred labels")
    p1 + p2
    ggsave(file = paste(output_dir,'/step7/figs7/',par_reference_name,'_UMAP_transferred_labels.pdf', sep=''), dpi = 300, height = 7, width = 14, unit = 'in')
 
##


## COPY OVER THE MERGED OBJECT TO SOME OTHER LOCATION OUTSIDE OF STEP 5 RNABOX
cp /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_limone/scrnabox/step7/objs7/seu_step7.rds /home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_SALS_BA4_merged.rds


####################################
####################################
#################################### Improve annotations annotations

## code All cells
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=300g

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    library(harmony)
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    
    
    ####################################
    ## All cells
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_SALS_BA4_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$annotation_cell_class) #"Exc_neuron" "Inh_neuron" "Non_neuron"

    ## Run Harmony
    length(unique(seu$Sample_ID))
    
    seu <- ScaleData(seu)
    seu <- RunPCA(seu, npcs = 20)
    seu <- RunHarmony(seu, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)
    seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
    seu <- FindClusters(seu, resolution = 1.5)
    
    ## By CellType subtype
    DimPlot(seu, reduction = "umap", group.by = "Pineda_BA4_predictions", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 

    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    ## Seurat clusters
    DimPlot(seu, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 

    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))


    ## modify annotations according to Seurat clusters
    meta_df <- data.frame(seu@meta.data)
    meta_df$cell_type_temp <- meta_df$Pineda_BA4_predictions
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 0] <- "Oligo"
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 1] <- "Oligo"
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 2] <- "Astro"
    #meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 3] <- 
    #meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 4] <- 
    #meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 5] <- 
    #meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 6] <- 
    #meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 7] <- 
    #meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 8] <- 
    #meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 9] <- 
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 10] <- "PV"
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 11] <- "5HT3aR"
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 12] <- "OPC"
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 13] <- "Micro"
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 14] <- "SOM"
    #meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 15] <- 
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 16] <- "REMOVE"
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 17] <- "Rosehip"
    #meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 18] <- 
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 19] <- "REMOVE"
    meta_df$cell_type_temp[meta_df$RNA_snn_res.0.8 == 20] <- "REMOVE"

    ## Reprint UMAP

    seu <- AddMetaData(seu, meta_df)
    seu_subset <- subset(seu, subset = cell_type_temp != c("REMOVE"))
    unique(seu_subset@meta.data$cell_type_temp)

    ## add colours for Ex to make it clearer
    DimPlot(seu_subset, reduction = "umap", group.by = "cell_type_temp", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    
    ######################################################
    #excitatory neurons
    ######################################################
    
    ## subset to only include the Exc neurons to re-cluster
    seu_subset_ex <- subset(seu, subset = cell_type_temp %in% c("L3_L5","L2_L3","L4_L6",  "L4_L5", "L5_L6" ,  "L5", "L6"  ))
    nrow(seu_subset_ex)
    ncol(seu_subset_ex)
    unique(seu_subset_ex@meta.data$predicted.id)

    ## recluster
    seu_subset_ex <- ScaleData(seu_subset_ex)
    seu_subset_ex <- RunPCA(seu_subset_ex, npcs = 50)
    seu_subset_ex <- RunHarmony(seu_subset_ex, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_subset_ex <- RunUMAP(seu_subset_ex, reduction = "harmony", dims = 1:50)
    seu_subset_ex <- FindNeighbors(seu_subset_ex, reduction = "harmony", dims = 1:50)
    seu_subset_ex <- FindClusters(seu_subset_ex, resolution = 0.8)
    
    ## By CellType subtype
    DimPlot(seu_subset_ex, reduction = "umap", group.by = "Pineda_BA4_predictions", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    DimPlot(seu_subset_ex, reduction = "umap", group.by = "seurat_clusters", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    ## modify annotations according to Seurat clusters
    meta_df <- data.frame(seu_subset_ex@meta.data)
    meta_df$cell_type_temp_ex <- meta_df$Pineda_BA4_predictions
    
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 0] <- "L2_L3"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 1] <- "Unknown"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 2] <- "L3_L5"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 3] <- "L4_L5"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 4] <- "L5_L6"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 5] <- "L6"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 6] <- "L3_L5"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 7] <- "L4_L5"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 8] <- "L6"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 9] <- "L2_L3"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 10] <- "L5_L6"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 11] <- "Unknown"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 12] <- "L4_L6"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 13] <- "L4_L5"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 14] <- "L5_L6"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 15] <- "L5"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 16] <- "L5"
    meta_df$cell_type_temp_ex[meta_df$seurat_clusters == 17] <- "Unknown"

    seu_subset_ex <- AddMetaData(seu_subset_ex, meta_df) ## ADD THIS TO THE SEAURAT OBJECT
    seu_subset_ex_lim <- subset(seu_subset_ex, subset = cell_type_temp_ex != c("Unknown"))
    unique(seu_subset_ex_lim@meta.data$cell_type_temp_ex)

    ## add colours for Ex to make it clearer
    DimPlot(seu_subset_ex_lim, reduction = "umap", group.by = "cell_type_temp_ex", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    ## recluster
    seu_subset_ex_lim <- ScaleData(seu_subset_ex_lim)
    seu_subset_ex_lim <- RunPCA(seu_subset_ex_lim, npcs = 50)
    seu_subset_ex_lim <- RunHarmony(seu_subset_ex_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_subset_ex_lim <- RunUMAP(seu_subset_ex_lim, reduction = "harmony", dims = 1:50)
    seu_subset_ex_lim <- FindNeighbors(seu_subset_ex_lim, reduction = "harmony", dims = 1:50)
    seu_subset_ex_lim <- FindClusters(seu_subset_ex_lim, resolution = 0.8)
    
    ## By CellType subtype
    DimPlot(seu_subset_ex_lim, reduction = "umap", group.by = "cell_type_temp_ex", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    ## WE ARE GOING TO PRINT THE MAIN FIGURE FOR EX RIGHT HERE. 
    DimPlot(seu_subset_ex_lim, reduction = "umap", group.by = "cell_type_temp_ex", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
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

    ## recluster at lower PC
    seu_subset_ex_lim <- ScaleData(seu_subset_ex_lim)
    seu_subset_ex_lim <- RunPCA(seu_subset_ex_lim, npcs = 35)
    seu_subset_ex_lim <- RunHarmony(seu_subset_ex_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_subset_ex_lim <- RunUMAP(seu_subset_ex_lim, reduction = "harmony", dims = 1:35)
    seu_subset_ex_lim <- FindNeighbors(seu_subset_ex_lim, reduction = "harmony", dims = 1:35)
    seu_subset_ex_lim <- FindClusters(seu_subset_ex_lim, resolution = 0.8)

    ## WE ARE GOING TO PRINT THE MAIN FIGURE FOR EX RIGHT HERE. --> USE THIS ONE
    DimPlot(seu_subset_ex_lim, reduction = "umap", group.by = "cell_type_temp_ex", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
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


    ## Save the object
    saveRDS(seu_subset_ex_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_SALS_BA4_Ex_temp.rds')

    ######################################################
    #Inhibitory neurons
    ######################################################
    
    ## subset to only include the Exc neurons to re-cluster
    seu_subset_in <- subset(seu, subset = cell_type_temp %in% c("SOM","PV","Rosehip",  "5HT3aR"))
    nrow(seu_subset_in)
    ncol(seu_subset_in)
    unique(seu_subset_in@meta.data$predicted.id)

    ## recluster
    seu_subset_in <- ScaleData(seu_subset_in)
    seu_subset_in <- RunPCA(seu_subset_in, npcs = 20)
    seu_subset_in <- RunHarmony(seu_subset_in, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_subset_in <- RunUMAP(seu_subset_in, reduction = "harmony", dims = 1:20)
    seu_subset_in <- FindNeighbors(seu_subset_in, reduction = "harmony", dims = 1:20)
    seu_subset_in <- FindClusters(seu_subset_in, resolution = 0.8)
    
    ## By CellType subtype
    DimPlot(seu_subset_in, reduction = "umap", group.by = "Pineda_BA4_predictions", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    
    DimPlot(seu_subset_in, reduction = "umap", group.by = "seurat_clusters", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    ## modify annotations according to Seurat clusters
    meta_df <- data.frame(seu_subset_in@meta.data)
    meta_df$cell_type_temp_in <- meta_df$Pineda_BA4_predictions
    
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 0] <- "PV"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 1] <- "SOM"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 2] <- "PV"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 3] <- "5HT3aR"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 4] <- "5HT3aR"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 5] <- "Rosehip"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 6] <- "5HT3aR"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 7] <- "SOM"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 8] <- "Rosehip"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 9] <- "PV"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 10] <- "5HT3aR"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 11] <- "5HT3aR"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 12] <- "PV"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 13] <- "Unknown"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 14] <- "Unknown"
    meta_df$cell_type_temp_in[meta_df$seurat_clusters == 15] <- "SOM"


    seu_subset_in <- AddMetaData(seu_subset_in, meta_df) ## ADD THIS TO THE SEAURAT OBJECT
    seu_subset_in_lim <- subset(seu_subset_in, subset = cell_type_temp_in != c("Unknown"))
    unique(seu_subset_in_lim@meta.data$cell_type_temp_in)

    ## add colours for Ex to make it clearer
    DimPlot(seu_subset_in_lim, reduction = "umap", group.by = "cell_type_temp_in", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    ## recluster
    seu_subset_in_lim <- ScaleData(seu_subset_in_lim)
    seu_subset_in_lim <- RunPCA(seu_subset_in_lim, npcs = 35)
    seu_subset_in_lim <- RunHarmony(seu_subset_in_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_subset_in_lim <- RunUMAP(seu_subset_in_lim, reduction = "harmony", dims = 1:35)
    seu_subset_in_lim <- FindNeighbors(seu_subset_in_lim, reduction = "harmony", dims = 1:35)
    seu_subset_in_lim <- FindClusters(seu_subset_in_lim, resolution = 0.8)
    
    ## By CellType subtype
    DimPlot(seu_subset_in_lim, reduction = "umap", group.by = "cell_type_temp_in", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    ## WE ARE GOING TO PRINT THE MAIN FIGURE FOR EX RIGHT HERE. 
    DimPlot(seu_subset_in_lim, reduction = "umap", group.by = "cell_type_temp_in", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
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

    ## Save the object
    saveRDS(seu_subset_in_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_SALS_BA4_In_temp.rds')


    ######################################################
    #Non neurons
    ######################################################
    
    ## subset to only include the Exc neurons to re-cluster
    seu_subset_non <- subset(seu, subset = cell_type_temp %in% c("Oligo", "Astro", "OPC", "Micro", "T_Cell", "Mural", "Endo", "Fibro"))
    nrow(seu_subset_non)
    ncol(seu_subset_non)
    unique(seu_subset_non@meta.data$predicted.id)

    ## recluster
    seu_subset_non <- ScaleData(seu_subset_non)
    seu_subset_non <- RunPCA(seu_subset_non, npcs = 20)
    seu_subset_non <- RunHarmony(seu_subset_non, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_subset_non <- RunUMAP(seu_subset_non, reduction = "harmony", dims = 1:20)
    seu_subset_non <- FindNeighbors(seu_subset_non, reduction = "harmony", dims = 1:20)
    seu_subset_non <- FindClusters(seu_subset_non, resolution = 0.8)
    
    ## By CellType subtype
    DimPlot(seu_subset_non, reduction = "umap", group.by = "Pineda_BA4_predictions", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    
    DimPlot(seu_subset_non, reduction = "umap", group.by = "seurat_clusters", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))


    ## modify annotations according to Seurat clusters
    meta_df <- data.frame(seu_subset_non@meta.data)
    meta_df$cell_type_temp_non <- meta_df$Pineda_BA4_predictions
    
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 0] <- "Oligo"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 1] <- "Oligo"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 2] <- "Astro"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 3] <- "Oligo"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 4] <- "Oligo"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 5] <- "Endo"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 6] <- "OPC"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 7] <- "Micro"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 8] <- "Astro"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 9] <- "Astro"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 10] <- "Mural"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 11] <- "Unknown"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 12] <- "Fibro"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 13] <- "Oligo"
    meta_df$cell_type_temp_non[meta_df$seurat_clusters == 14] <- "T_Cell"



    seu_subset_non <- AddMetaData(seu_subset_non, meta_df) ## ADD THIS TO THE SEAURAT OBJECT
    seu_subset_non_lim <- subset(seu_subset_non, subset = cell_type_temp_non != c("Unknown"))
    unique(seu_subset_non_lim@meta.data$cell_type_temp_non)

    ## add colours for Ex to make it clearer
    DimPlot(seu_subset_non_lim, reduction = "umap", group.by = "cell_type_temp_non", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    ## recluster
    seu_subset_non_lim <- ScaleData(seu_subset_non_lim)
    seu_subset_non_lim <- RunPCA(seu_subset_non_lim, npcs = 35)
    seu_subset_non_lim <- RunHarmony(seu_subset_non_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_subset_non_lim <- RunUMAP(seu_subset_non_lim, reduction = "harmony", dims = 1:35)
    seu_subset_non_lim <- FindNeighbors(seu_subset_non_lim, reduction = "harmony", dims = 1:35)
    seu_subset_non_lim <- FindClusters(seu_subset_non_lim, resolution = 0.8)
    
    ## By CellType subtype
    DimPlot(seu_subset_non_lim, reduction = "umap", group.by = "cell_type_temp_non", label = T, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2") 
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))

    ## WE ARE GOING TO PRINT THE MAIN FIGURE FOR EX RIGHT HERE. 
    DimPlot(seu_subset_non_lim, reduction = "umap", group.by = "cell_type_temp_non", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
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

    ## Save the object
    saveRDS(seu_subset_non_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_SALS_BA4_Non_temp.rds')

    ######################################################
    #Create final object and save
    ######################################################

    ncol(seu_subset_ex) #31115
    ncol(seu_subset_in) #8758
    ncol(seu_subset_non) #31386
    ncol(seu)

    ## create Ex df
    df_ex <- data.frame(seu_subset_ex@meta.data)
    nrow(df_ex)
    df_ex$rownames <- rownames(df_ex)
    df_ex <- df_ex %>% dplyr::select(rownames, cell_type_temp_ex)
    colnames(df_ex) <- c("rownames", "merged_celltype")
    table(df_ex$merged_celltype)

    ## create in df
    df_in <- data.frame(seu_subset_in@meta.data)
    nrow(df_in)
    df_in$rownames <- rownames(df_in)
    df_in <- df_in %>% dplyr::select(rownames, cell_type_temp_in)
    colnames(df_in) <- c("rownames", "merged_celltype")
    table(df_in$merged_celltype)

    ## create non df
    df_non <- data.frame(seu_subset_non@meta.data)
    nrow(df_non)
    df_non$rownames <- rownames(df_non)
    df_non <- df_non %>% dplyr::select(rownames, cell_type_temp_non)
    colnames(df_non) <- c("rownames", "merged_celltype")
    table(df_non$merged_celltype)

    ## all merge
    df_all_merge <- rbind(df_ex, df_in, df_non)
    nrow(df_all_merge) == nrow(df_ex) + nrow(df_in) + nrow(df_non)
    
    seu_test <- AddMetaData(seu, df_all_merge)
    unique(seu_test@meta.data$merged_celltype)

    seu_test$merged_celltype[is.na(seu_test$merged_celltype)] <- "Unknown"
    unique(seu_test@meta.data$merged_celltype)

    ## remove NA
    seu_subset_no_NA <- subset(seu_test, subset = merged_celltype != "Unknown")

    ncol(seu_subset_no_NA) == nrow(df_ex) + nrow(df_in) + nrow(df_non)

    saveRDS(seu_subset_no_NA, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_SALS_BA4_merged_annotated.rds')     
##


####################################
####################################
####################################
####################################
####################################
####################################
#################################### Clean and create separate objects for each brain region -- make sure that the metadata column names match Pineda. 
## code 

    salloc -A def-sfarhan --time=0-8 -c 1 --mem=300g

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    library(harmony)
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    
    ####################################
    ## BA4 merge
    ####################################
    
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_SALS_BA4_merged_annotated.rds')

    ## only select BA4
    str(seu@meta.data)
    
    ## Add Region column
    seu@meta.data$brain_region <- "BA4"


    
    ## fix column names
    df <- seu@meta.data
    df <- df %>% dplyr::select(merged_celltype, orig.ident, Group, brain_region, Sample_ID)
    colnames(df)
    colnames(df) <- c('CellType', 'Sample_ID', 'Group', 'Region', 'Donor')
    seu <- AddMetaData(seu, df)

    ## print
    saveRDS(seu, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Limone_All_conditions_BA4.rds')

##


### AFTER cell type annotation, need to do combat. 
### After combat transfer to Narval.