## THIS WAS PERFORMED IN BELUGA AND THEN TRANSFERED TO NARVAL AFTER THE DATA DOWNLOAD.

cd /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download

module load StdEnv/2020 
module load python/3.8.10
PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
            source $PYTHONENV0/bin/activate
            export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

module load StdEnv/2020 
module load r/4.2.2 
R

#################################################################################################################################### Data download with synapse
################################################################################################################################ 
## code
    module load StdEnv/2020 
    module load python/3.8.10
    PYTHONENV0=/lustre03/project/6070393/COMMON/samplepooling/software/ENVscRNA
                source $PYTHONENV0/bin/activate
                export PYTHONPATH=$PYTHONENV0/lib/python3.8/site-packages

        synapse get -r syn51105515

        mv 'Processed Data' processed_data

        ## we'll copy the previous downloads
        cp -r /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download
        cp -r '/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/Supplemental Materials' /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download
        cp -r /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/SYNAPSE_METADATA_MANIFEST.tsv /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download
##

#################################################################################################################################### Create Seurat object -- this was the code we used for the original download. 
################################################################################################################################  DONT USE THIS??

## code

    cd /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download

    nano data_prep.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=00-10:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=40g          # memory per cor
    #SBATCH --job-name=data_prep
    #SBATCH --error=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.err
    #SBATCH --output=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download

    Rscript /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/data_prep.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano data_prep.R


    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    ## Create a processing function
    list<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data",sep=""),full.names = TRUE)
    sample_name<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data",sep=""))

    sample_name <- sample_name[!(sample_name %in% "SYNAPSE_METADATA_MANIFEST.tsv")]  

    i = '220616_BL6_MCX_snRNA-H12'
    for(i in unique(sample_name)){
        print(i)
        ###################
        ## Barcodes
        ###################
        meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/col_metadata.tsv'), header = T, sep = "\t")
        meta_data <- meta_data[,1]
        writeLines(meta_data, con = gzfile(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/barcodes.tsv.gz')))

        ###################
        ## Features
        ###################

        meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/row_metadata.tsv'), header = T, sep = "\t")
        meta_data <- meta_data[,c(1,2)]
        meta_data$type <- "Gene Expression"

        meta_data$Gene <- ifelse(meta_data$Gene == "N/A", paste0("NA-", seq_along(meta_data$Gene)), meta_data$Gene)

        write.table(meta_data, file = gzfile(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/features.tsv.gz')), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )

        ###################
        ## Matrix
        ###################
        ## rename file
        # Specify the full paths for old and new files
        old_file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/counts_fil.mtx')
        new_file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/matrix.mtx')

        # Rename the file
        file.rename(from = old_file, to = new_file)

        # Set name
        mtx_file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/matrix.mtx')

        # Specify the path where you want to save the gzipped file
        gzipped_file <- paste0(mtx_file, ".gz")

        # Execute gzip command using system()
        system(paste("gzip -c", mtx_file, ">", gzipped_file))


        ###################
        ## Test create Seurat object
        ###################
        sparse_matrix <- Seurat::Read10X(data.dir = paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i))
        seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=1,min.features= 1, project = i) # only keep cells expressing atleast one gene

        ## import
        meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/col_metadata.tsv'), header = T, sep = "\t")
        seurat_object <- AddMetaData(seurat_object, meta_data)

        saveRDS(seurat_object, paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/',i,'.RDS'))
    }
##

#################################################################################################################################### Create Seurat object -- USE THIS ONE.
################################################################################################################################ 

## code

    cd /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download

    nano data_prep.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=00-10:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=40g          # memory per cor
    #SBATCH --job-name=data_prep
    #SBATCH --error=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.err
    #SBATCH --output=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 

    cd /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download

    Rscript /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/data_prep.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano data_prep.R


    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    ## Create a processing function
    list<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data",sep=""),full.names = TRUE)
    sample_name<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data",sep=""))

    sample_name <- sample_name[!(sample_name %in% "SYNAPSE_METADATA_MANIFEST.tsv")]  

    ## we are going to remove the mouse samples
    #sample_name <- sample_name[grepl("MCX|PFC", sample_name)]

    i = '191112_ALS_110_snRNA-B9'
    for(i in unique(sample_name)){
        print(i)
        ###################
        ## Barcodes
        ###################
        meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/col_metadata.tsv'), header = T, sep = "\t")
        meta_data <- meta_data[,1]
        writeLines(meta_data, con = gzfile(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/barcodes.tsv.gz')))

        ###################
        ## Features
        ###################

        meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/row_metadata.tsv'), header = T, sep = "\t")
        meta_data <- meta_data[,c(1,2)]
        meta_data$type <- "Gene Expression"

        meta_data$Gene[ meta_data$Gene == "N/A" ] <- meta_data$ENSEMBL[meta_data$Gene == "N/A" ]
        
        write.table(meta_data, file = gzfile(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/features.tsv.gz')), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE )

        ###################
        ## Matrix
        ###################
        ## rename file
        # Specify the full paths for old and new files
        #old_file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/counts_fil.mtx')
        #new_file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/matrix.mtx')

        # Rename the file
        #file.rename(from = old_file, to = new_file)

        # Set name
        #mtx_file <- paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/matrix.mtx')

        # Specify the path where you want to save the gzipped file
        #gzipped_file <- paste0(mtx_file, ".gz")

        # Execute gzip command using system()
        #system(paste("gzip -c", mtx_file, ">", gzipped_file))


        ###################
        ## Test create Seurat object
        ###################
        sparse_matrix <- Seurat::Read10X(data.dir = paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i))
        seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=1,min.features= 1, project = i) # only keep cells expressing atleast one gene

        ## import
        meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/col_metadata.tsv'), header = T, sep = "\t")
        
        ## merge
        colnames(seurat_object)
        meta_data$Barcode %in% colnames(seurat_object)
        rownames(meta_data) <- meta_data$Barcode

        seurat_object <- AddMetaData(seurat_object, meta_data)

        saveRDS(seurat_object, paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/',i,'.RDS'))
    }
##

#################################################################################################################################### identify different sample groups (sALS-control; sFTLD-control)
################################################################################################################################ 
## code
    list<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data",sep=""),full.names = TRUE)
    sample_name<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data",sep=""))

    sample_name <- sample_name[!(sample_name %in% "SYNAPSE_METADATA_MANIFEST.tsv")]

        
    Sample_ID <- "fill"
    Group <- "fill"
    Region <- "fill"
    fill <- data.frame(Sample_ID, Group, Region)

    i = "201019_ALS_102_snRNA-B2"
    for(i in unique(sample_name)){
        print(i)
        ###################
        ## Barcodes
        ###################
        meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/',i,'/col_metadata.tsv'), header = T, sep = "\t")
        meta_data <- meta_data %>% dplyr::select(Sample_ID, Group, Region)
        meta_data <- aggregate(. ~ Sample_ID, data = meta_data, FUN = function(x) x[1])
        fill <- rbind(fill, meta_data)
    }

    table(fill$Group, fill$Region)

    directory <- "/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data"

    # List all files ending with .RDS
    rds_files <- list.files(directory, pattern = "\\.RDS$", full.names = TRUE, recursive = TRUE)

    # Print the list of files
    print(rds_files)

    for(i in unique(rds_files)){
        source_file <- i
        dest_dir <- "/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/RDS_files"

        # Construct the destination path
        dest_file <- file.path(dest_dir, basename(source_file))

        # Copy the file
        file.copy(source_file, dest_file)
    }

    #sALS-PN BA4
    sALS_PN_BA4 <- subset(fill, Group == "SALS" | Group == "PN" | Group == "C9ALS")
    sALS_PN_BA4 <- subset(sALS_PN_BA4, Region == "BA4")
    sALS_PN_BA4 <- unique(sALS_PN_BA4$Sample_ID)

    #SFTLD-PN BA4
    sFTLD_PN_BA4 <- subset(fill, Group == "SFTLD" | Group == "PN" | Group == "C9FTLD")
    sFTLD_PN_BA4 <- subset(sFTLD_PN_BA4, Region == "BA4")
    sFTLD_PN_BA4 <- unique(sFTLD_PN_BA4$Sample_ID)

    #sALS-PN BA9
    sALS_PN_BA9 <- subset(fill, Group == "SALS" | Group == "PN" | Group == "C9ALS")
    sALS_PN_BA9 <- subset(sALS_PN_BA9, Region == "BA9")
    sALS_PN_BA9 <- unique(sALS_PN_BA9$Sample_ID)

    #SFTLD-PN BA9
    sFTLD_PN_BA9 <- subset(fill, Group == "SFTLD" | Group == "PN" | Group == "C9FTLD")
    sFTLD_PN_BA9 <- subset(sFTLD_PN_BA9, Region == "BA9")
    sFTLD_PN_BA9 <- unique(sFTLD_PN_BA9$Sample_ID)

    ## create directory for each and move RDS into them
    ##############
    #mkdir ALS_BA4
    ##############
    for(i in unique(sALS_PN_BA4)){
        # Define paths
        source_file <- paste0("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/RDS_files/",i,".RDS")
        dest_dir <- "/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/ALS_BA4"

        # Construct the destination path
        dest_file <- file.path(dest_dir, basename(source_file))

        # Copy the file
        file.copy(source_file, dest_file)
    }

    ##############
    #mkdir FTLD_BA4
    ##############
    for(i in unique(sFTLD_PN_BA4)){
        # Define paths
        source_file <- paste0("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/RDS_files/",i,".RDS")
        dest_dir <- "/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/FTLD_BA4"

        # Construct the destination path
        dest_file <- file.path(dest_dir, basename(source_file))

        # Copy the file
        file.copy(source_file, dest_file)
    }

    ##############
    #mkdir ALS_BA9
    ##############
    for(i in unique(sALS_PN_BA9)){
        # Define paths
        source_file <- paste0("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/RDS_files/",i,".RDS")
        dest_dir <- "/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/ALS_BA9"

        # Construct the destination path
        dest_file <- file.path(dest_dir, basename(source_file))

        # Copy the file
        file.copy(source_file, dest_file)
    }

    ##############
    #mkdir FTLD_BA9
    ##############
    for(i in unique(sFTLD_PN_BA9)){
        # Define paths
        source_file <- paste0("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/RDS_files/",i,".RDS")
        dest_dir <- "/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/FTLD_BA9"

        # Construct the destination path
        dest_file <- file.path(dest_dir, basename(source_file))

        # Copy the file
        file.copy(source_file, dest_file)
    }
##


cd /home/fiorini9/scratch/machine_learning_ALS
mkdir scRNAbox_pineda
cd scRNAbox_pineda

mkdir scrnabox_ALS_BA4
mkdir scrnabox_ALS_BA9
mkdir scrnabox_FTLD_BA4
mkdir scrnabox_FTLD_BA9

################################################################################################################################ scrnabox_ALS_BA4
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
export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 0 --method SCRNA --container TRUE

####################################
####################################
#################################### Step 1
## RDS files are here: /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/ALS_BA4

####################################
####################################
#################################### Step 2

## code -- here we  process the split objects with Step 2 of scRNAbox

    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
    mkdir step2
    cd step2
    mkdir objs

    nano step2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
    nano step2.sh  

    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=00-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=40g          # memory per cor
    #SBATCH --job-name=pineda_ALS_BA4_step2
    #SBATCH --error=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.err
    #SBATCH --output=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4

    Rscript  /step2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano step2.R

    ## load parameters
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    
    list<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/ALS_BA4",sep=""),full.names = TRUE)
    sample_name<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/ALS_BA4",sep=""))
     
    i = 1
    for(i in 1:length(sample_name)){
        
        seu <- readRDS(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/ALS_BA4/', sample_name[i] ) )
        
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
      saveRDS(seu,paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4/step2/objs/',sample_name[i], sep=''),compress=TRUE)
      
      ## print violin plot for QC metrics
      Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend()
      ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4/step2/',sample_name[i],".pdf", sep=""))
      
       ## print summary information
        sink(paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4/step2/summary_',sample_name[i],".txt", sep=""))
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
        write.csv(colnames(seu[[]]), file= paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4/step2/meta_info_',sample_name[i],".txt", sep=""))
    }

###

#cp /home/fiorini9/scratch/machine_learning_ALS/pineda_data/scrnabox_ALS_BA4/step2.R /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
#cp /home/fiorini9/scratch/machine_learning_ALS/pineda_data/scrnabox_ALS_BA4/step2.sh /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
#cp -r /home/fiorini9/scratch/machine_learning_ALS/pineda_data/scrnabox_ALS_BA4/step2 /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4

####################################
####################################
#################################### Step 3

/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4/step2/objs


export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 3 --method SCRNA --container TRUE

####################################
####################################
#################################### Step 4

/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4/step2/objs


export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 4 --method SCRNA --container TRUE

####################################
####################################
#################################### Step 5 --> merge

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 5 --method SCRNA --container TRUE

## COPY OVER THE MERGED OBJECT TO SOME OTHER LOCATION OUTSIDE OF STEP 5 RNABOX
cp /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA4/step5/objs5/seu_step5.rds /home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA4_merged.rds

####################################
####################################
#################################### Manual code to integrate by cell type class

## code Excitatory neurons
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

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
    ## Excitatory neurons
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA4_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("Ex")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_ALS_BA4_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_ALS_BA4_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_ALS_BA4_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Excitatory_neurons_ALS_BA4.rds')

    ####################################
    ## Inhibitory neurons
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA4_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("In")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_ALS_BA4_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_ALS_BA4_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_ALS_BA4_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Inhibitory_neurons_ALS_BA4.rds')

    ####################################
    ## non-neuronal
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA4_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("Glia", "Vasc")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_ALS_BA4_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_ALS_BA4_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_ALS_BA4_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Non_neuronal_ALS_BA4.rds')
##


################################################################################################################################ scrnabox_ALS_BA9
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
export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 0 --method SCRNA --container TRUE

####################################
####################################
#################################### Step 1
## RDS files are here: /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/ALS_BA9

####################################
####################################
#################################### Step 2

## code -- here we  process the split objects with Step 2 of scRNAbox
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
    mkdir step2
    cd step2
    mkdir objs

    nano step2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
    nano step2.sh  

    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=00-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=40g          # memory per cor
    #SBATCH --job-name=pineda_ALS_BA4_step2
    #SBATCH --error=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.err
    #SBATCH --output=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9

    Rscript  /step2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano step2.R

    ## load parameters
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    
    list<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/ALS_BA9",sep=""),full.names = TRUE)
    sample_name<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/ALS_BA9",sep=""))
     
    i = 1
    for(i in 1:length(sample_name)){
        
        seu <- readRDS(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/ALS_BA9/', sample_name[i] ) )
        
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
      saveRDS(seu,paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9/step2/objs/',sample_name[i], sep=''),compress=TRUE)
      
      ## print violin plot for QC metrics
      Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend()
      ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9/step2/',sample_name[i],".pdf", sep=""))
      
       ## print summary information
        sink(paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9/step2/summary_',sample_name[i],".txt", sep=""))
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
        write.csv(colnames(seu[[]]), file= paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9/step2/meta_info_',sample_name[i],".txt", sep=""))
    }

###

#cp /home/fiorini9/scratch/machine_learning_ALS/pineda_data/scrnabox_ALS_BA9/step2.R /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
#cp /home/fiorini9/scratch/machine_learning_ALS/pineda_data/scrnabox_ALS_BA9/step2.sh /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
#cp -r /home/fiorini9/scratch/machine_learning_ALS/pineda_data/scrnabox_ALS_BA9/step2 /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9

####################################
####################################
#################################### Step 3

/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9/step2/objs


export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 3 --method SCRNA --container TRUE

####################################
####################################
#################################### Step 4

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 4 --method SCRNA --container TRUE


####################################
####################################
#################################### Step 5 --> merge

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 5 --method SCRNA --container TRUE

## COPY OVER THE MERGED OBJECT TO SOME OTHER LOCATION OUTSIDE OF STEP 5 RNABOX
cp /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_ALS_BA9/step5/objs5/seu_step5.rds /home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA9_merged.rds

####################################
####################################
#################################### Manual code to integrate by cell type class

## code Excitatory neurons
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

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
    ## Excitatory neurons
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA9_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("Ex")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_ALS_BA9_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_ALS_BA9_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_ALS_BA9_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Excitatory_neurons_ALS_BA9.rds')

    ####################################
    ## Inhibitory neurons
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA9_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("In")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_ALS_BA9_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_ALS_BA9_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_ALS_BA9_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Inhibitory_neurons_ALS_BA9.rds')

    ####################################
    ## non-neuronal
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA9_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("Glia", "Vasc")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_ALS_BA9_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_ALS_BA9_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_ALS_BA9_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Non_neuronal_ALS_BA9.rds')
##




################################################################################################################################ scrnabox_FTLD_BA4
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
export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 0 --method SCRNA --container TRUE

####################################
####################################
#################################### Step 1
## RDS files are here: /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/FTLD_BA4

####################################
####################################
#################################### Step 2

## code -- here we  process the split objects with Step 2 of scRNAbox
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4
    mkdir step2
    cd step2
    mkdir objs

    nano step2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4
    nano step2.sh  

    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=00-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=40g          # memory per cor
    #SBATCH --job-name=pineda_ALS_BA4_step2
    #SBATCH --error=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.err
    #SBATCH --output=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4

    Rscript  /step2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano step2.R

    ## load parameters
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    
    list<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/FTLD_BA4",sep=""),full.names = TRUE)
    sample_name<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/FTLD_BA4",sep=""))
     
    i = 1
    for(i in 1:length(sample_name)){
        
        seu <- readRDS(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/FTLD_BA4/', sample_name[i] ) )
        
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
      saveRDS(seu,paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4/step2/objs/',sample_name[i], sep=''),compress=TRUE)
      
      ## print violin plot for QC metrics
      Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend()
      ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4/step2/',sample_name[i],".pdf", sep=""))
      
       ## print summary information
        sink(paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4/step2/summary_',sample_name[i],".txt", sep=""))
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
        write.csv(colnames(seu[[]]), file= paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4/step2/meta_info_',sample_name[i],".txt", sep=""))
    }

###


####################################
####################################
#################################### Step 3

/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4/step2/objs

############# [step3]
THREADS_ARRAY["step_3"]=4
MEM_ARRAY["step_3"]=40g
WALLTIME_ARRAY["step_3"]=00-10:00

############################################################################
# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
par_save_RNA= "yes"
par_save_metadata= "yes"


############################################################################
# If you already have a processed Seurat RDS object, and did not perform Step 2 of scRNAbox, 
# use parameter this to add the path to the directory containing you Seurat object(s). 
# Make sure that no other files/objects are present in the directory besides Seurat RDS objects. 
# Uncomment the line to activate the parameter
############################################################################
par_seurat_object= "/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4/step2/objs"


############################################################################
# Quality control parameters
# Uncomment the line to activate the parameter and add the desired value. Cells will be filtered out accordingly.
# L = lower bound threshold
# R = upper bound threshold
############################################################################
## Minimum number of unique RNA transcripts
par_nFeature_RNA_L= 700
## Maximum number of unique RNA transcripts
par_nFeature_RNA_U= 1000000
## Minimum number of total RNA transcripts
par_nCount_RNA_L= 0
## Maximum number of total RNA transcripts
par_nCount_RNA_U= 2000000
## Minimum mitochondrial RNA percentage
par_mitochondria_percent_L= 0
## Maximum mitochondrial RNA percentage
par_mitochondria_percent_U= 10
## Minimum ribosomal RNA percentage
par_ribosomal_percent_L= 0
## Maximum ribosomal RNA percentage
par_ribosomal_percent_U= 10


############################################################################
# Parameters to filter out genes
############################################################################
## If you want to filter out mitochondrial and ribosomal genes set the following parameters to "yes". If you do not want to remove them keep the default as "no".
par_remove_mitochondrial_genes= "no"
par_remove_ribosomal_genes= "no"

## If you have specific genes that you want to remove, enter a vector of the genes. Uncomment the line to activate the parameter.
#par_remove_genes= c("gene1", "gene2")


############################################################################
# Regress genes
############################################################################
## If you want to regress cell cycle genes, set the following parameters to "yes". If you do not want to regress them, keep the default as "no". Note: if you are using your own Seurat object (i.e. not from Step 2), you >
par_regress_cell_cycle_genes= "no"

## If you want to regress a custom list of genes, set the following parameters to "yes". If you do not want to regress a custom list, keep the default as "no". 
par_regress_custom_genes= "no"

## Enter the genes that you want to regress in the list below.
par_regress_genes= c("gene1", "gene2")

############################################################################
# Parameters for normalization and scaling after quality control 
############################################################################
## Normalization method
par_normalization.method= "LogNormalize"

## Scale factor
par_scale.factor= 10000

## Method for choosing the top variable features. vst, mean.var.plot (mvp), dispersion (disp). 
par_selection.method= "vst"

## Number of features to select as top variable features
par_nfeatures= 2500

## Number of most variable features to be reported in csv file
par_top= 10

## Total Number of PCs to compute and store for RunPCA
par_npcs_pca= 30


export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 3 --method SCRNA --container TRUE

####################################
####################################
#################################### Step 4
  THREADS_ARRAY["step_4"]=4
  MEM_ARRAY["step_4"]=40g
  WALLTIME_ARRAY["step_4"]=00-10:00


############################################################################
# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
par_save_RNA= "yes"
par_save_metadata= "yes"


############################################################################
# If you already have a processed Seurat RDS object, and did not perform Step 3 of scRNAbox use this to add the path to the directory containing you Seurat object(s). 
# Make sure that no other files/objects are present in the directory besides Seurat RDS objects. 
# Uncomment the line to activate the parameter
############################################################################
#par_seurat_object= "/path/to/directory/containing/seurat/object"


############################################################################
# Parameters for UMAP dimensional reduction
############################################################################
## Number of dimensions to use as input into UMAP
par_RunUMAP_dims= 25

## Number of neighbouring points to use in local approximation of manifold structure
par_RunUMAP_n.neighbors= 45


############################################################################
# Parameters for doublet detection and removal (optional)
############################################################################
## If you want to remove predicted doublets from downstream analyses set the following to "yes"
## If you want to keep predicted doublets for further analysis set the following to "no"
par_dropDN= "yes"

## Number of principal components to use as input doublet analysis. 
## This can be determined by the bend in the by elbow plot produced in Step 3
par_PCs= 25

## The number of artificial doublets to generate. DoubletFinder is largely invariant to this parameter. We suggest keeping 0.25
par_pN= 0.25

## Logical representing whether SCTransform was used during original Seurat object pre-processing
par_sct= FALSE

##rate_nExp: the doublet rate according to the number of cells
par_rate_nExp=0.05

par_sample_names = "temp"
par_expected_doublet_rate = "temp"



export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 4 --method SCRNA --container TRUE


####################################
####################################
#################################### Step 5 --> merge

############# [step5]
 THREADS_ARRAY["step_5"]=1
 MEM_ARRAY["step_5"]=150g
 WALLTIME_ARRAY["step_5"]=00-:00

############################################################################
# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
par_save_RNA= "yes"
par_save_metadata= "yes"


############################################################################
# If you already have a processed Seurat RDS object(s), and did not perform Step 4 of scRNAbox use this to add the path to the directory containing you Seurat object(s). 
# Make sure that no other files/objects are present in the directory besides Seurat RDS objects. 
# Uncomment the line to activate the parameter
############################################################################
#par_seurat_object= "/path/to/directory/containing/seurat/objects"


############################################################################
# If you only have one Seurat object and want to skip integration set the following to "yes"
############################################################################
par_one_seurat= "no"


############################################################################
# If you have multiple Seurat objects, choose whether you want to integrate or merge the objects
############################################################################
## Integrate Seurat objects
par_integrate_seurat= "no"

## Merge Seurat objects
par_merge_seurat= "yes"


############################################################################
# Parameters for normalization and scaling
# Even if you opt to skip integration, adjust the following parameters 
############################################################################
## Assay to perform normalization and scaling (prior to integration). For most use cases this will be RNA
par_DefaultAssay= "RNA"

## Normalization method
par_normalization.method= "LogNormalize"

## Scale factor
par_scale.factor= 10000


############################################################################
# Parameters for integration
############################################################################
## Method for detecting top variable features. vst, mean.var.plot (mvp), dispersion (disp)
par_selection.method= "vst"

## Number of features to select as top variable features for integration
par_nfeatures= 2500

## Which dimensions to use from the CCA to specify the neighbour search space
par_FindIntegrationAnchors_dim= 25

############################################################################
# Parameters for linear dimensional reduction
# even if you opt to skip integration, adjust the following parameters 
############################################################################
## Total Number of PCs to compute and store for RunPCA
par_RunPCA_npcs= 30

## Which dimensions to use as input features for RunUMAP
par_RunUMAP_dims= 25

## The number of neighbouring points used in local approximations of manifold structure.
par_RunUMAP_n.neighbors= 45

## Whether or not to perform JackStraw computation. This computation takes a long time.
par_compute_jackstraw= "no"

export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 5 --method SCRNA --container TRUE

## COPY OVER THE MERGED OBJECT TO SOME OTHER LOCATION OUTSIDE OF STEP 5 RNABOX
cp /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA4/step5/objs5/seu_step5.rds /home/fiorini9/scratch/machine_learning_ALS/base_objects/FTLD_BA4_merged.rds

####################################
####################################
#################################### Manual code to integrate by cell type class

## code Excitatory neurons
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

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
    ## Excitatory neurons
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/FTLD_BA4_merged.rds')
    length(unique(seu$Sample_ID))
    length(unique(seu$Donor))
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("Ex")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    table(seu_lim$Sample_ID, seu_lim$Donor)
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_FTLD_BA4_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_FTLD_BA4_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_FTLD_BA4_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Excitatory_neurons_FTLD_BA4.rds')

    ####################################
    ## Inhibitory neurons
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/FTLD_BA4_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("In")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_FTLD_BA4_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_FTLD_BA4_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_FTLD_BA4_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Inhibitory_neurons_FTLD_BA4.rds')

    ####################################
    ## non-neuronal
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/FTLD_BA4_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("Glia", "Vasc")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_FTLD_BA4_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_FTLD_BA4_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_FTLD_BA4_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Non_neuronal_FTLD_BA4.rds')
##




################################################################################################################################ scrnabox_FTLD_BA9
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
export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 0 --method SCRNA --container TRUE

####################################
####################################
#################################### Step 1
## RDS files are here: /home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/FTLD_BA9

####################################
####################################
#################################### Step 2

## code -- here we  process the split objects with Step 2 of scRNAbox
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
    mkdir step2
    cd step2
    mkdir objs

    nano step2.sh

    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
    nano step2.sh  

    #!/bin/bash  
    #SBATCH --account=def-sfarhan
    #SBATCH --time=00-12:00           # time (DD-HH:MM)
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=40g          # memory per cor
    #SBATCH --job-name=pineda_ALS_BA9_step2
    #SBATCH --error=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.err
    #SBATCH --output=/lustre03/project/6070393/COMMON/samplepooling/MF/PD_machine_learning/temp_error/job.%x-%j.out

    module load StdEnv/2020 
    module load r/4.2.2 
    R

    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9

    Rscript  /step2.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano step2.R

    ## load parameters
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    
    list<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/FTLD_BA9",sep=""),full.names = TRUE)
    sample_name<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/FTLD_BA9",sep=""))
     
    i = 1
    for(i in 1:length(sample_name)){
        
        seu <- readRDS(paste0('/home/fiorini9/scratch/machine_learning_ALS/pineda_data/data_download/processed_data/FTLD_BA9/', sample_name[i] ) )
        
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
      saveRDS(seu,paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9/step2/objs/',sample_name[i], sep=''),compress=TRUE)
      
      ## print violin plot for QC metrics
      Seurat::VlnPlot(seu, group.by= "orig.ident", features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo"), pt.size = 0.01,ncol = 3,raster = FALSE) + NoLegend()
      ggsave(paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9/step2/',sample_name[i],".pdf", sep=""))
      
       ## print summary information
        sink(paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9/step2/summary_',sample_name[i],".txt", sep=""))
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
        write.csv(colnames(seu[[]]), file= paste('/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9/step2/meta_info_',sample_name[i],".txt", sep=""))
    }

###

#cp /home/fiorini9/scratch/machine_learning_ALS/pineda_data/scrnabox_FTLD_BA9/step2.R /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
#cp /home/fiorini9/scratch/machine_learning_ALS/pineda_data/scrnabox_FTLD_BA9/step2.sh /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
#cp -r /home/fiorini9/scratch/machine_learning_ALS/pineda_data/scrnabox_FTLD_BA9/step2 /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9

####################################
####################################
#################################### Step 3

/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9/step2/objs

############# [step3]
THREADS_ARRAY["step_3"]=4
MEM_ARRAY["step_3"]=40g
WALLTIME_ARRAY["step_3"]=00-10:00

############################################################################
# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
par_save_RNA= "yes"
par_save_metadata= "yes"


############################################################################
# If you already have a processed Seurat RDS object, and did not perform Step 2 of scRNAbox, 
# use parameter this to add the path to the directory containing you Seurat object(s). 
# Make sure that no other files/objects are present in the directory besides Seurat RDS objects. 
# Uncomment the line to activate the parameter
############################################################################
par_seurat_object= "/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9/step2/objs"


############################################################################
# Quality control parameters
# Uncomment the line to activate the parameter and add the desired value. Cells will be filtered out accordingly.
# L = lower bound threshold
# R = upper bound threshold
############################################################################
## Minimum number of unique RNA transcripts
par_nFeature_RNA_L= 700
## Maximum number of unique RNA transcripts
par_nFeature_RNA_U= 1000000
## Minimum number of total RNA transcripts
par_nCount_RNA_L= 0
## Maximum number of total RNA transcripts
par_nCount_RNA_U= 2000000
## Minimum mitochondrial RNA percentage
par_mitochondria_percent_L= 0
## Maximum mitochondrial RNA percentage
par_mitochondria_percent_U= 10
## Minimum ribosomal RNA percentage
par_ribosomal_percent_L= 0
## Maximum ribosomal RNA percentage
par_ribosomal_percent_U= 10


############################################################################
# Parameters to filter out genes
############################################################################
## If you want to filter out mitochondrial and ribosomal genes set the following parameters to "yes". If you do not want to remove them keep the default as "no".
par_remove_mitochondrial_genes= "no"
par_remove_ribosomal_genes= "no"

## If you have specific genes that you want to remove, enter a vector of the genes. Uncomment the line to activate the parameter.
#par_remove_genes= c("gene1", "gene2")


############################################################################
# Regress genes
############################################################################
## If you want to regress cell cycle genes, set the following parameters to "yes". If you do not want to regress them, keep the default as "no". Note: if you are using your own Seurat object (i.e. not from Step 2), you >
par_regress_cell_cycle_genes= "no"

## If you want to regress a custom list of genes, set the following parameters to "yes". If you do not want to regress a custom list, keep the default as "no". 
par_regress_custom_genes= "no"

## Enter the genes that you want to regress in the list below.
par_regress_genes= c("gene1", "gene2")

############################################################################
# Parameters for normalization and scaling after quality control 
############################################################################
## Normalization method
par_normalization.method= "LogNormalize"

## Scale factor
par_scale.factor= 10000

## Method for choosing the top variable features. vst, mean.var.plot (mvp), dispersion (disp). 
par_selection.method= "vst"

## Number of features to select as top variable features
par_nfeatures= 2500

## Number of most variable features to be reported in csv file
par_top= 10

## Total Number of PCs to compute and store for RunPCA
par_npcs_pca= 30


export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 3 --method SCRNA --container TRUE

####################################
####################################
#################################### Step 4

THREADS_ARRAY["step_4"]=4
  MEM_ARRAY["step_4"]=40g
  WALLTIME_ARRAY["step_4"]=00-10:00


############################################################################
# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
par_save_RNA= "yes"
par_save_metadata= "yes"


############################################################################
# If you already have a processed Seurat RDS object, and did not perform Step 3 of scRNAbox use this to add the path to the directory containing you Seurat object(s). 
# Make sure that no other files/objects are present in the directory besides Seurat RDS objects. 
# Uncomment the line to activate the parameter
############################################################################
#par_seurat_object= "/path/to/directory/containing/seurat/object"


############################################################################
# Parameters for UMAP dimensional reduction
############################################################################
## Number of dimensions to use as input into UMAP
par_RunUMAP_dims= 25

## Number of neighbouring points to use in local approximation of manifold structure
par_RunUMAP_n.neighbors= 45


############################################################################
# Parameters for doublet detection and removal (optional)
############################################################################
## If you want to remove predicted doublets from downstream analyses set the following to "yes"
## If you want to keep predicted doublets for further analysis set the following to "no"
par_dropDN= "yes"

## Number of principal components to use as input doublet analysis. 
## This can be determined by the bend in the by elbow plot produced in Step 3
par_PCs= 25

## The number of artificial doublets to generate. DoubletFinder is largely invariant to this parameter. We suggest keeping 0.25
par_pN= 0.25

## Logical representing whether SCTransform was used during original Seurat object pre-processing
par_sct= FALSE

##rate_nExp: the doublet rate according to the number of cells
par_rate_nExp=0.05

par_sample_names = "temp"
par_expected_doublet_rate = "temp"



export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 4 --method SCRNA --container TRUE


####################################
####################################
#################################### Step 5 --> merge

############# [step5]
 THREADS_ARRAY["step_5"]=1
 MEM_ARRAY["step_5"]=150g
 WALLTIME_ARRAY["step_5"]=00-10:00
 

############################################################################
# If you want to save an RNA expression matrix and metadata data frame set the following to "yes"
############################################################################
par_save_RNA= "yes"
par_save_metadata= "yes"


############################################################################
# If you already have a processed Seurat RDS object(s), and did not perform Step 4 of scRNAbox use this to add the path to the directory containing you Seurat object(s). 
# Make sure that no other files/objects are present in the directory besides Seurat RDS objects. 
# Uncomment the line to activate the parameter
############################################################################
#par_seurat_object= "/path/to/directory/containing/seurat/objects"


############################################################################
# If you only have one Seurat object and want to skip integration set the following to "yes"
############################################################################
par_one_seurat= "no"


############################################################################
# If you have multiple Seurat objects, choose whether you want to integrate or merge the objects
############################################################################
## Integrate Seurat objects
par_integrate_seurat= "no"

## Merge Seurat objects
par_merge_seurat= "yes"


############################################################################
# Parameters for normalization and scaling
# Even if you opt to skip integration, adjust the following parameters 
############################################################################
## Assay to perform normalization and scaling (prior to integration). For most use cases this will be RNA
par_DefaultAssay= "RNA"

## Normalization method
par_normalization.method= "LogNormalize"

## Scale factor
par_scale.factor= 10000


############################################################################
# Parameters for integration
############################################################################
## Method for detecting top variable features. vst, mean.var.plot (mvp), dispersion (disp)
par_selection.method= "vst"

## Number of features to select as top variable features for integration
par_nfeatures= 2500

## Which dimensions to use from the CCA to specify the neighbour search space
par_FindIntegrationAnchors_dim= 25

############################################################################
# Parameters for linear dimensional reduction
# even if you opt to skip integration, adjust the following parameters 
############################################################################
## Total Number of PCs to compute and store for RunPCA
par_RunPCA_npcs= 30

## Which dimensions to use as input features for RunUMAP
par_RunUMAP_dims= 25

## The number of neighbouring points used in local approximations of manifold structure.
par_RunUMAP_n.neighbors= 45

## Whether or not to perform JackStraw computation. This computation takes a long time.
par_compute_jackstraw= "no"


export SCRNABOX_HOME=/lustre03/project/6070393/COMMON/Dark_Genome/pipeline_history/version0.1.53.00/scrnabox.slurm
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 5 --method SCRNA --container TRUE

## COPY OVER THE MERGED OBJECT TO SOME OTHER LOCATION OUTSIDE OF STEP 5 RNABOX
cp /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_pineda/scrnabox_FTLD_BA9/step5/objs5/seu_step5.rds /home/fiorini9/scratch/machine_learning_ALS/base_objects/FTLD_BA9_merged.rds

####################################
####################################
#################################### Manual code to integrate by cell type class

## code Excitatory neurons
    salloc -A def-sfarhan --time=0-8 -c 1 --mem=100g

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
    ## Excitatory neurons
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/FTLD_BA9_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("Ex")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_FTLD_BA9_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_FTLD_BA9_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_FTLD_BA9_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Excitatory_neurons_FTLD_BA9.rds')

    ####################################
    ## Inhibitory neurons
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/FTLD_BA9_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("In")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_FTLD_BA9_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_FTLD_BA9_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_FTLD_BA9_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Inhibitory_neurons_FTLD_BA9.rds')

    ####################################
    ## non-neuronal
    ####################################
    ## read in merged object
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/FTLD_BA9_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"

    ## subset to only include cell class of interest
    seu_celltype <- c("Glia", "Vasc")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$CellClass)
    xx <- xx[xx %in% c(seu_celltype)]
    Idents(seu) <- "CellClass"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$CellClass)
    dim(seu_lim)

    ## Run Harmony
    length(unique(seu_lim$Sample_ID))
    length(unique(seu_lim$Donor))
    seu_lim <- ScaleData(seu_lim)
    seu_lim <- RunPCA(seu_lim, npcs = 30)
    seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

    ## Plot
    seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
    seu_lim <- FindClusters(seu_lim)

    ## By Sample ID
    DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_FTLD_BA9_Sample_ID.pdf'))

    ## By CellType
    DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_FTLD_BA9_CellType.pdf'))

    ## By CellType subtype
    DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neuronal_FTLD_BA9_full_label.pdf'))

    ## Save RDS
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Non_neuronal_FTLD_BA9.rds')
##



################################################################################################################################ All samples merged
################################################################################################################################ 
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
## code


    screen -S merge_pineda all
    screen -r merge_pineda all

    salloc -A def-sfarhan --time=0-8 -c 1 --mem=300g

    module load StdEnv/2020 
    module load r/4.2.2 
    R



    ## code
    library(harmony)
    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(scCustomize, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    ####################################
    ####################################
    ####################################
    ####################################
    ####################################
    ####################################
    #################################### Create individual objects for all conditions (ALS, FTLD, PN), stratified by brain region (BA4, BA9)
    ## code
        ####################################
        ## ALS BA4 -- create own object
        ####################################
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA4_merged.rds')

        str(seu@meta.data)
        unique(seu@meta.data$Condition) #"ALS" "PN"

        ## subset to only include cell class of interest
        seu_condition <- c("ALS")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$Condition)
        xx <- xx[xx %in% c(seu_condition)]
        Idents(seu) <- "Condition"
        seu_ALS_BA4=subset(seu,idents=xx)
        unique(seu_ALS_BA4@meta.data$Condition)
        dim(seu_ALS_BA4)

        ####################################
        ## PN BA4 -- create own object
        ####################################
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA4_merged.rds')

        str(seu@meta.data)
        unique(seu@meta.data$Condition) #"ALS" "PN"

        ## subset to only include cell class of interest
        seu_condition <- c("PN")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$Condition)
        xx <- xx[xx %in% c(seu_condition)]
        Idents(seu) <- "Condition"
        seu_PN_BA4=subset(seu,idents=xx)
        unique(seu_PN_BA4@meta.data$Condition)
        dim(seu_PN_BA4)

        ####################################
        ## ALS BA9 -- create own object
        ####################################
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA9_merged.rds')

        str(seu@meta.data)
        unique(seu@meta.data$Condition) #"ALS" "PN"

        ## subset to only include cell class of interest
        seu_condition <- c("ALS")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$Condition)
        xx <- xx[xx %in% c(seu_condition)]
        Idents(seu) <- "Condition"
        seu_ALS_BA9=subset(seu,idents=xx)
        unique(seu_ALS_BA9@meta.data$Condition)
        dim(seu_ALS_BA9)

        ####################################
        ## PN BA9 -- create own object
        ####################################
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/ALS_BA9_merged.rds')

        str(seu@meta.data)
        unique(seu@meta.data$Condition) #"ALS" "PN"

        ## subset to only include cell class of interest
        seu_condition <- c("PN")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$Condition)
        xx <- xx[xx %in% c(seu_condition)]
        Idents(seu) <- "Condition"
        seu_PN_BA9=subset(seu,idents=xx)
        unique(seu_PN_BA9@meta.data$Condition)
        dim(seu_PN_BA9)


        ####################################
        ## FTLD BA4 -- create own object
        ####################################
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/FTLD_BA4_merged.rds')

        str(seu@meta.data)
        unique(seu@meta.data$Condition) #"FTLD" "PN"

        ## subset to only include cell class of interest
        seu_condition <- c("FTLD")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$Condition)
        xx <- xx[xx %in% c(seu_condition)]
        Idents(seu) <- "Condition"
        seu_FTLD_BA4=subset(seu,idents=xx)
        unique(seu_FTLD_BA4@meta.data$Condition)
        dim(seu_FTLD_BA4)

        ####################################
        ## FTLD BA9 -- create own object 
        ####################################
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/FTLD_BA9_merged.rds')

        str(seu@meta.data)
        unique(seu@meta.data$Condition) #"FTLD" "PN"

        ## subset to only include cell class of interest
        seu_condition <- c("FTLD")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$Condition)
        xx <- xx[xx %in% c(seu_condition)]
        Idents(seu) <- "Condition"
        seu_FTLD_BA9=subset(seu,idents=xx)
        unique(seu_FTLD_BA9@meta.data$Condition)
        dim(seu_FTLD_BA9)


        ## NOTE: We get the same number of cells for PN BA9 and BA4 when processed with FTLD and ALS. This is good. 

        rm(seu)
    ##

    ####################################
    ####################################
    ####################################
    ####################################
    ####################################
    ####################################
    #################################### Merge condition-specific objects for each brain region -- save base objects without harmony.
    ## code 
        ####################################
        ## BA4 merge
        ####################################

        seu_list <- list(seu_ALS_BA4 = seu_ALS_BA4, seu_PN_BA4 = seu_PN_BA4, seu_FTLD_BA4 = seu_FTLD_BA4)


        seu_merge_BA4 <- Merge_Seurat_List(
        list_seurat = seu_list,
        add.cell.ids = NULL,
        merge.data = TRUE,
        project = "MergeSeurat"
        )

        unique(seu_merge_BA4@meta.data$Region)
        unique(seu_merge_BA4@meta.data$Condition)

        rm(seu_ALS_BA4)
        rm(seu_PN_BA4)
        rm(seu_FTLD_BA4)

        saveRDS(seu_merge_BA4, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA4.rds')


        ####################################
        ## BA9
        ####################################

        seu_list <- list(seu_ALS_BA9 = seu_ALS_BA9, seu_PN_BA9 = seu_PN_BA9, seu_FTLD_BA9 = seu_FTLD_BA9)

        seu_merge_BA9 <- Merge_Seurat_List(
        list_seurat = seu_list,
        add.cell.ids = NULL,
        merge.data = TRUE,
        project = "MergeSeurat"
        )

        unique(seu_merge_BA9@meta.data$Region)
        unique(seu_merge_BA9@meta.data$Condition)

        rm(seu_ALS_BA9)
        rm(seu_PN_BA9)
        rm(seu_FTLD_BA9)
        
        saveRDS(seu_merge_BA9, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA9.rds')

    ##

    ####################################
    ####################################
    ####################################
    ####################################
    ####################################
    ####################################
    ####################################  BA4 -- Harmony and plot

    ####################################
    ## BA4 -- excitatory neurons 
    ####################################
    ## code
        ## NOTE: we save the object below, if want to redo for plotting purposes use the loaded object

        ## read in merged object
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA4.rds')
        str(seu@meta.data)
        unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu@meta.data$Region) #"Ex"   "In"   "Glia" "Vasc"

        ## subset to only include cell class of interest
        seu_celltype <- c("Ex")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu) <- "CellClass"
        seu_lim=subset(seu,idents=xx)
        unique(seu_lim@meta.data$CellClass)
        dim(seu_lim)

        ## Run Harmony
        length(unique(seu_lim$Sample_ID))
        length(unique(seu_lim$Donor))
        seu_lim <- ScaleData(seu_lim)    
        seu_lim<- FindVariableFeatures(seu_lim, selection.method = "vst", nfeatures = 2500)
        seu_lim <- RunPCA(seu_lim, npcs = 30)
        seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindClusters(seu_lim)

        ## By Sample ID
        DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "none",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_all_conditions_BA4_SampleID.pdf'))

        ## By CellType
        DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_all_conditions_BA4_CellType.pdf'), height = 5, width = 5)

        ## By CellType subtype
        DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_all_conditions_BA4_full_label.pdf'), height = 5, width = 7)


        ## SAVE RDS 
        saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Harmony_all_conditions_BA4_Excitatory_neurons.rds')
    ##

    ####################################
    ## BA4 -- inhibitory neurons 
    ####################################
    ## code
        ## NOTE: we save the object below, if want to redo for plotting purposes use the loaded object

        ## read in merged object
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA4.rds')
        str(seu@meta.data)
        unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu@meta.data$Region) #"Ex"   "In"   "Glia" "Vasc"

        ## subset to only include cell class of interest
        seu_celltype <- c("In")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu) <- "CellClass"
        seu_lim=subset(seu,idents=xx)
        unique(seu_lim@meta.data$CellClass)
        dim(seu_lim)

        ## Run Harmony
        length(unique(seu_lim$Sample_ID))
        length(unique(seu_lim$Donor))
        seu_lim <- ScaleData(seu_lim)    
        seu_lim<- FindVariableFeatures(seu_lim, selection.method = "vst", nfeatures = 2500)
        seu_lim <- RunPCA(seu_lim, npcs = 30)
        seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindClusters(seu_lim)

        ## By Sample ID
        DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "none",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_all_conditions_BA4_SampleID.pdf'))

        ## By CellType
        DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_all_conditions_BA4_CellType.pdf'), height = 5, width = 5)

        ## By CellType subtype
        DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_all_conditions_BA4_full_label.pdf'), height = 5, width = 7)


        ## SAVE RDS 
        saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Harmony_all_conditions_BA4_Inhibitory_neurons.rds')
    ##


    ####################################
    ## BA4 -- non-neuronal
    ####################################
    ## code
        ## NOTE: we save the object below, if want to redo for plotting purposes use the loaded object

        ## read in merged object
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA4.rds')
        str(seu@meta.data)
        unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu@meta.data$Region) #"Ex"   "In"   "Glia" "Vasc"

        ## subset to only include cell class of interest
        seu_celltype <- c("Glia", "Vasc")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu) <- "CellClass"
        seu_lim=subset(seu,idents=xx)
        unique(seu_lim@meta.data$CellClass)
        dim(seu_lim)

        ## Run Harmony
        length(unique(seu_lim$Sample_ID))
        length(unique(seu_lim$Donor))
        seu_lim <- ScaleData(seu_lim)    
        seu_lim<- FindVariableFeatures(seu_lim, selection.method = "vst", nfeatures = 2500)
        seu_lim <- RunPCA(seu_lim, npcs = 30)
        seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindClusters(seu_lim)

        ## By Sample ID
        DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "none",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neurons_all_conditions_BA4_SampleID.pdf'))

        ## By CellType
        DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neurons_all_conditions_BA4_CellType.pdf'), height = 5, width = 5)

        ## By CellType subtype
        DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neurons_all_conditions_BA4_full_label.pdf'), height = 5, width = 7)


        ## SAVE RDS 
        saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Harmony_all_conditions_BA4_Non_neurons.rds')
    ##


    ####################################
    ####################################
    ####################################
    ####################################
    ####################################
    ####################################
    ####################################  BA9 -- Harmony and plot

    ####################################
    ## BA9 -- excitatory neurons 
    ####################################
    ## code
        ## NOTE: we save the object below, if want to redo for plotting purposes use the loaded object

        ## read in merged object
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA9.rds')
        str(seu@meta.data)
        unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu@meta.data$Region) #"Ex"   "In"   "Glia" "Vasc"

        ## subset to only include cell class of interest
        seu_celltype <- c("Ex")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu) <- "CellClass"
        seu_lim=subset(seu,idents=xx)
        unique(seu_lim@meta.data$CellClass)
        dim(seu_lim)

        ## Run Harmony
        length(unique(seu_lim$Sample_ID))
        length(unique(seu_lim$Donor))
        seu_lim <- ScaleData(seu_lim)    
        seu_lim<- FindVariableFeatures(seu_lim, selection.method = "vst", nfeatures = 2500)
        seu_lim <- RunPCA(seu_lim, npcs = 30)
        seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindClusters(seu_lim)

        ## By Sample ID
        DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "none",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_all_conditions_BA9_SampleID.pdf'))

        ## By CellType
        DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_all_conditions_BA9_CellType.pdf'), height = 5, width = 5)

        ## By CellType subtype
        DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Excitatory_neurons_all_conditions_BA9_full_label.pdf'), height = 5, width = 7)


        ## SAVE RDS 
        saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Harmony_all_conditions_BA9_Excitatory_neurons.rds')
    ##

    ####################################
    ## BA9 -- inhibitory neurons 
    ####################################
    ## code
        ## NOTE: we save the object below, if want to redo for plotting purposes use the loaded object

        ## read in merged object
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA9.rds')
        str(seu@meta.data)
        unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu@meta.data$Region) #"Ex"   "In"   "Glia" "Vasc"

        ## subset to only include cell class of interest
        seu_celltype <- c("In")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu) <- "CellClass"
        seu_lim=subset(seu,idents=xx)
        unique(seu_lim@meta.data$CellClass)
        dim(seu_lim)

        ## Run Harmony
        length(unique(seu_lim$Sample_ID))
        length(unique(seu_lim$Donor))
        seu_lim <- ScaleData(seu_lim)    
        seu_lim<- FindVariableFeatures(seu_lim, selection.method = "vst", nfeatures = 2500)
        seu_lim <- RunPCA(seu_lim, npcs = 30)
        seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindClusters(seu_lim)

        ## By Sample ID
        DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "none",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_all_conditions_BA9_SampleID.pdf'))

        ## By CellType
        DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_all_conditions_BA9_CellType.pdf'), height = 5, width = 5)

        ## By CellType subtype
        DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Inhibitory_neurons_all_conditions_BA9_full_label.pdf'), height = 5, width = 7)


        ## SAVE RDS 
        saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Harmony_all_conditions_BA9_Inhibitory_neurons.rds')
    ##


    ####################################
    ## BA9 -- non-neuronal
    ####################################
    ## code
        ## NOTE: we save the object below, if want to redo for plotting purposes use the loaded object

        ## read in merged object
        seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA9.rds')
        str(seu@meta.data)
        unique(seu@meta.data$CellClass) #"Ex"   "In"   "Glia" "Vasc"
        unique(seu@meta.data$Region) #"Ex"   "In"   "Glia" "Vasc"

        ## subset to only include cell class of interest
        seu_celltype <- c("Glia", "Vasc")
            
        ## select only cell type of interest
        xx <- unique(seu@meta.data$CellClass)
        xx <- xx[xx %in% c(seu_celltype)]
        Idents(seu) <- "CellClass"
        seu_lim=subset(seu,idents=xx)
        unique(seu_lim@meta.data$CellClass)
        dim(seu_lim)

        ## Run Harmony
        length(unique(seu_lim$Sample_ID))
        length(unique(seu_lim$Donor))
        seu_lim <- ScaleData(seu_lim)    
        seu_lim<- FindVariableFeatures(seu_lim, selection.method = "vst", nfeatures = 2500)
        seu_lim <- RunPCA(seu_lim, npcs = 30)
        seu_lim <- RunHarmony(seu_lim, group.by.vars = "Sample_ID")  # Replace "batch" with your column name

        ## Plot
        seu_lim <- RunUMAP(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindNeighbors(seu_lim, reduction = "harmony", dims = 1:30)
        seu_lim <- FindClusters(seu_lim)

        ## By Sample ID
        DimPlot(seu_lim, reduction = "umap", group.by = "Sample_ID", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "none",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neurons_all_conditions_BA9_SampleID.pdf'))

        ## By CellType
        DimPlot(seu_lim, reduction = "umap", group.by = "CellType", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neurons_all_conditions_BA9_CellType.pdf'), height = 5, width = 5)

        ## By CellType subtype
        DimPlot(seu_lim, reduction = "umap", group.by = "full_label", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
        theme(axis.text.x = element_text(size =10),
        axis.text.y = element_text(size =10),
        axis.title = element_text(face="bold", size =10),
        legend.position = "right",
        plot.title = element_text( size =10)) + 
        xlab("UMAP1") + ylab("UMAP2")
        ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/figures/data_prep/Non_neurons_all_conditions_BA9_full_label.pdf'), height = 5, width = 7)


        ## SAVE RDS 
        saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Harmony_all_conditions_BA9_Non_neurons.rds')
    ##

###