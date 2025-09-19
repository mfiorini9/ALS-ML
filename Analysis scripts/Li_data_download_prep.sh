salloc -A def-sfarhan --time=0-2 -c 1 --mem=20g


# Notes: 
# All sequencing data generated in this study have been deposited in the Gene Expression Omnibus repository under the accession code GSE219281. 
# Curated snRNA-seq cell-by-gene tables are provided in the following Zenodo repository:123 https://doi.org/10.5281/zenodo.8190317. 
# Metadata for each nucleus in snRNA-seq and snATAC-seq are provided in Supplementary Dataset 2 and Supplementary Dataset 6. 
# An IGV browser session showing ChIP-Seq, snRNA-seq, and snATAC-seq tracks, and a UCSC single cell browser session showing the snRNA-seq data, are available at https://brainome.ucsd.edu/C9_ALS_FTD/. 
# Source data are provided with this paper.

#################################################################################################################################### Data download from GEO
################################################################################################################################ 

## code 
    #####################
    # snRNA, C9-ALS donor A1, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A1_FCx
    cd C9ALS_A1_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781906/suppl/GSM6781906%5FsnRNA%5FC9ALS%5FA1%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781906/suppl/GSM6781906%5FsnRNA%5FC9ALS%5FA1%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781906/suppl/GSM6781906%5FsnRNA%5FC9ALS%5FA1%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # C9-ALS donor A1, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A1_MCx
    cd C9ALS_A1_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781907/suppl/GSM6781907%5FsnRNA%5FC9ALS%5FA1%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781907/suppl/GSM6781907%5FsnRNA%5FC9ALS%5FA1%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781907/suppl/GSM6781907%5FsnRNA%5FC9ALS%5FA1%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-ALS donor A2, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A2_FCx
    cd C9ALS_A2_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781908/suppl/GSM6781908%5FsnRNA%5FC9ALS%5FA2%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781908/suppl/GSM6781908%5FsnRNA%5FC9ALS%5FA2%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781908/suppl/GSM6781908%5FsnRNA%5FC9ALS%5FA2%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-ALS donor A2, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A2_MCx
    cd C9ALS_A2_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781909/suppl/GSM6781909%5FsnRNA%5FC9ALS%5FA2%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781909/suppl/GSM6781909%5FsnRNA%5FC9ALS%5FA2%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781909/suppl/GSM6781909%5FsnRNA%5FC9ALS%5FA2%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-ALS donor A3, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A3_FCx
    cd C9ALS_A3_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781910/suppl/GSM6781910%5FsnRNA%5FC9ALS%5FA3%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781910/suppl/GSM6781910%5FsnRNA%5FC9ALS%5FA3%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781910/suppl/GSM6781910%5FsnRNA%5FC9ALS%5FA3%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-ALS donor A3, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A3_MCx
    cd C9ALS_A3_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781911/suppl/GSM6781911%5FsnRNA%5FC9ALS%5FA3%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781911/suppl/GSM6781911%5FsnRNA%5FC9ALS%5FA3%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781911/suppl/GSM6781911%5FsnRNA%5FC9ALS%5FA3%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-ALS donor A4, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A4_FCx
    cd C9ALS_A4_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781912/suppl/GSM6781912%5FsnRNA%5FC9ALS%5FA4%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781912/suppl/GSM6781912%5FsnRNA%5FC9ALS%5FA4%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781912/suppl/GSM6781912%5FsnRNA%5FC9ALS%5FA4%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-ALS donor A4, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A4_MCx
    cd C9ALS_A4_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781913/suppl/GSM6781913%5FsnRNA%5FC9ALS%5FA4%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781913/suppl/GSM6781913%5FsnRNA%5FC9ALS%5FA4%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781913/suppl/GSM6781913%5FsnRNA%5FC9ALS%5FA4%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-ALS donor A5, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A5_FCx
    cd C9ALS_A5_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781914/suppl/GSM6781914%5FsnRNA%5FC9ALS%5FA5%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781914/suppl/GSM6781914%5FsnRNA%5FC9ALS%5FA5%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781914/suppl/GSM6781914%5FsnRNA%5FC9ALS%5FA5%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-ALS donor A5, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A5_MCx
    cd C9ALS_A5_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781915/suppl/GSM6781915%5FsnRNA%5FC9ALS%5FA5%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781915/suppl/GSM6781915%5FsnRNA%5FC9ALS%5FA5%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781915/suppl/GSM6781915%5FsnRNA%5FC9ALS%5FA5%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-ALS donor A6, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A6_FCx
    cd C9ALS_A6_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781916/suppl/GSM6781916%5FsnRNA%5FC9ALS%5FA6%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781916/suppl/GSM6781916%5FsnRNA%5FC9ALS%5FA6%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781916/suppl/GSM6781916%5FsnRNA%5FC9ALS%5FA6%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-ALS donor A6, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9ALS_A6_MCx
    cd C9ALS_A6_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781917/suppl/GSM6781917%5FsnRNA%5FC9ALS%5FA6%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781917/suppl/GSM6781917%5FsnRNA%5FC9ALS%5FA6%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781917/suppl/GSM6781917%5FsnRNA%5FC9ALS%5FA6%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-FTD donor F1, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9FTD_F1_FCx
    cd C9FTD_F1_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781918/suppl/GSM6781918%5FsnRNA%5FC9FTD%5FF1%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781918/suppl/GSM6781918%5FsnRNA%5FC9FTD%5FF1%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781918/suppl/GSM6781918%5FsnRNA%5FC9FTD%5FF1%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-FTD donor F1, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9FTD_F1_MCx
    cd C9FTD_F1_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781919/suppl/GSM6781919%5FsnRNA%5FC9FTD%5FF1%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781919/suppl/GSM6781919%5FsnRNA%5FC9FTD%5FF1%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781919/suppl/GSM6781919%5FsnRNA%5FC9FTD%5FF1%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-FTD donor F2, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9FTD_F2_FCx
    cd C9FTD_F2_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781920/suppl/GSM6781920%5FsnRNA%5FC9FTD%5FF2%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781920/suppl/GSM6781920%5FsnRNA%5FC9FTD%5FF2%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781920/suppl/GSM6781920%5FsnRNA%5FC9FTD%5FF2%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-FTD donor F2, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9FTD_F2_MCx
    cd C9FTD_F2_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781921/suppl/GSM6781921%5FsnRNA%5FC9FTD%5FF2%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781921/suppl/GSM6781921%5FsnRNA%5FC9FTD%5FF2%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781921/suppl/GSM6781921%5FsnRNA%5FC9FTD%5FF2%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-FTD donor F3, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9FTD_F3_FCx
    cd C9FTD_F3_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781922/suppl/GSM6781922%5FsnRNA%5FC9FTD%5FF3%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781922/suppl/GSM6781922%5FsnRNA%5FC9FTD%5FF3%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781922/suppl/GSM6781922%5FsnRNA%5FC9FTD%5FF3%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-FTD donor F3, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9FTD_F3_MCx
    cd C9FTD_F3_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781923/suppl/GSM6781923%5FsnRNA%5FC9FTD%5FF3%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781923/suppl/GSM6781923%5FsnRNA%5FC9FTD%5FF3%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781923/suppl/GSM6781923%5FsnRNA%5FC9FTD%5FF3%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-FTD donor F4, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9FTD_F4_FCx
    cd C9FTD_F4_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781924/suppl/GSM6781924%5FsnRNA%5FC9FTD%5FF4%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781924/suppl/GSM6781924%5FsnRNA%5FC9FTD%5FF4%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781924/suppl/GSM6781924%5FsnRNA%5FC9FTD%5FF4%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-FTD donor F4, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9FTD_F4_MCx
    cd C9FTD_F4_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781925/suppl/GSM6781925%5FsnRNA%5FC9FTD%5FF4%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781925/suppl/GSM6781925%5FsnRNA%5FC9FTD%5FF4%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781925/suppl/GSM6781925%5FsnRNA%5FC9FTD%5FF4%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-FTD donor F5, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9FTD_F5_FCx
    cd C9FTD_F5_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781926/suppl/GSM6781926%5FsnRNA%5FC9FTD%5FF5%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781926/suppl/GSM6781926%5FsnRNA%5FC9FTD%5FF5%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781926/suppl/GSM6781926%5FsnRNA%5FC9FTD%5FF5%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, C9-FTD donor F5, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir C9FTD_F5_MCx
    cd C9FTD_F5_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781927/suppl/GSM6781927%5FsnRNA%5FC9FTD%5FF5%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781927/suppl/GSM6781927%5FsnRNA%5FC9FTD%5FF5%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781927/suppl/GSM6781927%5FsnRNA%5FC9FTD%5FF5%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, Control donor C1, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C1_FCx
    cd Control_C1_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781928/suppl/GSM6781928%5FsnRNA%5FControl%5FC1%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781928/suppl/GSM6781928%5FsnRNA%5FControl%5FC1%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781928/suppl/GSM6781928%5FsnRNA%5FControl%5FC1%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"


    #####################
    # snRNA, Control donor C1, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C1_MCx
    cd Control_C1_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781929/suppl/GSM6781929%5FsnRNA%5FControl%5FC1%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781929/suppl/GSM6781929%5FsnRNA%5FControl%5FC1%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781929/suppl/GSM6781929%5FsnRNA%5FControl%5FC1%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"
    
    
    #####################
    # snRNA, Control donor C2, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C2_FCx
    cd Control_C2_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781930/suppl/GSM6781930%5FsnRNA%5FControl%5FC2%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781930/suppl/GSM6781930%5FsnRNA%5FControl%5FC2%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781930/suppl/GSM6781930%5FsnRNA%5FControl%5FC2%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"
    
    
    #####################
    # snRNA, Control donor C2, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C2_MCx
    cd Control_C2_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781931/suppl/GSM6781931%5FsnRNA%5FControl%5FC2%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781931/suppl/GSM6781931%5FsnRNA%5FControl%5FC2%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781931/suppl/GSM6781931%5FsnRNA%5FControl%5FC2%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"
    
    
    #####################
    # snRNA, Control donor C3, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C3_FCx
    cd Control_C3_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781932/suppl/GSM6781932%5FsnRNA%5FControl%5FC3%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781932/suppl/GSM6781932%5FsnRNA%5FControl%5FC3%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781932/suppl/GSM6781932%5FsnRNA%5FControl%5FC3%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"
    
    
    #####################
    # snRNA, Control donor C3, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C3_MCx
    cd Control_C3_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781933/suppl/GSM6781933%5FsnRNA%5FControl%5FC3%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781933/suppl/GSM6781933%5FsnRNA%5FControl%5FC3%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781933/suppl/GSM6781933%5FsnRNA%5FControl%5FC3%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"
    
    
    #####################
    # snRNA, Control donor C4, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C4_FCx
    cd Control_C4_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781934/suppl/GSM6781934%5FsnRNA%5FControl%5FC4%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781934/suppl/GSM6781934%5FsnRNA%5FControl%5FC4%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781934/suppl/GSM6781934%5FsnRNA%5FControl%5FC4%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"
    
    
    #####################
    # snRNA, Control donor C4, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C4_MCx
    cd Control_C4_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781935/suppl/GSM6781935%5FsnRNA%5FControl%5FC4%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781935/suppl/GSM6781935%5FsnRNA%5FControl%5FC4%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781935/suppl/GSM6781935%5FsnRNA%5FControl%5FC4%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"
    
    
    #####################
    # snRNA, Control donor C5, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C5_FCx
    cd Control_C5_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781936/suppl/GSM6781936%5FsnRNA%5FControl%5FC5%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781936/suppl/GSM6781936%5FsnRNA%5FControl%5FC5%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781936/suppl/GSM6781936%5FsnRNA%5FControl%5FC5%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"
    
    
    #####################
    # snRNA, Control donor C5, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C5_MCx
    cd Control_C5_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781937/suppl/GSM6781937%5FsnRNA%5FControl%5FC5%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781937/suppl/GSM6781937%5FsnRNA%5FControl%5FC5%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781937/suppl/GSM6781937%5FsnRNA%5FControl%5FC5%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"
    
    
    #####################
    # snRNA, Control donor C6, Frontal cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C6_FCx
    cd Control_C6_FCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781938/suppl/GSM6781938%5FsnRNA%5FControl%5FC6%5FFrontalCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781938/suppl/GSM6781938%5FsnRNA%5FControl%5FC6%5FFrontalCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781938/suppl/GSM6781938%5FsnRNA%5FControl%5FC6%5FFrontalCortex%5Fraw%5Fmatrix.mtx.gz"
    
    
    #####################
    # snRNA, Control donor C6, Motor cortex
    #####################
    cd /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download
    mkdir Control_C6_MCx
    cd Control_C6_MCx

    ## barcodes
    wget -O barcodes.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781939/suppl/GSM6781939%5FsnRNA%5FControl%5FC6%5FMotorCortex%5Fraw%5Fbarcodes.tsv.gz"
    ## features
    wget -O features.tsv.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781939/suppl/GSM6781939%5FsnRNA%5FControl%5FC6%5FMotorCortex%5Fraw%5Ffeatures.tsv.gz"
    ## matrix
    wget -O matrix.mtx.gz "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6781nnn/GSM6781939/suppl/GSM6781939%5FsnRNA%5FControl%5FC6%5FMotorCortex%5Fraw%5Fmatrix.mtx.gz"
###


################################################################################################################################ scrnabox for Li data -- we are processing MCx and FCx together
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
export SCRNABOX_PWD=/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_li/scrnabox
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_li/scrnabox
bash $SCRNABOX_HOME/launch_scrnabox.sh -h 
bash $SCRNABOX_HOME/launch_scrnabox.sh -d ${SCRNABOX_PWD} --steps 0 --method SCRNA --container TRUE



####################################
####################################
#################################### Create initial Seurat objects -- because we are using the nuclei from the metadata file, we do not have to do quality contro. We can skip to the merge. 

## code

    cd /home/fiorini9/scratch/machine_learning_ALS_ALS/pineda_data/data_download

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

    cd /home/fiorini9/scratch/machine_learning_ALS_ALS/pineda_data/data_download

    Rscript /home/fiorini9/scratch/machine_learning_ALS_ALS/pineda_data/data_download/data_prep.R


    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    nano data_prep.R


    library(Seurat, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggplot2, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(dplyr)
    library(foreach, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(doParallel, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(Matrix, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")
    library(ggpubr, lib="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2")

    ##########################################
    ## Create a processing function for FCx
    ##########################################
    list<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/li_data/data_download",sep=""),full.names = TRUE)
    list <- list[grep("FCx", list)] ## CHANGE

    sample_name<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/li_data/data_download",sep=""))
    sample_name <- sample_name[grep("FCx", sample_name)] ## CHANGE

    ## no high quality cells for -- as per supplemental Table 2 from manuscript
    #'C9FTD_F2_FCx'

    ## remove from the loop
    sample_name <- sample_name[sample_name != "C9FTD_F2_FCx"]


    for(i in unique(sample_name)){
        print(i)
        
        ###################
        ## Test create Seurat object
        ###################
        sparse_matrix <- Seurat::Read10X(data.dir = paste0('/home/fiorini9/scratch/machine_learning_ALS/li_data/data_download/',i))
        seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=1,min.features= 1, project = i) # only keep cells expressing atleast one gene
        dim(seurat_object)

        ## import
        meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/li_data/Li_metadata.csv'), header = T, sep = ",")
        j = gsub("_FCx", "", i)
        sample_keep_j <- paste0('snRNA_', j, "_MedialFrontalCortex") ## CHANGE
        meta_data <- subset(meta_data, sample == sample_keep_j)
        
        ## check overlap
        colnames(seurat_object)
        meta_data$cell_barcode %in% colnames(seurat_object)
        
        ## only retain barcodes from Seurat that are present in the metdata
        seurat_subset <- subset(seurat_object, cells = meta_data$cell_barcode)
        dim(seurat_subset)

        ## merge
        colnames(seurat_subset)
        meta_data$cell_barcode %in% colnames(seurat_subset)
        rownames(meta_data) <- meta_data$cell_barcode
        seurat_subset <- AddMetaData(seurat_subset, meta_data)

        saveRDS(seurat_subset, paste0('/home/fiorini9/scratch/machine_learning_ALS/li_data/data_download/processed_data/',i,'.RDS'))
    }

    ##########################################
    ## Create a processing function for MCx
    ##########################################
    list<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/li_data/data_download",sep=""),full.names = TRUE)
    list <- list[grep("MCx", list)] ## CHANGE

    sample_name<-dir(path = paste("/home/fiorini9/scratch/machine_learning_ALS/li_data/data_download",sep=""))
    sample_name <- sample_name[grep("MCx", sample_name)] ## CHANGE

    ## no high quality cells for -- as per supplemental Table 2 from manuscript
    #'C9FTD_F2_FCx'

    ## remove from the loop
    #sample_name <- sample_name[sample_name != "C9FTD_F2_FCx"]


    for(i in unique(sample_name)){
        print(i)
        
        ###################
        ## Test create Seurat object
        ###################
        sparse_matrix <- Seurat::Read10X(data.dir = paste0('/home/fiorini9/scratch/machine_learning_ALS/li_data/data_download/',i))
        seurat_object <- Seurat::CreateSeuratObject(counts = sparse_matrix, min.cells=1,min.features= 1, project = i) # only keep cells expressing atleast one gene
        dim(seurat_object)

        ## import
        meta_data <- read.delim(paste0('/home/fiorini9/scratch/machine_learning_ALS/li_data/Li_metadata.csv'), header = T, sep = ",")
        j = gsub("_MCx", "", i)
        sample_keep_j <- paste0('snRNA_', j, "_MotorCortex") ## CHANGE
        meta_data <- subset(meta_data, sample == sample_keep_j)
        
        ## check overlap
        colnames(seurat_object)
        meta_data$cell_barcode %in% colnames(seurat_object)
        
        ## only retain barcodes from Seurat that are present in the metdata
        seurat_subset <- subset(seurat_object, cells = meta_data$cell_barcode)
        dim(seurat_subset)

        ## merge
        colnames(seurat_subset)
        meta_data$cell_barcode %in% colnames(seurat_subset)
        rownames(meta_data) <- meta_data$cell_barcode
        seurat_subset <- AddMetaData(seurat_subset, meta_data)

        saveRDS(seurat_subset, paste0('/home/fiorini9/scratch/machine_learning_ALS/li_data/data_download/processed_data/',i,'.RDS'))
    }
##

####################################
####################################
#################################### Create step 4 objs directory and transfer the objects so that we can merge directly with scRNAbox
cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_li/scrnabox
mkdir step4
cd step4
mkdir objs4
cd objs4

cp /home/fiorini9/scratch/machine_learning_ALS/li_data/data_download/processed_data/*.* /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_li/scrnabox/step4/objs4

## change .RDS to .rds

cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_li/scrnabox/step4/objs4
for f in *.RDS; do mv -- "$f" "${f%.RDS}.rds"; done


####################################
####################################
#################################### Step 5 --> merge

## we have to run this manually because the object does not have a "Sample_ID" metdata column
## R code for manual running of Step 5 -- merge.
    ## load parameters
    output_dir="/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_li/scrnabox"
    r_lib_path="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2"
    .libPaths(r_lib_path)

    ## load library
    packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix','scCustomize')
    invisible(lapply(packages, library, character.only = TRUE))


    ## load parameters
    source(paste(output_dir,'/job_info/parameters/step5_par.txt',sep=""))


    ## load Seurat objects
    if (exists("par_seurat_object")) {                                                  
        sample_name<-list.files(path = par_seurat_object)
        sample_nameb<-gsub(".rds","",sample_name)
        if(length(sample_name)<1) {
        print("You do not have any existing Seurat object")
        }
        for (i in 1:length(sample_name)) {
            if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
                print(c(sample_name[i],"is not R rds"))
            }
            }  
    } else {
        sample_name<-list.files(path = paste(output_dir, "/step4/objs4",sep=""),pattern = "*.rds")
        sample_nameb<-gsub(".rds","",sample_name)
        if(length(sample_name)<1) {
        print("You do not have any object from step 4 ")
        }
        for (i in 1:length(sample_name)) {
            if (!grepl(".rds",tolower(sample_name[i]), fixed = TRUE)){
                print(c(sample_name[i],"is not R rds"))
            }
            }
    }  

    ## create empty list to be populated by existing Seurat objects
    seu_list<-list()

    ################################################################################################################################################
    ## if users just want to merge their Seurat objects
    ################################################################################################################################################
    ## normalize exisiting Seurat objects
    seu_list<-foreach (i_s=1:length(sample_name)) %do% {  
        seu<-readRDS(paste(output_dir,'/step4/objs4/',sample_name[i_s], sep=""))
        DefaultAssay(seu) <- par_DefaultAssay
        seu <- Seurat::NormalizeData(seu,normalization.method = par_normalization.method,scale.factor =par_scale.factor)
    }  

    ## merge Seurat objects
    sample_name_2 <- gsub("\\..*","",sample_name)

    seu_int <- Merge_Seurat_List(
    list_seurat = seu_list,
    add.cell.ids = dput(as.character(sample_name_2)),
    merge.data = TRUE,
    project = "MergeSeurat"
    )

    ## set default assay to integrated
    Seurat::DefaultAssay(seu_int) <- "RNA"

    ## find variable features
    seu_int<- FindVariableFeatures(seu_int, selection.method = par_selection.method, nfeatures = par_nfeatures)

    ## scale integrated assay
    seu_int <- ScaleData(seu_int, verbose = FALSE)

    ## run PCA and UMAP on integrated Seurat object
    seu_int <- RunPCA(seu_int, npcs = par_RunPCA_npcs, verbose = FALSE)
    seu_int <- RunUMAP(seu_int, dims = 1:par_RunUMAP_dims, n.neighbors =par_RunUMAP_n.neighbors)

    ## print PCA
    DimPlot(seu_int, reduction = "pca", group.by="orig.ident", raster = FALSE )
    ggsave(paste(output_dir,'/step5/figs5',"/merge_DimPlot_pca.pdf", sep=""))

    ## print elbow plot
    ElbowPlot(seu_int, ndims = par_RunPCA_npcs)
    ggsave(paste(output_dir,'/step5/figs5',"/merge_elbow.pdf", sep=""))

    ## print UMAP
    DimPlot(seu_int, reduction = "umap", group.by="orig.ident", raster = FALSE)
    ggsave(paste(output_dir,'/step5/figs5',"/merge_DimPlot_umap.pdf", sep=""))

    ## save Seurat object as RDS
    saveRDS(seu_int, paste(output_dir,'/step5/objs5',"/seu_step5.rds", sep=""))
##

####################################
####################################
#################################### We need to transfer the cell type annotations from Pineda to Li. 
## code
    cd /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_li/scrnabox
    mkdir step7
    cd step7
    mkdir objs7
    mkdir figs7
    mkdir info7

    salloc -A def-sfarhan --time=0-8 -c 1 --mem=150g
    module load StdEnv/2020 
    module load r/4.2.2 

    output_dir="/home/fiorini9/scratch/machine_learning_ALS/scRNAbox_li/scrnabox"
    r_lib_path="/lustre03/project/6070393/COMMON/Dark_Genome/R/x86_64-pc-linux-gnu-library/4.2"
    .libPaths(r_lib_path)

    ## load library
    packages<-c('Seurat','ggplot2', 'dplyr', 'foreach', 'doParallel','Matrix','scCustomize')
    invisible(lapply(packages, library, character.only = TRUE))

    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_C9ALS_C9FTD_BA4_BA9_merged.rds')

    unique(seu@meta.data$annotation_cell_class)
    unique(seu@meta.data$annotation_major_cell_type)
    unique(seu@meta.data$annotation_cell_subtype)

    DimPlot(seu, reduction = "umap", group.by="annotation_major_cell_type", raster = FALSE)
    ggsave('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf')

    ##############################################################################
    ## use Pineda data to transfer labels. 
    ## BA4
    ##############################################################################
    ## load reference Seurat object
    par_reference_name = "Pineda_BA4"
    reference0 <-readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA4.rds')
    dim(reference0)
    length(unique(reference0@meta.data$CellType)) #19

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


    ##############################################################################
    ## use Pineda data to transfer labels. 
    ## BA9
    ##############################################################################
    ## load reference Seurat object
    par_reference_name = "Pineda_BA9"
    reference0 <-readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/All_conditions_BA9.rds')
    dim(reference0)
    length(unique(reference0@meta.data$CellType)) #19

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
cp /home/fiorini9/scratch/machine_learning_ALS/scRNAbox_li/scrnabox/step7/objs7/seu_step7.rds /home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_C9ALS_C9FTD_BA4_BA9_merged.rds


####################################
####################################
#################################### Check tranfer annotations

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
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_C9ALS_C9FTD_BA4_BA9_merged.rds')
    str(seu@meta.data)
    unique(seu@meta.data$annotation_cell_class) #"Exc_neuron" "Inh_neuron" "Non_neuron"

    ## Run Harmony
    length(unique(seu$sample))
    length(unique(seu$donor_id))
    seu <- ScaleData(seu)
    seu <- RunPCA(seu, npcs = 20)
    seu <- RunHarmony(seu, group.by.vars = "sample")  # Replace "batch" with your column name

    ## Plot
    seu <- RunUMAP(seu, reduction = "harmony", dims = 1:20)
    seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:20)
    seu <- FindClusters(seu)

    ## By CellType subtype
    DimPlot(seu, reduction = "umap", group.by = "Pineda_BA9_predictions", label = FALSE, label.size = 3,repel = TRUE, raster=FALSE) +
    theme(axis.text.x = element_text(size =10),
    axis.text.y = element_text(size =10),
    axis.title = element_text(face="bold", size =10),
    legend.position = "right",
    plot.title = element_text( size =10)) + 
    xlab("UMAP1") + ylab("UMAP2")
    ggsave(paste0('/home/fiorini9/scratch/machine_learning_ALS/temp_figures/temp.pdf'))
       
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
    
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_C9ALS_C9FTD_BA4_BA9_merged.rds')

    ## only select BA4
    str(seu@meta.data)
    unique(seu@meta.data$brain_region) #"ALS" "PN"

    ## subset to only include cell class of interest
    seu_region <- c("motor cortex")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$brain_region)
    xx <- xx[xx %in% c(seu_region)]
    Idents(seu) <- "brain_region"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$brain_region)
    dim(seu_lim)

    ## fix column names
    df <- seu_lim@meta.data
    df <- df %>% dplyr::select(Pineda_BA4_predictions, sample, diagnosis, brain_region, donor_id)
    colnames(df)
    colnames(df) <- c('CellType', 'Sample_ID', 'Group', 'Region', 'Donor')
    seu_lim <- AddMetaData(seu_lim, df)

    ## print
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_All_conditions_BA4.rds')


    ####################################
    ## BA9
    ####################################
    
    seu <- readRDS('/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_C9ALS_C9FTD_BA4_BA9_merged.rds')

    ## only select BA9
    str(seu@meta.data)
    unique(seu@meta.data$brain_region) #"ALS" "PN"

    ## subset to only include cell class of interest
    seu_region <- c("medial frontal cortex")
        
    ## select only cell type of interest
    xx <- unique(seu@meta.data$brain_region)
    xx <- xx[xx %in% c(seu_region)]
    Idents(seu) <- "brain_region"
    seu_lim=subset(seu,idents=xx)
    unique(seu_lim@meta.data$brain_region)
    dim(seu_lim)

    ## fix column names
    df <- seu_lim@meta.data
    df <- df %>% dplyr::select(Pineda_BA9_predictions, sample, diagnosis, brain_region, donor_id)
    colnames(df)
    colnames(df) <- c('CellType', 'Sample_ID', 'Group', 'Region', 'Donor')
    seu_lim <- AddMetaData(seu_lim, df)

    ## print
    saveRDS(seu_lim, '/home/fiorini9/scratch/machine_learning_ALS/base_objects/Li_All_conditions_BA9.rds')

##


### AFTER cell type annotation, need to do combat. 
### After combat transfer to Narval.

