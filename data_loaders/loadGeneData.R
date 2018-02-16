##########################################################################
# Load Gene Expression Data from GEO
#
# Available variables:
# gdata.d, gpheno.d, gdata.v, gpheno.v (where d and v denote discovery and validation)
# gdata, gpheno (merged d and v sets)
##########################################################################

if(!exists("data_cache")){
  data_cache <- "./data/"
}

cache_file <- paste(data_cache, "gdata.RData", sep="")

if(file.exists(cache_file)){
  message("Loading cached gene expression data")
  load(cache_file)
}else{
  message("Downloading gene expression data. This may take some time.")
  
  ################################################################################
  # Load CIS gene expression data
  #
  # Data is available from GEO as described in the manuscript.
  # Included in this repository is processed data in RData format.
  # Illumina data (discovery set) is quantile normalised.
  # Affymetrix data (validation set) is normalised using the rma method of the affy package.
  ################################################################################
  
  # This provides the following variables:
  # gdata, gpheno, gdata.d, gpheno.d, gdata.v, gpheno.v 
  load('./resources/gxnData.RData')
  
  ################################################################################
  # Add TCGA data
  ################################################################################
  
  # Load TCGA data
  source('data_loaders/downloadTcgaData.R')
  x <- downloadTcgaData() # Default options are for GXN RNAseq data.
  tcga.gdata  <- x[[1]]
  tcga.gpheno <- x[[2]]
  
  # Voom-transform RNAseq data for comparison with microarray data
  library(limma)
  tcga.gdata <- data.frame(voom(tcga.gdata))
  
  # Merge CIS and TCGA data together
  x <- runComBat(gdata, tcga.gdata)
  tcga.gdata.all  <- cbind(x[[1]], x[[2]])
  tcga.gpheno.all <- rbind(
    data.frame(
      name=gpheno$name,
      progression=gpheno$progression,
      source="Surveillance",
      dose=gpheno$progression + 1
    ),
    tcga.gpheno
  )
  
  # Cache results
  save(
    gdata, gpheno,
    gdata.d, gpheno.d, gdata.v, gpheno.v,
    tcga.gdata, tcga.gpheno, 
    tcga.gdata.all, tcga.gpheno.all,
    file=cache_file
  )
}

