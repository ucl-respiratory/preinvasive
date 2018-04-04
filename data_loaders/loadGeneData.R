##########################################################################
# Load Gene Expression Data from GEO
# Loads CIS expression data from GEO, and TCGA data using Genomic Data Commons. 
# Data is cached in RData format.
#
# Available variables:
# gdata.d, gpheno.d, gdata.v, gpheno.v (where d and v denote discovery and validation)
# gdata, gpheno (merged d and v sets)
##########################################################################

library(affycoretools)
library(oligo)
source('utility_functions/runComBat.R')
source('data_loaders/downloadTcgaData.R')
library(limma)

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
  
  cache_dir <- paste(data_cache, "gxn/geo", sep="")
  geo.file.d <- paste(cache_dir, "/discovery/gxn.discovery.txt", sep="")
  geo.file.v <- paste(cache_dir, "/validation/gxn.validation.tar", sep="")
  
  # Download Illumina GXN data (Discovery set)
  if(!file.exists(geo.file.d)){
    # Download non-normalised txt file directly from GEO
    url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE94611&format=file&file=GSE94611%5Fnon%5Fnormalized%2Etxt%2Egz"
    dir.create(cache_dir, recursive = T, showWarnings = F)
    x <- download.file(url, destfile = geo.file.d)
  }
  gdata.d <- read.table(geo.file.d, sep="\t", header=T, quote="", fill=T)
  rownames(gdata.d) <- gdata.d$ID_REF
  cols <- grep("Log", colnames(gdata.d))
  gdata.d <- gdata.d[,cols]
  
  # Add clinical data to this. This matches data available on GEO, and is included in RData format for convenience.
  load("resources/gxnPhenoDiscovery.RData")
  colnames(gdata.d) <- gpheno.d$name
  
  # Normalise using quantile normalisation
  gdata.d <- normalizeQuantiles(gdata.d)
  
  # Download Affymetrix GXN data (Discovery set)
  if(!file.exists(geo.file.v)){
    # Download IDAT files directly from GEO
    url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108082&format=file"
    dir.create(cache_dir, recursive = T, showWarnings = F)
    x <- download.file(url, destfile = geo.file.v)
  }
  
  # Untar the CEL files
  wd <- getwd()
  setwd(paste(cache_dir, "/validation", sep=""))
  cmd <- paste("tar -xvzf", paste(wd, geo.file.v, sep="/"))
  system(cmd)
  setwd(wd)
  
  # Now unzip them
  cmd <- paste("gunzip ", cache_dir, "/validation/*.gz", sep="")
  system(cmd)
  
  # Load the CEL files 
  cel.files <- list.files(paste(cache_dir, "validation", sep="/"), pattern = "CEL", full.names = T)
  gdata.v <- read.celfiles(cel.files)
  gdata.v <- rma(gdata.v)
  
  # Convert to matrix format with gene names as rownames and samples as columns
  # Throw away probes not mapped to a gene symbol
  gdata.v <- annotateEset(gdata.v, pd.clariom.d.human)
  sel <- which(!is.na(fData(gdata.v)$SYMBOL))
  gene.names <- fData(gdata.v)$SYMBOL[sel]
  gdata.v <- exprs(gdata.v)[sel,]
  rownames(gdata.v) <- gene.names
  
  # Convert column names to match existing data
  cnames <- list("GSM2889208"= "Sample_34 Affy",
                 "GSM2889209"=	"Sample_35 Affy",
                 "GSM2889210"=	"Sample_36 Affy",
                 "GSM2889211"=	"Sample_37 Affy",
                 "GSM2889212"=	"Sample_38 Affy",
                 "GSM2889213"=	"Sample_39 Affy",
                 "GSM2889214"=	"Sample_40 Affy",
                 "GSM2889215"=	"Sample_41 Affy",
                 "GSM2889216"=	"Sample_42 Affy",
                 "GSM2889217"=	"Sample_46 Affy",
                 "GSM2889218"=	"Sample_51 Affy",
                 "GSM2889219"=	"Sample_52 Affy",
                 "GSM2889220"=	"Sample_57 Affy",
                 "GSM2889221"=	"Sample_58 Affy",
                 "GSM2889222"=	"Sample_59 Affy",
                 "GSM2889223"=	"Sample_62 Affy",
                 "GSM2889224"=	"Sample_63 Affy",
                 "GSM2889225"=	"Sample_64 Affy"
  )
  colnames(gdata.v) <- unlist(lapply(colnames(gdata.v), function(x){
    unlist(strsplit(x, "_"))[[1]]
  }))
  colnames(gdata.v) <- as.character(cnames[colnames(gdata.v)])
  
  # Link to pheno data - again stored in resources directory for convenience
  load("resources/gxnPhenoValidation.RData")
  gdata.v <- gdata.v[,gpheno.v$name]
  
  
  # Lastly create a combined data set, gdata/gpheno, using ComBat to correct for batch effects
  # Note that we adjust the validation set to be comparable with the discovery set. The discovery set is not affected by this process.
  x <- runComBat(gdata.d, gdata.v)
  gdata <- cbind(x[[1]], x[[2]])
  gpheno <- rbind(gpheno.d, gpheno.v)
  
  ################################################################################
  # Add TCGA data
  ################################################################################
  
  # Load TCGA data
  x <- downloadTcgaData() # Default options are for GXN RNAseq data.
  tcga.gdata  <- x[[1]]
  tcga.gpheno <- x[[2]]
  
  # Voom-transform RNAseq data for comparison with microarray data
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

