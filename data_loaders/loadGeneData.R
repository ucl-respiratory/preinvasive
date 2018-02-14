##########################################################################
# Load Gene Expression Data from GEO
#
# Available variables:
# gdata.d, gpheno.d, gdata.v, gpheno.v (where d and v denote discovery and validation)
# gdata, gpheno (merged d and v sets)
##########################################################################

cache_file <- paste(data_cache, "gdata.RData", sep="")

if(file.exists(cache_file)){
  load(cache_file)
}else{
  
  source('data_loaders/loadFromGEO.R')
  x <- loadFromGEO(geo.gxn.d, paste(data_cache, "gxn/discovery", sep=""))
  gdata.d  <- x[[1]]
  gpheno.d <- x[[2]]
  
  x <- loadFromGEO(geo.gxn.v, paste(data_cache, "gxn/validation", sep=""))
  gdata.v  <- x[[1]]
  gpheno.v <- x[[2]]
  
  # Combine discovery and validation sets using ComBat
  source('utility_functions/runComBat.R')
  x <- runComBat(gdata.d, gdata.v)
  gdata  <- cbind(x[[1]], x[[2]])
  gpheno <- rbind(gpheno.d, gpheno.v)
  
  
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

