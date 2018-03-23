# Script to check for required packages and install if required

message("Checking for installed dependencies:")

pkgs <- c("gdata",
          "devtools",
          "ChAMP",
          "ggplot2",
          "ggsignif",
          "pheatmap",
          "pROC",
          "RColorBrewer",
          "limma",
          "stringr",
          "WriteXLS",
          "pamr",
          "gage",
          "gageData",
          "org.Hs.eg.db",
          "annotate",
          "knitr",
          "lme4",
          "sva",
          "preprocessCore",
          "DNAcopy",
          "httr",
          "GenomicDataCommons",
          "biomaRt",
          "GEOquery",
          "affy", 
          "oligo",
          "pd.clariom.d.human")

# By default the latest version of these packages will be installed. 
# Versions used to conduct the analysis can be found in the resources folder.
# The code used to create this dependency file from our development environment is as follows:
#
# version.file <- "resources/package.versions.csv"
# version.data <- installed.packages()[pkgs,]
# write.csv(version.data, file = version.file)


# Find missing packages:
to.install <- pkgs[which(!(pkgs %in% rownames(installed.packages())))]

if(length(to.install) > 0){
  message(paste("Packages to be installed:", paste(to.install, collapse=", ")))
  source("https://bioconductor.org/biocLite.R")
  for(pkg in to.install){
    if(pkg == "ChAMP"){
      library(devtools)
      install_github("adamp83/ChAMP")
    }else{
      biocLite(pkg, suppressUpdates=T,
               suppressAutoUpdate=T, ask=F)
    }
  }
}else{
  message("All required packages found.")
}
