# Script to check for required packages and install if required

message("Checking for installed dependencies:")

pkgs <- c("gdata",
          "devtools",
          "GenomeInfoDbData",
          "ChAMPdata",
          "Illumina450ProbeVariants.db",
          "IlluminaHumanMethylationEPICmanifest",
          "IlluminaHumanMethylation450kmanifest", 
          "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
          "IlluminaHumanMethylation450kanno.ilmn12.hg19",
          "DMRcate",
          "geneLenDataBase",
          "ChAMP",
          "ggplot2",
          "ggsignif",
          "pheatmap",
          "pROC",
          "PRROC",
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
          "pd.clariom.d.human",
          "Homo.sapiens",
          "BSgenome.Hsapiens.UCSC.hg19",
          "TxDb.Hsapiens.UCSC.hg19.knownGene",
          "biovizBase",
          "maftools",
          "dndscv",
          "affycoretools",
          "sciClone",
          "foreach",
          "doParallel", 
          "GenVisR",
          "MutationalPatterns")

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
    # The ChAMP package has been slightly modified to fix a bug in plotting CNA profiles
    # It is therefore downloaded from the authors' github page, which is publically accessible
    if(pkg == "ChAMP"){
      library(devtools)
      install_github("ucl-respiratory/ChAMP")
    }else{
      if(pkg == "dndscv"){
        library(devtools)
        install_github("im3sanger/dndscv")
      }else{
        if(pkg == "sciClone"){
          library(devtools)
          install_github("genome/bmm")
          install_github("genome/sciClone")
        }else{
          biocLite(pkg, suppressUpdates=T,
                   suppressAutoUpdate=T, ask=F)
        }
      }
      
    }
  }
}else{
  message("All required packages found.")
}
