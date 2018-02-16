# Script to check for required packages and install if required

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
          "GEOquery")

to.install <- pkgs[which(!(pkgs %in% rownames(installed.packages())))]

if(length(to.install) > 0){
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
}
