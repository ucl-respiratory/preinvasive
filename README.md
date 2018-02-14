## Preinvasive Study

This repository contains the code used in our publication "Deciphering the genomic, epigenomic and transcriptomic landscapes of pre-invasive lung cancer lesions". 

All R code used to analyse gene expression and microarray data is included.

Genomic analysis uses freely available, well established tools available at https://github.com/cancerit, as described in the Methods section.

## Data Downloads

This code downloads large volumes of data from NCBI GEO and the Cancer Genome Atlas (TCGA). It is best run on a cluster; downloads may take several hours. Data will be automatically cached as .RData files - contact us to request access to processed RData files.

*** GEO Downloads will not work currently ***

Microarray data is stored in private GEO repositories which will be made public on publication. For reviewer access please contact the authors.


## Dependencies

This code depends on many R packages which are publically available.

The ChAMP package contains issues in plotting CNA profiles. As such, please install ChAMP from the authors' account using the following:

library(devtools)
install_github("adamp83/ChAMP")