## Preinvasive Study

This repository contains the code used in our publication "Deciphering the genomic, epigenomic and transcriptomic landscapes of pre-invasive lung cancer lesions". 

All R code used to analyse gene expression and microarray data is included.

Genomic analysis uses freely available, well established tools available at https://github.com/cancerit, as described in the Methods section.

## Data Downloads

This code downloads large volumes of data from NCBI GEO and the Cancer Genome Atlas (TCGA). It is best run on a cluster; downloads may take several hours. Data will be automatically cached as .RData files - contact us to request access to processed RData files.

*** GEO Downloads will not work currently ***

Microarray data is stored in private GEO repositories which will be made public on publication. For reviewer access please contact the authors.

With reviewer access to relevant GEO datasets, files should be downloaded as follows:

```
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE94611&format=file&file=GSE94611%5Fnon%5Fnormalized%2Etxt%2Egz downloaded to data/gxn/discovery/gxn.discovery.txt

https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108082&format=file downloaded to data/gxn/validation/gxn.validation.tar

https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108123&format=file downloaded to data/meth/geo/meth.geo.data.tar
```

All processing functions and downloading of TCGA data are handled by the scripts in data_loaders.

## Dependencies

Analysis code is written in R. The relevant version numbers used in our analysis are as follows:

* R version 3.4.1
* Bioconductor version 3.6
* Package versions as defined in resources/package.versions.csv

This code depends on many R packages which are publically available. You should be able to install these automatically by running 

```r
Rscript install_dependencies.R
```
  
This installs the latest versions of all required packages. Should you experience any issues we recommend installing the specific versions detailed in resources/package.versions.csv.

The ChAMP package contains issues in plotting CNA profiles. As such, please install ChAMP from the authors' account using the following:

```r
library(devtools)
install_github("adamp83/ChAMP")
```