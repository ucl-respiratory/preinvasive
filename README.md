## Preinvasive Study

This repository contains the code used in our publication "Deciphering the genomic, epigenomic and transcriptomic landscapes of pre-invasive lung cancer lesions". 

All R code used to analyse genomic, gene expression and methylation data is included.

## Microarray Data

This code downloads large volumes of data from NCBI GEO and the Cancer Genome Atlas (TCGA). It is best run on a cluster; downloads may take several hours, and combining methylation data from these samples has high memory requirements. Data will be automatically cached as .RData files - contact us to request access to processed RData files.

*** GEO Downloads will not work until after publication ***

Microarray data is stored in private GEO repositories which will be made public on publication. For reviewer access please contact the authors.

With reviewer access to relevant GEO datasets, files should be downloaded as follows:

```
https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE94611&format=file&file=GSE94611%5Fnon%5Fnormalized%2Etxt%2Egz downloaded to data/gxn/discovery/gxn.discovery.txt

https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108082&format=file downloaded to data/gxn/validation/gxn.validation.tar

https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108123&format=file downloaded to data/meth/geo/meth.geo.data.tar
```

All processing functions and downloading of TCGA data are handled by the scripts in data_loaders.

## Genomic Data

Raw genomic data is stored in the European Genome Archive (EGA), accession number EGAD00001003883. Due to the potentially identifiable nature of the data, this is not openly accessible. To recreate this analysis, the reader must go through the EGA data access procedure. 

This code repository does not cover variant calling; this is performed using established tools freely available at https://github.com/cancerit, as described in the Methods section of our paper. Code in this repository works from downstream VCF files, which the reader must generate from BAM files downloaded from EGA.

To run the analysis code in this repository you will need a data/wgs directory containing the following directories:

* ascat - raw ASCAT output as tsv files
* ascat_summary - summary copy number profiles
* caveman - output of Caveman SNV calling in VCF format
* pindel - output of Pindel insertion-deletion calling in VCF format
* brass - output of Brass rearrangement calling in VCF format
* pileups - multi-sample pileup files for each patient

Please contact the EGA directly with any queries regarding data access.

## TCGA Data

Gene expression and methylation analyses use open-access microarray data downloaded directly from TCGA. Code to download these is included in the data_loaders directory.

Mutation data from TCGA is not open access. It was downloaded under an agreement between TCGA and the Sanger Institute. Variant calling was performed on raw data using the same methods described in the paper. We are not able to share these data directly with this paper. Should the reader wish to repeat our comparisons with TCGA data, similar results can be obtained using open-access masked TCGA data (removing germline mutations), but please be aware that minor differences may be present.

To create the figure 2 circos plot comparing CIS to TCGA data, and extended data figure 1 assessing TCGA mutational signatures, the user will need to add two additional files in the resources/ directory derived from TCGA data:

* TCGA_prop_samples_with_subs_or_indels_in_gene.txt - tab-delimited text file with two columns representing 1) gene name and 2) proportion of samples in which that gene is affected by a mutation
* TCGA_trinucleotide_counts.txt - tab-delimited text file containing a 96 x n matrix with sample names for n samples as rows, and base changes with trinucleotide contexts as columns e.g. C.A.in.ACA, C.A.in.ACC, C.A.in.ACG ... The cell entries should be the number of occurrences of each base change/trinucleotide context in each sample.

## Dependencies

Analysis code is written in R. The relevant version numbers used in our analysis are as follows:

* R version 3.5.0
* Bioconductor version 3.7
* Package versions as defined in resources/package.versions.csv

This code depends on many R packages which are publically available. You should be able to install these automatically by running 

```r
Rscript install_dependencies.R
```
  
This installs the latest versions of all required packages. Should you experience any issues we recommend installing the specific versions detailed in resources/package.versions.csv.

The ChAMP package contains issues in plotting CNA profiles. As such, please install ChAMP from the authors' account using the following (as per the install_dependencies.R file):

```r
library(devtools)
install_github("ucl-respiratory/ChAMP")
```