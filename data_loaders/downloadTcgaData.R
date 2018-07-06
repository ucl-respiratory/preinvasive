# Download data directly from the TCGA
# Data is exported with only required attribute of Progression (1/0) for cancer/control samples
library(GenomicDataCommons)
library(biomaRt)
library(httr)
library(maftools)

downloadTcgaData <- function(
  w_type = "HTSeq - Counts",
  d_type = "Gene Expression Quantification",
  pform = NULL#,
  #cache_dir="./data/gxn/tcga"
){
  cache_dir <- './data/gdc'
  gdc_set_cache(cache_dir)
  
  print(paste("Downloading", d_type, "data from TCGA"))
  source("utility_functions/parallel.setup.R")
  
  # Download tumour and controls matching our file specifications
  fileResults <- files() %>%
    GenomicDataCommons::filter(~ "data_type" == d_type & 
                                 "cases.project.project_id" == "TCGA-LUSC" &
                                 ("cases.samples.sample_type" == "Primary Tumor" | "cases.samples.sample_type" == "Solid Tissue Normal") & 
                                 "analysis.workflow_type" == w_type
                               ) %>%
    GenomicDataCommons::expand(c("cases.samples")) %>%
    results(size=50000)
  cases <- as.data.frame(fileResults)
  
  if(!is.null(pform) & "platform" %in% colnames(cases)){
    sel <- which(cases$platform == pform)
    cases <- cases[sel,]
  }
  
  print(paste("Found", dim(cases)[1], "samples, downloading to", cache_dir))
  dir.create(cache_dir, recursive = T, showWarnings = F)
  
  # Somatic mutations are a different file format so handle separately.
  if(d_type == "Masked Somatic Mutation"){
    # Download the files
    files <- gdcdata(unique(as.character(cases$file_id)), progress=F)
    
    # We expect only one file here - gzipped MAF
    if(length(files) != 1){
      message("WARNING: Multiple SNV MAF files found, using the first only")
    }
    file <- files[[1]]
    system(paste("gunzip", file))
    file <- gsub(".gz", "", file)
    
    # Read MAF file
    #s <- read.table(file, sep="\t", header=T)
    s <- read.maf(file)
    
    # Return the raw data file
    return(s@data)
  }
  
  # Ongoing code is for any type other than somatic mutations
  cases$sample_type <- unlist(lapply(1:dim(cases)[1], function(i){
    # cases$cases.samples[i][[1]]$sample_type
    fileResults$cases[[as.character(cases$file_id[i])]]$samples[[1]]$sample_type
  }))
  
  # Download the files - ignore cached files
  cached_files <- c()
  to.download <- c()
  for(i in 1:dim(cases)[1]){
    f <- list.files(cache_dir, pattern=as.character(cases$file_name[i]), full.names = T, recursive = T)
    if(length(f) > 0){
      cached_files <- c(cached_files, f)
    }else{
      to.download <- c(to.download, as.character(cases$file_name[i]))
    }
  }
  if(length(to.download) > 0){
    files <- gdcdata(as.character(cases$file_id), progress=F)
    # files <- foreach(f=cases$file_id) %dopar% {
    #   library(GenomicDataCommons)
    #   return( gdcdata(as.character(f), progress=F) )
    # }
    # files <- unlist(files)
  }
  #files <- list.files(cache_dir, pattern=paste(as.character(cases$file_name), collapse="|"), recursive = T, full.names = T)
  files <- unlist(lapply(1:dim(cases)[1], function(i){
    paste(cache_dir, as.character(cases$file_id[i]), as.character(cases$file_name[i]), sep="/")
  }))
  
  # Unzip - don't delete original
  for(file in files){
    if(length(grep(pattern = ".gz$", x=file)) == 0){ next }
    system(paste("gunzip -kf", file))
  }
  files <- gsub(".gz$", "", files)
  
  # Read into R
  # For Copy Number data map to bands. Otherwise merge by gene.
  if(d_type == "Copy Number Segment"){
    
    # Parse the files and add a pheno data frame
    x <- parseCnaFiles(files, is.TCGA = T)
    pheno <- data.frame(
      name=as.character(cases$file_name),
      progression=as.numeric(list("Solid Tissue Normal"=0,"Primary Tumor"=1)[cases$sample_type]),
      source="TCGA"
    )
    x[[4]] <- pheno
    return(x)
    
  }else{
    for(i in 1:length(files)){
      print(paste("Parsing file", i, "/", length(files)))
      if(d_type == "Methylation Beta Value"){
        s <- read.table(files[i], sep="\t", header=T)
        s <- s[,1:2]
      }else{
        s <- read.table(files[i], sep="\t", header=F)
      }
      colnames(s) <- c("gene", as.character(files[i]))
      if(i == 1){
        alldata <- s
      }else{
        alldata <- merge(alldata, s, by="gene")
      }
    }
    alldata.raw <- alldata # Just in case we make mistakes going forward...
  }
  
  
  pheno <- data.frame(
    name=files,
    progression=NA,
    source=rep("TCGA", length(files)),
    dose=NA
  )
  pheno$progression[which(cases$sample_type == "Primary Tumor")] <- 1
  pheno$progression[which(cases$sample_type == "Solid Tissue Normal")] <- 0
  pheno$dose[which(cases$sample_type == "Primary Tumor")] <- 3
  pheno$dose[which(cases$sample_type == "Solid Tissue Normal")] <- 0
  
  # For gene expression data, convert to gene symbols
  if(d_type == "Gene Expression Quantification"){
    # Convert rownames to gene name
    ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
    enids <- alldata.raw$gene
    sel <- grep("^ENS", enids)
    enids <- enids[sel]
    alldata <- alldata.raw[sel,]
    enids <- unlist(lapply(as.character(enids), function(x){unlist(strsplit(x, "[.]"))[[1]]}))
    alldata$gene <- enids
    bm <- getBM(
      attributes=c("ensembl_gene_id", "hgnc_symbol"),
      filters=c("ensembl_gene_id"),
      values=list(enids),
      mart=ensembl
    )
    
    alldata <- merge(alldata, bm, by.x="gene", by.y="ensembl_gene_id")
    alldata$gene <- NULL
    genes <- alldata$hgnc_symbol
    alldata$hgnc_symbol <- NULL
    alldata <- aggregate(alldata, by=list(genes), FUN=sum)
    sel <- which(alldata$Group.1 == "")
    if(length(sel) > 0){alldata <- alldata[-sel,]}
    rownames(alldata) <- alldata$Group.1
    alldata$Group.1 <- NULL
  
  }else{
    # Remove the gene column
    rownames(alldata) <- alldata$gene
    alldata$gene <- NULL
  }
  
  if(!all(colnames(alldata) == pheno$name)){
    message("ERROR: data and pheno data do not match")
  }
  
  return(list(alldata, pheno))
}
