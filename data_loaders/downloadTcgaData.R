# Download data directly from the TCGA
# Data is exported with only required attribute of Progression (1/0) for cancer/control samples
library(GenomicDataCommons)

downloadTcgaData <- function(
  w_type = "HTSeq - Counts",
  d_type = "Gene Expression Quantification",
  pform = NULL,
  cache_dir="./data/gxn/tcga"
){
  
  print(paste("Downloading", d_type, "data from TCGA"))
  
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
  
  # Extract extra info
  cases$sample_type <- unlist(lapply(1:dim(cases)[1], function(i){
    cases$cases.samples[i][[1]]$sample_type
  }))
  
  # Download the files
  dir.create(cache_dir, recursive = T, showWarnings = F)
  files <- gdcdata(as.character(cases$file_id), progress=F, destination_dir = cache_dir)
  
  # Unzip
  for(file in files){
    if(length(grep(pattern = ".gz$", x=file)) == 0){ next }
    system(paste("gunzip", file))
  }
  files <- gsub(".gz$", "", files)
  
  # Read into R
  for(i in 1:length(files)){
    print(paste("Parsing file", i, "/", length(files)))
    s <- read.table(files[i], sep="\t", header=F)
    colnames(s) <- c("gene", as.character(files[i]))
    if(i == 1){
      alldata <- s
    }else{
      alldata <- merge(alldata, s, by="gene")
    }
  }
  alldata.raw <- alldata # Just in case we make mistakes going forward...
  
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
  
  if(d_type == "Gene Expression Quantification"){
    # Convert rownames to gene name
    library(biomaRt)
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
  
  }
  
  if(!all(colnames(alldata) == pheno$name)){
    message("ERROR: data and pheno data do not match")
  }
  
  return(list(alldata, pheno))
}
