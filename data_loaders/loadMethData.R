##########################################################################
# Load Methylation Data from GEO and from TCGA
#
# Available variables:
# mdata.d, mpheno.d, mdata.v, mpheno.v (where d and v denote discovery and validation)
# mdata, mpheno (merged d and v sets)
# tcga.mdata, tcga.mpheno
# tcga.mdata.all, tcga.mpheno.all (merged CIS and TCGA data)
#
# Please note this script is best run on a cluster computer
# TCGA data is large (~60Gb) and parsing these data is memory-intensive.
##########################################################################

cache_file <- paste(data_cache, "mdata.RData", sep="")
cache_file_cna <- paste(data_cache, "cdata.RData", sep="")

if(file.exists(cache_file) & file.exists(cache_file_cna)){
  load(cache_file)
  load(cache_file_cna)
}else{
  
  cache_dir <- paste(data_cache, "meth/geo", sep="")
  geo.file <- paste(cache_dir, "/meth.geo.data.tar", sep="")
  
  if(!file.exists(geo.file)){
    # Download IDAT files directly from GEO
    url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108123&format=file"
    dir.create(cache_dir, recursive = T, showWarnings = F)
    x <- download.file(url, destfile = geo.file)
  }
  
  # Unzip IDAT tar file
  pwd <- getwd()
  setwd(cache_dir)
  cmd <- paste("tar -xvf ", pwd, "/", geo.file,sep = "")
  system(cmd)
  setwd(pwd)
  
  # Unzip individual idat files
  files <- list.files(cache_dir, pattern=".gz$", full.names = T)
  for(file in files){
    system(paste("gunzip", file))
  }
  # Remove CSV files as these confuse the ChAMP loader
  csv.files <- list.files(cache_dir, pattern=".csv$", full.names = T)
  file.remove(csv.files)
  
  # Copy sample sheet (stored in resources directory)
  file.copy(
    "./resources/methylation_sample_sheet.csv",
    cache_dir
  )
  
  # Run ChAMP - filter options explicitly defined here
  # QC steps are excluded as only samples passing QC are included in our uploaded data set
  myLoad <- champ.load(directory = cache_dir,
                       method="ChAMP",
                       methValue="B",
                       autoimpute=TRUE,
                       filterDetP=TRUE,
                       ProbeCutoff=0,
                       SampleCutoff=0.1,
                       detPcut=0.01,
                       filterBeads=TRUE,
                       beadCutoff=0.05,
                       filterNoCG=FALSE, 
                       filterSNPs=FALSE, 
                       population=NULL,
                       filterMultiHit=FALSE,
                       filterXY=TRUE,
                       force=FALSE,
                       arraytype="450K")
  mpheno   <- myLoad$pd
  myNorm   <- champ.norm(beta=myLoad$beta, rgSet=myLoad$rgSet, mset=myLoad$mset, cores=1)
  myCombat <- champ.runCombat(beta=myNorm, pd=mpheno)
  mdata    <- data.frame(myCombat)
  
  sel.d    <- which(mpheno$Cohort == "D")
  mdata.d  <- mdata[,sel.d]
  mpheno.d <- mpheno[sel.d,]
  sel.v    <- which(mpheno$Cohort == "V")
  mdata.v  <- mdata[,sel.v]
  mpheno.v <- mpheno[sel.v,]
  
  ##########################################################################
  # Load TCGA data
  ##########################################################################
  source('data_loaders/downloadTcgaData.R')
  x <- downloadTcgaData(
    w_type = "Liftover", 
    d_type = "Methylation Beta Value",
    pform="Illumina Human Methylation 450",
    cache_dir="./data/meth/tcga"
  )
  tcga.mdata  <- x[[1]]
  tcga.mpheno <- x[[2]]
  
  # Add Sample Group to TCGA data
  tcga.mpheno$Sample_Group <- NA
  tcga.mpheno$Sample_Group[which(tcga.mpheno$progression == 0)] <- "TCGA Control"
  tcga.mpheno$Sample_Group[which(tcga.mpheno$progression == 1)] <- "TCGA SqCC"
  
  # Merge CIS and TCGA data together
  x <- runComBat(mdata, tcga.mdata)
  tcga.mdata.all  <- cbind(x[[1]], x[[2]])
  tcga.mpheno.all <- rbind(
    data.frame(
      name=mpheno$Sample_Name,
      Sample_Group=mpheno$Sample_Group,
      progression=mpheno$progression,
      source="Surveillance",
      dose=mpheno$progression + 1
    ),
    tcga.mpheno
  )
  
  # Cache results
  save(
    mdata, mpheno,
    mdata.d, mpheno.d, mdata.v, mpheno.v,
    tcga.mdata, tcga.mpheno, 
    tcga.mdata.all, tcga.mpheno.all,
    file=cache_file
  )
  
  ##########################################################################
  # Generate methylation-derived CNA data
  ##########################################################################
  # Generate overall CNA profiles, including plot generation
  champ_cna_dir <- paste(data_cache, "meth_cna/CHAMP_CNA", sep="")
  dir.create(champ_cna_dir, recursive = T, showWarnings = F)
  cna <- champ.CNA(
    controlGroup="Control",
    resultsDir = champ_cna_dir
  )
  
  # Now we need to extract probe-level CNA data
  # To do this, we modify the champ.CNA function:
  ints <- myLoad$intensity
  controlGroup="Control"
  data(probe.features)
  pd <- myLoad$pd
  
  library(preprocessCore)
  names <- colnames(ints)
  intsqn <- normalize.quantiles(as.matrix(ints))
  colnames(intsqn) <- names
  intsqnlog <- log2(intsqn)
  
  controlSamples = pd[which(pd$Sample_Group == controlGroup),]
  caseSamples = pd[which(pd$Sample_Group != controlGroup),]
  case.intsqnlog <- intsqnlog[, which(colnames(intsqnlog) %in%
                                        caseSamples$Sample_Name)]
  control.intsqnlog <- intsqnlog[, which(colnames(intsqnlog) %in%
                                           controlSamples$Sample_Name)]
  control.intsqnlog <- rowMeans(control.intsqnlog)
  intsqnlogratio <- case.intsqnlog
  for (i in 1:ncol(case.intsqnlog)) {
    intsqnlogratio[, i] <- case.intsqnlog[, i] - control.intsqnlog
  }
  
  ints <- data.frame(ints, probe.features$MAPINFO[match(rownames(ints),
                                                        rownames(probe.features))])
  names(ints)[length(ints)] <- "MAPINFO"
  ints <- data.frame(ints, probe.features$CHR[match(rownames(ints),
                                                    rownames(probe.features))])
  names(ints)[length(ints)] <- "CHR"
  levels(ints$CHR)[levels(ints$CHR) == "X"] = "23"
  levels(ints$CHR)[levels(ints$CHR) == "Y"] = "24"
  CHR <- as.numeric(levels(ints$CHR))[ints$CHR]
  ints$MAPINFO <- as.numeric(ints$MAPINFO)
  MAPINFO = probe.features$MAPINFO[match(rownames(ints), rownames(probe.features))]
  MAPINFO <- as.numeric(MAPINFO)
  
  message("Saving Copy Number information for each Sample")
  library(DNAcopy)
  for (i in 1:ncol(case.intsqnlog)) {
    CNA.object <- CNA(cbind(intsqnlogratio[, i]), CHR,
                      MAPINFO, data.type = "logratio", sampleid = paste(colnames(case.intsqnlog)[i],
                                                                        "qn"))
    
    ####
    # This object has CNA information by location
    # Store for all samples in a large matrix
    smoothed.CNA.object <- smooth.CNA(CNA.object)
    
    # Combine these objects into a massive matrix
    if(i == 1){
      cnas <- smoothed.CNA.object
    }else{
      if(all(rownames(cnas) == rownames(smoothed.CNA.object))){
        s <- colnames(smoothed.CNA.object)[dim(smoothed.CNA.object)[2]]
        cnas[s] <- smoothed.CNA.object[s]
      }else{
        print(paste("FAILED SAMPLE", i))
      }
    }
  }
  
  # We now have a large matrix of CNA values by probe in cnas
  # Rownames are the indices of the correct row in the intensity dataset
  rownames(cnas) <- rownames(myLoad$intensity)[as.numeric(rownames(cnas))]
  
  cnas.probes <- cnas
  cnas.probes$chrom <- NULL
  cnas.probes$maploc <- NULL
  
  cnas.pheno <- pd
  cnas.pheno$Sample_Name <- make.names(paste(cnas.pheno$Sample_Name, ".qn", sep=""))
  # Remove Controls
  cnas.pheno <- cnas.pheno[-which(cnas.pheno$Sample_Group == "Control"),]
  
  if(!all(cnas.pheno$Sample_Name == colnames(cnas.probes))){
    message("ERROR: cna.probes colnames and pheno don't match")
  }
  
  # Add cytoband info to the dict reference
  # First get the location data from Ensembl:
  ensembl_url <- "http://rest.ensembl.org/info/assembly/homo_sapiens?content-type=application/json&bands=1"
  library(httr)
  req <- GET(url=ensembl_url)
  # JSON data stored in content(req)
  # Turn this into a data object with each band represented by chromosome, name, start, end
  bands <- c()
  chroms <- c()
  ids <- c()
  starts <- c()
  ends <- c()
  
  for(reg in content(req)[["top_level_region"]]){
    for(band in reg$bands){
      chroms <- c(chroms, band$seq_region_name)
      ids <- c(ids, band$id)
      starts <- c(starts, band$start)
      ends <- c(ends, band$end)
      bands <- c(bands, c(band$seq_region_name, band$id, band$start, band$end))
    }
  }
  
  bands <- data.frame(chrom=chroms, id=ids, start=starts, end=ends)
  rownames(bands) <- paste(bands$chrom, bands$id, sep="")
  
  cnas$band <- NA
  # Now go through all probes and label them with a band
  for(i in 1:dim(cnas)[1]){
    sel <- which(bands$chrom == cnas$chrom[i] & bands$start < cnas$maploc[i] & bands$end > cnas$maploc[i])
    if(length(sel) == 1){
      cnas$band[i] <- rownames(bands)[sel]
    }else{
      print(paste("Unable to match row", i))
    }
  }
  # Aggregate over bands
  cnas.band <- cnas[,3:(dim(cnas)[2])]
  cnas.band <- aggregate(cnas.band, by=list(cnas.band$band), FUN=mean)
  rownames(cnas.band) <- cnas.band$Group.1
  cnas.band$Group.1 <- NULL
  cnas.band$band <- NULL
  cnas.band <- cnas.band[,cnas.pheno$Sample_Name]
  
  # Included data set from Dutch group (van Boerdonk et al)
  source('./data_loaders/loadFromGEO.R')
  x <- loadFromGEO("GSE45287", destdir = paste(data_cache, "cna/dutch", sep=""))
  
  save(cnas.band, cnas.pheno, file=cache_file_cna)
}

