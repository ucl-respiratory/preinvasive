##########################################################################
# Figure 2E: genome-wide copy number plot
##########################################################################
plot.genomic.cna.genomewide <- function(filename="results/figures/cna.genomewide.png"){
  
  # GenVisR copy number visualisation
  
  ascat.dir <- "./data/wgs/ascat_summary/"
  library(stringr)
  library(GenVisR)
  
  # Function to read in a file from a sample name e.g. PD21884a
  a <- function(name){
    # Find the file
    x <- list.files(ascat.dir, pattern=name, full.names = T)
    # read data and set column names
    data <- read.csv(x, header=FALSE)
    data <- data[,c(2,3,4,7)]
    colnames(data) <- c("chromosome", "start", "end", "segmean")
    
    # Correct for ploidy
    ploidy <- wgs.pheno$ploidy[which(wgs.pheno$name == name)]
    data$segmean <- 2*data$segmean / ploidy
    
    # get the sample name from the file path
    sampleName <- str_extract(x, "PD[0-9]+[a-z]")
    data$sample <- sampleName
    
    # return the data
    return(data)
  }
  
  o <- order(wgs.pheno$progression)
  cnData <- lapply(wgs.pheno$name[o], a)
  cnData <- do.call("rbind", cnData)
  
  # construct genomic boundaries from cytoGeno
  genomeBoundaries <- aggregate(chromEnd ~ chrom, data=cytoGeno[cytoGeno$genome=="hg19",], max)
  genomeBoundaries$chromStart <- 0
  colnames(genomeBoundaries) <- c("chromosome", "end", "start")
  
  # Ignore sex chromosomes
  genomeBoundaries <- genomeBoundaries[which(genomeBoundaries$chromosome %in% paste0("chr", 1:22)),]
  cnData <- cnData[which(cnData$chromosome %in% 1:22),]
  
  # create the plot
  png(filename, width=2560, height=1440)
  cnSpec(cnData, y=genomeBoundaries)
  dev.off()
  
}
