#' loadFromGEO
#'
#' Load a data set from GEO and process it into data/pheno data frames

library(GEOquery)
loadFromGEO <- function(geo.num, destdir="./data/geo", gene.col="Gene Symbol", bandAdjust=F){
  dir.create(destdir, recursive = T, showWarnings = F)
  
  gse <- getGEO(geo.num, GSEMatrix = TRUE, destdir=destdir)
  
  geo.gdata <- data.frame(exprs(gse[[1]]))
  geo.pheno <- pData(gse[[1]])
  
  # Get gene names from the GSE annotation data
  annot <- pData(featureData(gse[[1]]))
  
  # Specific to CNA downloading - do band adjustment on cytobands 
  # This is only used for downloading van Boerdonk et al data
  if(bandAdjust){
    chr <- unlist(lapply(annot$CHROMOSOMAL_LOCATION, function(x){
      y <- unlist(strsplit(x, ":"))[1]
      y <- gsub("chr", "", y)
      return(y)
    }))
    annot$CYTOBAND <- gsub("[hs|]", "", annot$CYTOBAND)
    annot$CYTOBAND <- paste(chr, annot$CYTOBAND, sep="")
  }
  
  genes <- as.character(annot[,gene.col])
  # Discard multiple genes per probe
  genes <- unlist(lapply(genes, function(x){
    unlist(strsplit(x, " /// "))[1]
  }))
  if(dim(geo.gdata)[1] > 0){
    # Aggregate - use mean values for multiple probes per gene
    geo.gdata <- aggregate(geo.gdata, by=list(genes), FUN=mean)
    rownames(geo.gdata) <- geo.gdata$Group.1
    geo.gdata$Group.1 <- NULL
  }
  
  
  return(list(geo.gdata, geo.pheno))
}
