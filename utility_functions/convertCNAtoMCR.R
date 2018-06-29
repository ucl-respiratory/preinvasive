# Convert a table of CNA ranges to minimum consisntent regions
convertCNAtoMCR <- function(c){
  library(GenomicRanges)
  range <- GRanges(
    seqnames = Rle(c$Chr),
    ranges = IRanges(start=c$Start, end=c$End)
  )
  allranges <- disjoin(range)
  
  # Convert back to data frame with CN per sample
  samples <- unique(c$SampleID)
  df <- data.frame(chr=seqnames(allranges), start=start(allranges), end=end(allranges))
  for(i in 1:length(samples)){
    df[,samples[i]] <- NA
  }
  
  for(i in 1:dim(c)[1]){
    # For each segment, find all matching subsegments
    pt <- c$SampleID[i]
    sel <- which(df$start >= c$Start[i] & df$end <= c$End[i])
    df[sel,pt] <- c$cn[i]
  }
  
  return(df)
}