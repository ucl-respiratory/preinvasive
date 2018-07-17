# Convert a table of CNA ranges to minimum consisntent regions
source('utility_functions/parallel.setup.R')
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
  
  # Fast method - for each region, look up MCRs in that region
  # Return only MCRs completely covered by the region
  df2 <- foreach(samp=samples, .combine=cbind) %dopar% {
    s <- matrix(nrow=dim(df)[1], ncol=1)
    c.pt <- c[which(c$SampleID == samp),]
    for(i in 1:dim(c.pt)[1]){
      s[which(
        df$chr == c.pt$Chr[i] &
        df$start >= c.pt$Start[i] &
        df$end <= c.pt$End[i]
      ),1] <- c.pt$cn[i]
    }
    return(data.frame(s))
  }
  rownames(df2) <- rownames(df)
  colnames(df2) <- samples
  
  # There are some gaps - assume these have no CN change so fill in with 2s:
  sel <- which(is.na(df2), arr.ind = T)
  df2[sel] <- 2
  
  return(cbind(df, df2))
  
}