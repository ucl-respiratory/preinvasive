# Parse a group of CNA segment files
# Used for ASCAT summary output, and TCGA CNA files
# Should output:
#   A data frame of minimum consistent CNAs (cnas.segmented)
#   A matrix of copy number by band
#   A matrix of copy number by gene
#
# Should work on:
#   tcga data: filenames <- list.files('data/cna/tcga', full.names=T, pattern=".seg.txt$")
#   CIS data:  filenames <- list.files("data/wgs/ascat_summary/", full.names = T, pattern="summary.csv$")
library(GenomicRanges)
library(Homo.sapiens)
library(BSgenome.Hsapiens.UCSC.hg19)
parseCnaFiles <- function(filenames=c(), is.TCGA=F, my.names=NULL, cnadata=NULL){
  if(is.null(my.names) & is.null(cnadata)){
    my.names <- c()
    cnadata <- list()
    for(i in 1:length(filenames)){
      
      if(is.TCGA){
        c <- read.table(filenames[i], sep="\t", stringsAsFactors = F, header = T)
        pt <- c$Sample[1]
        my.names <- c(my.names, pt)
        c$Sample <- NULL
        c$Num_Probes <- NULL
        colnames(c) <- c('chr', 'start', 'end', 'cn')
      }else{
        c <- read.csv(filenames[i], stringsAsFactors = F, header = F, row.names = 1)
        ptindex <- regexpr("PD[0-9]+", filenames[i])
        pt <- substr(filenames[i], ptindex, ptindex + 7)
        my.names <- c(my.names, pt)
        colnames(c) <- c('chr', 'start', 'end', 'pl', 'pl.min', 'cn', 'cn.min')
        c <- c[,c("chr", "start", "end", "cn")]
      }
      
      
      
      
      # Remove X and Y chromosomes
      sel <- which(c$chr %in% c("X", "Y"))
      if(length(sel) > 0){c <- c[-sel,]}
      
      cnadata[[i]] <- c
      
    }
  }
  for(i in 1:length(cnadata)){
    c <- cnadata[[i]]
    range <- GRanges(
      seqnames = Rle(c$chr),
      ranges = IRanges(start=c$start, end=c$end)
    )
    
    if(i == 1){
      scnas <- c[,c('chr', 'start', 'end', 'cn')]
      
      allranges <- range
    }else{
      
      allranges <- c(allranges, range)
    }
    
  }
  
  # Now we have all ranges, populate them with sample data
  allranges <- disjoin(allranges)
  
  df <- data.frame(chr=seqnames(allranges), start=start(allranges), end=end(allranges))
  for(i in 1:length(cnadata)){
    pt <- my.names[i]
    df[,pt] <- NA
  }
  
  for(j in 1:length(cnadata)){
    pt <- my.names[j]
    print(paste("Parsing file", j, "/", length(cnadata), ":", pt))
    c <- cnadata[[j]]
    
    # For each segment, find all matching subsegments
    for(i in 1:dim(c)[1]){
      sel <- which(df$start >= c$start[i] & df$end <= c$end[i])
      df[sel,pt] <- c$cn[i]
    }
  }
  
  parsed.segmented <- df
  
  
  # Map to genes and bands
  load('resources/bands.dict.RData')
  bands <- bands.dict
  bands <- bands[which(bands$chrom %in% 1:22),]
  bands.cn <- data.frame(row.names = rownames(bands))
  
  # Combine files into a large matrix
  for(i in 1:length(cnadata)){
    if(i %% 10 == 0){print(paste("Mapping to bands: file", i, "/", length(cnadata)))}
    s <- cnadata[[i]]
    
    # For each band, find CN segments which overlap and take the mean Segment_Mean of those bands
    cns <- unlist(lapply(1:dim(bands)[1], function(x){
      segs <- s[which(s$chr == bands$chrom[x] & s$end > bands$start[x] & s$start < bands$end[x]), "cn"]
      if(length(segs) > 0){
        return(mean(segs, na.rm=T))
      }else{
        return(NA)
      }
    }))
    bands.cn[,my.names[i]] <- cns
  }
  
  # Map to genes
  
  refgenome <- BSgenome.Hsapiens.UCSC.hg19
  gene_coord <- genes(Homo.sapiens, columns=c("GENEID", "SYMBOL"))
  genes <- data.frame(gene_coord)
  genes <- genes[-which(is.na(genes$SYMBOL)),]
  genes$chr <- gsub("chr", "", as.character(genes$seqnames))
  
  genes.cn <- data.frame(row.names = genes$SYMBOL)
  for(i in 1:length(cnadata)){
    if(i %% 10 == 0){print(paste("Mapping to genes: file", i, "/", length(cnadata)))}
    s <- cnadata[[i]]
    
    genes.cn[,my.names[[i]]] <- NA
    # For each segment, find genes within it
    for(segid in 1:dim(s)[1]){
      sel <- which(genes$chr == s$chr[segid] & genes$end > s$start[segid] & genes$start < s$end[segid])
      genes.cn[sel,my.names[[i]]] <- s$cn[segid]
    }
  }
  
  
  return(list(parsed.segmented, bands.cn, genes.cn))
}
