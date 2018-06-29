##########################################################################
# Load Whole Genome Sequencing Data 
# Loads and processes SNV, indel and copy number data into usable RData files
# Works from variant call files - not from raw BAM files
#
# We want to analyse:
#   Overall mutational burden (subs and indels)
#   All mutations - need details for PvR comparison
#   Coding mutations - need details for comparison with TCGA
#   CNAs - by gene and by band
#   Rearrangements
#
#
# Available variables are listed at the bottom of this file.
##########################################################################


if(!exists("data_cache")){
  data_cache <- "./data/"
}

cache_file <- paste(data_cache, "wgsdata.RData", sep="")

if(file.exists(cache_file)){
  message("Loading cached WGS data")
  load(cache_file)
}else{
  message("Parsing sequencing data. This may take some time.")
  source('data_loaders/downloadTcgaData.R')
  source('utility_functions/parseCnaFiles.R')
  
  # Location of WGS data files
  # These files are output from Caveman, ASCAT, Pindel and Brass. They are not raw BAM files.
  # Output VCF files have been annotated using Vagrent.
  # These processed files are available from the authors on request.
  wgs.data.dir <- paste(data_cache, "wgs/", sep="")
  
  subs.dir <- paste(wgs.data.dir, "caveman/", sep="")
  cna.raw.dir <- paste(wgs.data.dir, "ascat/", sep="")
  cna.dir <- paste(wgs.data.dir, "ascat_summary/", sep="")
  indel.dir <- paste(wgs.data.dir, "pindel/", sep="")
  brass.dir <- paste(wgs.data.dir, "brass/", sep="")
  
  # Pheno file is shared in RData format in this repository.
  load("resources/wgsPheno.RData")
  
  # Add purity and ploidy data - output from ASCAT analysis
  ascat.output <- read.csv('resources/Ascat_ploidy_purity.csv', stringsAsFactors = F)
  wgs.pheno$purity <- ascat.output$ABBR_CELL_FRAC[match(wgs.pheno$name, ascat.output$SAMPLE)]
  #wgs.pheno$ploidy <- ascat.output$TUM_PLOIDY[match(wgs.pheno$name, ascat.output$SAMPLE)]

  # Additionally mark some samples as 'query regressive' - these samples regressed but subsequently showed evidence of new disease on longer follow up 
  # These were identified from our analysis, hence not included in the input pheno file
  wgs.pheno$query.reg <- 0
  wgs.pheno$query.reg[which(wgs.pheno$name %in% c("PD21884a", "PD21893a", "PD38326a"))] <- 1
  
  library(stringr)
  library(VariantAnnotation)
  
  ####################################################################################
  # Concatenate subs and indels
  #
  # This creates a large list of mutations with all patients in one file
  # All mutations passing filters are included
  ####################################################################################
  for(i in 1:dim(wgs.pheno)[1]){
    pt <- wgs.pheno$name[i]
    print(paste("Processing patient", pt))
    
    subs.file <- list.files(subs.dir, pattern=pt, full.names = T)
    indel.file <- list.files(indel.dir, pattern=pt, full.names = T)
    brass.file <- list.files(brass.dir, pattern=pt, full.names = T)
    
    # Read subs
    vcf <- readVcf(subs.file)
    filters <- fixed(vcf)$FILTER
    clpm <- as.numeric(readInfo(subs.file, x='CLPM'))
    asmd <- as.numeric(readInfo(subs.file, x='ASMD'))
    
    # Read other variables
    vd <- as.character(readInfo(subs.file, x='VD'))
    vc <- as.character(readInfo(subs.file, x='VC'))
    genes <- unlist(lapply(as.character(vd), function(x){
      unlist(strsplit(x, "[|]"))[1]
    }))
    rr <- rowRanges(vcf)
    chr <- as.character(seqnames(rr))
    start <- start(ranges(rr))
    end <- end(ranges(rr))
    ref <- (as.character(ref(vcf)))
    alt <- (as.character(unlist(alt(vcf))))
    
    # Record VAF = proportion of mutant allele in tumour sample
    vaf <- geno(vcf)$PM[,2]
    depth <- as.character(readInfo(subs.file, x="DP"))
    # Get read counts - need to do for each of A,C,T,G for efficiency
    ref.reads <- rep(NA, length(ref))
    alt.reads <- rep(NA, length(alt))
    for(base in unique(ref)){
      ref.reads[which(ref == base)] <- as.numeric(geno(vcf)[[paste("F", base, "Z", sep="")]][,2])[which(ref == base)] + as.numeric(geno(vcf)[[paste("R", base, "Z", sep="")]][,2])[which(ref == base)]
      alt.reads[which(alt == base)] <- as.numeric(geno(vcf)[[paste("F", base, "Z", sep="")]][,2])[which(alt == base)] + as.numeric(geno(vcf)[[paste("R", base, "Z", sep="")]][,2])[which(alt == base)]
    }
    tumour.reads <- as.numeric(geno(vcf)$FAZ[,2]) + as.numeric(geno(vcf)$FCZ[,2]) + as.numeric(geno(vcf)$FGZ[,2]) + as.numeric(geno(vcf)$FTZ[,2]) + 
      as.numeric(geno(vcf)$RAZ[,2]) + as.numeric(geno(vcf)$RCZ[,2]) + as.numeric(geno(vcf)$RGZ[,2]) + as.numeric(geno(vcf)$RTZ[,2])
    
    
    # Combine useful data
    subs.pt <- data.frame(
      patient=pt,
      gene=genes,
      class="SNV",
      type=vc,
      ref=ref,
      alt=alt,
      chr=chr,
      start=start,
      end=end,
      filters=fixed(vcf)$FILTER,
      asmd=asmd,
      clpm=clpm,
      vaf=vaf,
      ref.reads=ref.reads,
      alt.reads=alt.reads,
      tumour.reads=tumour.reads,
      depth=depth,
      exonic=grepl("exon", vd),
      protein.change=gsub("p.", "", str_extract(vd, "p.[A-Z][0-9]+[A-Z]"), fixed = T)
    )
    
    
    # Read indels
    vcf <- readVcf(indel.file)
    vd <- as.character(readInfo(indel.file, x='VD'))
    vc <- readInfo(indel.file, x='VC')
    m <- as.matrix(vc)
    vc <- as.character(m[,1])
    
    # Pindel call - is it insertion or deletion
    pc <- as.character(readInfo(indel.file, x='PC'))
    
    rr <- rowRanges(vcf)
    chr <- as.character(seqnames(rr))
    start <- start(ranges(rr))
    end <- end(ranges(rr))
    ref <- as.character(ref(vcf))
    alt <- as.character(alt(vcf)@unlistData)
    filters <- as.character(fixed(vcf)$FILTER)
    
    genes <- unlist(lapply(vd, function(x){
      unlist(strsplit(x, "[|]"))[1]
    }))
    
    
    indels.pt <- data.frame(
      patient=pt,
      gene=genes,
      class=pc,
      type=vc,
      ref=ref,
      alt=alt,
      chr=chr,
      start=start,
      end=end,
      filters=filters,
      asmd=1000, # Indels should automatically pass these filters
      clpm=0,
      vaf=NA, # Clonality based on subs only
      ref.reads=NA,
      alt.reads=NA,
      tumour.reads=NA,
      depth=NA,
      exonic=grepl("exon", vd),
      protein.change=gsub("p.", "", str_extract(vd, "p.[A-Z][0-9]+[A-Z]"), fixed = T)
    )
    
    # Read rearrangements
    # has.rearr <- length(brass.file) > 0
    # if(has.rearr){
    #   vcf <- readVcf(brass.file)
    #   rr <- rowRanges(vcf)
    #   rearr.pt <- data.frame(
    #     patient=pt,
    #     gene=readInfo(brass.file, x='GENE'),
    #     class="Rearrangement",
    #     type=readInfo(brass.file, x='SVTYPE'),
    #     ref=as.character(ref(vcf)),
    #     alt=as.character(alt(vcf)),
    #     chr=as.character(seqnames(rr)),
    #     start=start(ranges(rr)),
    #     end=end(ranges(rr)),
    #     filters=fixed(vcf)$FILTER,
    #     asmd=1000, # Rearrangements should automatically pass these filters
    #     clpm=0,
    #     vaf=NA, # Clonality based on subs only
    #     ref.reads=NA,
    #     alt.reads=NA,
    #     tumour.reads=NA,
    #     depth=NA,
    #     exonic=F, # We don't compare rearrangements to TCGA so can skip this
    #     protein.change=NA
    #   )
    # }
    
    
    # Merge subs and indels
    muts.pt <- rbind(subs.pt, indels.pt)
    # Order by location
    o <- order(muts.pt$chr, muts.pt$start)
    muts.pt <- muts.pt[o,]
    
    # Merge
    if(i == 1){
      subs.all <- subs.pt
      indels.all <- indels.pt
      if(has.rearr){ rearr.all <- rearr.pt }
      muts.all <- muts.pt
    }else{
      subs.all <- rbind(subs.all, subs.pt)
      indels.all <- rbind(indels.all, indels.pt)
      if(has.rearr){ rearr.all <- rbind(rearr.all, rearr.pt) }
      muts.all <- rbind(muts.all, muts.pt)
    }
  }
  
  # Save intermediary step
  save(muts.all, file="data/wgsTemp.RData")
  
  # Define our filters (but don't perform filtering yet):
  muts.all$filters.passed <- muts.all$filters == "PASS" & muts.all$asmd >= 140& muts.all$clpm == 0
  
  ####################################################################################
  # Check consistency between samples
  #
  # If a mutation is confidently called in one sample, and is present in another sample
  # from the same patient but failed filters, add it to the mutation list
  # Here we refer to pileup files directly so we do not miss mutations not called by CaveMan
  ####################################################################################
  pileup.dir <- paste0(wgs.data.dir, "pileups/")
  
  # Define a unique mutation string to compare between patients
  muts.all$mid <- paste(
    muts.all$chr, muts.all$start, muts.all$end, muts.all$ref, muts.all$alt, sep="-"
  )
  # Check each patient with multiple samples
  pts <- unique(wgs.pheno$Patient[which(duplicated(wgs.pheno$Patient))])
  for(pt in pts){
    samples <- wgs.pheno$name[which(wgs.pheno$Patient == pt)]
    # samples.data <- muts.all[which(muts.all$patient %in% samples),]
    # mt.ids <- samples.data$mid[which(samples.data$filters.passed)]
    print(paste("Checking patient", pt, "(", length(samples), "samples )"))
    # Find actual mutations. 
    # All of these mutations have been confirmed in at least one sample for this patient.
    # Therefore assume they are also present in other samples, if they have failed filters
    pileup.file <- list.files(pileup.dir, pattern=substr(samples[[1]], 1, 7), full.names=T)
    pileup <- read.table(pileup.file, sep="\t", header=T, row.names = NULL, skip = 54+length(samples))
    pileup$mid <- paste(pileup$Chrom, pileup$Pos, pileup$Pos, pileup$Ref, pileup$Alt, sep="-")
    for(sample in samples){
      # Mutations present in this sample from pileup data
      pileup.pt <- pileup[which(!is.na(pileup[,paste0(sample, "_OFS")])),]
      
      # Mutations not picked up by CaveMan:
      sel.missing <- which(!(pileup.pt$mid %in% muts.all$mid[which(muts.all$patient == sample)]))
      # Mutations picked up by CaveMan but failing filters:
      sel.failed  <- which(pileup.pt$mid %in% muts.all$mid[which(muts.all$patient == sample)] & !(pileup.pt$mid %in% muts.all$mid[which(muts.all$patient == sample & muts.all$filters.passed == T)]))
      
      if(length(sel.failed) > 0){
        print(paste("Rescued", length(sel.failed), "failed mutations from", sample))
        muts.all$filters.passed[which(muts.all$patient == sample & muts.all$mid %in% pileup.pt$mid[sel.failed])] <- TRUE
      }
      
      if(length(sel.missing) > 0){
        print(paste("Rescued", length(sel.missing), "missing mutations from", sample))
        df <- data.frame(
          patient=sample,
          gene=pileup$Gene[sel.missing],
          class="SNV",
          type=pileup$Effect[sel.missing],
          ref=pileup$Ref[sel.missing],
          alt=pileup$Alt[sel.missing],
          chr=pileup$Chrom[sel.missing],
          start=pileup$Start[sel.missing],
          end=pileup$End[sel.missing],
          filters=pileup[sel.missing, paste0(sample, "_OFS")],
          asmd=140,
          clpm=0,
          vaf=pileup[sel.missing, paste0(sample, "_MTR")] / pileup[sel.missing, paste0(sample, "_DEP")],
          ref.reads=pileup[sel.missing, paste0(sample, "_WTR")],
          alt.reads=pileup[sel.missing, paste0(sample, "_MTR")],
          tumour.reads=pileup[sel.missing, paste0(sample, "_DEP")],
          depth=pileup[sel.missing, paste0(sample, "_DEP")],
          exonic= pileup$VD[sel.missing],
          protein.change=pileup$Protein[sel.missing],
          filters.passed=T,
          mid=pileup$mid[sel.missing]
        )
        muts.all <- rbind(muts.all, df)
      }
    }

    # sel <- which(muts.all$mid %in% mt.ids & muts.all$patient %in% samples)
    # print(paste("Correcting", length(which(!muts.all$filters.passed[sel])), "mismatches, from total", dim(samples.data)[1], "mutations"))
    # muts.all$filters.passed[sel] <- TRUE
    
  }
  
  
  
  # Define a unique mutation string to compare between patients
  # muts.all$mid <- paste(
  #   muts.all$chr, muts.all$start, muts.all$end, muts.all$ref, muts.all$alt, sep="-"
  # )
  # # Check each patient with multiple samples
  # pts <- unique(wgs.pheno$Patient[which(duplicated(wgs.pheno$Patient))])
  # for(pt in pts){
  #   samples <- wgs.pheno$name[which(wgs.pheno$Patient == pt)]
  #   print(paste("Checking patient", pt, "(", length(samples), "samples )"))
  #   # Find actual mutations. 
  #   # All of these mutations have been confirmed in at least one sample for this patient.
  #   # Therefore assume they are also present in other samples, if they have failed filters
  #   samples.data <- muts.all[which(muts.all$patient %in% samples),]
  #   mt.ids <- samples.data$mid[which(samples.data$filters.passed)]
  #   sel <- which(muts.all$mid %in% mt.ids & muts.all$patient %in% samples)
  #   print(paste("Correcting", length(which(!muts.all$filters.passed[sel])), "mismatches, from total", dim(samples.data)[1], "mutations"))
  #   muts.all$filters.passed[sel] <- TRUE
  #   
  # }
  
  ####################################################################################
  # Filter subs and indels
  #
  # For comparison with TCGA data, we need to access only coding subs/indels
  ####################################################################################
  muts.unfiltered <- muts.all
  muts.all <- muts.unfiltered[which(muts.all$filters.passed),]
  
  # Extract only exonic mutations (for TCGA comparisons)
  sel <- which(muts.all$exonic)
  muts.coding <- muts.all[sel,]
  
  # Coerce into a data frame showing per-patient mutations
  # Muts shows mutation counts; muts.detail shows the type of mutation(s)
  for(i in 1:dim(wgs.pheno)[1]){
    pt <- wgs.pheno$name[i]
    t <- table(muts.coding$gene[which(muts.coding$patient == pt)])
    df <- data.frame(gene=names(t), count=as.numeric(t))
    colnames(df) <- c("gene", pt)
    if(i == 1){
      muts <- df
    }else{
      muts <- merge(muts, df, by="gene", all=T)
    }
  }
  rownames(muts) <- muts$gene
  muts$gene <- NULL
  
  # Process mutation frequencies for TCGA comparison. 
  # Ignore multiple mutations per sample for this, we want to know the percentage of samples with a particular gene mutated.
  muts.coding.freq <- apply(muts, 1, function(x){
    100 * length(which(x > 0)) / dim(muts)[2]
  })
  
  # Count coding mutations by sample
  muts.coding.counts <- apply(muts, 2, sum)
  
  ####################################################################################
  # Add Rearrangements
  #
  # These are pre-processed as described in the main text.
  # 
  ####################################################################################
  rearrangements.all <- read.table('resources/private/no_germ_reassembled_3_reads.txt', header=T, stringsAsFactors = F, sep="\t")
  
  rearrs <- data.frame(
    patient=rearrangements.all$sample,
    gene=rearrangements.all$gene1,
    class="R",
    type="Rearrangement",
    ref=NA,
    alt=NA,
    chr=rearrangements.all$chr1,
    start=rearrangements.all$start1,
    end=rearrangements.all$end1,
    filters="PASS",
    asmd=1000, # Should automatically pass these filters
    clpm=0,
    vaf=NA, # Clonality based on subs only
    ref.reads=NA,
    alt.reads=NA,
    tumour.reads=NA,
    depth=NA,
    exonic=grepl("exon", rearrangements.all$region1) | grepl("exon", rearrangements.all$region2),
    protein.change=NA,
    filters.passed=T,
    translocation.partner=rearrangements.all$gene2,
    mid=paste(rearrangements.all$chr1, rearrangements.all$start1, rearrangements.all$start2, rearrangements.all$chr2, rearrangements.all$start2, rearrangements.all$end2),
    chr2=rearrangements.all$chr2,
    start2=rearrangements.all$start2,
    end2=rearrangements.all$end2
  )
  
  muts.all$translocation.partner <- NA
  muts.all$chr2 <- NA
  muts.all$start2 <- NA
  muts.all$end2 <- NA
  
  muts.all <- rbind(muts.all, rearrs)
  
  ####################################################################################
  # Load Copy Number Alteration (CNA) summary profiles
  #
  # Load from summary files
  # Reduce to minimum consistent regions across all samples
  # These data are used for genome-wide plotting of copy number, as in figure 2E
  ####################################################################################
  cna.s.files <- list.files(cna.dir, full.names = T)
  
  for(i in 1:length(cna.s.files)){
    pt_index <- regexpr("PD[0-9]+", cna.s.files[i])
    pt <- substr(cna.s.files[i], pt_index, pt_index+7)
    print(paste("Reading CNAS for patient", pt))
    c <- read.table(cna.s.files[i], sep=",", row.names = 1, stringsAsFactors = F)
    colnames(c) <- c('chr', 'start', 'end', 'pl', 'pl.min', 'cn', 'cn.min')
    
    if(i == 1){
      scnas <- c[,c('chr', 'start', 'end', 'cn')]
      colnames(scnas)[4] <- pt
      
      allranges <- GRanges(
        seqnames = Rle(c$chr),
        ranges = IRanges(start=c$start, end=c$end)
      )
    }else{
      range <- GRanges(
        seqnames = Rle(c$chr),
        ranges = IRanges(start=c$start, end=c$end)
      )
      allranges <- c(allranges, range)
    }
  }
  
  # Now we have all ranges, populate them with sample data
  allranges <- disjoin(allranges)
  
  df <- data.frame(chr=seqnames(allranges), start=start(allranges), end=end(allranges))
  for(i in 1:length(cna.s.files)){
    pt_index <- regexpr("PD[0-9]+", cna.s.files[i])
    pt <- substr(cna.s.files[i], pt_index, pt_index+7)
    df[,pt] <- NA
  }
  
  
  for(j in 1:length(cna.s.files)){
    pt_index <- regexpr("PD[0-9]+", cna.s.files[j])
    pt <- substr(cna.s.files[j], pt_index, pt_index+7)
    c <- read.table(cna.s.files[j], sep=",", row.names = 1, stringsAsFactors = F)
    colnames(c) <- c('chr', 'start', 'end', 'pl', 'pl.min', 'cn', 'cn.min')
    
    for(i in 1:dim(df)[1]){
      
      sel <- which(
        c$chr == df$chr[i] & 
          c$start <= df$start[i] &
          c$end >= df$end[i]
      )
      if(length(sel) > 1){
        message(paste("WARNING: Multiple matching rows",pt, df$chr[i], df$start[i], df$end[i]))
        next
      }
      if(length(sel) == 0){
        message(paste("WARNING: No rows matched",pt, df$chr[i], df$start[i], df$end[i]))
        next
      }
      df[i, pt] <- c[sel, 'cn']
    }
  }
  
  
  # There are some gaps - assume these have no CN change so fill in with 2s:
  sel <- which(is.na(df), arr.ind = T)
  df[sel] <- 2
  
  cnas.segmented <- df
  
  # Include the mean CNA segment values across all samples
  cnas.segmented.mean <- cnas.segmented[,1:3]
  cnas.segmented.mean$cn <- apply(cnas.segmented[,4:dim(cnas.segmented)[2]], 1, mean)
  
  # Sometimes we can collapse this into a smaller data frame using Genomic Ranges if adjacent regions have the same mean:
  range <- GRanges(seqnames = Rle(cnas.segmented.mean$chr), 
                   ranges=IRanges(start=cnas.segmented.mean$start, end=cnas.segmented.mean$end),
                   cn=cnas.segmented.mean$cn)
  x <- sort(unlist(split(range, elementMetadata(range)$cn)))
  cnas.segmented.mean <- data.frame(
    chr=seqnames(x), start=start(x), end=end(x), cn=elementMetadata(x)$cn
  )
  
  ####################################################################################
  # Estimate ploidy for each sample and add to the pheno data frame
  ####################################################################################
  source('utility_functions/ploidyFunctions.R')
  names <- colnames(cnas.segmented)[4:dim(cnas.segmented)[2]]
  ploidys <- unlist(lapply(names, function(x){
    df <- data.frame(
      SampleID=x,
      Chr=cnas.segmented$chr,
      Start=cnas.segmented$start,
      End=cnas.segmented$end,
      n.probes=0,
      cn=cnas.segmented[,x]
    )
    ploidy <- fun.ploidy(x, df)
    return(as.numeric(ploidy))
  }))
  names(ploidys) <- names
  wgs.pheno$ploidy <- ploidys[wgs.pheno$name]
  
  
  
  ####################################################################################
  # Load CNAs by gene and by band
  #
  # The above ASCAT analysis gives genome-wide profile results.
  # To identify which genes are affected directly, we use the below:
  ####################################################################################
  read.ascat <- function(FILE.CN) {
    cv.data <- read.table(FILE.CN, header=FALSE, sep=',') # ASCAT
    cv.data <- cv.data[,1:8]
    cat( paste( dim(cv.data)[1], ' copy-number segments \n'))
    colnames(cv.data) <- c('seg_no', 'Chromosome', 'chromStart', 'chromEnd', 'total.copy.number.inNormal', 'minor.copy.number.inNormal', 'total.copy.number.inTumour', 'minor.copy.number.inTumour')
    cv.data$seg_no <- NULL
    cv.data$Chromosome <- as.character(cv.data$Chromosome)
    cv.data$Chromosome[cv.data$Chromosome=='23'] <- 'X'
    cv.data$Chromosome[cv.data$Chromosome=='24'] <- 'Y'
    
    if (nrow(cv.data)>0) {
      cv.data$Chromosome <- paste('chr', cv.data$Chromosome,sep='')
    }
    
    cv.data$major.copy.number.inTumour <- cv.data$total.copy.number.inTumour - cv.data$minor.copy.number.inTumour
    
    cv.data$major.copy.number.inTumour.temp <- pmax(cv.data$major.copy.number.inTumour, cv.data$minor.copy.number.inTumour)
    cv.data$minor.copy.number.inTumour.temp <- pmin(cv.data$major.copy.number.inTumour, cv.data$minor.copy.number.inTumour)
    
    cv.data$major.copy.number.inTumour <- cv.data$major.copy.number.inTumour.temp
    cv.data$minor.copy.number.inTumour <- cv.data$minor.copy.number.inTumour.temp
    cv.data$major.copy.number.inTumour.temp <- NULL
    cv.data$minor.copy.number.inTumour.temp <- NULL
    
    return(cv.data)
  }
  
  cna.summaries <- list()
  for (f in cna.s.files) {
    pt_index <- regexpr("PD[0-9]+", f)
    pt <- substr(f, pt_index, pt_index+7)
    ptdata <- read.ascat(f)
    ptdata$sample <- pt
    cna.summaries[[pt]] <- ptdata
  }
  library(GenomicDataCommons)
  cna.summary.list <- rbindlist2(cna.summaries)
  
  # FIND ALL AMPLIFIED AND DELETED GENES IN MULTIPLE SAMPLES 
  library(GenomicRanges)
  library(Homo.sapiens)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(AnnotationHub)
  library(biovizBase)
  
  # Get genes
  refgenome <- BSgenome.Hsapiens.UCSC.hg19
  gene_coord <- genes(Homo.sapiens, columns=c("GENEID", "SYMBOL"))
  
  # Get cytobands
  band_coord <- biovizBase::getIdeogram("hg19", cytobands=TRUE)

  amps <- list()
  dels <- list()
  amps.band <- list()
  dels.band <- list()
  for (i in 1:length(cna.summaries)) {
    pt <- names(cna.summaries)[i]
    samp <- cna.summaries[[i]]
    # Just get the regions of amplification - defined as 2*ploidy + 1
    amp.limit <- 2*wgs.pheno$ploidy[which(wgs.pheno$name == pt)] + 1
    amp <- samp[samp$total.copy.number.inTumour > amp.limit, ]
    amp_loci <- GRanges(seqnames=amp$Chromosome, ranges=IRanges(start=amp$chromStart, end=amp$chromEnd), seqinfo=seqinfo(refgenome))
    # Find genes
    amplified_genes <- subsetByOverlaps(gene_coord, amp_loci)
    amplified_df <- as.data.frame(amplified_genes)
    amplified_df <- amplified_df[order(amplified_df$seqnames, amplified_df$start),]
    amps[[pt]] <- amplified_df
    # Find bands
    amplified_bands <- subsetByOverlaps(band_coord, amp_loci)
    amplified_df <- as.data.frame(amplified_bands)
    amplified_df <- amplified_df[order(amplified_df$seqnames, amplified_df$start),]
    amplified_df$band <- paste(amplified_df$seqnames, amplified_df$name, sep="")
    amps.band[[pt]] <- amplified_df
    
    # Just get the regions of homozygous deletion
    homdel <- samp[samp$total.copy.number.inTumour==0, ]
    homdel_loci <- GRanges(seqnames=homdel$Chromosome, ranges=IRanges(start=homdel$chromStart, end=homdel$chromEnd), seqinfo=seqinfo(refgenome))
    # Find genes
    deleted_genes <- subsetByOverlaps(gene_coord, homdel_loci)
    deleted_df <- as.data.frame(deleted_genes)
    deleted_df <- deleted_df[deleted_df$seqnames != "chrY",] # Ignore Y chromosome
    deleted_df <- deleted_df[order(deleted_df$seqnames, deleted_df$start),]
    dels[[pt]] <- deleted_df
    # Find bands
    deleted_bands <- subsetByOverlaps(band_coord, homdel_loci)
    deleted_df <- as.data.frame(deleted_bands)
    deleted_df <- deleted_df[deleted_df$seqnames != "chrY",] # Ignore Y chromosome
    deleted_df <- deleted_df[order(deleted_df$seqnames, deleted_df$start),]
    deleted_df$band <- paste(deleted_df$seqnames, deleted_df$name, sep="")
    dels.band[[pt]] <- deleted_df
    
  }
  cnas.amps <- amps
  cnas.dels <- dels
  cnas.amps.band <- amps.band
  cnas.dels.band <- dels.band
  
  # Combine into a summary data frame with genes as rows, patients as columns, -1 for del, +1 for amps, 0 otherwise
  allgenes <- unique(c(
    as.character(unlist(lapply(cnas.amps, function(x){return(as.character(x$SYMBOL))}))),
    as.character(unlist(lapply(cnas.dels, function(x){return(as.character(x$SYMBOL))})))
  ))
  sel.na <- which(is.na(allgenes))
  if(length(sel.na) > 0){
    allgenes <- allgenes[-sel.na]
  }
    
  cnas.genes.summary <- data.frame(matrix(0, nrow=length(allgenes), ncol=length(cnas.amps)))
  rownames(cnas.genes.summary) <- allgenes
  colnames(cnas.genes.summary) <- names(cnas.amps)
  
  for(j in 1:dim(cnas.genes.summary)[2]){
    pt <- colnames(cnas.genes.summary)[j]
    print(pt)
    amps <- which(rownames(cnas.genes.summary) %in% as.character(cnas.amps[[pt]]$SYMBOL))
    dels <- which(rownames(cnas.genes.summary) %in% as.character(cnas.dels[[pt]]$SYMBOL))
    
    cnas.genes.summary[amps, j] <- 1
    cnas.genes.summary[dels, j] <- -1
  }
  
  # Combine band data
  allbands <- unique(c(
    as.character(unlist(lapply(cnas.amps.band, function(x){return(as.character(x$band))}))),
    as.character(unlist(lapply(cnas.dels.band, function(x){return(as.character(x$band))})))
  ))
  sel.na <- which(is.na(allbands))
  if(length(sel.na) > 0){
    allbands <- allbands[-sel.na]
  }
  
  cnas.bands.summary <- data.frame(matrix(0, nrow=length(allbands), ncol=length(cnas.amps.band)))
  rownames(cnas.bands.summary) <- allbands
  colnames(cnas.bands.summary) <- names(cnas.amps.band)
  
  for(j in 1:dim(cnas.bands.summary)[2]){
    pt <- colnames(cnas.bands.summary)[j]
    amps <- which(rownames(cnas.bands.summary) %in% as.character(cnas.amps.band[[pt]]$band))
    dels <- which(rownames(cnas.bands.summary) %in% as.character(cnas.dels.band[[pt]]$band))
    
    cnas.bands.summary[amps, j] <- 1
    cnas.bands.summary[dels, j] <- -1
  }
  
  # Order by the most affected bands first:
  counts <- apply(cnas.bands.summary, 1, function(x){
    return(length(which(x != 0)))
  })
  o <- order(counts, decreasing = T)
  cnas.bands.summary <- cnas.bands.summary[o,]
  
  ####################################################################################
  # Weighted Genome Integrity Index (WGII)
  #
  # Calculate this general measure of aneuploidy across the genome
  ####################################################################################
  
  getWgii <- function(data){
    # data should have columns chr, start, end, cn
    # For each chromosome (except X,Y) find the percentage for which CN !=2
    # Then take the mean across all chromosomes
    chrdata <- c()
    for(chr in unique(data$chr)){
      if(!(chr %in% 1:22)){ next }
      segs <- data[which(data$chr == chr),]
      seglengths <- as.numeric(segs$end) - as.numeric(segs$start)
      x <- sum(seglengths[which(segs$cn !=2)]) / sum(seglengths)
      chrdata <- c(chrdata, x)
    }
    return(mean(chrdata))
  }
  
  # Calculate and append to wgs.pheno
  wgs.pheno$wgii <- NA
  
  for(file in cna.s.files){
    pt_index <- regexpr("PD[0-9]+", file)
    pt <- substr(file, pt_index, pt_index + 7)
    if(!(pt %in% wgs.pheno$name)){ next }
    segs <- read.csv(file, header=F)
    ptdata <- data.frame(chr=segs$V2, start=segs$V3, end=segs$V4, cn=segs$V7)
    wgs.pheno$wgii[which(wgs.pheno$name == pt)] <- getWgii(ptdata)
  }
  
  
  ####################################################################################
  # Combine the above
  #
  # Combine SNV, indel, CNA data (by gene) to get master list of affected genes
  # Include a data frame muts.nocn that excludes copy number alterations for comparison with TCGA somatic mutations
  ####################################################################################
  # genes <- unique(c(rownames(subs), rownames(indels), rownames(cnas.genes.summary)))
  # muts <- data.frame(row.names = genes)
  # muts.nocn <- data.frame(row.names = genes)
  # muts.detail <- data.frame(row.names=genes)
  # 
  # for(i in 1:dim(wgs.pheno)[1]){
  #   pt <- wgs.pheno$name[i]
  #   # Simple 0/1 mutation DB:
  #   print(paste("Merging patient", pt))
  #   muts[,pt] <- unlist(lapply(genes, function(x){
  #     subs[x, pt] > 0 | indels[x, pt] > 0 | cnas.genes.summary[x, pt] != 0
  #   }))
  #   muts.nocn[,pt] <- unlist(lapply(genes, function(x){
  #     subs[x, pt] > 0 | indels[x, pt] > 0
  #   }))
  #   
  #   # More complex naming type of mutation:
  #   muts.detail[,pt] <- unlist(lapply(genes, function(x){
  #     gene.muts <- c()
  #     subs.sel <- which(subs.all$patient == pt & subs.all$gene == x)
  #     if(length(subs.sel) > 0){
  #       gene.muts <- c(gene.muts, unique(as.character(subs.all[subs.sel, 'type'])))
  #     }
  #     indels.sel <- which(indels.all$patient == pt & indels.all$gene == x)
  #     if(length(indels.sel) > 0){
  #       gene.muts <- c(gene.muts, unique(as.character(indels.all[indels.sel, 'type'])))
  #     }
  #     if(x %in% rownames(cnas.genes.summary)){
  #       if(cnas.genes.summary[x, pt] == 1){ gene.muts <- c(gene.muts, "CN amplification") }
  #       if(cnas.genes.summary[x, pt] == -1){ gene.muts <- c(gene.muts, "CN deletion") }
  #     }
  #     
  #     
  #     return(paste(gene.muts, collapse="|"))
  #     
  #   }))
  # }
  
  
  ####################################################################################
  # TCGA SNV data
  #
  # We compare our coding substitutions to TCGA data
  # Download that here
  # The file used is generated from the TCGA portal using our list of 200 driver genes as input
  ####################################################################################

  
  # tcga.snvs <- downloadTcgaData(
  #   d_type="Masked Somatic Mutation", 
  #   w_type = "MuTect2 Variant Aggregation and Masking", 
  #   cache_dir="./data/wgs/tcga"
  # )
  # 
  # tcga.mutect <- downloadTcgaData(
  #   d_type="Masked Somatic Mutation",
  #   w_type = "MuTect2 Variant Aggregation and Masking",
  #   cache_dir="./data/wgs/tcga"
  # )
  # tcga.varscan <- downloadTcgaData(
  #   d_type="Masked Somatic Mutation",
  #   w_type = "VarScan2 Variant Aggregation and Masking",
  #   cache_dir="./data/wgs/tcga"
  # )
  # tcga.ss <- downloadTcgaData(
  #   d_type="Masked Somatic Mutation",
  #   w_type = "SomaticSniper Variant Aggregation and Masking",
  #   cache_dir="./data/wgs/tcga"
  # )
  # tcga.muse <- downloadTcgaData(
  #   d_type="Masked Somatic Mutation",
  #   w_type = "MuSE Variant Aggregation and Masking",
  #   cache_dir="./data/wgs/tcga"
  # )
  # # Merge all four SNV calling datasets
  # # First add a UUID to identify the precise mutation
  # tcga.mutect$uuid <- paste(tcga.mutect$case_id, tcga.mutect$Chromosome, tcga.mutect$Start_Position, tcga.mutect$End_Position, tcga.mutect$Reference_Allele, tcga.mutect$Tumor_Seq_Allele1, tcga.mutect$Tumor_Seq_Allele2, sep="-")
  # tcga.varscan$uuid <- paste(tcga.varscan$case_id, tcga.varscan$Chromosome, tcga.varscan$Start_Position, tcga.varscan$End_Position, tcga.varscan$Reference_Allele, tcga.varscan$Tumor_Seq_Allele1, tcga.varscan$Tumor_Seq_Allele2, sep="-")
  # tcga.ss$uuid <- paste(tcga.ss$case_id, tcga.ss$Chromosome, tcga.ss$Start_Position, tcga.ss$End_Position, tcga.ss$Reference_Allele, tcga.ss$Tumor_Seq_Allele1, tcga.ss$Tumor_Seq_Allele2, sep="-")
  # tcga.muse$uuid <- paste(tcga.muse$case_id, tcga.muse$Chromosome, tcga.muse$Start_Position, tcga.muse$End_Position, tcga.muse$Reference_Allele, tcga.muse$Tumor_Seq_Allele1, tcga.muse$Tumor_Seq_Allele2, sep="-")
  # 
  # # Merge data frames and remove duplicates
  # tcga.snvs <- rbind(tcga.mutect, tcga.varscan, tcga.ss, tcga.muse)
  # tcga.snvs <- tcga.snvs[-which(duplicated(tcga.snvs$uuid)),]
  # 
  # # This data frame gives a good estimate of mutational burden in coding regions (TCGA uses WXS)
  # 
  # # Alternative driver method using TCGA portal downloads:
  # # tcga.drivers <- read.csv("resources/tcga.driver.genes.csv", header=T, stringsAsFactors=F)
  # # tcga.drivers$cases <- as.numeric(unlist(lapply(tcga.drivers$X..Affected.Cases.in.Cohort, function(x){
  # #   unlist(strsplit(x, " / "))[[1]]
  # # })))
  # # tcga.drivers$pc <- 100 * tcga.drivers$cases / 504
  # 
  # 
  # # Calculate mutation rates (i.e. the number of samples with 1+ mutation in a given gene)
  # # Remove multiple mutations per gene, per sample
  # tcga.snvs$pt.gene <- paste(tcga.snvs$Tumor_Sample_UUID, tcga.snvs$Hugo_Symbol)
  # tcga.snvs.rmdups <- tcga.snvs[-which(duplicated(tcga.snvs$pt.gene)),]
  # tcga.snvs.rates <- table(tcga.snvs.rmdups$Hugo_Symbol) / length(unique(tcga.snvs.rmdups$Tumor_Sample_UUID))
  
  
  # Load TCGA mutation rates
  # Data are generated from Caveman calls of TCGA data, using the same algorithms as for our CIS mutation calls
  # Data are similar to those presented on the TCGA portal
  tcga.snvs.rates.data <- read.table('resources/TCGA_prop_samples_with_subs_or_indels_in_gene.txt', header=T, stringsAsFactors = F)
  tcga.snvs.rates <- as.list(tcga.snvs.rates.data$V2)
  names(tcga.snvs.rates) <- tcga.snvs.rates.data$V1


  # TCGA CNA data
  # Relative GISTIC data used for comparison with methylation-derived CNA probes in prediction
  x <- downloadTcgaData(
    d_type="Copy Number Segment",
    w_type = "DNAcopy",
    cache_dir="./data/wgs/tcga_cnv"
  )
  # tcga.cnas.segmented <- x[[1]]
  tcga.cnas.bands     <- x[[2]]
  tcga.cnas.genes     <- x[[3]]
  tcga.cnas.pheno <- x[[4]]
  
  # TCGA CNA data - ASCAT
  # ASCAT CN data used for comparison with CIS ASCAT data (fig 2)
  load("resources/private/tcga.lusc.seg.hg19.rdata")
  tcga.cnas.segmented.data <- seg.mat.copy.list$segments
  tcga.cna.samples <- unique(tcga.cnas.segmented.data$SampleID)
  
  # Convert to minimum consistent regions format
  tcga.cnas.segmented <- convertCNAtoMCR(tcga.cnas.segmented.data)
  
  # Correct for ploidy
  tcga.ploidys <- lapply(tcga.cna.samples, function(x){
    fun.ploidy(x, tcga.cnas.segmented.data)
  })
  names(tcga.ploidys) <- tcga.cna.samples
  for(sample in tcga.cna.samples){
    tcga.cnas.segmented[,sample] <- tcga.cnas.segmented[,sample] / as.numeric(tcga.ploidys[[sample]])
  }
  # Add a mean value for tcga.cnas.segmented
  # Include the mean CNA segment values across all cancer samples, corrected for ploidy
  tcga.cnas.segmented.mean <- tcga.cnas.segmented[,1:3]
  tcga.cnas.segmented.mean$cn <- apply(tcga.cnas.segmented[,4:dim(tcga.cnas.segmented)[2]] , 1, mean )
  
  
  
  ####################################################################################
  # Save output in RData format
  ####################################################################################
  save(
    subs.all, indels.all, rearrangements.all,
    cna.summary.list,
    cnas.segmented, cnas.segmented.mean,
    cnas.genes.summary, cnas.amps, cnas.dels, 
    cnas.bands.summary, cnas.amps.band, cnas.dels.band,
    muts, muts.all, muts.unfiltered, muts.coding, muts.coding.counts,
    tcga.snvs.rates,
    tcga.cnas.segmented, tcga.cnas.segmented.mean,
    tcga.cnas.bands, tcga.cnas.genes, tcga.cnas.pheno, #tcga.cnas.ploidys,
    wgs.pheno, 
    file=cache_file
  )
  
}