# Clonality analysis uses the sciClone package
# By default sciClone ignores mutations outside CN-neutral, LOH-free areas.
# As our samples have extensive CN change we correct for local CN first, then
# cluster by Cancer Cell Fraction (CCF)
# Correction is based on forumla from Mcgranahan et al (http://science.sciencemag.org/content/sci/suppl/2016/03/02/science.aaf1490.DC1/McGranahan-SM.pdf)
# ccf = VAF * (1/p) * (p*CNt + CNn(1-p))
# Where p = purity, CNt = local tumour copy number, CNn = local normal copy number

if(!exists("muts.all")){
  load("data/wgsdata.RData")
}

source('utility_functions/parallel.setup.R')

# Make gene annotation data frame for sciClone - used to label mutations with genes only
annot <- muts.all[,c("chr", "start", "gene")]
colnames(annot) <- c("chr", "pos", "gene")
sel.rm <- which(duplicated(annot$pos) | is.na(annot$gene))
if(length(sel.rm) > 0){ annot <- annot[-sel.rm,] }

muts.cache.file <- paste0(results_dir, "clonality/muts.with.ccf.RData")
clones.cache.file <- paste0(results_dir, "clonality/clusterdata.RData")

if(file.exists(muts.cache.file)){
  load(muts.cache.file)
}else{
  muts.with.ccf <- foreach(i = 1:dim(wgs.pheno)[1]) %dopar% {
    library(GenomicDataCommons)
    source('utility_functions/absolute.cancer.cell.fraction.R')
    sample <- wgs.pheno$name[i]
    
    # Find the mutation and CNA data
    sample.muts <- muts.all[which(muts.all$patient == sample & muts.all$class == "SNV"),]
    sample.cnas <- cna.summary.list[which(cna.summary.list$sample == sample),]
    
    # Assign local tumour and normal CNs to sample.muts (from ASCAT output)
    sample.muts$CNt <- NA
    sample.muts$CNn <- NA
    for(j in 1:dim(sample.cnas)[1]){
      sel <- which(paste0("chr", as.character(sample.muts$chr)) == sample.cnas$Chromosome[j] & sample.muts$start <= sample.cnas$chromEnd[j] & sample.muts$end >= sample.cnas$chromStart[j])
      sample.muts$CNt[sel] <- sample.cnas$total.copy.number.inTumour[j]
      sample.muts$CNn[sel] <- sample.cnas$total.copy.number.inNormal[j]
    }
    sample.muts <- sample.muts[which(!is.na(sample.muts$CNt)),]
    
    # Calculate Cancer Cell Fraction (CCF) for each mutation
    ccf <- lapply(1:dim(sample.muts)[1], function(x){
      absolute.cancer.cell.fraction(
        n.alt=sample.muts$alt.reads[x], 
        depth=sample.muts$tumour.reads[x], 
        purity=wgs.pheno$purity[i], 
        local.copy.number=sample.muts$CNt[x]
      )
    })
    ccf <- rbindlist2(ccf)
    sample.muts <- cbind(sample.muts, ccf)
    return(sample.muts)
  }
  names(muts.with.ccf) <- wgs.pheno$name
  save(muts.with.ccf, file=muts.cache.file)
}


if(file.exists(clones.cache.file)){
  load(clones.cache.file)
}else{
  # Create cluster data for each sample from CCFs
  clusterdata <- foreach(i = 1:dim(wgs.pheno)[1]) %dopar% {
    library(sciClone)
    sample <- wgs.pheno$name[i]
    sample.muts <- muts.with.ccf[[sample]]
    
    # Create "VAFs" data frame for sciClone clustering (actually CCF data)
    # Do not use copy numbers here as we have already corrected for this
    # Our input "VAFs" are 100 * CCF / 2. sciClone expects clonality at VAF=0.5 -> correct by dividing CCF by 2
    # Also we expect alt/(ref+alt) to be the VAF, so set ref such that alt/(ref+alt)=CCF/2
    # -> ref = 2alt/CCF - alt
    sc.vafs <- sample.muts[,c("chr", "start", "ref.reads", "alt.reads", "ccf.est")]
    sc.vafs$ccf.est <- 100*sc.vafs$ccf.est / 2
    sc.vafs$ref.reads <- (2*sc.vafs$alt.reads / sc.vafs$ccf.est) - sc.vafs$alt.reads
    sc = sciClone(
      vafs=sc.vafs,
      sampleNames=sample,
      minimumDepth = 10,
      annotation = annot
    )
    return(sc)
  }
  
  names(clusterdata) <- wgs.pheno$name
  save(clusterdata, file=clones.cache.file)
}

