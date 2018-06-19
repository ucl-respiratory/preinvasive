plot.genomic.clonality.matrix <- function(filename){
  pdf(filename)
  # For patients with multiple samples per patient, assess the proportion of mutations shared between each sample
  # Plot as a large heatmap-matrix
  pts.uniq <- unique(wgs.pheno$Patient[which(duplicated(wgs.pheno$Patient))])
  samps.uniq <- wgs.pheno[which(wgs.pheno$Patient %in% pts.uniq),] 
  o <- order(samps.uniq$progression, samps.uniq$Patient)
  samps.uniq <- samps.uniq[o,]
  pt.samples <- samps.uniq$name
  
  m <- matrix(ncol = length(pt.samples), nrow = length(pt.samples))
  for(i in 1:length(pt.samples)){
    for(j in 1:length(pt.samples)){
      # Only fill in if samples are from the same patient
      pt.i <- wgs.pheno$Patient[which(wgs.pheno$name == pt.samples[i])]
      pt.j <- wgs.pheno$Patient[which(wgs.pheno$name == pt.samples[j])]
      if(!(pt.i == pt.j)){ next }
      m[i,j] <- 
        length(intersect(
          muts.all$mid[which(muts.all$patient == pt.samples[i])],
          muts.all$mid[which(muts.all$patient == pt.samples[j])]
        )) / length(union(
          muts.all$mid[which(muts.all$patient == pt.samples[i])],
          muts.all$mid[which(muts.all$patient == pt.samples[j])]
        ))
    }
  }
  pt.samples <- paste(pt.samples, ifelse(wgs.pheno$progression[match(pt.samples, wgs.pheno$name)] == 1, "P", "R"), sep="-")
  rownames(m) <- pt.samples
  colnames(m) <- pt.samples
  breaksList = seq(0, 1, by = 0.01)
  pheatmap(m, cluster_rows = F, cluster_cols = F, display_numbers = T, breaks = breaksList)
  
  
  dev.off()
}
