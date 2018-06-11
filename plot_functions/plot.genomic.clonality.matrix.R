plot.genomic.clonality.matrix <- function(filename){
  pdf(filename)
  # For patients with multiple samples per patient, assess the proportion of mutations shared between each sample
  pts.uniq <- unique(wgs.pheno$Patient)
  for(pt in pts.uniq){
    pt.samples <- wgs.pheno[which(wgs.pheno$Patient == pt),]$name
    if(length(pt.samples) == 1){next}
    m <- matrix(ncol = length(pt.samples), nrow = length(pt.samples))
    for(i in 1:length(pt.samples)){
      for(j in 1:length(pt.samples)){
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
  }
  dev.off()
}
