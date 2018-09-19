# Plot prog vs reg comparisons alongside TCGA cancer vs control comparisons to show similarities
# Plot differentially methylated regions (DMRs) and copy number changes
# Use plot functions from the ggbio package
library(GenomicRanges)
library(ggbio)
data(ideoCyto, package = "biovizBase")

plot.overlaps <- function(filename1, filename2){
  
  # Create a genomic ranges object from our DMR data
  dn.dmrs <- GRanges(
    seqnames = dmrs$ProbeLassoDMR$seqnames,
    ranges = IRanges(
      start = dmrs$ProbeLassoDMR$start, end = dmrs$ProbeLassoDMR$end
    ),
    col = dmrs$ProbeLassoDMR$betaAv_Progressive - dmrs$ProbeLassoDMR$betaAv_Regressive,
    source = "CIS"
  )
  
  # Create a simliar object of TCGA DMR data
  dn.dmrs.tcga <- GRanges(
    seqnames = dmrs.tcga$ProbeLassoDMR$seqnames,
    ranges = IRanges(
      start = dmrs.tcga$ProbeLassoDMR$start, end = dmrs.tcga$ProbeLassoDMR$end
    ),
    col = dmrs.tcga$ProbeLassoDMR$betaAv_TCGA.SqCC - dmrs.tcga$ProbeLassoDMR$betaAv_TCGA.Control,
    source = "TCGA"
  )
  
  # Combine the two
  dn.dmrs <- c(dn.dmrs, dn.dmrs.tcga)
  dn.dmrs$levels <- as.numeric(factor(dn.dmrs$source, levels=c("TCGA", "CIS")))
  seqlengths(dn.dmrs) <- seqlengths(ideoCyto$hg19)[names(seqlengths(dn.dmrs))]
  
  # Plot as a karyogram, splitting TCGA and CIS data sets on the y-axis
  p.ylim <- autoplot(dn.dmrs, layout = "karyogram", aes(color=col, fill = col, 
                                                        ymin = ifelse(source == "CIS", 5.5, 1), 
                                                        ymax = ifelse(source == "CIS", 9, 4.5)
  ))
  # Use the same colour scale as for methylation heatmaps
  p.ylim + scale_colour_distiller(palette = 'RdYlBu')
  ggsave(filename1, scale=2)
  # pdf(filename1)
  # dmr.plot
  # dev.off()
  
  ##################################################################
  # Create an analagous plot for copy number
  ##################################################################
  cnas.segmented.mean.p <- cnas.segmented[,1:3]
  sel <- wgs.pheno$name[which(wgs.pheno$progression == 1)]
  cnas.segmented.mean.p$cn <- apply(cnas.segmented[,sel] / wgs.pheno$ploidy[match(sel, wgs.pheno$name)], 1, mean)
  
  cnas.segmented.mean.r <- cnas.segmented[,1:3]
  sel <- wgs.pheno$name[which(wgs.pheno$progression == 0)]
  cnas.segmented.mean.r$cn <- apply(cnas.segmented[,sel] / wgs.pheno$ploidy[match(sel, wgs.pheno$name)], 1, mean)
  
  dn.cnas.p <- GRanges(
    seqnames = paste0("chr", cnas.segmented.mean.p$chr),
    ranges = IRanges(
      start = cnas.segmented.mean.p$start, end = cnas.segmented.mean.p$end
    ),
    cn = cnas.segmented.mean.p$cn,
    source = "Prog"
  )
  dn.cnas.r <- GRanges(
    seqnames = paste0("chr", cnas.segmented.mean.r$chr),
    ranges = IRanges(
      start = cnas.segmented.mean.r$start, end = cnas.segmented.mean.r$end
    ),
    cn = cnas.segmented.mean.r$cn,
    source = "Reg"
  )
  sel <- which(tcga.cnas.segmented$chr %in% 1:22) # This step keeps factor names consistent with hg19
  dn.cnas.tcga <- GRanges(
    seqnames = paste0("chr", tcga.cnas.segmented.mean$chr[sel]),
    ranges = IRanges(
      start = tcga.cnas.segmented.mean$start[sel], end = tcga.cnas.segmented.mean$end[sel]
    ),
    cn = tcga.cnas.segmented.mean$cn[sel],
    source = "TCGA"
  )
  
  dn.cnas <- c(dn.cnas.p, dn.cnas.r, dn.cnas.tcga)
  # Remove sex chromosomes. 
  sel <- which(as.character(seqnames(dn.cnas)) %in% paste0("chr", 1:22))
  dn.cnas <- dn.cnas[sel]
  
  dn.cnas$levels <- as.numeric(factor(dn.cnas$source, levels=c("TCGA", "Prog", "Reg")))
  seqlengths(dn.cnas) <- seqlengths(ideoCyto$hg19)[names(seqlengths(dn.cnas))]
  
  # Specify colours explicitly to match ext. data fig 4
  dn.cnas$col <- "#E9EDF8"
  dn.cnas$col[which(dn.cnas$cn < 0.75)] <- "#90BEDA"
  dn.cnas$col[which(dn.cnas$cn < 0.5)] <- "darkblue"
  dn.cnas$col[which(dn.cnas$cn > 1.25)] <- "#EB7C64"
  dn.cnas$col[which(dn.cnas$cn > 2)] <- "#B81321"
  
  p.ylim <- autoplot(dn.cnas, layout = "karyogram", aes(color=col, fill = col, 
                                                        ymin = (levels - 1) * 10/3 + 0.5,
                                                        ymax = levels * 10 /3 - 0.5)
  )
  cols <- unique(dn.cnas$col)
  names(cols) <- cols
  
  # Do the plot 
  p.ylim + scale_color_manual(values=cols) + scale_fill_manual(values=cols)
  ggsave(filename2, scale=2)
}