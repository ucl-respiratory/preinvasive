##########################################################################
# Figure 2E: genome-wide copy number plot
##########################################################################
plot.genomic.cna.genomewide <- function(filename="results/figures/cna.genomewide.pdf"){
  
  if(!exists("cnas.segmented")){
    stop("ERROR: cnas.segmented not defined. Please run loadWgsData.")
  }
  
  pdf(filename)
  segs <- cnas.segmented[which(as.numeric(as.character(cnas.segmented$chr)) %in% 1:22),]
  
  # Sort segs by chr and position
  o <- order(as.numeric(as.character(segs$chr)), as.numeric(segs$start))
  segs <- segs[o,]
  
  seglength <- segs$end - segs$start
  # Round to the nearest 60000 and divide
  seglength <- round(seglength/60000)
  # Ignore small segments (<60,000bp)
  segs <- segs[-which(seglength == 0),]
  chrs <- as.numeric(as.character(segs$chr))
  segs <- segs[,4:dim(segs)[2]]
  seglength <- seglength[-which(seglength == 0)]
  
  segs.plot <- (matrix(rep(segs[1,], seglength[1]), nrow = length(segs[1,])))
  chr <- rep(chrs[1], seglength[1])
  for(i in 2:dim(segs)[1]){
    m <- (matrix(rep(segs[i,], seglength[i]), nrow = length(segs[i,])))
    segs.plot <- cbind(segs.plot, m)
    chr <- c(chr, rep(chrs[i], seglength[i]))
  }
  
  # Correct for ploidy:
  segs.plot <- apply(segs.plot, c(1,2), as.numeric)
  segs.plot <- segs.plot / ploidys
  
  # Cap outliers - we just want to show loss, normal, gain here:
  segs.plot <- apply(segs.plot, c(1,2), function(x){
    x <- as.numeric(x)
    if(x > 3){
      return(3)
    }else{
      return(x)
    }
  })
  
  sel <- match(colnames(cnas.segmented)[4:dim(cnas.segmented)[2]], wgs.pheno$name)
  annot <- data.frame(
    row.names = wgs.pheno[sel,]$name,
    ploidy = as.factor(ploidys[wgs.pheno[sel,]$name]),
    patient=as.factor(make.names(wgs.pheno[sel,]$Patient)),
    status=as.factor(c("Regressive", "Progressive")[wgs.pheno[sel,]$progression + 1])
  )
  rownames(segs.plot) <- colnames(cnas.segmented)[4:dim(cnas.segmented)[2]]
  # Sort by progression status:
  o <- order(as.numeric(annot$status))
  # Annotate the columns with chromosome number
  annot_cols <- data.frame(
    chr=as.factor(chr)
  )
  colnames(segs.plot) <- make.names(1:dim(segs.plot)[2])
  rownames(annot_cols) <- colnames(segs.plot)
  g.cols <- as.character(pt_cols[1:length(unique(levels(annot$patient)))])
  names(g.cols) <- levels(annot$patient)
  annot_colors <- list(
    status=c(Progressive="red", Regressive="green"),
    ploidy=c("2"="darkgreen", "3"="blue","4"="red","6"="yellow"),
    patient=g.cols
  )
  annot_colors$chr = pt_cols[1:length(levels(factor(chr)))]
  names(annot_colors$chr) <- levels(factor(chr))
  breaks <- c(-0.01, 0.01, 0.99,1.01,1.99, 4)
  break_labels <- c("0 (deletion)", "<1 (loss)", "1", "1-2 (gain)", "2+ (amplification)", "")
  # Remove patient from annot
  annot$patient <- NULL
  pheatmap(
    segs.plot[o,],
    color=c(cols.cn, 'firebrick'),
    cluster_rows=F, cluster_cols = F,
    annotation_row = annot,
    annotation_col = annot_cols,
    annotation_colors = annot_colors,
    legend=F,
    breaks=breaks,
    legend_breaks=breaks,
    legend_labels=break_labels,
    border_color=NA,
    show_rownames=T, show_colnames=F
    
  )
  dev.off()
}
