##########################################################################
# Figure 2D: Driver mutation comparison
##########################################################################
plot.genomic.drivers <- function(filename){
  if(!exists("muts") | !exists("driver.genes") | !exists("tcga.drivers")){
    stop("ERROR: please run loadWgsData before plot.genomic.drivers")
  }
  
  # Find mutation rates in CIS - all, prog and reg - for driver genes
  dgenes <- as.character(driver.genes[which(driver.genes %in% rownames(muts))])
  sel.p <- which(wgs.pheno$progression == 1)
  sel.r <- which(wgs.pheno$progression == 0)
  sel.tr <- which(wgs.pheno$name %in% c("PD21884c", "PD21885a", "PD21885c", "PD21904d", "PD21908a"))
  driver.plot <- data.frame(
    gene=dgenes,
    # tcga=driver.genes$tcga.pc[match(dgenes, driver.genes$Gene)],
    #tcga=100*as.numeric(tcga.snvs.rates[dgenes]),
    tcga=tcga.drivers$pc[match(dgenes, tcga.drivers$Symbol)],
    cis.all=apply(muts[dgenes, ], 1, function(x){
      100 * length(which(x == 1)) / length(x)
    }),
    cis.prog=apply(muts[dgenes, sel.p], 1, function(x){
      100 * length(which(x == 1)) / length(x)
    }),
    cis.reg=apply(muts[dgenes, sel.r], 1, function(x){
      100 * length(which(x == 1)) / length(x)
    }),
    cis.truereg=apply(muts[dgenes, sel.tr], 1, function(x){
      100 * length(which(x == 1)) / length(x)
    })
  )
  driver.plot <- driver.plot[order(driver.plot$cis.all, decreasing = T),]
  
  sel.plot <- match(c("TP53", "CDKN2A", "PTEN", "PIK3CA", "KEAP1", "KMT2B", "HLA-A", "NFE2L2", "NOTCH1", "RB1"), driver.plot$gene)
  # Barplot these data:
  plot.cols <- c("tcga", "cis.all", "cis.prog", "cis.reg")
  
  pdf(filename)
  barplot(
    as.matrix(t(driver.plot[sel.plot,plot.cols])), beside = T,
    legend=plot.cols,
    col=c("black", "orange", "red", "green")
  )
  dev.off()
}
