# Function to plot a gene over different 'dose' groups
plotGeneByGroup <- function(gene, legend.pos='topright'){
  d0 <- density(as.numeric(tcga.gdata.all[gene,which(tcga.gpheno.all$dose == 0)]))
  d1 <- density(as.numeric(tcga.gdata.all[gene,which(tcga.gpheno.all$dose == 1)]))
  d2 <- density(as.numeric(tcga.gdata.all[gene,which(tcga.gpheno.all$dose == 2)]))
  d3 <- density(as.numeric(tcga.gdata.all[gene,which(tcga.gpheno.all$dose == 3)]))
  #d4 <- density(as.numeric(tcga.gdata.all[gene,which(tcga.gpheno.all$dose == -1)]))
  plot(
    d0, 
    col="darkgreen",
    xlim=c(
      min(as.numeric(tcga.gdata.all[gene,])), max(as.numeric(tcga.gdata.all[gene,]))
    ),
    xlab="Gene Expression Raw Value",
    ylim=c(0,max(d0$y, d1$y, d2$y, d3$y)),
    main=gene
  )
  lines(d1, col='green')
  lines(d2, col='red')
  lines(d3, col='orange')
  #lines(d4, col='black')
  if(show.legends){
    legend(legend.pos, 
           c("TCGA.cont", "CIS.reg", "CIS.prog", "TCGA.LUSC"),
           col=c("darkgreen", "green", "red", "orange"), lty=1)
  }
  
}