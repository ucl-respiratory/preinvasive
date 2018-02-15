plotMethylationByGroup <- function(mvp, gene="", legend.pos="topright"){
  d0 <- density(as.numeric(tcga.mdata.all[mvp,which(tcga.mpheno.all$dose == 0 & tcga.mpheno.all$Sample_Group != "Control")]))
  d1 <- density(as.numeric(tcga.mdata.all[mvp,which(tcga.mpheno.all$dose == 1)]))
  d2 <- density(as.numeric(tcga.mdata.all[mvp,which(tcga.mpheno.all$dose == 2)]))
  d3 <- density(as.numeric(tcga.mdata.all[mvp,which(tcga.mpheno.all$dose == 3)]))
  plot(
    d0, 
    col="darkgreen",
    xlim=c(0,1
           #min(as.numeric(tcga.mdata.all[mvp,])), max(as.numeric(tcga.mdata.all[mvp,]))
    ),
    ylim=c(0,max(d0$y, d1$y, d2$y, d3$y)),
    main=paste(mvp, gene, sep=':'),
    xlab="Methylation beta value"
  )
  lines(d1, col='green')
  lines(d2, col='red')
  lines(d3, col='orange')
  if(show.legends){
    legend(legend.pos, 
           c("TCGA.cont", "CIS.reg", "CIS.prog", "TCGA.LUSC"),
           col=c("darkgreen", "green", "red", "orange"), lty=1)
  }
}