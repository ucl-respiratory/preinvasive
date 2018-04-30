########################################################################################################
# Methylation Distributions
########################################################################################################
plot.meth.distribution <- function(filename){
  # Density plot of MVPs in cancer vs control:
  pdf(filename)
  
  plot(density(apply(tcga.mdata.all[,which(tcga.mpheno.all$Sample_Group == "TCGA Control")], 1, mean)), col='darkgreen',
       main="Distribution of beta values across all probes", xlab="Beta value")
  lines(density(apply(tcga.mdata.all[,which(tcga.mpheno.all$Sample_Group == "Regressive")], 1, mean)), col='green')
  lines(density(apply(tcga.mdata.all[,which(tcga.mpheno.all$Sample_Group == "Progressive")], 1, mean)), col='red')
  lines(density(apply(tcga.mdata.all[,which(tcga.mpheno.all$Sample_Group == "TCGA SqCC")], 1, mean)), col='orange')
  if(show.legends){
    legend('topright', c("TCGA Control", "CIS Regressive", "CIS Progressive", "TCGA Cancer"), col=c('darkgreen', 'green', 'red', 'orange'), lty=1)
  }
  dev.off()
}