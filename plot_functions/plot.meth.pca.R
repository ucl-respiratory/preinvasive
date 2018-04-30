##########################################################################
# Figure 3D: Methylation PCA
##########################################################################
plot.meth.pca <- function(filename){
  if(!exists("mdata")){
    stop("ERROR: Please run loadMethData")
  }
  
  pdf(filename)
  plotPCA(t(mdata), as.numeric(factor(mpheno$Sample_Group, levels=c("Regressive", "Progressive", "Control"))), title = "Methylation PCA")
  if(show.legends){
    legend('topleft', c('Regressive', 'Progressive', 'Control'), col=c('green', 'red', 'blue'), pch=1)
  }
  dev.off()
}
