##########################################################################
# Figure 3C: GXN PCA
##########################################################################
plot.gxn.pca <- function(filename){
  if(!exists("gdata")){
    stop("ERROR: Please run loadGxnData")
  }
  
  pdf(filename)
  plotPCA(t(gdata), gpheno$progression+1, title = "Gene Expression PCA")
  if(show.legends){
    legend('topleft', c('Regressive', 'Progressive'), col=c('green', 'red'), pch=1)
  }
  dev.off()
}
