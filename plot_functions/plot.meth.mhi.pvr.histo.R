########################################################################################################
# Methylation Heterogeneity Index for PvR using Random Sampling
########################################################################################################
plot.meth.mhi.pvr.histo <- function(filename){
  
  # Plot for Prog vs Reg (with AUCs:)
  # As for fig 4F use pre-saved data for reproducibility
  pdf(filename)
  for(size in c(2000)){
    auc.mean <- mean(aucs.all$pvr[which(aucs.all$size == size)])
    auc.sd <- sd(aucs.all$pvr[which(aucs.all$size == size)])
    hist(aucs.all$pvr[which(aucs.all$size == size)], 
         main=paste("Prediction with random samples of", size, "probes"),
         xlab=paste("AUC for Progressive vs Regressive CIS. Mean ", signif(auc.mean, 2), "(CI ", signif(auc.mean-2*auc.sd, 2), "-", signif(auc.mean+2*auc.sd, 2),")", sep=""))
  }
  dev.off()
}