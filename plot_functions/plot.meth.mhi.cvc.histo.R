########################################################################################################
# Histogram of MHI AUC values using random probe samples for cancer vs control
# Based on data generated in mhi_analysis.R
########################################################################################################
plot.meth.mhi.cvc.histo <- function(filename){
  # Plot MHI for cancer vs control (with AUCs:)
  pdf(filename)
  for(size in c(2000)){
    auc.mean <- mean(aucs.all$cvc[which(aucs.all$size == size)])
    auc.sd <- sd(aucs.all$cvc[which(aucs.all$size == size)])
    hist(aucs.all$cvc[which(aucs.all$size == size)], 
         main=paste("Prediction with random samples of", size, "probes"),
         xlab=paste("AUC for Cancer vs Control. Mean ", signif(auc.mean, 2), "(CI ", signif(auc.mean-2*auc.sd, 2), "-", signif(auc.mean+2*auc.sd, 2),")", sep=""))
  }
  dev.off()
}