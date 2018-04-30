########################################################################################################
# Methylation Heterogeneity Index plotted by sample
########################################################################################################
plot.meth.mhi <- function(filename){
  
  if(!exists("mhi")){
    stop("ERROR: Please calculate MHI before plotting, see full_analysis.R")
  }
  
  pdf(filename)
  
  # AUC for cancer vs control:
  sel <- which(mpheno.mhi$dose %in% c(0,3))
  auc.cvc.all <- auc(roc(predictor=mhi[sel], response=mpheno.mhi$progression[sel]))
  # AUC for PvR - using validation set only
  sel <- which(mpheno.mhi$dose %in% c(1,2) & mpheno.mhi$name %in% mpheno.v$Sample_Name)
  auc.pvr.all <- auc(roc(predictor=mhi[sel], response=mpheno.mhi$progression[sel]))
  plot(mhi, col=tcga.cols[as.numeric(factor(mpheno.mhi$Sample_Group, levels=c("TCGA Control", "Regressive", "Progressive", "TCGA SqCC", "Control")))], 
       main=paste("Probes with", thresh.low, "< beta <", thresh.up), ylab="Number of probes", xlab="Sample")
  if(show.legends){
    legend('bottomleft', c(
      paste("Cancer vs Control AUC=", signif(auc.cvc.all, 3)),
      paste("Prog vs Reg AUC=", signif(auc.pvr.all, 3))
    ))
  }
  
  
  dev.off()
}