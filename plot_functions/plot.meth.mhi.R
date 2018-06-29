########################################################################################################
# Methylation Heterogeneity Index plotted by sample
########################################################################################################
plot.meth.mhi <- function(filename){
  
  if(!exists("mhi")){
    stop("ERROR: Please calculate MHI before plotting, see full_analysis.R")
  }
  
  pdf(filename)
  
  # AUC for cancer vs control:
  sel.cvc <- which(mpheno.mhi$dose %in% c(0,3))
  auc.cvc.all <- auc(roc(predictor=mhi[sel.cvc], response=mpheno.mhi$progression[sel.cvc]))
  
  # AUC for PvR - using validation set only
  sel.pvr <- which(mpheno.mhi$dose %in% c(1,2) & mpheno.mhi$name %in% mpheno.v$Sample_Name)
  auc.pvr.all <- auc(roc(predictor=mhi[sel.pvr], response=mpheno.mhi$progression[sel.pvr]))
  plot(mhi, col=tcga.cols[as.numeric(factor(mpheno.mhi$Sample_Group, levels=c("TCGA Control", "Regressive", "Progressive", "TCGA SqCC", "Control")))], 
       main=paste("Probes with", thresh.low, "< beta <", thresh.up), ylab="Number of probes", xlab="Sample")
  
  # Plot ROC and PR curves
  plotRocs(predictor = mhi[sel.cvc], response=mpheno.mhi$progression[sel.cvc], title="C vs C")
  plotRocs(predictor = mhi[sel.pvr], response=mpheno.mhi$progression[sel.pvr], title="P vs R")
  
  if(show.legends){
    legend('bottomleft', c(
      paste("Cancer vs Control AUC=", signif(auc.cvc.all, 3)),
      paste("Prog vs Reg AUC=", signif(auc.pvr.all, 3))
    ))
  }
  
  
  
  
  dev.off()
}