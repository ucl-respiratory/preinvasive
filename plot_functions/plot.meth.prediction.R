plot.meth.prediction <- function(filename){
  
  if(!exists("meth.pamr.mycv")){
    stop("Predictive model code not yet run, see full_analysis.R")
  }
  
  pdf(filename)
  
  # Plot discovery set probabilities from cross-validated model
  meth.prob.prog <- meth.pamr.mycv$prob[,2,meth.threshold.id]
  plot(meth.prob.prog, col=c("green", "red", "blue")[factor(meth.pamr.traindata$group, levels=c("Regressive", "Progressive", "Control"))], ylab="Progression Score", main='Cross-validated prediction model - discovery set')
  abline(v=length(which(meth.pamr.traindata$y == meth.pamr.traindata$y[1]))+0.5, col='grey')
  
  if(show.legends){
    legend("topleft", paste("Features used:", dim(pamr.meth.features)[1]))
  }
  plotRocs(meth.prob.prog, meth.pamr.traindata$y)
  
  # Plot validation set
  meth.pamr.pred.v <- pamr.predict(meth.pamr.trainfit, newx=meth.pamr.testdata$x, type='posterior', threshold=meth.threshold)
  plot(meth.pamr.pred.v[,2], col=c("green", "red", "blue")[factor(meth.pamr.testdata$group, levels=c("Regressive", "Progressive", "Control"))], main="Validation Set", ylab="Progression Score", ylim=c(0,1))
  abline(v=length(which(meth.pamr.testdata$y == meth.pamr.testdata$y[1]))+0.5, col='grey')
  meth.roc.v <- roc(predictor=meth.pamr.pred.v[,2], response=meth.pamr.testdata$y)
  if(show.legends){
    legend('bottomright', legend=paste("AUC", signif(auc(meth.roc.v), 3), sep="="))
  }
  plotRocs(meth.pamr.pred.v[,2], meth.pamr.testdata$y)
  
  # Repeat validation using TCGA data
  meth.pamr.pred.t <- pamr.predict(meth.pamr.trainfit, newx=meth.pamr.tcgadata$x, type='posterior', threshold=meth.threshold)
  plot(meth.pamr.pred.t[,2], col=c("darkgreen", "orange")[as.numeric(as.character(meth.pamr.tcgadata$y)) + 1], main="TCGA Data", ylab="Progression Score", ylim=c(0,1))
  abline(v=length(which(meth.pamr.tcgadata$y == meth.pamr.tcgadata$y[1]))+0.5, col='grey')
  meth.roc.t <- roc(predictor=meth.pamr.pred.t[,2], response=meth.pamr.tcgadata$y)
  if(show.legends){
    legend('right', legend=paste("AUC", signif(auc(meth.roc.t), 3), sep="="))
  }
  plotRocs(meth.pamr.pred.t[,2], meth.pamr.tcgadata$y)
  
  
  dev.off()
}