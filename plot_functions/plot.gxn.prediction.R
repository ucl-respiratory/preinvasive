plot.gxn.prediction <- function(filename){
  
  if(!exists("gxn.pamr.mycv")){
    stop("Predictive model code not yet run, see full_analysis.R")
  }
  
  pdf(filename)
  
  # Plot discovery set probabilities from cross-validated model
  gxn.prob.prog <- gxn.pamr.mycv$prob[,2,gxn.threshold.id]
  plot(gxn.prob.prog, col=c("green", "red")[gxn.pamr.traindata$y + 1], ylab="Progression Score", main='Cross-validated prediction model - discovery set')
  abline(v=length(which(gxn.pamr.traindata$y == gxn.pamr.traindata$y[1]))+0.5, col='grey')
  
  if(show.legends){
    legend("topleft", paste("Features used:", dim(pamr.gxn.features)[1]))
  }
  plotRocs(gxn.prob.prog, gxn.pamr.traindata$y)
  
  # Plot validation set
  gxn.pamr.pred.v <- pamr.predict(gxn.pamr.trainfit, newx=gxn.pamr.testdata$x, type='posterior', threshold=gxn.threshold)
  plot(gxn.pamr.pred.v[,2], col=c("green", "red")[as.numeric(as.character(gxn.pamr.testdata$y)) + 1], main="Validation Set", ylab="Progression Score", ylim=c(0,1))
  abline(v=length(which(gxn.pamr.testdata$y == gxn.pamr.testdata$y[1]))+0.5, col='grey')
  gxn.roc.v <- roc(predictor=gxn.pamr.pred.v[,2], response=gxn.pamr.testdata$y)
  if(show.legends){
    legend('bottomright', legend=paste("AUC", signif(auc(gxn.roc.v), 3), sep="="))
  }
  plotRocs(gxn.pamr.pred.v[,2], gxn.pamr.testdata$y)
  
  # Repeat validation using TCGA data
  gxn.pamr.pred.t <- pamr.predict(gxn.pamr.trainfit, newx=gxn.pamr.tcgadata$x, type='posterior', threshold=gxn.threshold)
  plot(gxn.pamr.pred.t[,2], col=c("darkgreen", "orange")[as.numeric(as.character(gxn.pamr.tcgadata$y)) + 1], main="TCGA Data", ylab="Progression Score", ylim=c(0,1))
  abline(v=length(which(gxn.pamr.tcgadata$y == gxn.pamr.tcgadata$y[1]))+0.5, col='grey')
  gxn.roc.t <- roc(predictor=gxn.pamr.pred.t[,2], response=gxn.pamr.tcgadata$y)
  if(show.legends){
    legend('right', legend=paste("AUC", signif(auc(gxn.roc.t), 3), sep="="))
  }
  plotRocs(gxn.pamr.pred.t[,2], gxn.pamr.tcgadata$y)
  
  dev.off()
}