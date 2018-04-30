plot.mcna.prediction <- function(filename){
  
  if(!exists("cna.pamr.mycv")){
    stop("Predictive model code not yet run, see full_analysis.R")
  }
  
  pdf(filename)
  
  cna.prob.prog <- cna.pamr.mycv$prob[,2,cna.threshold.id]
  
  plot(cna.prob.prog, col=c("green", "red")[cna.pamr.traindata$y + 1], ylab="Progression Score", main='Cross-validated prediction model')
  abline(v=length(which(cna.pamr.traindata$y == cna.pamr.traindata$y[1]))+0.5, col='grey')
  
  if(show.legends){
    legend("topleft", paste("Features used:", dim(pamr.cna.features)[1]))
  }
  plotRocs(cna.prob.prog, cna.pamr.traindata$y)
  
  # Plot validation with Van Boerdonk data
  cna.pamr.pred.v <- pamr.predict(cna.pamr.trainfit, newx=cna.pamr.testdata$x, type='posterior', threshold=cna.threshold)
  plot(cna.pamr.pred.v[,2], col=c("green", "red")[as.numeric(as.character(cna.pamr.testdata$y)) + 1], main="CNA Validation with Van Boerdonk data", ylab="Progression Score", ylim=c(0,1))
  abline(v=length(which(cna.pamr.testdata$y == cna.pamr.testdata$y[1]))+0.5, col='grey')
  gxn.roc.v <- roc(predictor=cna.pamr.pred.v[,2], response=cna.pamr.testdata$y)
  if(show.legends){
    legend('topleft', legend=paste("AUC", signif(auc(gxn.roc.v), 3), sep="="))
  }
  plotRocs(cna.pamr.pred.v[,2], cna.pamr.testdata$y)
  
  # Plot validation with TCGA data
  cna.pamr.pred.t <- pamr.predict(cna.pamr.trainfit, newx=cna.pamr.testdata2$x, type='posterior', threshold=cna.threshold)
  plot(cna.pamr.pred.t[,2], col=c("darkgreen", "orange")[as.numeric(as.character(cna.pamr.testdata2$y)) + 1], main="CNA Validation with TCGA data", ylab="Progression Score", ylim=c(0,1))
  abline(v=length(which(cna.pamr.testdata2$y == cna.pamr.testdata2$y[1]))+0.5, col='grey')
  gxn.roc.t <- roc(predictor=cna.pamr.pred.t[,2], response=cna.pamr.testdata2$y)
  if(show.legends){
    legend('topleft', legend=paste("AUC", signif(auc(gxn.roc.t), 3), sep="="))
  }
  plotRocs(cna.pamr.pred.t[,2], cna.pamr.testdata2$y)
  
  dev.off()
  
}