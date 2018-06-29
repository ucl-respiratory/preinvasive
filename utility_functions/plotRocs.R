# Simple function to plot ROC and precision-recall plots based on our predictors
# Uses the PRROC package
library(PRROC)
plotRocs <- function(predictor, response, title=""){

  roc<-roc.curve(
    scores.class0 = predictor[which(response == 1)], 
    scores.class1 = predictor[which(response == 0)],
    curve=T
  )
  plot(roc, color=F, main=paste(title, "ROC curve"))
  pr<-pr.curve(
    scores.class0 = predictor[which(response == 1)], 
    scores.class1 = predictor[which(response == 0)],
    curve=T
  )
  plot(pr, color=F, main=paste(title, "PR curve"))
}