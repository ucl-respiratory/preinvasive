########################################################################################################
# Methylation PCAs with p-values
########################################################################################################
plot.meth.pcas <- function(filename){
  if(!exists("mdata")){
    stop("Error: Please run loadMethylData first")
  }
  
  pcfit <- prcomp(t(mdata))
  # Calculate p-values using multivariate ANOVA (all predictors are categorical)
  aovdata <- mpheno
  aovdata$pca <- pcfit$x[,1]
  myaov <- aov(pca ~ progression + smoking_group + COPD + Previous.Lung.CA + age_group + Gender + Slide, data=aovdata)
  pdf(filename)
  par(mfrow=c(1,1), xpd=T)
  inset <- c(0,-0.15)
  # Progression/Regression
  plot(pcfit$x[,1:2], col=c('blue', 'red', 'green')[as.numeric(factor(mpheno$Sample_Group, levels=c("Control", "Progressive", "Regressive")))],
       xlab=NA, ylab=NA,
       main=paste("Progression Status"))
  legend('bottom', col=c('red', 'green', 'blue'), 
         legend=c('Prog.','Reg.', 'Cont.'), pch=1, horiz=T, 
         bty="n", inset=inset)
  p <- summary(myaov)[[1]][["Pr(>F)"]][[1]]
  if(show.legends){
    legend('bottomright', paste("ANOVA p=", signif(p, 2), sep=""))
  }
  
  # Smoking
  plot(pcfit$x[,1:2], col=as.character(smoking_group_names[mpheno$smoking_group]),
       xlab=NA, ylab=NA,
       main=paste("Smoking (pack years)"))
  legend('bottom', col=as.character(smoking_group_names), 
         legend=names(smoking_group_names), pch=1, horiz=T, 
         bty="n", inset=inset)
  p <- summary(myaov)[[1]][["Pr(>F)"]][[2]]
  if(show.legends){
    legend('bottomright', paste("ANOVA p=", signif(p, 2), sep=""))
  }
  
  # COPD
  plot(pcfit$x[,1:2], col=c("blue", "red")[as.numeric(factor(mpheno$COPD, levels=c("NO", "YES")))],
       xlab=NA, ylab=NA,
       main=paste("COPD"))
  legend('bottom', col=c("blue", "red"), 
         legend=c("No", "Yes"), pch=1, horiz=T, 
         bty="n", inset=inset)
  p <- summary(myaov)[[1]][["Pr(>F)"]][[3]]
  if(show.legends){
    legend('bottomright', paste("ANOVA p=", signif(p, 2), sep=""))
  }
  
  # Lung cancer history
  plot(pcfit$x[,1:2], col=c("blue", "red")[as.numeric(as.character(mpheno$Previous.Lung.CA))+1],
       xlab=NA, ylab=NA,
       main=paste("Previous lung cancer history"))
  legend('bottom', col=c("blue", "red"), 
         legend=c("No", "Yes"), pch=1, horiz=T, 
         bty="n", inset=inset)
  p <- summary(myaov)[[1]][["Pr(>F)"]][[4]]
  if(show.legends){
    legend('bottomright', paste("ANOVA p=", signif(p, 2), sep=""))
  }
  
  # Age at Bronchoscopy
  plot(pcfit$x[,1:2], col=as.character(age_group_names[mpheno$age_group]),
       xlab=NA, ylab=NA,
       main=paste("Age at Bronchoscopy"))
  legend('bottom', col=as.character(age_group_names), 
         legend=names(age_group_names), pch=1, horiz=T, 
         bty="n", inset=inset)
  p <- summary(myaov)[[1]][["Pr(>F)"]][[5]]
  if(show.legends){
    legend('bottomright', paste("ANOVA p=", signif(p, 2), sep=""))
  }
  
  # Gender
  plot(pcfit$x[,1:2], col=c("red", "blue")[as.numeric(factor(mpheno$Gender))],
       xlab=NA, ylab=NA,
       main=paste("Gender"))
  legend('bottom', col=c("red", "blue"), 
         legend=c("Female", "Male"), pch=1, horiz=T, 
         bty="n", inset=inset)
  p <- summary(myaov)[[1]][["Pr(>F)"]][[6]]
  if(show.legends){
    legend('bottomright', paste("ANOVA p=", signif(p, 2), sep=""))
  }
  
  # Sentix ID
  plot(pcfit$x[,1:2], col=myPalette[as.numeric(factor(mpheno$Slide))],
       xlab=NA, ylab=NA,
       main=paste("Sentix ID"))
  legend('bottom', col=myPalette, 
         legend=1:length(levels(factor(mpheno$Slide))), pch=1, horiz=T, 
         bty="n", inset=inset)
  p <- summary(myaov)[[1]][["Pr(>F)"]][[7]]
  if(show.legends){
    legend('bottomright', paste("ANOVA p=", signif(p, 2), sep=""))
  }
  
  dev.off()
}
