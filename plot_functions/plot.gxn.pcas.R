########################################################################################################
# Gene Expression PCAs with p-values
########################################################################################################
plot.gxn.pcas <- function(filename){
  
  if(!exists("gdata")){
    stop("GXN data not loaded. Please run loadGxnData first.")
  }
  
  pcfit <- prcomp(t(gdata))
  aovdata <- gpheno
  aovdata$pca <- pcfit$x[,1]
  myaov <- aov(pca ~ progression + smoking_group + COPD + Prev.History.of.LC + age_group + Gender, data=aovdata)
  pdf(filename)
  par(mfrow=c(1,1), xpd=T)
  inset <- c(0,-0.15)
  
  # Progression/Regression
  plot(pcfit$x[,1:2], col=c('red', 'green')[as.numeric(gpheno$progression) + 1],
       xlab=NA, ylab=NA,
       main=paste("Progression Status"))
  legend('bottom', col=c('red', 'green'), 
         legend=c('Progressive','Regressive'), pch=1, horiz=T, 
         bty="n", inset=inset)
  p <- summary(myaov)[[1]][["Pr(>F)"]][[1]]
  if(show.legends){
    legend('bottomright', paste("ANOVA p=", signif(p, 2), sep=""))
  }
  
  # Smoking
  plot(pcfit$x[,1:2], col=as.character(smoking_group_names[gpheno$smoking_group]),
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
  plot(pcfit$x[,1:2], col=c("blue", "red")[as.numeric(gpheno$COPD)],
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
  plot(pcfit$x[,1:2], col=c("blue", "red")[as.numeric(gpheno$Prev.History.of.LC)+1],
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
  plot(pcfit$x[,1:2], col=as.character(age_group_names[gpheno$age_group]),
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
  plot(pcfit$x[,1:2], col=c("red", "blue")[as.numeric(factor(gpheno$Gender))],
       xlab=NA, ylab=NA,
       main=paste("Gender"))
  legend('bottom', col=c("red", "blue"), 
         legend=c("Female", "Male"), pch=1, horiz=T, 
         bty="n", inset=inset)
  p <- summary(myaov)[[1]][["Pr(>F)"]][[6]]
  if(show.legends){
    legend('bottomright', paste("ANOVA p=", signif(p, 2), sep=""))
  }
  
  dev.off()
}