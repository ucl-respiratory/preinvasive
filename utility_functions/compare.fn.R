# Mixed effects model for comparing two groups
library("lme4")
compare.fn <- function(modelinfo, dependent_variable, compared_variable, fixed_effects, random_effects, offset=NULL, strip.method="jitter", title=NULL) {
  
  # scale the numeric predictors 
  # first need to work out which are the numeric ones.
  nums <- sapply(modelinfo, is.numeric)
  fixedtoscale <- intersect(fixed_effects, names(nums[which(nums)]))
  fixednottoscale <- fixed_effects[!fixed_effects %in% names(nums[which(nums)])]
  
  # now do the scaling.
  scaleddat <- as.data.frame(cbind(modelinfo[,c(dependent_variable, compared_variable,random_effects, fixednottoscale)], scale(modelinfo[,c(fixedtoscale)])))
  scaleddat <- scaleddat[complete.cases(scaleddat),]
  
  # if the dependent variable is not an integer, round it.
  scaleddat[,dependent_variable] <- round(scaleddat[,dependent_variable], digits=0)
  
  # fit the full model
  fullform <- paste0(dependent_variable,  " ~ ", compared_variable, " + ", paste(fixed_effects, collapse=" + "), " + (1|", paste(random_effects, collapse=") + (1|"), ")" )
  full.model <- glmer(fullform, data=scaleddat, family="poisson", offset=offset)
  
  print("Full model")
  print(full.model)
  cat("\n")
  
  # fit the model removing the variable that is being assessed.
  reducedform <- paste0(dependent_variable,  " ~ ", paste(fixed_effects, collapse=" + "), " + (1|", paste(random_effects, collapse=") + (1|"), ")" )
  reduced.model <- glmer(reducedform, data=scaleddat, family="poisson", offset=offset)
  print("Reduced model")
  print(reduced.model)
  cat("\n")
  
  # compare the models using anova
  compare_anova <- anova(full.model, reduced.model)
  print("Comparing models")
  print(compare_anova)
  cat("\n")
  
  # plot the differences between the two groups
  # only plot if there is no offset. 
  if (is.null(offset)){
    par(mfrow=c(1,1))
    if(is.null(title)){
      title <- paste0(dependent_variable, " ~ ", compared_variable)
    }
    
    
    boxplot(modelinfo[,(colnames(modelinfo)==dependent_variable)] ~ modelinfo[,(colnames(modelinfo)==compared_variable)],
            add=F, notch=T, at=c(1.15,1.85), range=0, main=title, boxwex=1*1*0.5)
    stripchart(modelinfo[,(colnames(modelinfo)==dependent_variable)] ~ modelinfo[,(colnames(modelinfo)==compared_variable)],
               vertical=T, pch=16, col=c("green", "red"), ylab=dependant_variable,
               method=strip.method, at=c(1.15,1.85), add=T, cex=2)
    if(exists("show.legends") & show.legends){
      legend("top", legend=c(paste0('LRT p.val = ', round(compare_anova$`Pr(>Chisq)`[2], digits=3))), bty="n")
    }
    
    # Mark the 'query regressive' samples in orange
    sel.qr <- which(modelinfo$query.reg == 1)
    points(
      rep(0.75,length(sel.qr)),
      modelinfo[sel.qr, dependent_variable],
      col='orange', cex=2, pch=8
    )
    # Mark the other odd sample in blue (PD21908a)
    sel.pd21908a <- which(modelinfo$name == "PD21908a")
    points(
      rep(0.75,length(sel.pd21908a)),
      modelinfo[sel.pd21908a, dependent_variable],
      col='blue', cex=2, pch=8
    )
  }
}