########################################################################################################
# Methylation Heterogeneity Index Analysis
#
# This file details the derivation of the MHI classified.
# This classifier simply counts probes with intermediate methylation, defined as lo < beta < hi
# In this file we determine optimum thresholds for lo and hi based on CIS data
# We also perform random sampling to demonstrate that randomly selected probes perform well as a predictor
#
# This file depends on the variable tcga.mdata.all - please run data_loaders/loadMethData.R if you have not already done so
########################################################################################################

if(!exists("tcga.mdata.all") | !exists("tcga.mpheno.all")){
  message("WARNING: tcga.mdata.all not found. Running data_loaders/loadMethData.R")
  source('data_loaders/loadMethData.R')
}

library(pROC)

# For plots:
cols <- c("darkgreen", "green", "red", "orange")

o <- order(tcga.mpheno.all$dose)
mdata.mhi  <- tcga.mdata.all[,o]
mpheno.mhi <- tcga.mpheno.all[o,]
# Exclude all Control samples
sel.control <- which(mpheno.mhi$Sample_Group == "Control")
mdata.mhi <- mdata.mhi[,-sel.control]
mpheno.mhi <- mpheno.mhi[-sel.control,]


# Define a function that calculates AUC for given probes and samples, and thresholds
mhiAucs <- function(probes=rownames(mdata.mhi), samples=1:dim(mdata.mhi)[2], thresh.up=0.9, thresh.low=0.1, doplot=F){
  rundata <- mdata.mhi[probes, samples]
  runpheno <- mpheno.mhi[samples,]
  mhi <- apply(rundata, 2, function(x){length(which(x > thresh.low & x < thresh.up))})  
  
  sel <- which(runpheno$dose %in% c(0,3))
  if(length(sel) > 0){
    auc.cvc.all2 <- auc(roc(predictor=mhi[sel], response=runpheno$progression[sel]))
  }else{
    auc.cvc.all2 <- NA
  }
  
  sel <- which(runpheno$dose %in% c(1,2))
  if(length(sel) > 0){
    auc.pvr.all2 <- auc(roc(predictor=mhi[sel], response=runpheno$progression[sel]))
  }else{
    auc.pvr.all2 <- NA
  }
  
  
  if(doplot){
    o <- order(runpheno$dose)
    plot(mhi[o], col=cols[runpheno$dose[o]+1], main=paste("Probes with", thresh.low, "< beta <", thresh.up), ylab="Number of probes", xlab="Sample")
    legend('bottomleft', c(
      paste("Cancer vs Control AUC=", signif(auc.cvc.all2, 3)),
      paste("Prog vs Reg AUC=", signif(auc.pvr.all2, 3))
    ))
  }
  return(c(auc.cvc.all2, auc.pvr.all2))
}


# Optimise thresholds based on discovery samples.
# Use all probes for this analysis.
# Note that CVC AUCs are NA because we only use CIS discovery samples here.
sel.d <- which(as.character(mpheno.mhi$name) %in% mpheno$Sample_Name[which(mpheno$Cohort == "D")])
sel.v <- which(as.character(mpheno.mhi$name) %in% mpheno$Sample_Name[which(mpheno$Cohort == "V")])
for(low in seq(from=0.07, to=0.4, by=0.01)){
  for(high in seq(from=0.75, to=0.95, by=0.01)){
    print(paste("Testing", low, "-", high))
    myauc <- mhiAucs(probes=rownames(mdata.mhi), samples=sel.d, doplot=F, thresh.up = high, thresh.low = low)
    df <- data.frame(low=low, high=high, cvc=myauc[[1]], pvr=myauc[[2]])
    if(low==0.07 & high==0.75){
      aucs <- df
    }else{
      aucs <- rbind(aucs, df)
    }
  }
}
save(aucs, file="./results/MHIoptimisation.RData")

# Choose the best thresholds
o <- order(aucs$pvr, decreasing = T)
row <- aucs[o,][1,]
thresh.low <- row$low
thresh.up  <- row$high

message(paste("Found optimised thresholds: Low", thresh.low, "High", thresh.up))

# Run on validation set + TCGA data
sel.vtcga <- c(sel.v, which(mpheno.mhi$source == "TCGA"))
mhiAucs(rownames(mdata.mhi), samples=sel.vtcga, doplot=T, thresh.up = thresh.up, thresh.low = thresh.low)

# Check for correlations with other variables
# Here we use mpheno rather than mpheno.mhi as this variable has full clinical data
sel.noctl <- which(mpheno$Sample_Group != "Control")
m.clin <- mpheno[sel.noctl,]
m.clin$mhi <- apply(mdata[,sel.noctl], 2, function(x){length(which(x > thresh.low & x < thresh.up))})
myLm <- glm(progression ~ mhi + Pack.years  + COPD + Gender + Age.at.specimen.collected + Previous.Lung.CA.........No..0..Yes..1., 
            data=m.clin, family = binomial)
summary(myLm)


# Using identified thresholds, use randomly sampled probes to demonstrate that this signature is robust with random sampling
# In the paper we use 10,000 simulations - for speed we use only 100 in this script. The user can increase n.runs to replicate this.
# We use CIS validation set and TCGA data for this analysis.
n.runs <- 100
samp.sizes <- c(100, 500, 1000, 2000, 4000)
for(samp.size in samp.sizes){
  for(i in 1:n.runs){
    myauc <- mhiAucs(sample(rownames(mdata.mhi), samp.size), samples=sel.vtcga, thresh.up=thresh.up, thresh.low = thresh.low)
    df <- data.frame(size=samp.size, cvc=myauc[[1]], pvr=myauc[[2]])
    if(samp.size == 100 & i == 1){
      aucs.samp <- df
    }else{
      aucs.samp <- rbind(aucs.samp, df)
    }
  }
}
save(aucs.samp, file="results/mhi_sampling_aucs.RData")

# Histograms of these AUC values are plotted in figure 4F and extended data figure 8