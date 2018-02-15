##########################################################################
# Setup
#
# Load gene expression and methylation data from NCBI GEO
##########################################################################

# Optionally change data_dir to set where data is downloaded to
# This variable is required before sourcing the data_loader scripts
data_cache <- "./data/"
results_dir <- "./results/"
for(dir in c(data_cache, results_dir)){
  dir.create(dir, recursive = T, showWarnings = F)
}

# Define GEO datasets used in this study:
geo.gxn.d <- "GSE94611" # GXN Discovery
geo.gxn.v <- "GSE108082" # GXN Validation
geo.meth  <- "GSE108123"

# Load data
# These functions load CIS and TCGA data, and make many variables available
source('data_loaders/loadGeneData.R')
source('data_loaders/loadMethData.R')

# Load CIN genes - both the CIN70 signature, and CIN70 with cell-cycle genes removed
load('resources/cin_genes.RData')

# Load libraries in a specific order to avoid overwriting functions
library(gdata)
library(ChAMP)
library(ggplot2)
library(ggsignif)
library(pheatmap)
library(pROC)
library(RColorBrewer)
library(limma)
library(stringr)
library(WriteXLS)
library(pamr)

# Load utility functions
for(file in list.files("./utility_functions/", full.names = T)){
  source(file)
}

# Define colour palettes for heatmaps: hmcol is green/red, hmcol2 is yellow/blue
hmcol <- colorRampPalette(c("Green","Black","Red"))(256)
hmcol2 <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(256)
# Colours for differentiating patients:
pt_cols <- c("#ff0000", "#cc0000", "#594343", "#7f4840", "#4c2213", "#f2c6b6", "#f26100", "#7f3300", "#f2aa79", "#593a16", "#8c7c69", "#bf8000", "#665200", "#d9c36c", "#e5d600", "#4a592d", "#e1ffbf", "#65b359", "#00cc1b", "#60806c", "#00331b", "#26332d", "#79f2ca", "#008c83", "#005c73", "#40d9ff", "#002233", "#0088ff", "#003059", "#668fcc", "#000e66", "#4059ff", "#5e53a6", "#3c394d", "#4700b3", "#290033", "#e63df2", "#967399", "#8c2377", "#ffbff2", "#331a27", "#660029", "#cc3370", "#ff4059", "#e6acb4")
# Copy Number Colours:
cols.cn <- c("#2b69ca", "#68aeff", "#ffffff", "#ff5468", "#ff0825")
# Colours for PCA plots
myPalette <- c("green", "red", "blue", "yellow", "black", "magenta", "cyan", "orange")
# Colours for TCGA plots
tcga.cols <- c('darkgreen', 'green', 'red', 'orange', 'blue')

# Define colours used for clinical variables
smoking_group_names <-c("<20"="green", "20-39"="yellow", "40-69"="orange", "70+"="red")
age_group_names <- c("<50"="blue","50-59"="green", "60-69"="orange", "70+"="red")

# Choose whether to include legends - turned off for production plots
show.legends <- T

##########################################################################
# Start of analysis
##########################################################################

##########################################################################
# Differential expression analysis of GXN and methylation
#
# Uses limma with FDR cutoff set to 1 (significant genes are identified downstream)
# For methylation, we use the ChAMP package (which is itself built on limma)
# Outputs are stored in the gdiff and mdiff variables
##########################################################################
source('utility_functions/limmaCompare.R')
# Gene Expression
gdiff <- limmaCompare(data=gdata.d, pheno=gpheno.d, fdr_limit = 1)

# Methylation
myDMP <- champ.DMP(
  beta=mdata.d,
  pheno=mpheno.d$Sample_Group,
  adjPVal = 1
)
mdiff <- myDMP$Progressive_to_Regressive

# Methylation-derived copy number (using logR aggregated by cytogenetic band)
cdiff <- limmaCompare(cnas.band, cnas.pheno, fdr_limit = 0.01)

# Repeat differential analysis using a continuous variable of 'dose'
# This is defined as 0=TCGA control, 1=regressive CIS, 2=progressive CIS, 3=TCGA cancer
AAprog <- as.numeric(tcga.gpheno.all$dose)
design <- model.matrix(~AAprog)
fit <- lmFit(tcga.gdata.all, design)
fit2 <- eBayes(fit)
p <- fit2$p.value[,"AAprog"]
fdr <- p.adjust(p, method="BH")
fc <- logratio2foldchange(fit2$coef[,"AAprog"])
t <- fit2$t[,"AAprog"]
uvv <- data.frame(
  row.names=rownames(tcga.gdata.all),
  fc=fc,
  fdr=fdr,
  t=t
)
o <- order(abs(t), decreasing = T)
gdiff.dose <- uvv[o,]

# Repeat for methylation
mdiff.dose <- champ.DMP(
  beta=tcga.mdata.all,
  pheno=as.numeric(tcga.mpheno.all$dose)
)
mdiff.dose <- mdiff.dose$NumericVariable


##########################################################################
# Start of plots
##########################################################################

##########################################################################
# Figure 1: study design, not generated with R
# Figure 2: based on genomic data, not included here
# Figure 3: heatmaps and PCAs, included
##########################################################################

##########################################################################
# Figure 3. Altered methylation and gene expression in progressive lung carcinoma-in-situ (CIS) lesions
##########################################################################
# Figure 3A: GXN heatmap
##########################################################################
sig_genes <- gdata.d[rownames(gdiff[which(gdiff$fdr < 0.01),]),]
g.annot <- data.frame(
  pack.years=gpheno.d$smoking_group,
  age.group=gpheno.d$age_group,
  gender=gpheno.d$Gender,
  COPD=gpheno.d$COPD,
  status=c("Regressive", "Progressive")[gpheno.d$progression+1]
)
g.annot_colors <- list(
  status=c(Progressive="red", Regressive="green"),
  pack.years=smoking_group_names,
  age.group=age_group_names,
  gender=c("F"="pink", "M"="blue"),
  COPD=c("N"="green", "Y"="red")
)
rownames(g.annot) <- colnames(sig_genes)

pdf(paste(results_dir, 'Fig3A_gxn_heatmap.pdf', sep=""))
pheatmap(removeOutliers(sig_genes), cluster_rows=T, cluster_cols=T, scale="row", main=paste("Gene Expression (",dim(sig_genes)[1]," genes)", sep=""),
         annotation_col=g.annot, treeheight_row=0, treeheight_col=0, show_rownames=F, show_colnames=F,
         annotation_colors=g.annot_colors,
         color=hmcol, legend=F)
dev.off()

##########################################################################
# Figure 3B: Methylation heatmap
##########################################################################
mdata.sig <- mdata[rownames(mdiff)[1:1000],]

m.annot <- data.frame(
  pack.years=mpheno$smoking_group,
  age.group=mpheno$age_group,
  gender=mpheno$Gender,
  COPD=substr(mpheno$COPD, 1,1),
  status=mpheno$Sample_Group
)
m.annot_colors <- list(
  status=c(Progressive="red", Regressive="green", Control="blue"),
  pack.years=smoking_group_names,
  age.group=age_group_names,
  gender=c("F"="pink", "M"="blue"),
  COPD=c("N"="green", "Y"="red")
)
rownames(m.annot) <- colnames(mdata.sig)

pdf(paste(results_dir, 'Fig3B_meth_heatmap.pdf', sep=""))
pheatmap(removeOutliers(mdata.sig), cluster_rows=T, cluster_cols=T, scale="row", main=paste("Methylation (Top ",dim(mdata.sig)[1]," MVPs)", sep=""),
         annotation_col=m.annot, treeheight_row=0, treeheight_col=0, show_rownames=F, show_colnames=F,
         annotation_colors=m.annot_colors,
         color=hmcol2, legend=F)
dev.off()

##########################################################################
# Figure 3C: GXN PCA
##########################################################################
pdf(paste(results_dir, 'Fig3C_gxn_pca.pdf', sep=""))
plotPCA(t(gdata), gpheno$progression+1, title = "Gene Expression PCA")
if(show.legends){
  legend('topleft', c('Regressive', 'Progressive'), col=c('green', 'red'), pch=1)
}
dev.off()

##########################################################################
# Figure 3D: Methylation PCA
##########################################################################
pdf(paste(results_dir, 'Fig3D_meth_pca.pdf', sep=""))
plotPCA(t(mdata), as.numeric(factor(mpheno$Sample_Group, levels=c("Regressive", "Progressive", "Control"))), title = "Methylation PCA")
if(show.legends){
  legend('topleft', c('Regressive', 'Progressive', 'Control'), col=c('green', 'red', 'blue'), pch=1)
}
dev.off()


########################################################################################################
# Figure 4A-C - GXN Predictive analyses
########################################################################################################

# GXN PAM prediction
set.seed(2)

# Use differentially expressed genes which are present in TCGA data
genes.shared <- intersect(rownames(tcga.gdata.all), rownames(gdiff)[which(gdiff$fdr < 0.01)])

pamr.pheno <- tcga.gpheno.all
pamr.pheno$train <- 0
pamr.pheno$train[which(pamr.pheno$name %in% as.character(gpheno$name)[which(gpheno$training == 1)])] <- 1

o <- order(pamr.pheno$dose)
pamr.pheno <- pamr.pheno[o,]
pamr.gdata <- tcga.gdata.all[,o]

sel <- which(pamr.pheno$train == 1 & pamr.pheno$source == "Surveillance")
gxn.pamr.traindata <- list(
  x=as.matrix(pamr.gdata[genes.shared,sel]),
  y=pamr.pheno$progression[sel],
  geneid = genes.shared
)
sel <- which(pamr.pheno$train == 0 & pamr.pheno$source == "Surveillance")
gxn.pamr.testdata <- list(
  x=as.matrix(pamr.gdata[genes.shared,sel]),
  y=pamr.pheno$progression[sel],
  geneid = genes.shared
)
sel <- which(pamr.pheno$source == "TCGA")
gxn.pamr.tcgadata <- list(
  x=as.matrix(pamr.gdata[genes.shared,sel]),
  y=pamr.pheno$progression[sel],
  geneid = genes.shared
)
gxn.pamr.trainfit <- pamr.train(gxn.pamr.traindata)
gxn.pamr.mycv <- pamr.cv(gxn.pamr.trainfit, gxn.pamr.traindata)
# Manually choose threshold = 2.5 from experimentation:
threshold=2.5

pdf(paste(results_dir, 'Fig4A-C_gxn_prediction.pdf', sep=""))
gxn.pamr.pred.d <- pamr.predict(gxn.pamr.trainfit, newx=gxn.pamr.traindata$x, type='posterior', threshold=threshold)
plot(gxn.pamr.pred.d[,2], col=c("green", "red")[as.numeric(as.character(gxn.pamr.traindata$y)) + 1], main="Discovery Set", ylab="Progression Score", ylim=c(0,1))
abline(v=length(which(gxn.pamr.traindata$y == gxn.pamr.traindata$y[1]))+0.5, col='grey')

features <- pamr.listgenes(gxn.pamr.trainfit, data=gxn.pamr.traindata, threshold=threshold)
if(show.legends){
  legend("topleft", paste("Features used:", dim(features)[1]))
}

gxn.pamr.pred.v <- pamr.predict(gxn.pamr.trainfit, newx=gxn.pamr.testdata$x, type='posterior', threshold=threshold)
plot(gxn.pamr.pred.v[,2], col=c("green", "red")[as.numeric(as.character(gxn.pamr.testdata$y)) + 1], main="Validation Set", ylab="Progression Score", ylim=c(0,1))
abline(v=length(which(gxn.pamr.testdata$y == gxn.pamr.testdata$y[1]))+0.5, col='grey')
gxn.roc.v <- roc(predictor=gxn.pamr.pred.v[,2], response=gxn.pamr.testdata$y)
if(show.legends){
  legend('bottomright', legend=paste("AUC", signif(auc(gxn.roc.v), 3), sep="="))
}

gxn.pamr.pred.t <- pamr.predict(gxn.pamr.trainfit, newx=gxn.pamr.tcgadata$x, type='posterior', threshold=threshold)
plot(gxn.pamr.pred.t[,2], col=c("darkgreen", "orange")[as.numeric(as.character(gxn.pamr.tcgadata$y)) + 1], main="TCGA Data", ylab="Progression Score", ylim=c(0,1))
abline(v=length(which(gxn.pamr.tcgadata$y == gxn.pamr.tcgadata$y[1]))+0.5, col='grey')
gxn.roc.t <- roc(predictor=gxn.pamr.pred.t[,2], response=gxn.pamr.tcgadata$y)
if(show.legends){
  legend('right', legend=paste("AUC", signif(auc(gxn.roc.t), 3), sep="="))
}
dev.off()


########################################################################################################
# Figure 4D - Methylation Distributions
########################################################################################################

# Density plot of MVPs in cancer vs control:
pdf(paste(results_dir, 'Fig4D_mvp_distribution.pdf', sep=""))

plot(density(apply(tcga.mdata.all[,which(tcga.mpheno.all$Sample_Group == "TCGA Control")], 1, mean)), col='darkgreen',
     main="Distribution of beta values across all probes", xlab="Beta value")
lines(density(apply(tcga.mdata.all[,which(tcga.mpheno.all$Sample_Group == "Regressive")], 1, mean)), col='green')
lines(density(apply(tcga.mdata.all[,which(tcga.mpheno.all$Sample_Group == "Progressive")], 1, mean)), col='red')
lines(density(apply(tcga.mdata.all[,which(tcga.mpheno.all$Sample_Group == "TCGA SqCC")], 1, mean)), col='orange')
if(show.legends){
  legend('topright', c("TCGA Control", "CIS Regressive", "CIS Progressive", "TCGA Cancer"), col=c('darkgreen', 'green', 'red', 'orange'), lty=1)
}
dev.off()

########################################################################################################
# Figure 4E - Methylation Heterogeneity Index
########################################################################################################

o <- order(tcga.mpheno.all$dose)
mdata.mhi <- tcga.mdata.all[,o]
mpheno.mhi <- tcga.mpheno.all[o,]
# Remove Controls
sel.control <- which(mpheno.mhi$Sample_Group == "Control")
if(length(sel.control) > 0){
  mdata.mhi <- mdata.mhi[,-sel.control]
  mpheno.mhi <- mpheno.mhi[-sel.control,]
}
# Thresholds are defined on our discovery set - see mhi_analysis.R
thresh.up <- 0.88
thresh.low <- 0.26

pdf(paste(results_dir, 'Fig4E_mvp_probe_counts.pdf', sep=""))

mhi <- apply(mdata.mhi, 2, function(x){length(which(x > thresh.low & x < thresh.up))})  
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


########################################################################################################
# Figure 4F - Methylation Heterogeneity Index based on 2000 probe samples
########################################################################################################
# To plot consistently we store data based on 10000 sample runs.
# Code to generate this is in mhi_analysis.R
load("resources/randomMvpV2AUCs10k.RData")

# Plot for cancer vs control (with AUCs:)
pdf(paste(results_dir, 'Fig4F_mhi_sample_histograms.pdf', sep=""))
for(size in c(2000)){
  auc.mean <- mean(aucs.all$cvc[which(aucs.all$size == size)])
  auc.sd <- sd(aucs.all$cvc[which(aucs.all$size == size)])
  hist(aucs.all$cvc[which(aucs.all$size == size)], 
       main=paste("Prediction with random samples of", size, "probes"),
       xlab=paste("AUC for Cancer vs Control. Mean ", signif(auc.mean, 2), "(CI ", signif(auc.mean-2*auc.sd, 2), "-", signif(auc.mean+2*auc.sd, 2),")", sep=""))
}
dev.off()

# Plot for Prog vs Reg (with AUCs:)
pdf(paste(results_dir, 'Sup_figure_8_mhi_sample_histograms_pvr.pdf', sep=""))
for(size in c(2000)){
  auc.mean <- mean(aucs.all$pvr[which(aucs.all$size == size)])
  auc.sd <- sd(aucs.all$pvr[which(aucs.all$size == size)])
  hist(aucs.all$pvr[which(aucs.all$size == size)], 
       main=paste("Prediction with random samples of", size, "probes"),
       xlab=paste("AUC for Progressive vs Regressive CIS. Mean ", signif(auc.mean, 2), "(CI ", signif(auc.mean-2*auc.sd, 2), "-", signif(auc.mean+2*auc.sd, 2),")", sep=""))
}
dev.off()



########################################################################################################
# Figure 5A - GXN Pathway analysis
# Figure plotted in Prism. Data output to Excel file.
########################################################################################################

sel <- which(tcga.gpheno.all$dose == 0 | tcga.gpheno.all$dose == 1)
gage_0_1 <- gage_analysis(tcga.gdata.all[,sel], tcga.gpheno.all[sel,]$dose, source='msig.kegg')
sel <- which(tcga.gpheno.all$dose == 1 | tcga.gpheno.all$dose == 2)
gage_1_2 <- gage_analysis(tcga.gdata.all[,sel], tcga.gpheno.all[sel,]$dose-1, source='msig.kegg')
sel <- which(tcga.gpheno.all$dose == 2 | tcga.gpheno.all$dose == 3)
gage_2_3 <- gage_analysis(tcga.gdata.all[,sel], tcga.gpheno.all[sel,]$dose-2, source='msig.kegg')
# Save
cont.vs.reg.up <- data.frame(gage_0_1$greater)
cont.vs.reg.down <- data.frame(gage_0_1$less)
reg.vs.prog.up <- data.frame(gage_1_2$greater)
reg.vs.prog.down <- data.frame(gage_1_2$less)
prog.vs.ca.up <- data.frame(gage_2_3$greater)
prog.vs.ca.down <- data.frame(gage_2_3$less)

# Make summary sheets 
summary.up <- merge(cont.vs.reg.up, reg.vs.prog.up,by=0)
rownames(summary.up) <- summary.up$Row.names
summary.up$Row.names <- NULL
summary.up <- merge(summary.up, prog.vs.ca.up, by=0)
rownames(summary.up) <- summary.up$Row.names
summary.up <- summary.up[,c("q.val.x", "q.val.y", "q.val")]
colnames(summary.up) <- c("q.val.cont.vs.reg", "q.val.reg.vs.prog", "q.val.prog.vs.ca")
o <- order(apply(summary.up, 1, sum))
summary.up <- summary.up[o,]

summary.down <- merge(cont.vs.reg.down, reg.vs.prog.down,by=0)
rownames(summary.down) <- summary.down$Row.names
summary.down$Row.names <- NULL
summary.down <- merge(summary.down, prog.vs.ca.down, by=0)
rownames(summary.down) <- summary.down$Row.names
summary.down <- summary.down[,c("q.val.x", "q.val.y", "q.val")]
colnames(summary.down) <- c("q.val.cont.vs.reg", "q.val.reg.vs.prog", "q.val.prog.vs.ca")
o <- order(apply(summary.down, 1, sum))
summary.down <- summary.down[o,]

WriteXLS(
  c("summary.up", "summary.down", "cont.vs.reg.up", "cont.vs.reg.down", "reg.vs.prog.up", "reg.vs.prog.down", "prog.vs.ca.up", "prog.vs.ca.down"),
  ExcelFileName = paste(results_dir,"Fig_5A_data_gxn_pathways.xlsx", sep=""), row.names = T, AdjWidth = T
)

########################################################################################################
# Figure 5B - Expression of CIN-related genes
########################################################################################################
pdf(paste(results_dir, 'Fig5B_CIN_mean_expression.pdf', sep=""))
plotGeneWithAuc(cin_genes[which(cin_genes %in% rownames(uvv)[which(uvv$fdr < 0.01)])], 
                title = "Mean expression of significant CIN genes")
dev.off()

########################################################################################################
# Figure 5C - NEK2 expression by group
########################################################################################################
pdf(paste(results_dir, 'Fig5C_NEK2_by_group.pdf', sep=""))
plotGeneByGroup("NEK2", legend.pos='topleft')
dev.off()

########################################################################################################
# Figure 5D - Expression of CIN-related genes
########################################################################################################
pdf(paste(results_dir, 'Fig5D_NEK2_expression_with_AUC.pdf', sep=""))
plotGeneWithAuc("NEK2")
dev.off()


########################################################################################################
# Extended Data Figures 1-3
# Genomic data not included here - based on Sanger pipeline
########################################################################################################

########################################################################################################
# Extended Data Figure 4A-B - Methylation-derived CNA Plots
########################################################################################################
# Plots are derived from ChAMP pipeline:
# TODO


########################################################################################################
# Extended Data Figure 4C - Methylation-derived CNA Heatmap
########################################################################################################
sig_bands <- cnas.band[rownames(cdiff[which(cdiff$fdr < 0.01),]),]

c.annot <- data.frame(
  pack.years=cnas.pheno$smoking_group,
  age.group=cnas.pheno$age_group,
  gender=cnas.pheno$Gender,
  COPD=cnas.pheno$COPD,
  status=c("Regressive", "Progressive")[cnas.pheno$progression+1]
)
c.annot_colors <- list(
  status=c(Progressive="red", Regressive="green"),
  pack.years=smoking_group_names,
  age.group=age_group_names,
  gender=c("F"="pink", "M"="blue"),
  COPD=c("NO"="green", "YES"="red")
)
rownames(c.annot) <- colnames(sig_bands)

pdf(paste(opdir, 'Ext_Data4C_cna_heatmap.pdf', sep=""))
pheatmap(removeOutliers(sig_bands), cluster_rows=T, cluster_cols=T, scale="row", main=paste("Copy number top (",dim(sig_bands)[1]," bands)", sep=""),
         annotation_col=c.annot, treeheight_row=0, treeheight_col = F, show_rownames=F, show_colnames=F,
         annotation_colors=c.annot_colors,
         color=hmcol, legend=T)
dev.off()


########################################################################################################
# Extended Data Figure 5A-F - Methylation PCAs
########################################################################################################
pcfit <- prcomp(t(mdata))
# Calculate p-values using multivariate ANOVA (all predictors are categorical)
aovdata <- mpheno
aovdata$pca <- pcfit$x[,1]
myaov <- aov(pca ~ progression + smoking_group + COPD + Previous.Lung.CA.........No..0..Yes..1. + age_group + Gender + Slide, data=aovdata)
pdf(paste(results_dir, 'Ext_Data5A-F_methylation_PCAs.pdf', sep=""))
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
plot(pcfit$x[,1:2], col=c("blue", "red")[as.numeric(as.character(mpheno$Previous.Lung.CA.........No..0..Yes..1.))+1],
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

########################################################################################################
# Extended Data Figure 5G-K - Gene Expression PCAs
########################################################################################################
pcfit <- prcomp(t(gdata))
aovdata <- gpheno
aovdata$pca <- pcfit$x[,1]
myaov <- aov(pca ~ progression + smoking_group + COPD + Prev.History.of.LC + age_group + Gender, data=aovdata)
pdf(paste(results_dir, 'Ext_Data5G-K_gxn_PCAs.pdf', sep=""))
par(mfrow=c(1,1), xpd=T)

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

########################################################################################################
# Extended Data Figure 6 - misclassified lesions
# Data not generated from R.
########################################################################################################

########################################################################################################
# Extended Data Figure 7A-C - Methylation Predictive analyses
########################################################################################################
# Base our model on significant MVPs
mvps.shared <- mdiff$probeID[which(mdiff$adj.P.Val < 0.05 & abs(mdiff$deltaBeta) > 0.3)]

pamr.mdata <- tcga.mdata.all[mvps.shared,]
sel.na <- which(apply(pamr.mdata, 1, function(x){any(is.na(x))}))
if(length(sel.na) > 0){pamr.mdata <- pamr.mdata[-sel.na,]}
mvps.shared <- rownames(pamr.mdata)

pamr.mpheno <- tcga.mpheno.all
pamr.mpheno$train <- 0
pamr.mpheno$train[which(pamr.mpheno$Sample_Name %in% mpheno.d$Sample_Name)] <- 1

o <- order(pamr.mpheno$progression)
pamr.mdata <- pamr.mdata[,o]
pamr.mpheno <- pamr.mpheno[o,]

sel <- which(pamr.mpheno$train == 1 & pamr.mpheno$source == "Surveillance")
meth.pamr.traindata <- list(
  x=as.matrix(pamr.mdata[mvps.shared,sel]),
  y=pamr.mpheno$progression[sel],
  geneid = mvps.shared
)
sel <- which(pamr.mpheno$train == 0 & pamr.mpheno$source == "Surveillance")
meth.pamr.testdata <- list(
  x=as.matrix(pamr.mdata[mvps.shared,sel]),
  y=pamr.mpheno$progression[sel],
  geneid = mvps.shared
)
sel <- which(pamr.mpheno$source == "TCGA")
meth.pamr.tcgadata <- list(
  x=as.matrix(pamr.mdata[mvps.shared,sel]),
  y=pamr.mpheno$progression[sel],
  geneid = mvps.shared
)
meth.pamr.trainfit <- pamr.train(meth.pamr.traindata)
meth.pamr.mycv <- pamr.cv(meth.pamr.trainfit, meth.pamr.traindata)
threshold <- max(meth.pamr.mycv$threshold[which(meth.pamr.mycv$error == min(meth.pamr.mycv$error))])
# Manually choose threshold 
threshold=8

pdf(paste(results_dir, 'Sup_Fig_07A-C_meth_prediction.pdf', sep=""))
meth.pamr.pred.d <- pamr.predict(meth.pamr.trainfit, newx=meth.pamr.traindata$x, type='posterior', threshold=threshold)
plot(meth.pamr.pred.d[,2], col=c("green", "red", "blue")[as.numeric(factor(pamr.mpheno[match(colnames(meth.pamr.traindata$x), pamr.mpheno$Sample_Name),]$Sample_Group, levels=c("Regressive", "Progressive", "Control")))], main="Discovery Set", ylab="Progression Score", ylim=c(0,1))
abline(v=length(which(meth.pamr.traindata$y == meth.pamr.traindata$y[1]))+0.5, col='grey')

# State how many features used:
features <- pamr.listgenes(meth.pamr.trainfit, data=meth.pamr.traindata, threshold=threshold)
if(show.legends){
  legend("topleft", paste("Features used:", dim(features)[1]))
}

meth.pamr.pred.v <- pamr.predict(meth.pamr.trainfit, newx=meth.pamr.testdata$x, type='posterior', threshold=threshold)
plot(meth.pamr.pred.v[,2], col=c("green", "red", "blue")[as.numeric(factor(pamr.mpheno[match(colnames(meth.pamr.testdata$x), pamr.mpheno$Sample_Name),]$Sample_Group, levels=c("Regressive", "Progressive", "Control")))], main="Validation Set", ylab="Progression Score", ylim=c(0,1))
abline(v=length(which(meth.pamr.testdata$y == meth.pamr.testdata$y[1]))+0.5, col='grey')
meth.roc.v <- roc(predictor=meth.pamr.pred.v[,2], response=meth.pamr.testdata$y)
if(show.legends){
  legend('bottomright', legend=paste("AUC", signif(auc(meth.roc.v), 3), sep="="))
}

meth.pamr.pred.t <- pamr.predict(meth.pamr.trainfit, newx=meth.pamr.tcgadata$x, type='posterior', threshold=threshold)
plot(meth.pamr.pred.t[,2], col=c("darkgreen", "orange")[as.numeric(as.character(meth.pamr.tcgadata$y)) + 1], main="TCGA Data", ylab="Progression Score", ylim=c(0,1))
abline(v=length(which(meth.pamr.tcgadata$y == meth.pamr.tcgadata$y[1]))+0.5, col='grey')
meth.roc.t <- roc(predictor=meth.pamr.pred.t[,2], response=meth.pamr.tcgadata$y)
if(show.legends){
  legend('left', legend=paste("AUC", signif(auc(meth.roc.t), 3), sep="="))
}
dev.off()


########################################################################################################
# Extended Data Figure 7D-F - Methylation-derived Copy Number Predictive analyses
########################################################################################################

# Include data from Dutch group (van Boerdonk et al)
pamr.cdata <- runComBat(cnas.band, dutch.bands)
pamr.cdata <- cbind(pamr.cdata[[1]], pamr.cdata[[2]])
pamr.cpheno <- rbind(
  data.frame(name=cnas.pheno$Sample_Name, progression=cnas.pheno$progression, train=1, source="Surveillance"),
  data.frame(name=rownames(dutch.pheno), progression=dutch.pheno$progression, train=0, source="Dutch")
)

# Reduce to differentially expressed bands:
pamr.cdata <- pamr.cdata[which(rownames(pamr.cdata) %in% rownames(cdiff)),]

# Include TCGA relative data
pamr.cdata <- runComBat(pamr.cdata, tcga.cnas.bands)
pamr.cdata <- cbind(pamr.cdata[[1]], pamr.cdata[[2]])
pamr.cpheno <- rbind(
  pamr.cpheno,
  data.frame(name=tcga.cnas.pheno$name, progression=tcga.cnas.pheno$progression, train=0, source="TCGA")
)

o <- order(pamr.cpheno$progression)
pamr.cdata <- pamr.cdata[,o]
pamr.cpheno <- pamr.cpheno[o,]

sel <- which(pamr.cpheno$source == "Surveillance")
cna.pamr.traindata <- list(
  x=as.matrix(pamr.cdata[,sel]),
  y=pamr.cpheno$progression[sel],
  geneid = rownames(pamr.cdata)
)
sel <- which(pamr.cpheno$source == "Dutch")
cna.pamr.testdata <- list(
  x=as.matrix(pamr.cdata[,sel]),
  y=pamr.cpheno$progression[sel],
  geneid = rownames(pamr.cdata)
)
sel <- which(pamr.cpheno$source == "TCGA")
cna.pamr.testdata2 <- list(
  x=as.matrix(pamr.cdata[,sel]),
  y=pamr.cpheno$progression[sel],
  geneid = rownames(pamr.cdata)
)

set.seed(2)
cna.pamr.trainfit <- pamr.train(cna.pamr.traindata)
cna.pamr.mycv <- pamr.cv(cna.pamr.trainfit, cna.pamr.traindata)
threshold <- max(cna.pamr.mycv$threshold[which(cna.pamr.mycv$error == min(cna.pamr.mycv$error))])
# Manually choose threshold 
threshold=1.5

pdf(paste(results_dir, 'Ext_Data7D-F_cna_prediction.pdf', sep=""))

plot(cna.pamr.mycv$prob[,2,1], 
     col=c("green", "red")[as.numeric(as.character(cna.pamr.traindata$y)) + 1], 
     main=paste("Cross-validated probabilities"), 
     ylab="Progression Score", ylim=c(0,1))
abline(v=length(which(cna.pamr.traindata$y == cna.pamr.traindata$y[1]))+0.5, col='grey')
myauc=auc(roc(predictor=cna.pamr.mycv$prob[,2,1], response=as.numeric(as.character(cna.pamr.traindata$y))))
if(show.legends){
  legend('right', paste("AUC=", signif(myauc, 3), sep=""))
  features <- pamr.listgenes(cna.pamr.trainfit, cna.pamr.traindata, threshold = threshold)
  legend('left', paste("Features used:", dim(features)[1]))
}


cna.pamr.pred.v <- pamr.predict(cna.pamr.trainfit, newx=cna.pamr.testdata$x, type='posterior', threshold=threshold)
plot(cna.pamr.pred.v[,2], col=c("green", "red")[as.numeric(as.character(cna.pamr.testdata$y)) + 1], main="CNA Validation with Van Boerdonk data", ylab="Progression Score", ylim=c(0,1))
abline(v=length(which(cna.pamr.testdata$y == cna.pamr.testdata$y[1]))+0.5, col='grey')
gxn.roc.v <- roc(predictor=cna.pamr.pred.v[,2], response=cna.pamr.testdata$y)
if(show.legends){
  legend('topleft', legend=paste("AUC", signif(auc(gxn.roc.v), 3), sep="="))
}

cna.pamr.pred.v <- pamr.predict(cna.pamr.trainfit, newx=cna.pamr.testdata2$x, type='posterior', threshold=threshold)
plot(cna.pamr.pred.v[,2], col=c("darkgreen", "orange")[as.numeric(as.character(cna.pamr.testdata2$y)) + 1], main="CNA Validation with TCGA data", ylab="Progression Score", ylim=c(0,1))
abline(v=length(which(cna.pamr.testdata2$y == cna.pamr.testdata2$y[1]))+0.5, col='grey')
gxn.roc.v <- roc(predictor=cna.pamr.pred.v[,2], response=cna.pamr.testdata2$y)
if(show.legends){
  legend('topleft', legend=paste("AUC", signif(auc(gxn.roc.v), 3), sep="="))
}

dev.off()


########################################################################################################
# Extended Data Figure 8 - Methylation Heterogeneity Index for PvR using Random Sampling
########################################################################################################
# TODO

########################################################################################################
# Extended Data Figure 9A - Methylation Pathway Analysis
# Figure plotted in Prism. Data output to Excel file.
########################################################################################################
# Map genes to probes by taking the mean beta value of all genes over a probe
# Aggregate by gene with mean function
data("probe.features")
tcga.mdata.all.genes <- tcga.mdata.all
gene.map <- probe.features$gene[match(rownames(tcga.mdata.all.genes), rowames(probe.features))]
tcga.mdata.all.genes <- aggregate(tcga.mdata.all.genes, by=list(gene.map), FUN=mean)
rownames(tcga.mdata.all.genes) <- tcga.mdata.all.genes$Group.1
tcga.mdata.all.genes$Group.1 <- NULL

sel <- which(tcga.mpheno.all$Sample_Group == "TCGA Control" | tcga.mpheno.all$Sample_Group == "Regressive")
gage_0_1 <- gage_analysis(tcga.mdata.all.genes[,sel], tcga.mpheno.all[sel,]$dose, source='msig.kegg')
sel <- which(tcga.mpheno.all$Sample_Group == "Regressive" | tcga.mpheno.all$Sample_Group == "Progressive")
gage_1_2 <- gage_analysis(tcga.mdata.all.genes[,sel], tcga.mpheno.all[sel,]$dose-1, source='msig.kegg')
sel <- which(tcga.mpheno.all$Sample_Group == "Progressive" | tcga.mpheno.all$Sample_Group == "TCGA SqCC")
gage_2_3 <- gage_analysis(tcga.mdata.all.genes[,sel], tcga.mpheno.all[sel,]$dose-2, source='msig.kegg')
# Save
cont.vs.reg.up <- data.frame(gage_0_1$greater)
cont.vs.reg.down <- data.frame(gage_0_1$less)
reg.vs.prog.up <- data.frame(gage_1_2$greater)
reg.vs.prog.down <- data.frame(gage_1_2$less)
prog.vs.ca.up <- data.frame(gage_2_3$greater)
prog.vs.ca.down <- data.frame(gage_2_3$less)

# Make summary sheets 
summary.up <- merge(cont.vs.reg.up, reg.vs.prog.up,by=0)
rownames(summary.up) <- summary.up$Row.names
summary.up$Row.names <- NULL
summary.up <- merge(summary.up, prog.vs.ca.up, by=0)
rownames(summary.up) <- summary.up$Row.names
summary.up <- summary.up[,c("q.val.x", "q.val.y", "q.val")]
colnames(summary.up) <- c("q.val.cont.vs.reg", "q.val.reg.vs.prog", "q.val.prog.vs.ca")
o <- order(summary.up$q.val.cont.vs.reg)
summary.up <- summary.up[o,]

summary.down <- merge(cont.vs.reg.down, reg.vs.prog.down,by=0)
rownames(summary.down) <- summary.down$Row.names
summary.down$Row.names <- NULL
summary.down <- merge(summary.down, prog.vs.ca.down, by=0)
rownames(summary.down) <- summary.down$Row.names
summary.down <- summary.down[,c("q.val.x", "q.val.y", "q.val")]
colnames(summary.down) <- c("q.val.cont.vs.reg", "q.val.reg.vs.prog", "q.val.prog.vs.ca")
o <- order(summary.down$q.val.cont.vs.reg)
summary.down <- summary.down[o,]

WriteXLS(
  c("summary.up", "summary.down", "cont.vs.reg.up", "cont.vs.reg.down", "reg.vs.prog.up", "reg.vs.prog.down", "prog.vs.ca.up", "prog.vs.ca.down"),
  ExcelFileName = paste(results_dir,"ExtData_Fig9A_methyl_pathways.xlsx", sep=""), row.names = T, AdjWidth = T
)


########################################################################################################
# Extended Data Figure 9B - Methylation of NEK2-associated probe by group
########################################################################################################
pdf(paste(results_dir, 'Ext_Data9B_methylation_of_NEK2_probe_by_group.pdf', sep=""))
plotMethylationByGroup("cg17931972", "NEK2", legend.pos="topleft")
dev.off()

########################################################################################################
# Extended Data Figure 10 - Correlation of WGII with CIN gene expression
########################################################################################################
# TODO