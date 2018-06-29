##########################################################################
# Setup
#
# Load gene expression and methylation data from NCBI GEO
##########################################################################

# Auto-install dependencies
source('install_dependencies.R')

# Optionally change data_dir to set where data is downloaded to
# This variable is required before sourcing the data_loader scripts
data_cache <- "./data/"
results_dir <- "./results/"
for(dir in c(data_cache, results_dir)){
  dir.create(dir, recursive = T, showWarnings = F)
}


# Load libraries in a specific order to avoid overwriting functions
library(gdata)
library(ChAMP)
library(ggplot2)
library(ggsignif)
library(pheatmap)
library(pROC)
library(PRROC)
library(RColorBrewer)
library(limma)
library(stringr)
library(WriteXLS)
library(pamr)
library(affycoretools)
library(foreach)
library(doParallel)
library(GenomicDataCommons)

# Load utility functions
for(file in list.files("./utility_functions/", full.names = T)){
  source(file)
}

# Load figure plotting functions
for(file in list.files("./plot_functions/", full.names = T)){
  source(file)
}

# For some analyses we use parallel processing to take advantage of multiple cores
cores=detectCores()
if(cores > 1){
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  registerDoParallel(cl)
}

##########################################################################
# Load Data
##########################################################################
# These functions load CIS and TCGA data, and make many variables available
# See the individual files for details of data pre-processing
source('data_loaders/loadGeneData.R')
source('data_loaders/loadMethData.R')
source('data_loaders/loadWgsData.R')

# Load CIN genes - both the CIN70 signature, and CIN70 with cell-cycle genes removed
load('resources/cin_genes.RData')

# Load a list of overlapping samples
overlap.pheno <- read.xls('resources/overlap.samples.xlsx')

# Read in a list of genes previously associated with lung cancer as defined in the text
#driver.genes <- read.csv('resources/driver_mutations.csv')
# Filter for pan-cancer or lusc-specific
#driver.genes.info <- read.table('resources/driver_genes.tsv', stringsAsFactors = F, sep="\t", header=T)
#driver.genes <- unique(driver.genes.info$Gene)
#load("resources/cgc.genes.RData")
driver.genes.info <- read.xls('resources/driver.genes.xls', stringsAsFactors=F)
driver.genes <- driver.genes.info$Gene.Symbol


##########################################################################
# Define additional variables 
#   e.g. colour palettes
##########################################################################
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

# Gene Expression
gdiff <- limmaCompare(data=gdata.d, pheno=gpheno.d, fdr_limit = 1)

# Methylation
myDMP <- champ.DMP(
  beta=mdata.d,
  pheno=mpheno.d$Sample_Group,
  adjPVal = 1
)
mdiff <- myDMP$Progressive_to_Regressive

# Also calculated differentially methylated regions (DMRs):
dmrs <- champ.DMR(beta=as.matrix(mdata.d), pheno=mpheno.d$Sample_Group, compare.group = c("Progressive", "Regressive"), method="ProbeLasso")
# Additionally calculate TCGA DMRs for comparison:
tcga.mpheno.tmp <- tcga.mpheno
tcga.mpheno.tmp$Sample_Group <- make.names(tcga.mpheno.tmp$Sample_Group)
# Impute to remove NAs
tcga.mdata.imputed <- champ.impute(beta=as.matrix(tcga.mdata), SampleCutoff = 0.5, ProbeCutoff = 0.5, pd=tcga.mpheno.tmp)
tcga.mpheno.tmp <- tcga.mdata.imputed$pd
tcga.mdata.imputed <- tcga.mdata.imputed$beta
# Strange bug - only use probes in package data(illumina450Gr) for DM (removes very few probes)
data(illumina450Gr)
tcga.mdata.imputed <- tcga.mdata.imputed[which(rownames(tcga.mdata.imputed) %in% names(illumina450Gr)),]
# Find DMRs
dmrs.tcga <- champ.DMR(beta=tcga.mdata.imputed, pheno=tcga.mpheno.tmp$Sample_Group, compare.group=c("TCGA.SqCC", "TCGA.Control"), method="ProbeLasso")

# Methylation-derived copy number (using logR aggregated by cytogenetic band)
cdiff <- limmaCompare(mcnas.band, mcnas.pheno, fdr_limit = 0.01)

# Repeat differential analysis using a continuous variable of 'dose'
# This is defined as 0=TCGA control, 1=regressive CIS, 2=progressive CIS, 3=TCGA cancer
AAprog <- as.numeric(tcga.gpheno.all$dose)
design <- model.matrix(~AAprog)
fit <- lmFit(tcga.gdata.all[,sel], design)
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

# Repeat for methylation - use ChAMP package which uses limma internally, and gives identical results plus additional annotation
mdiff.dose <- champ.DMP(
  beta=tcga.mdata.all,
  pheno=as.numeric(tcga.mpheno.all$dose)
)
mdiff.dose <- mdiff.dose$NumericVariable


##########################################################################
# Differential expression analysis of genomic data
##########################################################################
# Look for genes which perfectly differentiate P vs R
sel <- which(
  (all(muts[,which(wgs.pheno$progression == 1)] > 0) & all(muts[,which(wgs.pheno$progression == 0)] == 0)) |
  (all(muts[,which(wgs.pheno$progression == 1)] == 0) & all(muts[,which(wgs.pheno$progression > 0)] == 0))
)
if(length(sel) > 0){
  print(paste("Genes perfectly separating progressive and regressive:", paste(rownames(muts)[sel], collapse=", ")))
}else{
  print("No genes perfectly separated progressive and regressive.")
}
# We find no genes here.

# Now try applying a statistical test (Fisher's) for every gene and correcting for multiple testing
# Here we are comparing the presence of any subs or indels PvR
sig.genomic <- data.frame(
  gene=rownames(muts),
  p=unlist(lapply(rownames(muts), function(x){
    if(all(as.numeric(muts[x,]) == as.numeric(muts[x,1]))){
      return(1)
    }
    return(fisher.test(as.numeric(muts[x,]) > 0, wgs.pheno$progression)$p.value)
  }))
)
sig.genomic$p.adj <- p.adjust(sig.genomic$p)
print(paste("Genes with significantly different mutations P vs R:", length(which(sig.genomic$p.adj < 0.05))))
# Again we find no genes (not surprising given our sample numbers).


##########################################################################
# Predictive modelling of gene expression
#
# Uses PAM
# Data is divided into discovery and validation sets. Additional external validation sets are used (TCGA).
# k-fold cross-validation is applied to the discovery set
##########################################################################

# GXN PAM prediction
set.seed(2)

# Use differentially expressed genes which are present in TCGA data
genes.shared <- intersect(rownames(tcga.gdata.all), rownames(gdiff)[which(gdiff$fdr < 0.01)])

pamr.gpheno <- tcga.gpheno.all
pamr.gpheno$train <- 0
pamr.gpheno$train[which(pamr.gpheno$name %in% as.character(gpheno$name)[which(gpheno$training == 1)])] <- 1

o <- order(pamr.gpheno$dose)
pamr.gpheno <- pamr.gpheno[o,]
pamr.gdata <- tcga.gdata.all[,o]

# Training set. Cross-validation is performed on this, and a threshold selected.
sel <- which(pamr.gpheno$train == 1 & pamr.gpheno$source == "Surveillance")
gxn.pamr.traindata <- list(
  x=as.matrix(pamr.gdata[genes.shared,sel]),
  y=pamr.gpheno$progression[sel],
  geneid = genes.shared
)
# Independent validation set, of CIS data, using an orthogonal platform (Affymetrix)
sel <- which(pamr.gpheno$train == 0 & pamr.gpheno$source == "Surveillance")
gxn.pamr.testdata <- list(
  x=as.matrix(pamr.gdata[genes.shared,sel]),
  y=pamr.gpheno$progression[sel],
  geneid = genes.shared
)
# External validation set from the TCGA (here we apply our model to cancer/control samples)
sel <- which(pamr.gpheno$source == "TCGA")
gxn.pamr.tcgadata <- list(
  x=as.matrix(pamr.gdata[genes.shared,sel]),
  y=pamr.gpheno$progression[sel],
  geneid = genes.shared
)

# Train the model and apply k-fold cross-validation
gxn.pamr.trainfit <- pamr.train(gxn.pamr.traindata)

# For reproducibility, manually define folds (these were generated using a call to pamr.cv)
folds <- list(c(14, 12, 19, 26), c(4, 7, 24, 27), c(11, 32, 29), c(6, 2, 17), c(16, 8, 31, 23), c(13, 22), c(5, 15, 25, 20), c(1, 33, 28), c(10, 9, 30), c(3, 21, 18))
gxn.pamr.mycv <- pamr.cv(gxn.pamr.trainfit, gxn.pamr.traindata, folds = folds, nfold = length(folds))
# Manually choose threshold from experimentation - use pamr.plotcv(gxn.pamr.mycv) to help choose a threshold:
gxn.threshold.id <- 16
gxn.threshold <- gxn.pamr.mycv$threshold[[gxn.threshold.id]]

# Identify the genes used in this model
pamr.gxn.features <- pamr.listgenes(gxn.pamr.trainfit, data=gxn.pamr.traindata, threshold=gxn.threshold)


##########################################################################
# Predictive modelling of methylation
#
# Uses PAM, same method as for gene expression
# Data is divided into discovery and validation sets. Additional external validation sets are used (TCGA).
# k-fold cross-validation is applied to the discovery set
##########################################################################
set.seed(2)
# Base our model on significant MVPs
mvps.shared <- rownames(mdiff)[which(mdiff$adj.P.Val < 0.05 & abs(mdiff$deltaBeta) > 0.3)]

pamr.mdata <- tcga.mdata.all[mvps.shared,]
sel.na <- which(apply(pamr.mdata, 1, function(x){any(is.na(x))}))
if(length(sel.na) > 0){pamr.mdata <- pamr.mdata[-sel.na,]}
mvps.shared <- rownames(pamr.mdata)

pamr.mpheno <- tcga.mpheno.all
pamr.mpheno$train <- 0
pamr.mpheno$train[which(pamr.mpheno$name %in% mpheno.d$Sample_Name)] <- 1

o <- order(pamr.mpheno$progression)
pamr.mdata <- pamr.mdata[,o]
pamr.mpheno <- pamr.mpheno[o,]

sel <- which(pamr.mpheno$train == 1 & pamr.mpheno$source == "Surveillance")
# Note that for methylation we compare progressive vs (regressive + control) but include group data here for plots
meth.pamr.traindata <- list(
  x=as.matrix(pamr.mdata[mvps.shared,sel]),
  y=pamr.mpheno$progression[sel],
  geneid = mvps.shared,
  group = pamr.mpheno$Sample_Group[sel]
)
sel <- which(pamr.mpheno$train == 0 & pamr.mpheno$source == "Surveillance")
meth.pamr.testdata <- list(
  x=as.matrix(pamr.mdata[mvps.shared,sel]),
  y=pamr.mpheno$progression[sel],
  geneid = mvps.shared,
  group = pamr.mpheno$Sample_Group[sel]
)
sel <- which(pamr.mpheno$source == "TCGA")
meth.pamr.tcgadata <- list(
  x=as.matrix(pamr.mdata[mvps.shared,sel]),
  y=pamr.mpheno$progression[sel],
  geneid = mvps.shared,
  group = pamr.mpheno$Sample_Group[sel]
)

# Train the model and apply k-fold cross-validation
meth.pamr.trainfit <- pamr.train(meth.pamr.traindata)
# For reproducibility, manually define folds (these were generated using a call to pamr.cv)
folds <- list(c(14, 17, 11, 52, 59, 36),c(13, 9, 32, 37, 44, 42),c(6, 30, 3, 48, 60, 47),c(4, 25, 1, 16, 49, 56),c(28, 18, 10, 33, 50, 39),c(23, 31, 15, 40, 46),c(19, 8, 27, 38, 43),c(7, 2, 26, 41, 57, 45),c(29, 22, 5, 34, 54, 35, 51),c(24, 21, 12, 20, 55, 53, 58))
meth.pamr.mycv <- pamr.cv(meth.pamr.trainfit, meth.pamr.traindata, folds = folds, nfold=length(folds))

# Manually choose threshold from experimentation - use pamr.plotcv(meth.pamr.mycv) to help choose a threshold:
#threshold <- max(meth.pamr.mycv$threshold[which(meth.pamr.mycv$error == min(meth.pamr.mycv$error))])
#threshold=8
meth.threshold.id <- 23
meth.threshold <- meth.pamr.mycv$threshold[meth.threshold.id]

pamr.meth.features <- pamr.listgenes(meth.pamr.trainfit, data=meth.pamr.traindata, threshold=meth.threshold)


##########################################################################
# Predictive modelling of methylation-derived copy number data
#
# Uses PAM, as above
# Due to lower sample numbers, data is not divided into discovery/validation.
# Model is trained on CIS data, internally cross-validated, then applied to two external datasets:
# van Boerdonk et al (comparable pre-invasive lung data generated from arrayCGH) and TCGA (cancer/control data)
##########################################################################

# Include data from Dutch group (van Boerdonk et al)
pamr.cdata <- runComBat(mcnas.band, dutch.bands)
pamr.cdata <- cbind(pamr.cdata[[1]], pamr.cdata[[2]])
pamr.cpheno <- rbind(
  data.frame(name=mcnas.pheno$Sample_Name, progression=mcnas.pheno$progression, train=1, source="Surveillance"),
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

# Train PAMR model
# As above, use pre-defined folds for reproducibility
set.seed(2)
cna.pamr.trainfit <- pamr.train(cna.pamr.traindata)
folds <- list(c(11, 31, 48, 50),c(4, 15, 45, 37, 52),c(14, 21, 36, 26, 51),c(5, 7, 29, 44, 39, 24),c(13, 18, 41, 34, 53, 25),c(2, 6, 33, 43, 38, 42),c(10, 8, 23, 46, 32, 40),c(12, 17, 35, 22, 47, 54),c(16, 1, 28, 30, 20),c(3, 9, 49, 19, 27))
cna.pamr.mycv <- pamr.cv(cna.pamr.trainfit, cna.pamr.traindata, folds=folds, nfold = length(folds))
# Manually choose threshold 
cna.threshold.id <- 13
cna.threshold=cna.pamr.mycv$threshold[cna.threshold.id]

pamr.cna.features <- pamr.listgenes(cna.pamr.trainfit, data=cna.pamr.traindata, threshold=cna.threshold)


########################################################################################################
# Calculation of Methylation Heterogeneity Index
########################################################################################################
# This measure of methylation heterogeneity is defined in the main text.
# See mhi_analysis.R for further information and threshold calculation.

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

# Calculate MHI - number of intermediate value probes divided by total probe count
mhi <- apply(mdata.mhi, 2, function(x){length(which(x > thresh.low & x < thresh.up)) / dim(mdata.mhi)[1]})  

# To plot MHI consistently we store data based on 10000 sample runs.
# Code to generate this is in mhi_analysis.R - it is stored as an RData file for reproducible plots.
load("resources/MHIsimulation10k.RData")


########################################################################################################
# GXN Pathway analysis
# Uses Gage package to do pairwise comparison of control/reg, reg/prog, prog/cancer
########################################################################################################

gage.pvr <- gage_analysis(gdata.d, gpheno.d$progression, source = 'msig.kegg')
all.pathways <- unique(c(rownames(gage.pvr$greater), rownames(gage.pvr$less)))
gxn.gage.summary <- data.frame(row.names = all.pathways, 
                               q.val.up=gage.pvr$greater[all.pathways,'q.val'], 
                               q.val.down=gage.pvr$less[all.pathways,'q.val']
                    )
gxn.gage.summary <- gxn.gage.summary[order(gxn.gage.summary$q.val.up + gxn.gage.summary$q.val.down),]

########################################################################################################
# Methylation Pathway Analysis
########################################################################################################
# Map genes to probes by taking the mean beta value of all genes over a probe
# Aggregate by gene with mean function
data("probe.features")
mdata.all.genes <- mdata
gene.map <- probe.features$gene[match(rownames(mdata.all.genes), rownames(probe.features))]
mdata.all.genes <- aggregate(mdata.all.genes, by=list(gene.map), FUN=mean)
rownames(mdata.all.genes) <- mdata.all.genes$Group.1
mdata.all.genes$Group.1 <- NULL

gage.pvr <- gage_analysis(mdata.all.genes, mpheno$progression, source = 'msig.kegg')
all.pathways <- unique(c(rownames(gage.pvr$greater), rownames(gage.pvr$less)))
meth.gage.summary <- data.frame(row.names = all.pathways, 
                               q.val.up=gage.pvr$greater[all.pathways,'q.val'], 
                               q.val.down=gage.pvr$less[all.pathways,'q.val']
)
meth.gage.summary <- meth.gage.summary[order(meth.gage.summary$q.val.up + meth.gage.summary$q.val.down),]


##########################################################################
# Mutational Signature Analysis
##########################################################################

library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = T)
library(VariantAnnotation)

# Get colors for plots
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Read SNV VCF files (unfiltered)
# Read the filters also (ASMD >= 140 & CLPM == 0) - these aren't collected by the read_vcfs_as_granges function
vcf.files <- c()
passed.filters <- list()
is.exonic <- list()
for(i in 1:dim(wgs.pheno)[1]){
  print(paste("Checking filters for", wgs.pheno$name[i]))
  vcf <- list.files("data/wgs/caveman", pattern=wgs.pheno$name[i], full.names=T)
  filters <- fixed(readVcf(vcf))$FILTER
  clpm <- as.numeric(readInfo(vcf, x='CLPM'))
  asmd <- as.numeric(readInfo(vcf, x='ASMD'))
  vcf.files <- c(vcf.files, vcf)
  passed.filters[[i]] <- which(filters == "PASS" & clpm == 0 & asmd >= 140)
  vd <- as.character(readInfo(vcf, x='VD'))
  is.exonic[[i]] <- grepl("exon", vd)
}
# Read VCFs
vcfs <- read_vcfs_as_granges(vcf.files, sample_names = wgs.pheno$name, ref_genome, check_alleles = F)
vcfs.exonic <- list()
# Apply filters:
for(i in 1:length(vcfs)){
  sel <- passed.filters[[i]]
  sel.exonic <- which(is.exonic[[i]])
  vcfs.exonic[[i]] <- vcfs[[i]][intersect(sel, sel.exonic)]
  vcfs[[i]] <- vcfs[[i]][sel]
}
names(vcfs.exonic) <- names(vcfs)

# Calculate spectrum across all samples
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
#plot_spectrum(type_occurrences)

# Map to the 96 count system used by Alexandrov et al
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
#plot_96_profile(mut_mat, condensed = T)

# Download known signature data from COSMIC
sp_url <- "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
cancer_signatures = as.matrix(cancer_signatures[,4:33])

#plot_96_profile(cancer_signatures[,1:21], condensed = TRUE, ymax = 0.3)

# Cluster similar signatures together
hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
#plot(hclust_cosmic)

# Compare with our data
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)


# Fit COSMIC signatures to our data
fit_res <- fit_to_signatures(mut_mat, cancer_signatures)

# Find the top signatures by adding up relative contributions
fit_res.rel <- apply(fit_res$contribution, 2, function(x){x/sum(x)})
ordered.sigs <- names(sort(rowSums(fit_res.rel), decreasing = T))

# From this we pick top signatures and force all mutations to align to these.
# We are also influenced here by previous data in our choice of signatures.
sigs.to.analyse <- c(
  "Signature.1", "Signature.2", "Signature.4", "Signature.5", "Signature.13"
)
cancer_signatures.used <- cancer_signatures[,sigs.to.analyse]
fit_res <- fit_to_signatures(mut_mat, cancer_signatures.used)


# Add TCGA signatures - based on Caveman calls of TCGA data
tcga_mut_mat <- read.table("resources/TCGA_trinucleotide_counts.txt", header=T, stringsAsFactors = F)
# Coerce to correct format
tcga_mut_mat <- t(tcga_mut_mat)
strs <- strsplit(rownames(tcga_mut_mat), "[.]")
rnames <- unlist(lapply(strs, function(x){
  paste0(
    substr(x[4],1,1),
    "[",x[1],">",x[2],"]",
    substr(x[4],3,3)
  )
}))
rownames(tcga_mut_mat) <- rnames

fit_res_tcga <- fit_to_signatures(tcga_mut_mat, cancer_signatures.used)

# We can also make a mut_type_occurrences matrix
type_occurrences_tcga <- aggregate(tcga_mut_mat, by=list(substr(rownames(tcga_mut_mat), 3, 5)), FUN=sum)
rownames(type_occurrences_tcga) <- type_occurrences_tcga$Group.1
type_occurrences_tcga$Group.1 <- NULL
type_occurrences_tcga <- t(type_occurrences_tcga)
# Compare to exonic muts only
mut_mat_exonic <- mut_matrix(vcf_list = vcfs.exonic, ref_genome = ref_genome)
fit_res_exonic <- fit_to_signatures(mut_mat_exonic, cancer_signatures.used)
type_occurrences_exonic <- mut_type_occurrences(vcfs.exonic, ref_genome)


##########################################################################
# Clonality analysis
##########################################################################
# Due to complexity and long processing times, clonality analysis is in a separate file.
# Data is cached where available.
source("analysis_functions/clonality.analysis.ccf.R")


# Now analyse the cluster data
# Create output plots and calculate key attributes per sample - number of clusters and number of clonal/subclonal mutations
opdir <- paste0(results_dir, "clonality/")
dir.create(opdir, recursive = T, showWarnings = F)
wgs.pheno$nclusters <- NA
wgs.pheno$clonal.muts <- NA
wgs.pheno$subclonal.muts <- NA
wgs.pheno$muts.in.clonal.cluster <- NA
for(i in 1:dim(wgs.pheno)[1]){
  sample <- wgs.pheno$name[i]
  sc <- clusterdata[[sample]]
  sample.muts <- muts.with.ccf[[sample]]
  
  if(is.null(sc)){
    next
  }
  
  # Cluster plots
  writeClusterTable(sc, paste(opdir, "clusters_", sample, ".tsv", sep=""))
  sc.plot1d(sc, paste(opdir, "clusters_",sample,".1d.pdf", sep=""))
  
  # Cluster counts
  clusters <- unique(sc@vafs.merged$cluster)
  clusters <- clusters[which(!is.na(clusters) & clusters != 0)]
  cluster.means <- unlist(lapply(clusters, function(x){
    mean(sc@vafs.merged[which(sc@vafs.merged$cluster == x),5], na.rm = T)
  }))
  names(cluster.means) <- clusters
  clonal.cluster <- names(cluster.means)[which(abs(cluster.means - 50) == min(abs(cluster.means - 50)))]
  wgs.pheno$nclusters[i] <- length(clusters)
  wgs.pheno$clonal.muts[i] <- length(which(sample.muts$is.clonal))
  wgs.pheno$subclonal.muts[i] <- length(which(!sample.muts$is.clonal))
  wgs.pheno$muts.in.clonal.cluster[i] <- length(which(sc@vafs.merged$cluster == clonal.cluster))
}



##########################################################################
# DMR analysis
#
# Check for overlaps between CIS DMRs and TCGA cancer/control DMRs
##########################################################################
dmrs.range <- GRanges(
  seqnames = dmrs$ProbeLassoDMR$seqnames,
  ranges = IRanges(
    start = dmrs$ProbeLassoDMR$start,
    end = dmrs$ProbeLassoDMR$end
  ),
  strand = dmrs$ProbeLassoDMR$strand,
  deltaBeta = dmrs$ProbeLassoDMR$betaAv_Progressive - dmrs$ProbeLassoDMR$betaAv_Regressive
)
dmrs.tcga.range <- GRanges(
  seqnames = dmrs.tcga$ProbeLassoDMR$seqnames,
  ranges = IRanges(
    start = dmrs.tcga$ProbeLassoDMR$start,
    end = dmrs.tcga$ProbeLassoDMR$end
  ),
  strand = dmrs.tcga$ProbeLassoDMR$strand,
  deltaBeta = dmrs.tcga$ProbeLassoDMR$betaAv_TCGA.SqCC - dmrs.tcga$ProbeLassoDMR$betaAv_TCGA.Control
)
dmr.overlaps <- countOverlaps(dmrs.range, dmrs.tcga.range)
print(paste("Of", dim(dmrs$ProbeLassoDMR)[1], 'identified DMRs,', length(which(dmr.overlaps == 1)), 'are identified in TCGA cancer vs control -', 100*length(which(dmr.overlaps == 1))/length(dmr.overlaps), "%"))


##########################################################################
# Driver analysis
#
# Extract potential driver mutations for supplementary file 1
# Add CN amplifications and deletions to muts.all:
##########################################################################
cnas.amps2 <- lapply(names(cnas.amps), function(x){
  if(dim(cnas.amps[[x]])[1] == 0){return(NA)}
  cnas.amps[[x]]$patient <- x
  return(cnas.amps[[x]])
})
cnas.amps2 <- cnas.amps2[which(!is.na(cnas.amps2))]
cnas.amps2 <- rbindlist2(cnas.amps2)
cnas.amps2$class="AMP"
cnas.amps2$type="CN amplification"

cnas.dels2 <- lapply(names(cnas.dels), function(x){
  if(dim(cnas.dels[[x]])[1] == 0){return(NA)}
  cnas.dels[[x]]$patient <- x
  return(cnas.dels[[x]])
})
cnas.dels2 <- cnas.dels2[which(!is.na(cnas.dels2))]
cnas.dels2 <- rbindlist2(cnas.dels2)
cnas.dels2$class="DEL"
cnas.dels2$type="CN deletion"

cnas.all <- rbind(cnas.amps2, cnas.dels2)
df <- data.frame(
  patient=cnas.all$patient,
  gene=as.character(cnas.all$SYMBOL),
  class=cnas.all$class,
  type=cnas.all$type,
  ref=NA,
  alt=NA,
  chr=gsub("chr", "", cnas.all$seqnames),
  start=cnas.all$start,
  end=cnas.all$end,
  filters="PASS",
  asmd=140,
  clpm=0,
  vaf=NA,ref.reads=NA,alt.reads=NA,tumour.reads=NA,depth=NA,exonic=NA,protein.change=NA,filters.passed=T,mid=NA,translocation.partner=NA,chr2=NA,start2=NA,end2=NA
)
muts.all.cn <- rbind(
  muts.all[which(muts.all$type %in% c("missense", "nonsense", "start_lost", "stop_lost", "ess_splice", "splice_region", "frameshift", "inframe", "SO:0000010:protein_coding", "CN amplification", "CN deletion", "Rearrangement")),]
  , df)

# Select for driver genes:
t <- table(as.character(muts.all.cn$gene)[which(as.character(muts.all.cn$gene) %in% driver.genes)])
driver.genes.sorted <- names(t)[order(as.numeric(t), decreasing = T)]
driver.muts <- lapply(driver.genes.sorted, function(x){
  m <- muts.all.cn[which(muts.all.cn$gene == x),]
  if(dim(m)[1] == 0){ return(NA) }
  df <- data.frame(
    Sample=as.character(m$patient),
    Gene=x,
    Mutation.Type=m$type,
    Chromosome=m$chr,
    Start.Position=m$start,
    End.Position=m$end,
    Reference.Allele=m$ref,
    Variant.Allele=m$alt,
    Protein.Change=m$protein.change,
    Transloc.Chr=m$chr2,
    Transloc.Start=m$start2,
    Transloc.End=m$end2,
    Transloc.Gene=m$translocation.partner,
    Outcome=c("Regression", "Progression")[wgs.pheno$progression[match(m$patient, wgs.pheno$name)]+1]
  )
  return(df[order(df$Sample),])
})
names(driver.muts) <- driver.genes.sorted
driver.muts <- driver.muts[which(!is.na(driver.muts))]

# Find driver counts per sample for prog vs reg analysis
driver.muts.all <- rbindlist2(driver.muts)

# Annotate with gene function from driver genes list
driver.muts.all$Role.in.Cancer <- cgc.genes$Role.in.Cancer[match(driver.muts.all$Gene, cgc.genes$Gene.Symbol)]
driver.muts.all$CGC.Tier <- cgc.genes$Tier[match(driver.muts.all$Gene, cgc.genes$Gene.Symbol)]
driver.muts.all$Validated.Mut.Types <- cgc.genes$Mutation.Types[match(driver.muts.all$Gene, cgc.genes$Gene.Symbol)]
driver.muts.all$Validated.Transloc.Partner <- cgc.genes$Translocation.Partner[match(driver.muts.all$Gene, cgc.genes$Gene.Symbol)]
# Add VEP impact data (exclude CN changes from this)
sel <- which(!is.na(driver.muts.all$Reference.Allele))
vep.muts <- vep.annotate(driver.muts.all[sel,])
driver.muts.all$VEP.impact <- NA
driver.muts.all$VEP.impact[sel] <- vep.muts$vep.impact

# Add a flag for genuine drivers:
driver.muts.all$genuine.driver <- F

# Check the mutation type against CGC validated mutation types
# Types are: A, D, Mis, N, F, T, S, O
driver.muts.all$genuine.driver[intersect(grep("A", driver.muts.all$Validated.Mut.Types), which(driver.muts.all$Mutation.Type == "CN amplification"))] <- T
driver.muts.all$genuine.driver[intersect(grep("D", driver.muts.all$Validated.Mut.Types), which(driver.muts.all$Mutation.Type == "CN deletion"))] <- T
driver.muts.all$genuine.driver[intersect(grep("Mis", driver.muts.all$Validated.Mut.Types), which(driver.muts.all$Mutation.Type == "missense"))] <- T
driver.muts.all$genuine.driver[intersect(grep("N", driver.muts.all$Validated.Mut.Types), which(driver.muts.all$Mutation.Type == "nonsense"))] <- T
driver.muts.all$genuine.driver[intersect(grep("F", driver.muts.all$Validated.Mut.Types), which(driver.muts.all$Mutation.Type == "frameshift"))] <- T
driver.muts.all$genuine.driver[intersect(grep("S", driver.muts.all$Validated.Mut.Types), grep("splice", driver.muts.all$Mutation.Type))] <- T
# For tranlocations also check the translocated gene is the right gene
rearr.filt <- which(
  driver.muts.all$Mutation.Type == "Rearrangement" &
  grepl("T", driver.muts.all$Validated.Mut.Types) &
  driver.muts.all$Transloc.Gene == driver.muts.all$Validated.Transloc.Partner
)
if(length(rearr.filt) > 0){
  driver.muts.all$genuine.driver[rearr.filt] <- T
}

# As a further check, filter out any mutation with a VEP score less than moderate/high
vep.filter <- !is.na(driver.muts.all$VEP.impact) & !grepl("HIGH|MODERATE", driver.muts.all$VEP.impact)
driver.muts.all$genuine.driver[vep.filter] <- F

# For mutations with VEP predictions, mark as genuine driver if impact is high or moderate
# driver.muts.all$genuine.driver[which(!is.na(driver.muts.all$VEP.impact))] <- FALSE
# driver.muts.all$genuine.driver[grep("HIGH|MODERATE", driver.muts.all$VEP.impact)] <- TRUE
# # Oncogenes with CN gain are true drivers, with loss are not:
# driver.muts.all$genuine.driver[intersect(grep("oncogene", driver.muts.all$Role.in.cancer), which(driver.muts.all$Mutation.Type == "CN amplification"))] <- TRUE
# driver.muts.all$genuine.driver[which(!is.na(driver.muts.all$Role.in.cancer) & !grepl("oncogene", driver.muts.all$Role.in.cancer) & driver.muts.all$Mutation.Type == "CN amplification")] <- FALSE
# # Tumour suppressor genes with CN loss are true drivers, with gain are not:
# driver.muts.all$genuine.driver[intersect(grep("TSG", driver.muts.all$Role.in.cancer), which(driver.muts.all$Mutation.Type == "CN deletion"))] <- TRUE
# driver.muts.all$genuine.driver[which(!is.na(driver.muts.all$Role.in.cancer) & !grepl("TSG", driver.muts.all$Role.in.cancer) & driver.muts.all$Mutation.Type == "CN deletion")] <- FALSE
#WriteXLS(driver.muts.all, ExcelFileName = "results/drivers.for.henry.xls", row.names = T)

wgs.pheno$driver.count <- unlist(lapply(wgs.pheno$name, function(x){
  length(which(driver.muts.all$Sample == x & driver.muts.all$genuine.driver))
}))

# Additionally run dndscv to identify novel drivers
# (Code in separate file - no plots created therefore not included here)

##########################################################################
# Start of plots
##########################################################################

##########################################################################
# Figure 1: study design, not generated with R
##########################################################################
# Figure 2: Genomic analysis 
##########################################################################
plot.genomic.circos(paste(results_dir, 'Fig2_circos.png', sep=''))

##########################################################################
# Figure 3: Heatmaps and PCAs comparing gene expression and methylation data
##########################################################################
plot.gxn.heatmap(paste(results_dir, 'Fig3A_gxn_heatmap.pdf', sep=""))
plot.meth.heatmap(paste(results_dir, 'Fig3B_meth_heatmap.pdf', sep=""))
plot.gxn.pca(paste(results_dir, 'Fig3C_gxn_pca.pdf', sep=""))
plot.meth.pca(paste(results_dir, 'Fig3D_meth_pca.pdf', sep=""))

##########################################################################
# Figure 4: Prediction plots for GXN and MHI
##########################################################################
plot.gxn.prediction(paste(results_dir, 'Fig4A-C_gxn_prediction.pdf', sep=""))
plot.meth.distribution(paste(results_dir, 'Fig4D_mvp_distribution.pdf', sep=""))
plot.meth.mhi(paste(results_dir, 'Fig4E_mvp_probe_counts.pdf', sep=""))
plot.meth.mhi.cvc.histo(paste(results_dir, 'Fig4F_mhi_sample_histograms.pdf', sep=""))

##########################################################################
# Figure 5: GXN Pathway analysis and CIN gene plots. Methylation output also done here.
# Note we do not make plots in this file, just Excel output for plotting in Prism.
##########################################################################

plot.gxn.cin.expression(paste(results_dir, 'Fig5A-B_CIN_mean_expression.pdf', sep=""))
WriteXLS(
  "gxn.gage.summary",
  ExcelFileName = paste(results_dir,"Fig5C_data_gxn_pathways.xlsx", sep=""), row.names = T, AdjWidth = T
)
WriteXLS(
  "meth.gage.summary",
  ExcelFileName = paste(results_dir,"Fig5D_methyl_pathways.xlsx", sep=""), row.names = T, AdjWidth = T
)

##########################################################################
# Extended data Figure 1: Mutational Signature Analysis
##########################################################################
plot.genomic.signatures(paste0(results_dir, "Ext_Fig1_MutationalSignatures"))

##########################################################################
# Extended data Figure 2: Genome-wide copy number plots
##########################################################################
plot.genomic.cna.genomewide(paste0(results_dir, 'Ext_Fig2_cna_genomewide.png'))

##########################################################################
# Extended data Figure 3
# 
# Clinical follow up data, not plotted in R
##########################################################################

##########################################################################
# Extended data Figure 4: Genomic PvR boxplots
##########################################################################
plot.genomic.pvr(paste(results_dir, "Ext_Fig4_Genomic_PVR_plots.pdf", sep=""))

##########################################################################
# Extended data Figure 5: Clonality
#
# Plots are stored in results_dir/clonality by the above code.
# Here we plot an additional matrix to show shared mutations between samples from the same patient
##########################################################################
plot.genomic.clonality.matrix(paste(results_dir, 'Ext_Fig5_multisample_clonality_matrix.pdf', sep=''))

##########################################################################
# Extended data Figure 6: PvR Circos plot
##########################################################################
plot.pvr.circos(paste(results_dir, 'Ext_Fig6_pvr_circos.png', sep=''))

##########################################################################
# Extended data Figure 7: Extended PCA plots
##########################################################################
plot.meth.pcas(paste(results_dir, 'Ext_Fig7A-F_methylation_PCAs.pdf', sep=""))
plot.gxn.pcas(paste(results_dir, 'Ext_Fig7G-K_gxn_PCAs.pdf', sep=""))


##########################################################################
# Extended data Figure 8: Methylation and CNA predictive models
##########################################################################

plot.meth.prediction(paste(results_dir, 'Ext_Fig8A-C_meth_prediction.pdf', sep=""))
plot.mcna.prediction(paste(results_dir, 'Ext_Fig8D-F_cna_prediction.pdf', sep=""))

##########################################################################
# Extended data Figure 9: MHI Heatmap for PvR data
##########################################################################

plot.meth.mhi.pvr.histo(paste(results_dir, 'Ext_Fig9_mhi_sample_histograms_pvr.pdf', sep=""))

##########################################################################
# Extended data Figure 10: Correlation of wGII with CIN gene expression
##########################################################################

plot.cin.gxn.cor(paste0(results_dir, "Ext_Fig10_cin_gxn_cor.pdf"))


##########################################################################
# Table 1: Demographic data
##########################################################################
plot.demographic.table(paste0(results_dir, "Table1_demographics.xls"))

##########################################################################
# Supplementary data file 1: list of all driver mutations
##########################################################################
WriteXLS(driver.muts, ExcelFileName = paste0(results_dir,"Sup_Data1.driver_mutations.xlsx"), AdjWidth = T)

##########################################################################
# Supplementary data file 2: Differential expression of GXN, methylation and CNAs
##########################################################################
gdiff.sig <- gdiff[which(gdiff$fdr < 0.01),]
# Include both DMPs and DMRs
mdiff.sig <- mdiff[which(mdiff$adj.P.Val < 0.01 & abs(mdiff$deltaBeta) > 0.3),]
dmrs.sig <- dmrs$ProbeLassoDMR[which(dmrs$ProbeLassoDMR$dmrP < 0.01),]
cdiff.sig <- cdiff[which(cdiff$fdr < 0.01),]
WriteXLS(c("gdiff.sig", "mdiff.sig", "dmrs.sig", "cdiff.sig"), ExcelFileName = paste0(results_dir, "Sup_Data2_Differentially_Expressed_Genes.xlsx"), AdjWidth = T,
         SheetNames = c("DE Genes", "DMPs", "DMRs", "CN bands"), row.names = T)

# Identify consistent genes - increased expression and hypomethylation or decreased expression and hypermethylation
genes.consistent <- c(
  intersect(rownames(gdiff.sig)[gdiff.sig$fc > 1], mdiff.sig$gene[mdiff.sig$deltaBeta < 0]),
  intersect(rownames(gdiff.sig)[gdiff.sig$fc < 1], mdiff.sig$gene[mdiff.sig$deltaBeta > 0])
)

##########################################################################
# Supplementary data file 3: Biological and experimental details of all samples
# Not produced in R
##########################################################################

##########################################################################
# Supplementary data file 4: Pathways implicated in progressive vs regressive analysis
##########################################################################
WriteXLS(
  c("gxn.gage.summary", "meth.gage.summary"),
  ExcelFileName = paste(results_dir,"Sup_Data4_pathways.xlsx", sep=""), row.names = T, AdjWidth = T
)