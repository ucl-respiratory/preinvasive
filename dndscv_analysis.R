# dndscv analysis for potential drivers in CIS data
# Following instructions from here: http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/dNdScv.html
library("seqinr")
library("Biostrings")
library("MASS")
library("GenomicRanges")
library("dndscv")

source('data_loaders/loadWgsData.R')

# Use all mutations - subs and indels
dndscv.input <- muts.all[,c("patient", "chr", "start", "ref", "alt")]
colnames(dndscv.input) <- c("sampleID", "chr", "pos", "ref", "mut")

# Make sure we remove mutations duplicated in samples from the same patient
dndscv.input$uuid <- paste(substr(dndscv.input$sampleID, 1, 7), dndscv.input$chr, dndscv.input$pos, dndscv.input$ref, dndscv.input$mut, sep="-")
dndscv.input <- dndscv.input[-which(duplicated(dndscv.input$uuid)),]
dndscv.input$uuid <- NULL

# Run the analysis
dndsout = dndscv(dndscv.input)

# Check global dN/dS - should be around 1
print(dndsout$globaldnds)

# Find significant genes
sel_cv = dndsout$sel_cv
signif_genes = sel_cv[sel_cv$qglobal_cv<0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes) = NULL
print(signif_genes)

# Result - only TP53 is identified.