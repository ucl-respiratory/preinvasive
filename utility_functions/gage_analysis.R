#' Run gage analysis on a differential expression analysis
#'
#' @import gage
#' @import gageData
#' @import org.Hs.eg.db
#' @import annotate
#' @param data Matrix of genes as rows and samples as columns
#' @param pheno Vector of 1/0 values determining the case/control groups
#' @param source One of 'GO' or 'KEGG' or MSigDB pathways: 'hallmark', 'onc', 'msig.go', 'msig.kegg'
#' @param compare Paired or unpaired data. See ?gage for details.
#' @export

library(gage)
library(gageData)
library(org.Hs.eg.db)
library(annotate)

gage_analysis <- function(genedata, pheno, source="GO", compare='unpaired'){
  
  # Use custom gene sets - all KEGG gene sets + CIN genes
  load("resources/c2.kegg.gsets.RData")
  gsets <- c2.kegg.gsets
  
  mygage <- gage(
    genedata, gsets=gsets,
    ref=which(pheno == 0),
    samp=which(pheno == 1),
    compare=compare
  )
  
  return(mygage)
}
