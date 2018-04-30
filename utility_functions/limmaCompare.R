#' limmaCompare: Compare two data sets using limma
#' Assumes pheno has a column of 0 and 1 values to identify groups
#' @param data Data matrix with genes as rows and samples as columns
#' @param pheno Pheno - assumes there is a column with 0/1 values to identify two groups. By default this is progression.
#' @param colname Default is progression. Identifies the relevant pheno column.
#' @param fdr_limit Default 0.01. Cutoff for 'significance'.
#' @import limma
#' @import gtools
#' @export

limmaCompare <- function(data, pheno, colname='progression', fdr_limit=0.01){
  
  # Quick checks:
  # Data and pheno should have the same dimensions
  if(!(dim(data)[2] == dim(pheno)[1])){
    print("ERROR: different dimensions in data and pheno")
  }
  # pheno should have a 'colname' column
  if(!(colname %in% colnames(pheno))){
    print(paste("ERROR: supplied pheno data frame has no column '", colname, "'", sep=""))
  }
  
  AAprog <- factor(pheno[,colname])
  design <- model.matrix(~AAprog )
  
  # Fit a linear model
  fit <- lmFit(data, design)
  fit2 <- eBayes(fit)
  # Select the last column of the P values -> p value for progression group
  colId <- length(colnames(fit2))
  p<-fit2$p.value[, colId]
  # Adjust for multiple testing
  fdr<-p.adjust(p,method="BH")
  # Calculate fold change data
  fc<-logratio2foldchange(fit2$coef[,colId])
  
  # Identify significant genes
  sel <- which(fdr < fdr_limit)
  sig_genes <- data[sel,]
  uvv <- data.frame(row.names=rownames(sig_genes), fc=fc[sel], fdr=fdr[sel], t=fit2$t[sel,colId])
  # Sort by t
  uvv <- uvv[order(-abs(uvv$t)),]
  
  return(uvv)
}

# Depends on this function:
logratio2foldchange<-
  function (logratio, base = 2)
  {
    retval <- base^(logratio)
    retval <- ifelse(retval < 1, -1/retval, retval)
    retval
  }
