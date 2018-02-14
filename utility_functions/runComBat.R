#' runComBat: utility function to run ComBat on array data.
#'
#' Run ComBat on a new data set to bring it into line with a reference data set
#' @param refdata Your reference data. Will be reduced to only values which overlap with the new data, but otherwise unchanged.
#' @param newdata New data to be brought into line with the reference data
#' @param doPCA Default FALSE. If TRUE will plot a PCA of the old and new data.
#' @param par.prior Default TRUE. Do parametric adjustment.
#' @importFrom sva ComBat
#' @return Array of: [reduced set of refdata, corrected newdata]
#' @export

library(sva)
runComBat <- function(refdata, newdata, doPCA=F, par.prior=T, na.rm=T){
  if(na.rm){
    # Only include complete rows
    sel.na <- which(apply(newdata, 1, function(x){any(is.na(x))}))
    if(length(sel.na) > 0){ newdata <- newdata[-sel.na,] }
  }
  
  
  # Reduce the data sets to the shared probes and combine
  shared <- intersect(rownames(refdata), rownames(newdata))
  refdata.shared <- refdata[shared,]
  newdata.shared <- newdata[shared,]
  alldata <- cbind(refdata.shared, newdata.shared)
  
  # Make a pheno data frame
  pheno <- matrix(data='old', ncol=1, nrow=dim(alldata)[2])
  pheno[(dim(refdata)[2]+1):(dim(pheno)[1]),] <- 'new'
  pheno <- data.frame(pheno)
  colnames(pheno) <- c("group")
  
  # Run ComBat
  modcombat = model.matrix(~1, data=pheno)
  combat_data = ComBat(dat=as.matrix(alldata), batch=as.factor(pheno$group), mod=modcombat, ref.batch="old", par.prior = par.prior)
  
  # Check with a PCA (optional)
  if(doPCA){
    plotPCA(t(combat_data), as.factor(pheno$group))
  }
  
  # Split out the adjusted new data set and return
  newdata.corrected <- combat_data[,which(pheno$group == 'new')]
  
  # Returns the reduced set of reference data and the corrected new data
  return(list(refdata.shared, newdata.corrected))
}
