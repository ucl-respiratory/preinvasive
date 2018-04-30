#' removeOutliersMatrix
#' 
#' Remove outlier values from a matrix
#' @param X Matrix to remove outliers from
#' @export
removeOutliersMatrix<-function(X){
  
  num<-X
  for(i in 1:dim(num)[1]) {
    x<-num[i,]
    trim = 0.05
    lo = quantile(x, trim)
    hi = quantile(x, 1 - trim)
    x[x < lo] = lo
    x[x > hi] = hi
    num[i,]<-x
  }
  num;
}

removeOutliersVector<-function(X){
  
  num<-X
  trim = 0.05
  lo = quantile(X, trim)
  hi = quantile(X, 1 - trim)
  num[num < lo] = lo
  num[num > hi] = hi
  num;
}

#' removeOutliersDF
#' 
#' Remove outlier values from a data frame
#' @param X Data frame to remove outliers from
#' @export
removeOutliersDF<-function(X){
  
  num<-X
  for(i in 1:dim(num)[1]) {
    x<-num[i,]
    trim = 0.05
    lo = quantile(x, trim, na.rm=T)
    hi = quantile(x, 1 - trim, na.rm=T)
    x[which(x < as.numeric(lo))] = lo
    x[which(x > as.numeric(hi))] = hi
    num[i,]<-x
  }
  num;
}

#' removeOutliers
#' 
#' Utility functions to remove outliers from matrices and data frames. This is useful for heatmaps -> improves contrast
#' @param X Matrix or data frame to remove outliers from
#' @export
removeOutliers <- function(X){
  if(is.data.frame(X)){
    removeOutliersDF(X)
  }else{
    if(is.matrix(X)){
      removeOutliersMatrix(X)
    }else{
      if(is.double(X)){
        removeOutliersVector(X)
      }
    }
  }
}