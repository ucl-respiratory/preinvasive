#' plotPCA: Utility function to plot a PCA from progression/regression data with our default options
#'
#' @param data Matrix with rows as genes, columns as samples
#' @param progression Numerical array of group indices
#' @param labels If required use a legend
#' @export
#' @examples
#' plotPCA(t(gdata), gpheno$progression+1, labs=c('regressive', 'progressive'))

plotPCA <- function(data, group, cols=c('green', 'red', 'blue', 'black', 'yellow'), labs=c(), labpos = 'bottom', labhoriz=T, title="PCA plot"){
  pc <- prcomp(data)
  gcol <- cols[group]
  plot(pc$x[, 1], pc$x[, 2], col = gcol, main = title, xlab = "PC1", ylab = "PC2")
  
  if(length(labs) > 0){
    legend(labpos, col=cols[1:length(labs)],
           legend=labs, pch=1, horiz=labhoriz,
           bty="n")
  }
}
