##########################################################################
# Figure 3A: GXN heatmap
##########################################################################
plot.gxn.heatmap <- function(filename){
  if(!exists('gdata') | !exists('gdiff')){
    stop("ERROR: Please run loadGxnData first and make sure gdiff has been calculated")
  }
  
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
  
  pdf(filename)
  pheatmap(removeOutliers(sig_genes), cluster_rows=T, cluster_cols=T, scale="row", main=paste("Gene Expression (",dim(sig_genes)[1]," genes)", sep=""),
           annotation_col=g.annot, treeheight_row=0, treeheight_col=0, show_rownames=F, show_colnames=F,
           annotation_colors=g.annot_colors,
           color=hmcol, legend=F)
  dev.off()
}
