plot.mcna.heatmap <- function(filename){
  
  if(!exists("cdiff")){
    stop("ERROR: cdiff not found. Please run as part of full_analysis.R.")
  }
  
  sig_bands <- mcnas.band[rownames(cdiff[which(cdiff$fdr < 0.01),]),]
  
  c.annot <- data.frame(
    pack.years=mcnas.pheno$smoking_group,
    age.group=mcnas.pheno$age_group,
    gender=mcnas.pheno$Gender,
    COPD=mcnas.pheno$COPD,
    status=c("Regressive", "Progressive")[mcnas.pheno$progression+1]
  )
  c.annot_colors <- list(
    status=c(Progressive="red", Regressive="green"),
    pack.years=smoking_group_names,
    age.group=age_group_names,
    gender=c("F"="pink", "M"="blue"),
    COPD=c("NO"="green", "YES"="red")
  )
  rownames(c.annot) <- colnames(sig_bands)
  
  pdf(filename)
  pheatmap(removeOutliers(sig_bands), cluster_rows=T, cluster_cols=T, scale="row", main=paste("Copy number top (",dim(sig_bands)[1]," bands)", sep=""),
           annotation_col=c.annot, treeheight_row=0, treeheight_col = F, show_rownames=F, show_colnames=F,
           annotation_colors=c.annot_colors,
           color=hmcol, legend=T)
  dev.off()
}