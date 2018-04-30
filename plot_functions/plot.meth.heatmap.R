##########################################################################
# Methylation heatmap
##########################################################################
plot.meth.heatmap <- function(filename){
  if(!exists("mdata") | !exists("mdiff")){
    stop("ERROR: data not loaded. Please run loadMethData and make sure mdiff exists.")
  }
  
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
  
  pdf(filename)
  pheatmap(removeOutliers(mdata.sig), cluster_rows=T, cluster_cols=T, scale="row", main=paste("Methylation (Top ",dim(mdata.sig)[1]," MVPs)", sep=""),
           annotation_col=m.annot, treeheight_row=0, treeheight_col=0, show_rownames=F, show_colnames=F,
           annotation_colors=m.annot_colors,
           color=hmcol2, legend=F)
  dev.off()
}

