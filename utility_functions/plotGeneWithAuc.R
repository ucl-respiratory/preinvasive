# Function to plot a gene with AUC. If multiple genes are given, uses the mean.
plotGeneWithAuc <- function(genes, title=NA, legend.pos="topleft"){
  o <- order(tcga.gpheno.all$dose)
  if(is.na(title)){
    title <- paste("Mean expression of", paste(genes, collapse=", "))
  }
  sel.na <- which(!(genes %in% rownames(tcga.gdata.all)))
  if(length(sel.na) > 0){ genes <- genes[-sel.na] }
  
  index <- 1:length(o)
  plot(index, apply(tcga.gdata.all[genes,o], 2, mean),
       main=title,
       xlab="Sample",
       ylab="Gene Expression")
  cols.stage <- c('black', "darkgreen", "green", "red", "orange")
  for(i in -1:3){
    sel <- which(tcga.gpheno.all[o,]$dose == i)
    points(index[sel], apply(tcga.gdata.all[genes,o][,sel], 2, mean), col=cols.stage[i+2])
  }
  
  roc.cvc <- roc(
    predictor=apply(tcga.gdata.all[genes,which(tcga.gpheno.all$dose == 0 | tcga.gpheno.all$dose == 3)], 2, mean),
    response=tcga.gpheno.all$dose[which(tcga.gpheno.all$dose == 0 | tcga.gpheno.all$dose == 3)]
  )
  roc.pvr <- roc(
    predictor=apply(tcga.gdata.all[genes,which(tcga.gpheno.all$dose == 1 | tcga.gpheno.all$dose == 2)], 2, mean),
    response=tcga.gpheno.all$dose[which(tcga.gpheno.all$dose == 1 | tcga.gpheno.all$dose == 2)]
  )
  if(show.legends){
    legend(
      legend.pos,
      c(paste("Ca vs Control AUC", signif(auc(roc.cvc), 3)), paste("CIS Prog vs Reg AUC", signif(auc(roc.pvr), 3)))
    )
  }
}