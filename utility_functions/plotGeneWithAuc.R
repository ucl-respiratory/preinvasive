# Function to plot a gene with AUC. If multiple genes are given, uses the mean.
plotGeneWithAuc <- function(genes, title=NA, legend.pos="topleft"){
  plot.data <- gdata
  plot.pheno <- gpheno
  o <- order(plot.pheno$progression)
  
  if(is.na(title)){
    title <- paste("Mean expression of", paste(genes, collapse=", "))
  }
  sel.na <- which(!(genes %in% rownames(plot.data)))
  if(length(sel.na) > 0){ genes <- genes[-sel.na] }
  
  index <- 1:length(o)
  plot(
    index, apply(plot.data[genes, o], 2, mean),
    main=title,
    xlab="Sample",
    ylab="Gene Expression",
    col=c('green', 'red')[plot.pheno$progression[o]+1]
  )
  
  roc.pvr <- roc(
    predictor=apply(plot.data[genes,o], 2, mean),
    response=plot.pheno$progression[o]
  )
  if(show.legends){
    legend(
      legend.pos,
      paste("CIS Prog vs Reg AUC", signif(auc(roc.pvr), 3))
    )
  }
  
  # sel.noctl <- which(tcga.gpheno.all$dose > 0)
  # plot.data <- tcga.gdata.all[,sel.noctl]
  # plot.pheno <- tcga.gpheno.all[sel.noctl,]
  # o <- order(plot.pheno$dose[sel.noctl])
  # 
  # if(is.na(title)){
  #   title <- paste("Mean expression of", paste(genes, collapse=", "))
  # }
  # sel.na <- which(!(genes %in% rownames(plot.data)))
  # if(length(sel.na) > 0){ genes <- genes[-sel.na] }
  # 
  # index <- 1:length(o)
  # plot(index, apply(plot.data[genes,o], 2, mean),
  #      main=title,
  #      xlab="Sample",
  #      ylab="Gene Expression")
  # cols.stage <- c('black', "darkgreen", "green", "red", "orange")
  # for(i in -1:3){
  #   sel <- which(plot.pheno[o,]$dose == i)
  #   points(index[sel], apply(plot.data[genes,o][,sel], 2, mean), col=cols.stage[i+2])
  # }
  # 
  # # roc.cvc <- roc(
  # #   predictor=apply(plot.data[genes,which(plot.pheno$dose == 0 | plot.pheno$dose == 3)], 2, mean),
  # #   response=plot.pheno$dose[which(plot.pheno$dose == 0 | plot.pheno$dose == 3)]
  # # )
  # roc.pvr <- roc(
  #   predictor=apply(plot.data[genes,which(plot.pheno$dose == 1 | plot.pheno$dose == 2)], 2, mean),
  #   response=plot.pheno$dose[which(plot.pheno$dose == 1 | plot.pheno$dose == 2)]
  # )
  # if(show.legends){
  #   legend(
  #     legend.pos,
  #     paste("CIS Prog vs Reg AUC", signif(auc(roc.pvr), 3))
  #   )
  # }
}