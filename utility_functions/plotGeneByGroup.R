# Function to plot a gene over different 'dose' groups
plotGeneByGroup <- function(gene, legend.pos='topright'){
  
  colorder <- c( "darkgreen", "green", "red", "orange")
  stages <- c("TCGA control", "Regressive CIS", "Progressive CIS", "TCGA SqCC")
  
  bpdata <- tcga.gpheno.all
  bpdata$gene <- as.numeric(tcga.gdata.all[gene,])
  bpdata$stage <- factor(stages[tcga.gpheno.all$dose + 1], levels=stages)
  bplot <- ggplot(bpdata, aes(x=factor(stage),y=gene, fill=factor(stage))) +
    geom_boxplot() +
    geom_signif(comparisons=list(
      c(stages[1], stages[2]), c(stages[2], stages[3]), c(stages[3], stages[4])
    ), map_signif_level = TRUE) +
    ggtitle(paste(gene, "Expression")) + 
    scale_fill_manual(values=colorder) + 
    theme(axis.title.x=element_blank()) +
    theme(legend.position="none") +
    scale_y_continuous(name=gene)
  print(bplot)
  
}