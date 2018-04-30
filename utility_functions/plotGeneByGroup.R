# Function to plot a gene over different 'dose' groups
plotGeneByGroup <- function(gene, legend.pos='topright'){
  
  colorder <- c("green", "red", "orange")
  stages <- c("Regressive CIS", "Progressive CIS", "TCGA SqCC")
  
  sel.noctl <- which(tcga.gpheno.all$dose > 0)
  bpdata <- tcga.gpheno.all[sel.noctl,]
  bpdata$gene <- as.numeric(tcga.gdata.all[gene,sel.noctl])
  bpdata$stage <- factor(stages[tcga.gpheno.all$dose[sel.noctl]], levels=stages)
  bplot <- ggplot(bpdata, aes(x=factor(stage),y=gene, fill=factor(stage))) +
    geom_boxplot() +
    geom_signif(comparisons=list(
      c(stages[1], stages[2]), c(stages[2], stages[3])
    ), map_signif_level = TRUE) +
    ggtitle(paste(gene, "Expression")) + 
    scale_fill_manual(values=colorder) + 
    theme(axis.title.x=element_blank()) +
    theme(legend.position="none") +
    scale_y_continuous(name=gene)
  print(bplot)
  
}