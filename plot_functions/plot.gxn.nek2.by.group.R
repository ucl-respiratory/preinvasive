########################################################################################################
# Figure 5C - NEK2 expression by group
########################################################################################################
plot.gxn.nek2.by.group <- function(filename){
  pdf(filename)
  plotGeneByGroup("NEK2", legend.pos='topleft')
  plotGeneWithAuc("NEK2")
  dev.off()
}
