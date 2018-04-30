########################################################################################################
# Methylation of NEK2-associated probe by group
########################################################################################################
plot.meth.nek2 <- function(filename){
  
  pdf(filename)
  plotMethylationByGroup("cg17931972", "NEK2", legend.pos="topleft")
  dev.off()
  
}