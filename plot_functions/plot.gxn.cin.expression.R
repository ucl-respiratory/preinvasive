########################################################################################################
# CIN-related gene expression
########################################################################################################
plot.gxn.cin.expression <- function(filename){
  
  pdf(filename)
  plotGeneWithAuc(cin_genes[which(cin_genes %in% rownames(gdiff)[which(gdiff$fdr < 0.01)])], 
                  title = "Mean expression of significant CIN genes")
  # Plot NEK2 alone
  plotGeneWithAuc("NEK2")
  dev.off()
}
