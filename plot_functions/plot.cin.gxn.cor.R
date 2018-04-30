plot.cin.gxn.cor <- function(filename){
  
  if(!exists("overlap.pheno") | !exists("gdata") | !exists("wgs.pheno")){
    stop("ERROR: function depends on overlap.pheno, gdata, wgs.pheno as defined in full_analysis.R")
  }
  
  sel <- which(overlap.pheno$Gene.expression == "YES" & overlap.pheno$Whole.Genome.Sequencing == "YES")
  cin.expr <- apply(gdata[cin_genes[which(cin_genes %in% rownames(gdata))],as.character(overlap.pheno$Sample.Number..GXN.[sel])], 2, mean)
  cin.wgii <- wgs.pheno$wgii[match(overlap.pheno$Sample.Number..WGS.[sel], wgs.pheno$name)]
  
  pdf(filename)
  plot(cin.expr, cin.wgii, main="Correlation of CIN gene expression with wGII",
       xlab="Mean CIN gene expression", ylab="wGII",
       col=c("green", "red")[wgs.pheno$progression[match(overlap.pheno$Sample.Number..WGS.[sel], wgs.pheno$name)]+1])
  c <- cor(cin.expr, cin.wgii)
  if(show.legends){
    legend('bottomright', paste0("r^2=", signif(c, 3)))
  }
  dev.off()
}