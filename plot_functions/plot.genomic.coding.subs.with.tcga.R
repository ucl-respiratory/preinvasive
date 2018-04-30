plot.genomic.coding.subs.with.tcga <- function(filename){
  
  if(!exists("muts.coding.counts")){
    stop("ERROR: Please run loadWgsPheno")
  }
  
  # Plot mutation burden TCGA vs CIS
  o <- order(muts.coding.counts)
  df <- rbind(
    data.frame(
      source="TCGA", col='black', subs=sort(as.numeric(table(tcga.snvs$Tumor_Sample_UUID))), name=names(table(tcga.snvs$Tumor_Sample_UUID))
    ),
    data.frame(
      source="CIS", col=c('green','red')[wgs.pheno$progression[o]+1], subs=muts.coding.counts[o], name=wgs.pheno$name[o]
    )
  )
  # Colour the query regressives orange
  df$col <- as.character(df$col)
  df$col[match(wgs.pheno$name[which(wgs.pheno$query.reg==1)], as.character(df$name))] <- 'orange'
  
  pdf(filename)
  plot(df$subs, col=df$col, log="y", ylab="Coding Substitutions")
  dev.off()
  
}