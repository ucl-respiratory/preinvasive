########################################################################################################
# New genomic analysis follows:
########################################################################################################
# Comparison of progressive vs regressive samples for genomic data:
# Subs, incels, CNAs, mutation burden (subs+indels), rearrangements, coding substitutions
plot.genomic.pvr <- function(filename){
  
  if(!exists("wgs.pheno")){
    stop("Error: Please run loadWgsData first")
  }
  source('utility_functions/compare.fn.R')
  
  plotdata <- wgs.pheno
  
  # Total mutations (subs+indels+rearrangements):
  plotdata$burden <- as.numeric(table(muts.all$patient)[wgs.pheno$name])
  if(length(which(is.na(plotdata$burden))) > 0){ plotdata$burden[which(is.na(plotdata$burden))] <- 0 }
  
  # Coding mutations
  plotdata$burden.coding <- as.numeric(table(muts.coding$patient)[wgs.pheno$name])
  if(length(which(is.na(plotdata$burden.coding))) > 0){ plotdata$burden.coding[which(is.na(plotdata$burden.coding))] <- 0 }
  
  # Driver mutations
  #plotdata$driver.count <- as.numeric(table(muts.coding$patient[which(muts.coding$gene %in% driver.genes)])[wgs.pheno$name])
  plotdata$driver.count <- as.numeric(table(driver.muts.all$Sample[which(driver.muts.all$genuine.driver)])[wgs.pheno$name])
  if(length(which(is.na(plotdata$driver.count))) > 0){ plotdata$driver.count[which(is.na(plotdata$driver.count))] <- 0 }
  
  # Number of CN changes:
  plotdata$cna.counts <- unlist(lapply(plotdata$name, function(x){
    length(which(cna.summary.list$sample == x))
  }))
  # Number of genes affected by a CN change:
  plotdata$cna.counts <- as.numeric(apply(cnas.genes.summary[,wgs.pheno$name], 2, function(x){ sum(abs(x), na.rm=T)}))
  
  # Subs only:
  plotdata$subs.count <- as.numeric(table(muts.all$patient[which(muts.all$class == "SNV")])[wgs.pheno$name])
  if(length(which(is.na(plotdata$subs.count))) > 0){ plotdata$subs.count[which(is.na(plotdata$subs.count))] <- 0 }
  # Indels only:
  plotdata$indels.count <- as.numeric(table(muts.all$patient[which(muts.all$class %in% c("D", "I", "DI"))])[wgs.pheno$name])
  if(length(which(is.na(plotdata$indels.count))) > 0){ plotdata$indels.count[which(is.na(plotdata$indels.count))] <- 0 }
  
  # Rearrangements
  plotdata$rearr.count <- as.numeric(table(muts.all$patient[which(muts.all$class %in% c("R"))])[wgs.pheno$name])
  if(length(which(is.na(plotdata$rearr.count))) > 0){ plotdata$rearr.count[which(is.na(plotdata$rearr.count))] <- 0 }
  
  # Telomere lengths
  telomeres <- rbind(
    read.csv('resources/telomere_lengths_batch1.csv')[,1:2],
    read.csv('resources/telomere_lengths_batch2.csv')[,1:2]
  )
  plotdata$telomeres <- as.numeric(telomeres$Length)[match(plotdata$name, telomeres$Sample)]
  
  # # Number of genes affected by a CN change:
  plotdata$cna.gene.counts <- as.numeric(apply(cnas.genes.summary[,wgs.pheno$name], 2, function(x){ sum(abs(x), na.rm=T)}))
  
  # Proportion of mutations which are clonal
  plotdata$prop.clonal <- 100*plotdata$clonal.muts / (plotdata$clonal.muts + plotdata$subclonal.muts)
  
  
  # For this analysis, remove the 'query' regressives which later turned out to progress
  plotdata <- plotdata[-which(plotdata$query.reg == 1),]
  
  pdf(filename)
  
  compare.fn(dependent_variable = "subs.count", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"),
             random_effects = "Patient", modelinfo=plotdata, title="Substitutions")
  compare.fn(dependent_variable = "indels.count", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"),
             random_effects = "Patient", modelinfo=plotdata, title="Small Insertions and Deletions")
  compare.fn(dependent_variable = "rearr.count", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"),
             random_effects = "Patient", modelinfo=plotdata, title="Rearrangements")
  compare.fn(dependent_variable = "cna.counts", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"), 
             random_effects = "Patient", modelinfo=plotdata, title="Copy number segments")
  # Clonality data - generated in full_analysis.R
  compare.fn(dependent_variable = "clonal.muts", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"), 
             random_effects = "Patient", modelinfo=plotdata, strip.method="jitter", title="Clonal Substitutions")
  compare.fn(dependent_variable = "prop.clonal", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"), 
             random_effects = "Patient", modelinfo=plotdata, strip.method="jitter", title="Proportion of clonal substitutions")
  compare.fn(dependent_variable = "nclusters", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"), 
             random_effects = "Patient", modelinfo=plotdata, strip.method="jitter", title="Number of clones")
  compare.fn(dependent_variable = "driver.count", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"),
             random_effects = "Patient", modelinfo=plotdata, title="Putative driver mutations")
  compare.fn(dependent_variable = "telomeres", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"), 
             random_effects = "Patient", modelinfo=plotdata, strip.method="jitter", title="Telomere length")
  
  
  
  # Other comparisons not used in final paper
  compare.fn(dependent_variable = "burden", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"),
             random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "burden.coding", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"),
             random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "cna.gene.counts", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"), 
             random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "wgii", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years", "purity"), 
             random_effects = "Patient", modelinfo=plotdata, strip.method="jitter")
 dev.off() 
  
}