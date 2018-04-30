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
  # Total mutations (subs+indels):
  plotdata$burden <- as.numeric(table(muts.all$patient)[wgs.pheno$name])
  if(length(which(is.na(plotdata$burden))) > 0){ plotdata$burden[which(is.na(plotdata$burden))] <- 0 }
  # Unfiltered burden
  #plotdata$burden.unfiltered <- as.numeric(table(muts.unfiltered$patient)[wgs.pheno$name])
  #if(length(which(is.na(plotdata$burden.unfiltered))) > 0){ plotdata$burden.unfiltered[sel.na] <- 0 }
  # Coding mutations
  plotdata$burden.coding <- as.numeric(table(muts.coding$patient)[wgs.pheno$name])
  if(length(which(is.na(plotdata$burden.coding))) > 0){ plotdata$burden.coding[which(is.na(plotdata$burden.coding))] <- 0 }
  # Number of genes affected by a CN change:
  plotdata$cna.gene.counts <- as.numeric(apply(cnas.genes.summary[,wgs.pheno$name], 2, function(x){ sum(abs(x), na.rm=T)}))
  
  # Subs only:
  plotdata$subs.count <- as.numeric(table(muts.all$patient[which(muts.all$class == "SNV")])[wgs.pheno$name])
  if(length(which(is.na(plotdata$subs.count))) > 0){ plotdata$subs.count[which(is.na(plotdata$subs.count))] <- 0 }
  # Indels only:
  plotdata$indels.count <- as.numeric(table(muts.all$patient[which(muts.all$class %in% c("D", "I", "DI"))])[wgs.pheno$name])
  if(length(which(is.na(plotdata$indels.count))) > 0){ plotdata$indels.count[which(is.na(plotdata$indels.count))] <- 0 }
  # Rearrangements - load annotated rearrangements
  rearr <- read.table('resources/rearrangements_annotated', sep="\t", header = T)
  
  #plotdata$rearr.count <- as.numeric(table(rearr.all$patient)[wgs.pheno$name])
  plotdata$rearr.count <- as.numeric(table(rearr$SAMPLE)[wgs.pheno$name])
  if(length(which(is.na(plotdata$rearr.count))) > 0){ plotdata$rearr.count[which(is.na(plotdata$rearr.count))] <- 0 }
  
  # Telomere lengths
  telomeres <- read.csv('resources/telomere_lengths.csv')
  plotdata$telomeres <- as.numeric(telomeres$Length)[match(plotdata$name, telomeres$Sample)]
  
  # # Total subs:
  # plotdata$subcounts <- as.numeric(table(subs.all$patient)[wgs.pheno$name])
  # sel.na <- which(is.na(plotdata$subcounts))
  # if(length(sel.na) > 0){ plotdata$subcounts[sel.na] <- 0 }
  # #plotdata$subcounts <- as.numeric(subs.burden[plotdata$name])
  # # Total indels:
  # plotdata$indelcounts <- as.numeric(table(indels.all$patient)[wgs.pheno$name])
  # sel.na <- which(is.na(plotdata$indelcounts))
  # if(length(sel.na) > 0){ plotdata$indelcounts[sel.na] <- 0 }
  # #plotdata$indelcounts <- as.numeric(indels.burden[plotdata$name])
  # # Number of genes affected by a CN change:
  # plotdata$cna.gene.counts <- as.numeric(apply(cnas.genes.summary[,wgs.pheno$name], 2, function(x){ sum(abs(x), na.rm=T)}))
  # 
  
  
  
  pdf(filename)
  # compare.fn(dependent_variable = "subcounts", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"), 
  #            random_effects = "Patient", modelinfo=plotdata)
  # compare.fn(dependent_variable = "indelcounts", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"), 
  #            random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "burden", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"),
             random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "burden.coding", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"),
             random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "subs.count", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"),
             random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "indels.count", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"),
             random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "rearr.count", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"),
             random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "driver.count", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"),
             random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "cna.gene.counts", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"), 
             random_effects = "Patient", modelinfo=plotdata)
  compare.fn(dependent_variable = "wgii", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"), 
             random_effects = "Patient", modelinfo=plotdata, strip.method="jitter")
  compare.fn(dependent_variable = "telomeres", compared_variable = "progression", fixed_effects = c("Age.at.specimen.profiled", "Pack.years"), 
             random_effects = "Patient", modelinfo=plotdata, strip.method="jitter")
 dev.off() 
  
}