annotate.drivers <- function(driver.muts.all){
  # Annotate with CGC data
  driver.muts.all$Role.in.Cancer <- driver.genes.info$Role.in.Cancer[match(driver.muts.all$Gene, driver.genes.info$Gene.Symbol)]
  driver.muts.all$CGC.Tier <- driver.genes.info$Tier[match(driver.muts.all$Gene, driver.genes.info$Gene.Symbol)]
  driver.muts.all$Validated.Mut.Types <- driver.genes.info$Mutation.Types[match(driver.muts.all$Gene, driver.genes.info$Gene.Symbol)]
  driver.muts.all$Validated.Transloc.Partner <- driver.genes.info$Translocation.Partner[match(driver.muts.all$Gene, driver.genes.info$Gene.Symbol)]
  
  # Add VEP impact data (exclude CN changes from this)
  sel <- which(!is.na(driver.muts.all$Reference.Allele))
  vep.muts <- vep.annotate(driver.muts.all[sel,])
  driver.muts.all$VEP.impact <- NA
  driver.muts.all$VEP.impact[sel] <- vep.muts$vep.impact
  
  # Check against COSMIC data, count occurrences of each mutation
  # We use all lung carcinoma here as classifications are not perfect (e.g. LUSC may be categorised as either NSCLC or LUSC)
  cosmic.muts <- read.csv("resources/cosmic_alllungca_2018-07-18.csv", stringsAsFactors = F)
  # Treat all transcripts for a given gene as mutations in that gene
  cosmic.muts$Gene.Name <- unlist(lapply(cosmic.muts$Gene.Name, function(x){
    unlist(strsplit(x, "_"))[[1]]
  }))
  driver.muts.all$cosmic.count <- unlist(lapply(1:dim(driver.muts.all)[1], function(i){
    cds <- as.character(driver.muts.all$CDS.Mutation[i])
    if(is.na(cds)){
      return(NA)
    }
    gene <- as.character(driver.muts.all$Gene[i])
    return(length(which(cosmic.muts$Gene.Name == gene & cosmic.muts$CDS.Mutation == cds)))
  }))
  
  # Add a flag for genuine drivers:
  driver.muts.all$genuine.driver <- NA
  
  # Check each gene in turn:
  map <- list(
    "missense"="Mis",
    "nonsense"="N",
    "splice_region"="S",
    "frameshift"="F",
    "Rearrangement"="T",
    "CN amplification"="A",
    "ess_splice"="S",
    "inframe"="F",
    "CN deletion"="D",
    "start_lost"="N"
  )
  for(i in 1:dim(driver.muts.all)[1]){
    # Check the mutation type against CGC validated mutation types
    # Types are: A, D, Mis, N, F, T, S, O
    mtype = as.character(driver.muts.all$Mutation.Type[i])
    mtypes = unlist(strsplit(driver.muts.all$Validated.Mut.Types[i], ", "))
    if(!map[[mtype]] %in% mtypes){
      driver.muts.all$genuine.driver[i] <- F
      next
    }
    # For amplifications and deletions, if these are in the list of validated types then we are done
    if(mtype %in% c("CN amplification", "CN deletion")){
      driver.muts.all$genuine.driver[i] <- T
      next
    }
    
    roles <- unlist(strsplit(driver.muts.all$Role.in.Cancer[i], ", "))
    # For tumour suppressors, check for deleterious mutations using VEP
    if("TSG" %in% roles & grepl("HIGH|MODERATE", driver.muts.all$VEP.impact[i])){
      driver.muts.all$genuine.driver[i] <- T
      next
    }
    
    # For oncogenes, look for the specific mutation in COSMIC
    ccount <- driver.muts.all$cosmic.count[i]
    if("oncogene" %in% roles & !is.na(ccount) & ccount >= 3){
      driver.muts.all$genuine.driver[i] <- T
      next
    }
    
    # For fusion proteins, make sure the translocation partner is validated
    if("fusion" %in% roles & mtype == "Rearrangement" & driver.muts.all$Transloc.Gene[i] %in% unlist(strsplit(driver.muts.all$Validated.Transloc.Partner[i], ", "))){
      driver.muts.all$genuine.driver[i] <- T
      next
    }
    
    # If all above fail, the gene is not called as a genuine driver
    driver.muts.all$genuine.driver[i] <- F
  }
  
  return(driver.muts.all)
}