# VEP predictions for a given mutation
# A mutation list in the format of muts.all
# Script from here: http://grch37.rest.ensembl.org/documentation/info/vep_region_post

# Get driver gene variants in correct format
#vep.muts <- driver.muts.all[which(!is.na(driver.muts.all$Reference.Allele))]

vep.annotate <- function(vep.muts){
  
  vep.input <- unlist(lapply(1:dim(vep.muts)[1], function(i){
    paste(vep.muts$Chromosome[i], vep.muts$Start.Position[i], ".", vep.muts$Reference.Allele[i], vep.muts$Variant.Allele[i], ".", ".", ".", sep=" ")
  }))
  
  library(httr)
  library(jsonlite)
  library(xml2)
  
  server <- "http://grch37.rest.ensembl.org"
  ext <- "/vep/homo_sapiens/region"
  r <- POST(
    paste(server, ext, sep = ""), 
    content_type("application/json"), 
    accept("application/json"), 
    body = paste0('{ "variants" : ["',paste(vep.input, collapse='", "'),'"] }')
  )
  
  stop_for_status(r)
  
  # use this if you get a simple nested list back, otherwise inspect its structure
  x <- data.frame(t(sapply(httr::content(r),c)))
  #head(fromJSON(toJSON(content(r))))
  
  
  # Parse this list and add to input data frame. For each we want to know:
  #   Gene name and consequence (to check the same as ours)
  #   Impact of transcripts affecting the named gene (+ SIFT and polyphen scores)
  #vep.muts$vep.gene <- NA
  #vep.muts$vep.consequence <- NA
  vep.muts$vep.impact <- NA
  #vep.muts$sift.impact <- NA
  #vep.muts$polyphen.impact <- NA
  for(i in 1:dim(vep.muts)[1]){
    y <- x[,i][[1]]
    #vep.muts$vep.gene[i] <- paste(unique(unlist(lapply(y$transcript_consequences, function(z){z$gene_symbol}))), collapse=",")
    #vep.muts$vep.consequence[i] <- y$most_severe_consequence
    vep.muts$vep.impact[i] <- paste(unique(unlist(lapply(y$transcript_consequences, function(z){z$impact}))), collapse=",")
    #vep.muts$sift.impact[i] <- paste(unique(unlist(lapply(y$transcript_consequences, function(z){z$sift_prediction}))), collapse=",")
    #vep.muts$polyphen.impact[i] <- paste(unique(unlist(lapply(y$transcript_consequences, function(z){z$polyphen_prediction}))), collapse=",")
  }
  
  return(vep.muts)
}
