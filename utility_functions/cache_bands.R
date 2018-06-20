# Cache bands - create the bands.dict object in resources/ based on Ensembl data

library(httr)
ensembl_url <- "http://grch37.rest.ensembl.org/info/assembly/homo_sapiens?content-type=application/json&bands=1"
req <- httr::GET(url=ensembl_url)
# JSON data stored in content(req)
# Turn this into a data object with each band represented by chromosome, name, start, end
bands <- c()
chroms <- c()
ids <- c()
starts <- c()
ends <- c()
for(reg in httr::content(req)[["top_level_region"]]){
  for(band in reg$bands){
    chroms <- c(chroms, band$seq_region_name)
    ids <- c(ids, band$id)
    starts <- c(starts, band$start)
    ends <- c(ends, band$end)
    bands <- c(bands, c(band$seq_region_name, band$id, band$start, band$end))
  }
}

bands <- data.frame(chrom=as.character(chroms), id=ids, start=starts, end=ends)
rownames(bands) <- paste(bands$chrom, bands$id, sep="")