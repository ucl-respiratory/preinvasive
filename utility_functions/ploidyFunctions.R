
#######################################################################
fun.ploidy<- function(x, seg){
  # This function estimates the mode of the median ploidy of each chromosome
  # The weighted median ploidy is calculated on a chromosome by chromosome basis. 
  # The weights are equal to the segment length.
  # The variabel seg needs to me a segmentation matrix with columns in the follwing order:
  # Sample Name, Chromosome, Start Region, End Region, Number Probes, Copy Number Value.
  # The matrix for seg must be a character matrix and mustn't be a dataframe.
  # The variable x is should in include the unique sample names in the first column of seg .
  # The function returns the  ploidy for each sample.
  sub = subset(seg, seg$SampleID == x)
  require(limma)
  
  w.meds <- c()
  for (chr in unique(seg$Chr))
  {
    sub.chr <- subset(sub, sub$Chr==chr)
    w.med   <- weighted.median(as.numeric(sub.chr$cn), w = as.numeric(sub.chr$End) - as.numeric(sub.chr$Start)    )
    w.meds  <- c(w.meds,w.med)  
  }
  
  return(fun.mode(w.meds))
  
}

###########################################################################
fun.mode <- function(x, show_all_modes = FALSE)
{
  x_freq <- table(x)
  mode_location <- which.max(x_freq)
  The_mode <- names(x_freq)[mode_location]
  Number_of_modes <- length(mode_location)
  #
  if(show_all_modes) {
    if(Number_of_modes >1) {
      warning(paste("Multiple modes exist - returning all",Number_of_modes,"of them"))}
    return(The_mode)
  } 
  
  else {
    if(Number_of_modes >1) {
      warning(paste("Multiple modes exist - returning only the first one out of", Number_of_modes))}
    return(The_mode[1])
  }
}


