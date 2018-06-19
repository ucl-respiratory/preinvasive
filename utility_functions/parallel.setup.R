# This analysis is computationally intensive so do in parallel if possible
library(foreach)
library(doParallel)
cores=detectCores()
if(cores > 1){
  cl <- makeCluster(min(cores[1]-1, 20))
  registerDoParallel(cl)
}