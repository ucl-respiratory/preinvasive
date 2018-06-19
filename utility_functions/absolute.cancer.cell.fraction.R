# Function to calculate cancer cell fraction
# Adapted from Mcgranahan et al, Science
absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number)
{
  f.function <- function (c,purity,local.copy.number)
  {
    
    return(min(c((purity*c) / (2*(1-purity) + purity*local.copy.number),1)))
    
  }
  x              <- dbinom(n.alt,depth, prob=sapply(seq(0.01,1,length.out=100),f.function,purity,local.copy.number))
  if(min(x)==0)
  {
    x[length(x)] <- 1
  }
  
  names(x)       <- seq(0.01,1,length.out=100)
  sub.cint <- function(x, prob = 0.95,n.alt,depth) {
    xnorm   <- x/sum(x)
    xsort   <- sort(xnorm, decreasing = TRUE)
    xcumLik <- cumsum(xsort)
    n = sum(xcumLik < prob) + 1
    LikThresh <- xsort[n]
    cint  <- x[xnorm >= LikThresh]
    all   <- as.numeric(names(x))
    cellu <- as.numeric(names(cint))
    l.t   <- cellu[1]
    r.t   <- cellu[length(cellu)]
    m     <- cellu[which.max(cint)]
    
    prob.subclonal <- sum(xnorm[1:90])# 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative=‘less’)$p.val
    prob.clonal    <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative=‘greater’)$p.val
    
    data.frame(
      ccf.lower.ci = l.t, ccf.est = m, ccf.upper.ci = r.t,
      prob.subclonal=prob.subclonal,prob.clonal=prob.clonal, 
      # Define clonality strictly as a confidence interval overlapping 1
      is.clonal=(l.t <= 1 & r.t >= 1),
      # Define clonality differently as prob.clonal > prob.subclonal (only for comparison purposes)
      is.clonal.byprob=(prob.clonal > prob.subclonal)
    )
  }
  
  
  return(sub.cint(x,n.alt=n.alt,depth=depth))
  
  
}
