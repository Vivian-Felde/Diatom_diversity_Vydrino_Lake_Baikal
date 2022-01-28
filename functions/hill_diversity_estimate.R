#########################################################################
## Estimation of interpolation of individual-based Hill number
## Note: The following R function is code modified from Chao et al. 2014. Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. Ecological Monographs 84:45-67. See also R-package iNEXT (https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html)
## Relevant modfications are made to extract estimates at the minimum sample size 
## x is a matrix of samples and taxa
## Samplesize is the minimum sample size
## The output is a matrix with Hills number N0, N1, and N2, related evenness ratios, and error estimates


hill_diversity_estimate <- function(x, samplesize, MARGIN = 1){
  
  x <- as.matrix(x)
  
  if(missing(samplesize)){
    samplesize = min(apply(x, MARGIN, sum))
  }
  
  Hill0 <- function(x, samplesize){
    x <- x[x > 0]
    n <- sum(x)
    if(samplesize <= n){
      Sub <- sum(1-exp(lchoose(n - x, samplesize) - lchoose(n, samplesize)))
    }
    
    else {0}
  }  
  
  fk.hat <- function(x, samplesize){
    x <- x[x > 0]
    n <- sum(x)
    if(samplesize <= n){
      Sub <- function(k)	sum(exp(lchoose(x, k) + lchoose(n - x, samplesize - k) - lchoose(n, samplesize)))
      sapply(1:samplesize, Sub)
    }
    else {0}
  } 
  
  Hill1 <- function(x, samplesize){
    x <- x[x > 0]
    n <- sum(x)
    if(samplesize <= n){
      k <- 1:samplesize
      Sub <- exp(-sum(k/samplesize * log(k/ samplesize) * fk.hat(x, samplesize)))
    }      
    else{0}
  }
  
  Hill2 <- function(x, samplesize){
    x <- x[x > 0]
    n <- sum(x)
    Sub <- 1 / (1 / samplesize + (1 - 1 / samplesize) * sum(x * (x - 1) / n / (n - 1)))
  }
  
  # Estimates of Hill numbers
  estN0 <- sapply(samplesize, function(n) apply(x, MARGIN, Hill0, samplesize = n))
  estN1 <- sapply(samplesize, function(n) apply(x, MARGIN, Hill1, samplesize = n))
  estN2 <- sapply(samplesize, function(n) apply(x, MARGIN, Hill2, samplesize = n))
  hilldiv <- cbind("n0" = estN0, "n1" = estN1, "n2" = estN2)
  colnames(hilldiv) <- c("N0", "N1", "N2")
  rownames(hilldiv) <- rownames(x)    
  
  # Evenness ratios
  hillsratios <- cbind(estN1-estN2, estN2/estN1, estN1/estN0, (estN2-1)/(estN1-1), (estN1-1)/(estN0-1))  
  colnames(hillsratios) <- c("N1-N2", "N2/N1", "N1/N0", "modN2/N1", "modN1/N0")
  hilldiv <- cbind(hilldiv, hillsratios)
  
  # Confidence intervals
  conf.est <- function(x, samplesize){
    x <- x[x>0]
    est.prob <- EstiBootComm.Ind(x) #function from Chao et al 2014 below
    Abun.Mat <- rmultinom(200, samplesize, est.prob) 
    error.0 <- qnorm(0.975) * sd(apply(Abun.Mat, 2, function(x) Hill0(x, samplesize)))
    error.1 <- qnorm(0.975) * sd(apply(Abun.Mat, 2, function(x) exp(-sum(x[x>0]/sum(x[x>0]) * log(x[x>0] / sum(x[x>0]))))))
    error.2 <- qnorm(0.975) * sd(apply(Abun.Mat, 2, function(x) Hill2(x, samplesize)))
    err.est <- c("error.0" = error.0, "error.1" = error.1, "error.2" = error.2)
    return(err.est)
    
  }
  
  err.est <- lapply(nrow(x), function(z) apply(x, MARGIN, conf.est, samplesize = z))
  err.est <- t(err.est[[1]])
  hilldiv <- cbind(hilldiv, err.est)
  return(hilldiv)
}

#########################################################################
# Function to estimate confidence intervals from Chao et al. 2014
#########################################################################
EstiBootComm.Ind <- function(Spec)
{
  Sobs <- sum(Spec > 0) 	#observed species
  n <- sum(Spec)		  	#sample size
  f1 <- sum(Spec == 1) 	#singleton 
  f2 <- sum(Spec == 2) 	#doubleton
  a <- ifelse(f1 == 0, 0, (n - 1) * f1 / ((n - 1) * f1 + 2 * f2) * f1 / n)
  b <- sum(Spec / n * (1 - Spec / n) ^ n)
  w <- a / b  			#adjusted factor for rare species in the sample
  f0.hat <- ceiling(ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2))	#estimation of unseen species via Chao1
  Prob.hat <- Spec / n * (1 - w * (1 - Spec / n) ^ n)					#estimation of relative abundance of observed species in the sample
  Prob.hat.Unse <- rep(2 * f2/((n - 1) * f1 + 2 * f2), f0.hat)		#estimation of relative abundance of unseen species in the sample
  return(c(Prob.hat, Prob.hat.Unse))									#Output: a vector of estimated relative abundance
}



