

logmixcomp <- function(output){
  param <- output$variable$param
  x <- output$variable$data
  logtik <- matrix(0, length(x$z_class$completed), output$mixture$nbCluster)
  for (k in 1:output$mixture$nbCluster)  logtik[,k] <- log(param$z_class$pi$stat[k,1])

  for (j in 2:length(x)){
    if (output$variable$type[[j]] == "Gaussian_sjk"){
      for (k in 1:ncol(logtik)){
        tmp <- matrix(param[[j]]$NumericalParam$stat[,1], 2, output$mixture$nbCluster)
        logtik[,k] <- logtik[,k] + dnorm(x[[j]]$completed, tmp[1,k], tmp[2,k], log=TRUE)
      }
    }
    if (output$variable$type[[j]] == "Categorical_pjk"){
      for (k in 1:ncol(logtik)){
        tmp <- log(matrix(param[[j]]$NumericalParam$stat[,1], ncol=output$mixture$nbCluster))
        logtik[,k] <- logtik[,k] + tmp[x[[j]]$completed, k]
      }
    }

  }
  logtik <- sweep(logtik, 1, apply(logtik, 1, max), "-")
  logtik <- sweep(logtik, 1, log(rowSums(exp(logtik))), "-")
  logtik
}

dlogsinglemixcomp <- function(x, p, type, g){
  dlog <- rep(0, g)
  if (type=="Gaussian_sjk"){
    p <- matrix(p, g, 2)
    dlog <- dnorm(x, p[,1], p[,2], log = TRUE)
  }else if (type == "Categorical_pjk"){
    m <- length(p) / g
    p <- matrix(p, m, g)
    dlog <- log(p[x,])
  }else{
    stop("type de sortie non programmée")
  }
  dlog
}

rsinglemixcomp <- function(p, type, z, g, Tseq=NULL){
  x <- 0
  if (type=="Gaussian_sjk"){
    x <- rnorm(1, p[1+(z-1)*2], p[2+(z-1)*2])
  }else if (type == "Categorical_pjk"){
    m <- length(p) / g
    x <- sample(1:m, 1, prob=p[m*(z-1) + (1:m)])
  }else {
    stop("type de sortie non programmée")


  }
  x
}

rmixcomp <- function(output){
  z <- sample(1:output$mixture$nbCluster, 1, prob=output$variable$param$z_class$pi$stat[,1])
  x <- sapply(1:(length(output$variable$param)-1),
              function(j, output, z) rsinglemixcomp(output$variable$param[[j+1]]$NumericalParam$stat[,1], output$variable$type[[j+1]], z, output$mixture$nbCluster),
              output=output,
              z=z)
  names(x) <- names(output$variable$type)[-1]
  x
}

rlogtikmixcomp <- function(output){
  x=rmixcomp(output)
  dlog <- rowSums(sapply(1:(length(output$variable$param)-1),
                         function(j, output, x) dlogsinglemixcomp(x[j], output$variable$param[[j+1]]$NumericalParam$stat[,1], output$variable$type[[j+1]], output$mixture$nbCluster),
                         output=output,
                         x=x)) + log(output$variable$param$z_class$pi$stat[,1])
  dlog <- dlog - max(dlog)
  dlog <- dlog - log(sum(exp(dlog)))
  return(dlog)
}



###################################################################################
##' This function estimates the parameters used for visualization of model-based clustering performs with R package Rmixmod. To achieve the parameter infernece, it automatically samples probabilities of classification from the model parameters
##'
##'
##' @param resmixcomp MixmodCluster. It is an instance of class MixmodCluster returned by function mixmodCluster of R package Rmixmod.
##' @param sample.size numeric. Number of probabilities of classification sampled for parameter inference.
##' @param maxit numeric. It limits the number of iterations for the Quasi-Newton algorithm (default 1000)!!!!!!
##' @param nbrandomInit numeric. It defines the number of random initialization of the Quasi-Newton algorithm.
##' @param nbcpu numeric. It specifies the number of CPU (only for linux)
##'
##' @return Returns a list
##' @examples
##' # Data loading (categorical data)
##' data(birds)
##'
##' # Model-based clustering with 3 components
##' resmixmod <- mixmodCluster(birds, 3)
##'
##' # Inference of the parameters used for results visualization (specific for Rmixmod results)
##' resvisu <- clusvisMixmod(resmixmod)
##'
##' # Component interpretation graph
##' plotDensityClusVisu(resvisu)
##'
##' # Scatter-plot of the observation memberships
##' plotDensityClusVisu(resvisu,  add.obs = TRUE)
##' @export
##'
##'
clusvisMixcomp <- function(resmixcomp, sample.size=5000, maxit=10**3, nbrandomInit=12, nbcpu=3){
  logtik.estim <- t(replicate(sample.size, rlogtikmixcomp(resmixcomp)))
  out <- clusvis(logtik.estim, prop=resmixcomp$variable$param$z_class$pi$stat[,1], logtik.obs=logmixcomp(resmixcomp), maxit, nbrandomInit, nbcpu)
  return(out)
}

