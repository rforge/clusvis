
##' Model-based visualization of model-based clustering.
##'
##'  
##'
##' \tabular{ll}{
##'   Package: \tab ClusVis\cr
##'   Type: \tab Package\cr
##'   Version: \tab 1.0.0\cr
##'   Date: \tab 2017-07-10\cr
##'   License: \tab GPL-2\cr
##'   LazyLoad: \tab yes\cr
##' }
##'
##' @description 
##' The main function for parameter inference is \link{clusvis}.
##' However, specific functions for parameter inference  \link{clusvisMixmod} are implemented to deal with model-based clustering done with R packages Rmixmod and Rmixcomp respectively.
##' After parameter inference, visualization is done with function \link{plotDensityClusVisu}.
##'
##'
##' @name ClusVis-package
##' @aliases ClusVis
##' @rdname ClusVis-package
##' @docType package
##' @keywords package
##' @useDynLib ClusVis
##' @author Biernacki, C. and Marbac, M. and Vandewalle, V.
##' @import Rcpp
##' @import MASS
##' @import parallel
##' @import mgcv
##' @import mvtnorm
##' @import Rmixmod
##' @importFrom graphics contour image legend points text 
##' @importFrom stats dist dnorm optim optimize rnorm runif
##' @examples
##' \dontrun{
##' ### Categorical data clustering
##' # Package loading
##' require(Rmixmod)
##'  
##' # Data loading (categorical data)
##' data(birds)
##'
##' # Model-based clustering with 3 components
##' resmixmod <- mixmodCluster(birds, 3)
##'
##' # Inference of the parameters used for results visualization (general approach)
##' # Probabilities of classification are not sampled from the model parameter,
##' # but observed probabilities of classification are used for parameter estimation
##' resvisu <- clusvis(log(resmixmod@bestResult@proba),
##'                    resmixmod@bestResult@parameters@proportions)
##'
##' # Inference of the parameters used for results visualization
##' # (specific for Rmixmod results)
##' # It is better because probabilities of classification are generated
##' # by using the model parameters
##' resvisu <- clusvisMixmod(resmixmod)
##'
##' # Component interpretation graph
##' plotDensityClusVisu(resvisu)
##'
##' # Scatter-plot of the observation memberships
##' plotDensityClusVisu(resvisu,  add.obs = TRUE)
##' }
NULL


###################################################################################
##' This function estimates the parameters used for visualization
##'
##'
##' @param logtik.estim matrix. It contains the probabilities of classification used for parameter inference (should be sampled from the model parameter or computed from the observations).
##' @param prop vector. It contains the class proportions (by default, classes have same proportion).
##' @param logtik.obs   matrix. It contains the probabilities of classification of the clustered sample. If missing, logtik.estim is used.
##' @param maxit numeric. It limits the number of iterations for the Quasi-Newton algorithm (default 1000).
##' @param nbrandomInit numeric. It defines the number of random initialization of the Quasi-Newton algorithm.
##' @param nbcpu numeric. It specifies the number of CPU (only for linux)
##'
##' @return Returns a list
##' @export
##'  
##' @examples
##' \dontrun{
##' ### Package loading
##' require(Rmixmod)
##'  
##' # Data loading (categorical data)
##' data(birds)
##'
##' # Model-based clustering with 3 components
##' resmixmod <- mixmodCluster(birds, 3)
##'
##' # Inference of the parameters used for results visualization (general approach)
##' # Probabilities of classification are not sampled from the model parameter,
##' # but observed probabilities of classification are used for parameter estimation
##' resvisu <- clusvis(log(resmixmod@bestResult@proba),
##'                    resmixmod@bestResult@parameters@proportions)
##  
##'
##'
##'}
clusvis <- function(logtik.estim,
                    prop=rep(1/ncol(logtik.estim), ncol(logtik.estim)),
                    logtik.obs=NULL,
                    maxit=10**3,
                    nbrandomInit=12,
                    nbcpu=1){
  if (Sys.info()["sysname"] == "Linux"){
    nbcpu <- min(nbcpu, detectCores())
  }else{
    nbcpu <- 1
  }
  out <- list()
  if (length(prop)>=3){
    out$error <- FALSE
    if (any(logtik.estim == -Inf)) logtik.estim <- logtik.estim[which(rowSums(logtik.estim == -Inf) ==0),]
    if (is.null(logtik.obs))   logtik.obs <- logtik.estim
    
    out$prop <- prop
    out$startvec <- smartInit(logtik.estim, out$prop)
    out$startmat <- convertmuVecToMat(out$startvec,  ncol(logtik.estim) - 1 )
    loguref <- sweep(logtik.estim, 1, logtik.estim[,ncol(logtik.estim)], "-")[, -ncol(logtik.estim), drop=FALSE]
    logurefweighted <- sweep(loguref, 2, log(out$prop[length(out$prop)]) - log(out$prop[-length(out$prop)]), "+")
    if (length(prop)==2){
      out$centers <- out$startvec
    }else{
      allinit <- list(out$startvec)
      if (nbrandomInit>1)
        allinit <- c(allinit, lapply(1:(nbrandomInit-1), function(i) runif(length(out$startvec), max = 25)) )
      allresults <- mclapply(allinit,
                             function(start)
                               optim(par = out$startvec,
                                     fn = computeLikelihoodCPP,
                                     gr = computeGradientCPP,
                                     control = list(maxit=maxit, fnscale=-1),
                                     Rprop=prop,
                                     Rlogu=logurefweighted,
                                     Rtik=exp(logtik.estim),
                                     method="BFGS"), mc.cores = nbcpu)
      values <- sapply(allresults, function(i) i$value)
      out$optim <- allresults[[which.max(values)]]
      values <- round(values, 2)
      out$nbbest <- sum(values == max(values))
      out$centers <- convertmuVecToMat(out$optim$par, ncol(logtik.estim) - 1 )
      out$allresults <- allresults
    }
    out$logtik.obs <- logtik.obs
    out$dist <- as.matrix(dist(out$centers))
    out$partition <- list(hard= apply(logtik.obs, 1, which.max), fuzzy= exp(logtik.obs))
    loguobs <- sweep(logtik.obs, 1, logtik.obs[,ncol(logtik.obs)], "-")[, -ncol(logtik.obs), drop=FALSE]
    loguobsweighted <- sweep(loguobs, 2, log(out$prop[length(out$prop)]) - log(out$prop[-length(out$prop)]), "+")
    out$y <- phiinvall(loguobsweighted, out$centers)
    g <- as.numeric(t(out$centers) %*% out$prop)
    out$centers <- sweep(out$centers, 2, g, "-")
    tmp <- eigen(t(out$centers) %*% diag(out$prop) %*% out$centers)
    out$centers <- out$centers %*% tmp$vectors
    out$y <- sweep(out$y, 2, g, "-") %*% tmp$vectors
    out$inertia <- tmp$values
    out$logtik.obs <- logtik.obs
    # out$modes <- modesSearch(out$prop, out$centers)
    tmp <- rlogtikvisu(out, nrow(logtik.estim))
    out$EM <- -sum(exp(logtik.estim) * logtik.estim)/(log(ncol(tmp)) * nrow(tmp))
    out$EV <- -sum(exp(tmp) * tmp) / (log(ncol(tmp)) * nrow(logtik.estim))
    tmp2 <- rlogtikvisu(out, nrow(logtik.estim), 1:(length(out$prop)-1))
    out$EVtot <- -sum(exp(tmp2) * tmp2) / (log(ncol(tmp2)) * nrow(logtik.estim))
  }else{
    out$error <- TRUE
  }
    return(out)
}



rlogtikvisu <- function(out, sample.size, dim=c(1,2)){
  z <- sample(1:length(out$prop), sample.size, replace=TRUE, prob = out$prop)
  x <- sapply(dim, function(d) rnorm(sample.size, out$centers[z,d], 1))
  logprob <- sapply(1:length(out$prop),
                    function(z) rowSums(sapply(dim, function(d) dnorm(x[,d], out$centers[z,d], 1, log=TRUE))) + log(out$prop[z]))
  logprob <- sweep(logprob, 1, apply(logprob, 1, max), "-")
  sweep(logprob, 1, log(rowSums(exp(logprob))), "-")
}
