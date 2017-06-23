rmixmod <- function(resmixmod){
  x <- NA
  if (resmixmod@dataType == "quantitative"){
    z <- sample(1:resmixmod@bestResult@nbCluster, 1, prob = resmixmod@bestResult@parameters@proportions)
    x <- rmvnorm(1, resmixmod@bestResult@parameters@mean[z,], resmixmod@bestResult@parameters@variance[[z]])
  }else if (resmixmod@dataType == "qualitative"){
    z <- sample(1:resmixmod@bestResult@nbCluster, 1, prob = resmixmod@bestResult@parameters@proportions)
    x <- sapply(1:ncol(resmixmod@bestResult@parameters@center),
                function(j, resmixmod, z){
                  p <- resmixmod@bestResult@parameters@scatter[[z]][j,]
                  p[resmixmod@bestResult@parameters@center[z,j]] <- 1 -  p[resmixmod@bestResult@parameters@center[z,j]]
                  sample(1:length(p), 1, prob = p)
                },
                resmixmod=resmixmod,
                z=z)
  }else if (resmixmod@dataType == "composite"){
    z <- sample(1:resmixmod@bestResult@nbCluster, 1, prob = resmixmod@bestResult@parameters@proportions)
    xcont <- rmvnorm(1, resmixmod@bestResult@parameters@g_parameter@mean[z,], resmixmod@bestResult@parameters@g_parameter@variance[[z]])
    xcat <- sapply(1:ncol(resmixmod@bestResult@parameters@m_parameter@center),
                   function(j, resmixmod, z){
                     p <- resmixmod@bestResult@parameters@m_parameter@scatter[[z]][j,]
                     p[resmixmod@bestResult@parameters@m_parameter@center[z,j]] <- 1 -  p[resmixmod@bestResult@parameters@m_parameter@center[z,j]]
                     sample(1:length(p), 1, prob = p)
                   },
                   resmixmod=resmixmod,
                   z=z)
    x <- list(xcont=xcont, xcat=xcat)
  }else{
    stop("type non programe")
  }
  x
}

dlogsinglequalimixmod <- function(x, j, resmixmod, z){
  p <- resmixmod@bestResult@parameters@scatter[[z]][j,]
  p[resmixmod@bestResult@parameters@center[z,j]] <- 1 -  p[resmixmod@bestResult@parameters@center[z,j]]
  log(p[x])
}

dlogsinglequalimixmodcompo <- function(x, j, resmixmod, z){
  p <- resmixmod@bestResult@parameters@m_parameter@scatter[[z]][j,]
  p[resmixmod@bestResult@parameters@m_parameter@center[z,j]] <- 1 -  p[resmixmod@bestResult@parameters@m_parameter@center[z,j]]
  log(p[x])
}

dlogtik <- function(x, resmixmod){
  dlog <- NA
  if (resmixmod@dataType == "quantitative"){
    dlog <-  sapply(1:resmixmod@bestResult@nbCluster,
                    function(z, x, resmixmod)
                      dmvnorm(x, resmixmod@bestResult@parameters@mean[z,], resmixmod@bestResult@parameters@variance[[z]], log=TRUE) + log(resmixmod@bestResult@parameters@proportions[z]),
                    x=x,
                    resmixmod=resmixmod)
  }else if (resmixmod@dataType == "qualitative"){
    dlog <- rowSums(sapply(1:length(x),
                           function(j, resmixmod, x)
                             sapply(1:resmixmod@bestResult@nbCluster, dlogsinglequalimixmod, x=x[j], j=j, resmixmod=resmixmod),
                           resmixmod=resmixmod,
                           x=x)) + log(resmixmod@bestResult@parameters@proportions)
  }else if (resmixmod@dataType == "composite"){
    dlog <- log(resmixmod@bestResult@parameters@proportions) +
      sapply(1:resmixmod@bestResult@nbCluster,
             function(z, x, resmixmod)
               dmvnorm(x, resmixmod@bestResult@parameters@g_parameter@mean[z,], resmixmod@bestResult@parameters@g_parameter@variance[[z]], log=TRUE),
             x=x$xcont,
             resmixmod=resmixmod) +
      rowSums(sapply(1:length(x$xcat),
                     function(j, resmixmod, x)
                       sapply(1:resmixmod@bestResult@nbCluster, dlogsinglequalimixmodcompo, x=x[j], j=j, resmixmod=resmixmod),
                     resmixmod=resmixmod,
                     x=x$xcat))
  }else{
    stop("type de sortie non programmée")
  }
  dlog <- dlog - max(dlog)
  dlog <- dlog - log(sum(exp(dlog)))
  return(dlog)
}

rlogtikmixmod <- function(resmixmod){
  x <- rmixmod(resmixmod)
  dlogtik(x, resmixmod)
}



###################################################################################
##' This function performs the variable selection and the maximum likelihood estimation of the Latent Class Model
##'
##' @param logtik.estim matrix. It contains the probabilities of classification used for parameter inference.
##' @param prop vector. It contains the class proportions (by default, classes have same proportion).
##' @param logtik.obs   matrix. It contains the probabilities of classification of the clustered sample. If missing, logtik.estim is used.
##' @param maxit numeric. It limits the number of iterations for the Quasi-Newton algorithm (default 1000)
##' @param nbrandomInit numeric. It defines the number of random initialization of the Quasi-Newton algorithm.
##'
##' @return Returns a list
##' @examples
##' set.seed(123)
##' @export
##'
##'
clusvisMixmod <- function(mixmodResult, sample.size=5000, maxit=10**3, nbrandomInit=4*mixmodResult@bestResult@nbCluster, nbcpu=3,loccont=NULL){
  if (mixmodResult@dataType != "composite"){
    tmp <- t(apply(mixmodResult@data, 1, dlogtik, resmixmod=mixmodResult))
  }else{
    loccat <- (1:ncol(mixmodResult@data))[-loccont]
    tmp <- t(sapply(1:nrow(mixmodResult@data), function(i) dlogtik(list(xcont=mixmodResult@data[i,loccont], xcat=mixmodResult@data[i,loccat]), mixmodResult)))
  }
  if (sample.size>0){
    logtik.estim <-  t(replicate(sample.size, rlogtikmixmod(mixmodResult)))
  }else{
    logtik.estim <- tmp
  }


  out <- clusvis(logtik.estim, prop=mixmodResult@bestResult@parameters@proportions, logtik.obs=tmp, maxit, nbrandomInit, nbcpu)
  out$EM <- sum(exp(logtik.estim) * logtik.estim) / (log(length(mixmodResult@bestResult@parameters@proportions)) * sample.size)

  return(out)
}
