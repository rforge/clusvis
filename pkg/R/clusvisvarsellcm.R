logvarsellcm <- function(output){
  logtik <- matrix(log(output@param@pi), output@data@n, output@model@g, byrow=TRUE)
  if (output@data@withContinuous){
    for (j in 1:output@data@dataContinuous@d){
      who <- which(output@data@dataContinuous@notNA[,j]==1)
      for (k in 1:output@model@g){
        logtik[who,k] <- logtik[who,k] + dnorm(output@data@dataContinuous@data[who,j], output@param@paramContinuous@mu[j,k], output@param@paramContinuous@sd[j,k], log=TRUE)
      }
    }
  }
  if (output@data@withCategorical){
    for (j in 1:output@data@dataCategorical@d){
      for (k in 1:output@model@g){
        for (h in 1:length(output@data@dataCategorical@modalitynames[[j]])){
          who <- which(output@data@dataCategorical@data=="h")
          logtik[who,k] <- logtik[who,k] + log(output@param@paramCategorical@alpha[[j]][k,h])

        }
      }
    }
  }
  logtik <- sweep(logtik, 1, apply(logtik, 1, max), "-")
  logtik <- sweep(logtik, 1, log(rowSums(exp(logtik))), "-")
  logtik
}

dlogsinglevarsellcm <- function(x, j, output){
  dlog <- rep(0, output@model@g)
  if (j<=output@data@dataContinuous@d){
    dlog <- dnorm(x, output@param@paramContinuous@mu[j,], output@param@paramContinuous@sd[j,], log=TRUE)
  }else{
    j <- j - output@data@dataContinuous@d
    dlog <- log(output@param@paramCategorical@alpha[[j]][,x])
  }
  dlog
}

rsinglevarsellcm <- function(j, output, z){
  x <- 0
  if (j<=output@data@dataContinuous@d){
    x <- rnorm(1, output@param@paramContinuous@mu[j,z], output@param@paramContinuous@sd[j,z])
  }else{
    j <- j - output@data@dataContinuous@d
    x <- sample(1:length(output@param@paramCategorical@alpha[[j]][z,]), 1, prob=output@param@paramCategorical@alpha[[j]][z,])
  }
  x
}

rvarsellcm <- function(output){
  z <- sample(1:output@model@g, 1, prob=output@param@pi)
  x <- sapply(1:output@data@d,
              function(j, output, z) rsinglevarsellcm(j, output, z),
              output=output,
              z=z)
  x
}

rlogtikvarsellcm <- function(output){
  x=rvarsellcm(output)
  dlog <- rowSums(sapply(1:output@data@d,
                         function(j, output, x) dlogsinglevarsellcm(x[j],j, output),
                         output=output,
                         x=x)) + log(output@param@pi)
  dlog <- dlog - max(dlog)
  dlog <- dlog - log(sum(exp(dlog)))
  return(dlog)
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
clusvisVarselclm <- function(resvarsellcm, sample.size=2000, maxit=10**3, nbrandomInit=4*resvarsellcm@model@g, nbcpu=3){
  logtik.estim <- t(replicate(sample.size, rlogtikvarsellcm(resvarsellcm)))
  out <- clusvis(logtik.estim, prop=resvarsellcm@param@pi, logtik.obs=logvarsellcm(resvarsellcm), maxit, nbrandomInit, nbcpu)

  return(out)
}

