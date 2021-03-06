% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plotDensityClusVisu}
\alias{plotDensityClusVisu}
\title{Function for visualizing the clustering results}
\usage{
plotDensityClusVisu(res, dim = c(1, 2), threshold = 0.95, add.obs = FALSE,
  positionlegend = "topright", xlim = NULL, ylim = NULL)
}
\arguments{
\item{res}{object return by function \link{clusvis} or \link{clusvis}}

\item{dim}{numeric. This vector of size two choose the axes to represent.}

\item{threshold}{numeric. It contains the thersholds used for computing the level curves.}

\item{add.obs}{boolean. If TRUE, coordinnates of the observations are plotted.}

\item{positionlegend}{character. It specifies the legend location.}

\item{xlim}{numeric. It specifies the range of x-axis.}

\item{ylim}{numeric. It specifies the range of y-axis.}
}
\description{
Function for visualizing the clustering results
}
\examples{
\dontrun{
# Package loading
require(Rmixmod)

# Data loading (categorical data)
data(birds)

# Model-based clustering with 3 components
resmixmod <- mixmodCluster(birds, 3)

# Inference of the parameters used for results visualization
# (specific for Rmixmod results)
resvisu <- clusvisMixmod(resmixmod)

# Component interpretation graph
plotDensityClusVisu(resvisu)

# Scatter-plot of the observation memberships
plotDensityClusVisu(resvisu,  add.obs = TRUE)
}
}
