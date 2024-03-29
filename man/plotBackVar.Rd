% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VarID_functions.R
\name{plotBackVar}
\alias{plotBackVar}
\title{Function for plottinhg the background model of gene expression variability}
\usage{
plotBackVar(x)
}
\arguments{
\item{x}{List object returned by function \code{fitBackVar} or list object returned by function \code{pruneKnn} (if it was run with \code{FSelect=TRUE}).}
}
\value{
None
}
\description{
This function plots the variance against mean expression across all genes and a second order polynomial to the variance-mean dependence in log space. It also plots a local regression.
}
\examples{
bg <- fitBackVar(intestinalDataSmall)
plotBackVar(bg)
}
