% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VarID_functions.R
\name{createKnnMatrix}
\alias{createKnnMatrix}
\title{Function to create a knn matrix}
\usage{
createKnnMatrix(res, pvalue = 0.01)
}
\arguments{
\item{res}{List object with k nearest neighbour information returned by \code{pruneKnn} function.}

\item{pvalue}{Positive real number between 0 and 1. All nearest neighbours with link probability \code{< pvalue} are discarded. Default is 0.01.}
}
\value{
Adjacency matrix in sparse matrix format (see package \pkg{Matrix}) with positive non-zero entries only for k nearest neighours with link probability \code{>= pvalue}. The value of these entries equals the link probability.
}
\description{
This creates an adjacency matrix, keeping only nearest neighbour with a link probability above a minimum probability
}
\examples{
res <- pruneKnn(intestinalDataSmall,metric="pearson",knn=10,alpha=1,no_cores=1,FSelect=FALSE)
y <- createKnnMatrix(res,pvalue=0.01)
}
