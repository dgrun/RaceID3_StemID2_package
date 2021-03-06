% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RaceID.R
\name{comptsne}
\alias{comptsne}
\title{Computation of a two dimensional t-SNE representation}
\usage{
comptsne(
  object,
  dimRed = FALSE,
  initial_cmd = TRUE,
  perplexity = 30,
  rseed = 15555
)
}
\arguments{
\item{object}{\code{SCseq} class object.}

\item{dimRed}{logical. If \code{TRUE} then the t-SNE is computed from the feature matrix in slot \code{dimRed$x} (if not equal to \code{NULL}).
Default is \code{FALSE} and the t-SNE is computed from the distance matrix stored in slot \code{distances}. If slot \code{distances} equals \code{NULL}
\code{dimRed} is automatially set to \code{TRUE}.}

\item{initial_cmd}{logical. If \code{TRUE}, then the t-SNE map computation is initialized with a configuration obtained by classical
multidimensional scaling. Default is \code{TRUE}.}

\item{perplexity}{Positive number. Perplexity of the t-SNE map. Default is \code{30}.}

\item{rseed}{Integer number. Random seed to enforce reproducible t-SNE map.}
}
\value{
\code{SCseq} object with t-SNE coordinates stored in slot \code{tsne}.
}
\description{
This functions performs the computation of a t-SNE map from the distance object in slot \code{distances} using the \pkg{Rtsne} package.
}
\examples{
sc <- SCseq(intestinalDataSmall)
sc <- filterdata(sc)
sc <- compdist(sc)
sc <- clustexp(sc)
sc <- findoutliers(sc)
sc <- comptsne(sc)
}
