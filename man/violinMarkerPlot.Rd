% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VarID_functions.R
\name{violinMarkerPlot}
\alias{violinMarkerPlot}
\title{Violin plot of marker gene expression or noise}
\usage{
violinMarkerPlot(g, object, noise = NULL, set = NULL, ti = NULL)
}
\arguments{
\item{g}{Valid gene ID corresponding to a (set of) rownames of \code{object@ndata} or \code{noise}.}

\item{object}{\pkg{RaceID} \code{SCseq} object.}

\item{noise}{List of noise parameters returned by \code{compTBNoise}. If this argument is given, then the distribution of noise (epsilon) is plotted. Default is NULL and normalized gene expression (normalized by median count across all clusters in \code{set}) is plotted.}

\item{set}{Postive integer number or vector of integers corresponding to valid cluster numbers. Violin plots are shown for all clusters in \code{set}. Default is NULL and data are shown for all clusters in \code{object@cpart}.}

\item{ti}{String of characters representing the title of the plot. Default is \code{NULL} and the first element of \code{g} is chosen.}
}
\value{
None
}
\description{
Displaying violin plots of gene expression or gene expression noise (epsilon) across (a set of) clusters
}
