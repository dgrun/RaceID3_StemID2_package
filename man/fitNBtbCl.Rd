% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VarID_functions.R
\name{fitNBtbCl}
\alias{fitNBtbCl}
\title{Function for fitting a negative binomial noise model of technical and biological variability}
\usage{
fitNBtbCl(z, mu, rt, gamma = 2, x0 = 0.1, lower = 0, upper = 100)
}
\arguments{
\item{z}{Transcript count matrix with cells as columns and genes as rows.}

\item{mu}{Vector of mean expression values across cells in \code{z}.}

\item{rt}{Vector of dispersion parameters explaining global cell-to-cell variability of transcript counts across cells in \code{z}.}

\item{gamma}{Positive real number. Scale paramter of the cauchy prior. Default is 2.}

\item{x0}{Real number greater or equal to zero. Location parameter of the cauchy prior.}

\item{lower}{Real number greater or equal to zero. Lower bound for the maximum a posterior inference of the biological noise. Default is 0.}

\item{upper}{Real number greater or equal to zero. Upper bound for the maximum a posterior inference of the biological noise. Default is 100.}
}
\value{
Vector of biological noise parameters across  cells in \code{z}.
}
\description{
This function fits a negative binomial model to transcript counts of a group of cells thereby deconvoluting variability into sampling noise, global cell-to-cell variability of transcript counts, and residual variability, which corresponds to biological noise. Local mean and and global cell-to-cell variability of transcript counts are pre-computed arguments.
}
