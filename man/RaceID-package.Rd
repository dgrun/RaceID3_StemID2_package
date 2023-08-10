\name{RaceID-package}
\alias{RaceID-package}
\alias{RaceID}
\docType{package}
\title{
Identification of Cell Types, Inference of Lineage Trees, and Prediction of Noise Dynamics from Single-Cell RNA-Seq Data
}
\description{
RaceID is a clustering algorithm for the identification of cell types from single-cell RNA-sequencing data. It was specifically designed for the detection of rare cells which correspond to outliers in conventional clustering methods. The package contains RaceID3, the most recently published version of this algorithm, and StemID2, an algorithm for the identification of lineage trees based on RaceID3 analysis. RaceID3 utilizes single cell expression data, and was designed to work well with quantitative single-cell RNA-seq data incorporating unique molecular identifiers. It requires a gene-by-cell expression matrix as input and produces a clustering partition representing cell types. StemID2 assembles these cell types into a lineage tree.
The RaceID package (>= v0.1.4) also contains functions for a VarID analysis. VarID comprises a sensitive clustering method utilizing pruned k-nearest neighbor networks, connecting only cells with links supported by a background model of gene expression. These pruned k-nearest neighbor networks further enable the definition of homogenous neighborhoods for the quantification of local gene expression
variability in cell state space.
}
\details{
For details please see vignette.}
\author{
Dominic Grun, dominic.gruen@gmail.com.

Maintainer: Dominic Grun <dominic.gruen@gmail.com>
}
\references{
Herman, J.S., Sagar, Grun D. (2018) <DOI:10.1038/nmeth.4662>
Rosales-Alvarez, R.E., Rettkowski, J., Herman, J.S., Dumbovic, G., Cabezas-Wallscheid, N., Grun, D. (2023) <DOI:10.1186/s13059-023-02974-1>}
\keyword{ package }
