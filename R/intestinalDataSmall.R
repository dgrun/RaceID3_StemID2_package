#' Single-cell transcriptome data of intestinal epithelial cells
#'
#' This dataset is a smaller subset of the original dataset, which contains gene expression values, i. e. transcript counts, of 278 intestinal epithelial cells.
#' The dataset is included for quick testing and examples. Only cells with >10,000 transcripts per cell and only genes with >20 transcript counts in >10 cells were retained.
#'
#' @format A sparse matrix (using the \pkg{Matrix}) with cells as columns and genes as rows. Entries are raw transcript counts.
#'
#' @return None
#' @references Grün et al. (2016) Cell Stem Cell 19(2): 266-77 <DOI:10.1016/j.stem.2016.05.010>
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/27345837}{PubMed})
#'
"intestinalDataSmall"
