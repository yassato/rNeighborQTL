#' @title rNeighbor QTL: quantitative trait loci mapping for neighbor effects
#'
#' @description
#' The "rNeighborQTL" package provides a set of functions to perform a genome scan of neighbor QTL effects.
#' Theoretical background is described in Sato, Takeda & Nagano (2021).
#' See \code{vignette("rNeighborQTL")} for how to use this package.
#' @docType package
#' @name rNeighborQTL-package
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' Neighbor QTL is an extended method of single-marker regression.
#' Taking conditional genotype probabilities from \code{r/qtl} package, the \code{rNeighborQTL} package performs an interval mapping of neighbor genotypic identity.
#' Alternative to QTL mapping, \code{rNeighborGWAS} package is also available for genome-wide association studies (GWAS) via CRAN at \url{https://cran.r-project.org/package=rNeighborGWAS}.
#' @references
#' * Sato Y, Yamamoto E, Shimizu KK, Nagano AJ (2021) Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory. Heredity \url{https://doi.org/10.1038/s41437-020-00401-w}
#' * Sato Y, Takeda K, Nagano AJ (2021) Neighbor QTL: an interval mapping method for quantitative trait loci underlying plant neighborhood effects. G3; Genes|Genomes|Genetics \url{https://doi.org/10.1093/g3journal/jkab017}
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
