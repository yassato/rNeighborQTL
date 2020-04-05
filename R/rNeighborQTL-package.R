#' @title rNeighbor QTL: quantitative trait loci mapping for neighbor effects
#'
#' @description
#' The "rNeighborQTL" package provides a set of functions to perform a genome scan of neighbor QTL effects.
#' See \code{vignette("rNeighborQTL")} for how to use this package.
#' @docType package
#' @name rNeighborQTL-package
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' Neighbor QTL is an extended method of single-marker regression as proposed by Sato et al. (2019).
#' Taking conditional genotype probabilities from \code{r/qtl2} package, the \code{rNeighborQTL} package performs an interval mapping of neighbor genotypic identity.
#' @references
#' Sato Y, Yamamoto E, Shimizu KK, Nagano AJ (2019) Neighbor GWAS: incorporating neighbor genotypic identity into genome-wide association studies of field herbivory on *Arabidopsis thaliana*. bioRxiv https://doi.org/10.1101/845735
#' @keywords internal
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
