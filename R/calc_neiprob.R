#' Calculating a set of neighbor QTL effects from conditional genotype probabilities
#'
#' A function to calculate self QTL effects for all individuals, with given deviation coefficients and conditional genotype probabilities.
#' @param genoprobs Conditional genotype probabilities as taken from \code{qtl::calc.genoprob()}.
#' @param contrasts A vector composed of three TRUE/FALSE values, which represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial positions along an x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param a2 A numeric scalar indicating additive deviation.
#' @param d2 A numeric scalar indicating dominance deviation.
#' @param grouping An integer vector assigning each individual to a group. This argument can be used when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param d2sq0 An option to make AB/AB interaction effects zero.
#' @return A numeric matrix containing individuals x marker elements for neighbor QTL effects.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
calc_neiprob = function(genoprobs, a2, d2, contrasts=c(TRUE,TRUE,TRUE), smap, scale, grouping=rep(1,nrow(smap)), d2sq0=FALSE) {
  p <- dim(genoprobs$geno[[1]]$prob)[1]
  geno <- decompose_genoprobs(genoprobs=genoprobs, contrasts=contrasts)

  neiprob_i <- function(i) {
    id <- c(1:p)[grouping == grouping[i]]

    d_i <- mapply(function(x) { return(sqrt((smap[x,1]-smap[i,1])^2 + (smap[x,2]-smap[i,2])^2)) },id)
    prob_i <- 0
    j_id <- id[(d_i>0)&(d_i<=scale)]
    if(length(j_id)==0) {
      return(rep(0,ncol(geno$AA)))
    } else {
      for(j in j_id){
        prob_ij <- neiprob(i=i, j=j, a2=a2, d2=d2, AA=geno$AA, AB=geno$AB, BB=geno$BB, d2sq0=d2sq0)
        prob_i <- prob_i + prob_ij
      }
      prob_i <- prob_i/length(j_id)
      return(prob_i)
    }
  }
  neiList <- mapply(neiprob_i, 1:p)
  neiList <- t(neiList)

  marker_info <- get_markers(genoprobs=genoprobs)
  colnames(neiList) <- rownames(marker_info)
  return(neiList)
}
