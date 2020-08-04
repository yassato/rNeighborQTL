#' Calculating neighbor QTL effects
#'
#' A function to calculate neighbor QTL effects between two individuals, with given deviation coefficients and conditional genotype probabilities.
#' @param i ID of a target individual.
#' @param j ID of an interacting neighbor.
#' @param a2 A numeric scalar indicating additive deviation.
#' @param d2 A numeric scalar indicating dominance deviation.
#' @param AA An individual x marker matrix of conditional probabilities for AA genotype.
#' @param AB An individual x marker matrix of conditional probabilities for AB genotype. Input NA if heterozygotes are absent.
#' @param BB An individual x marker matrix of conditional probabilities for BB genotype. Input NA for backcross lines.
#' @param d2sq0 An option to make AB/AB interaction effects zero.
#' @return A numeric vector containing each marker effect for individual i.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
neiprob = function(i, j, a2, d2, AA, AB, BB, d2sq0=FALSE) {
  if(is.na(AB)[1]==TRUE) {
    AAi <- (a2^2)*AA[i,]*AA[j,]-(a2^2)*AA[i,]*BB[j,]
    BBi <- (a2^2)*BB[i,]*BB[j,]-(a2^2)*BB[i,]*AA[j,]
    return(AAi+BBi)
  } else if(is.na(BB)[1]==TRUE) {
    AAi <- (a2^2)*AA[i,]*AA[j,]+a2*d2*AA[i,]*AB[j,]
    ABi <- a2*d2*AB[i,]*AA[j,]+(d2^2)*AB[i,]*AB[j,]
    return(AAi+ABi)
  } else {
    AAi <- (a2^2)*AA[i,]*AA[j,]+a2*d2*AA[i,]*AB[j,]-(a2^2)*AA[i,]*BB[j,]
    if(d2sq0==TRUE) {
      ABi <- a2*d2*AB[i,]*AA[j,]+0*AB[i,]*AB[j,]-a2*d2*AB[i,]*BB[j,]
    } else {
      ABi <- a2*d2*AB[i,]*AA[j,]+(d2^2)*AB[i,]*AB[j,]-a2*d2*AB[i,]*BB[j,]
    }
    BBi <- (a2^2)*BB[i,]*BB[j,]-a2*d2*BB[i,]*AB[j,]-(a2^2)*BB[i,]*AA[j,]
    return(AAi+ABi+BBi)
  }
}
