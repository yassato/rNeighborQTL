#' Calculating self QTL effects
#'
#' A function to calculate self QTL effects for an individual, with given deviation coefficients and conditional genotype probabilities.
#' @param i ID of a target individual.
#' @param a1 A numeric scalar indicating additive deviation.
#' @param d1 A numeric scalar indicating dominance deviation.
#' @param AA An individual x marker matrix of conditional probabilities for AA genotype.
#' @param AB An individual x marker matrix of conditional probabilities for AB genotype. Input NA if heterozygotes are absent.
#' @param BB An individual x marker matrix of conditional probabilities for BB genotype. Input NA for backcross lines.
#' @return A numeric vector containing each marker effect for individual i.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
selfprob = function(i, a1, d1, AA, AB, BB) {
  if(is.na(AB)[1]==TRUE) {
    return(t(a1*t(BB[i,,drop=F])-a1*t(AA[i,,drop=F])))
  } else if(is.na(BB)[1]==TRUE) {
    return(t(d1*t(AB[i,,drop=F])-a1*t(AA[i,,drop=F])))
  } else {
    return(t(a1*t(BB[i,,drop=F])+d1*t(AB[i,,drop=F])-a1*t(AA[i,,drop=F])))
  }
}
