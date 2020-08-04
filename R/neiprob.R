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
    AAi <- (a2^2)*AA[i,]*t(AA[j,,drop=F])-(a2^2)*AA[i,]*t(BB[j,,drop=F])
    ABi <- 0
    BBi <- (a2^2)*BB[i,]*t(BB[j,,drop=F])-(a2^2)*BB[i,]*t(AA[j,,drop=F])
  } else if(is.na(BB)[1]==TRUE) {
    AAi <- (a2^2)*AA[i,]*t(AA[j,,drop=F])+a2*d2*AA[i,]*t(AB[j,,drop=F])
    ABi <- a2*d2*AB[i,]*t(AA[j,,drop=F])+(d2^2)*AB[i,]*t(AB[j,,drop=F])
    BBi <- 0
  } else {
    AAi <- (a2^2)*AA[i,]*t(AA[j,,drop=F])+a2*d2*AA[i,]*t(AB[j,,drop=F])-(a2^2)*AA[i,]*t(BB[j,,drop=F])
    if(d2sq0==TRUE) {
      ABi <- a2*d2*AB[i,]*t(AA[j,,drop=F])+0*AB[i,]*t(AB[j,,drop=F])-a2*d2*AB[i,]*t(BB[j,,drop=F])
    } else {
      ABi <- a2*d2*AB[i,]*t(AA[j,,drop=F])+(d2^2)*AB[i,]*t(AB[j,,drop=F])-a2*d2*AB[i,]*t(BB[j,,drop=F])
    }
    BBi <- (a2^2)*BB[i,]*t(BB[j,,drop=F])-a2*d2*BB[i,]*t(AB[j,,drop=F])-(a2^2)*BB[i,]*t(AA[j,,drop=F])
  }
  Z <- AAi+ABi+BBi
  s <- 0
  for (k in 1L:ncol(Z)) {
    s <- s + Z[,k]
  }
  return(s)
}
