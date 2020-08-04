#' Decomposition of conditional genotype probabilities
#'
#' A function to decompose \code{qtl}'s object of conditional genotype probabilities.
#' @param genoprobs Conditional genotype probabilities as taken from \code{qtl::calc.genoprob()}.
#' @param contrasts A vector composed of three TRUE/FALSE values, which represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB.
#' @return A list of three numeric matrices for genotype probabilities AA, AB, and BB. Each contains elements of individuals x markers.
#' \itemize{
#'  \item{\code{AA}} {Homozygote AA probabilities.}
#'  \item{\code{AB}} {Heterozygote AB probabilities for. NA if inbred lines}
#'  \item{\code{BB}} {Homozygote BB probabilities. NA if backcross lines}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
decompose_genoprobs = function(genoprobs, contrasts=NULL) {
  p <- dim(genoprobs$geno[[1]]$prob)[1]

  if (is.null(contrasts)) {
    if (!all(dimnames(genoprobs$geno[[1]]$prob)[[3]] %in% c("AA", "AB", "BB"))) {
      stop("error: genoprobs type error!")
    }
    contrasts = c("AA", "AB", "BB") %in% dimnames(genoprobs$geno[[1]]$prob)[[3]]
  } else {
    r <- dim(genoprobs$geno[[1]]$prob)[3]
    if(r!=sum(contrasts)) {
      stop("error: allele dimension does not match!")
    }
  }

  geno <- c()
  AA <- c(); AB <- c(); BB <- c()

  if(contrasts[1]==TRUE) {
    for(chr in 1:length(genoprobs$geno)) {
      q <- dim(genoprobs$geno[[chr]]$prob)[2]
      AA_chr <- matrix(genoprobs$geno[[chr]]$prob[1:(p*q)],p,q)
      AA <- cbind(AA, AA_chr)
    }
  } else { AA <- NA }

  if(contrasts[2]==TRUE) {
    for(chr in 1:length(genoprobs$geno)) {
      q <- dim(genoprobs$geno[[chr]]$prob)[2]
      AB_chr <- matrix(genoprobs$geno[[chr]]$prob[(p*q+1):(2*p*q)],p,q)
      AB <- cbind(AB, AB_chr)
    }

    if(contrasts[3]==TRUE) {
      for(chr in 1:length(genoprobs$geno)) {
        q <- dim(genoprobs$geno[[chr]]$prob)[2]
        BB_chr <- matrix(genoprobs$geno[[chr]]$prob[(2*p*q+1):(3*p*q)],p,q)
        BB <- cbind(BB, BB_chr)
      }
    } else { BB <- NA }

  } else if(contrasts[3]==TRUE) {
    AB <- NA
    for(chr in 1:length(genoprobs$geno)) {
      q <- dim(genoprobs$geno[[chr]]$prob)[2]
      BB_chr <- matrix(genoprobs$geno[[chr]]$prob[(p*q+1):(2*p*q)],p,q)
      BB <- cbind(BB, BB_chr)
    }
  } else { BB <- NA }

  geno[[1]] <- AA
  geno[[2]] <- AB
  geno[[3]] <- BB

  names(geno) <- c("AA","AB","BB")
  attr(geno, "contrasts") <- contrasts
  attr(geno, "marker_info") <- get_markers(genoprobs=genoprobs)
  return(geno)
}
