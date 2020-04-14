#' Decomposition of conditional genotype probabilities
#'
#' A function to reshape \code{qtl2}'s object of conditional genotype probabilities.
#' @param genoprobs Conditional genotype probabilities as obtained from \code{qtl2::calc_genoprob()}.
#' @param contrasts A vector composed of three TRUE/FALSE values. Depending on crossing design, it represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB.
#' @return A list of three numeric matrices for genotype probabilities AA, AB, and BB. Each contains elements of individuals x markers.
#' \itemize{
#'  \item{\code{AA}} {Homozygote AA probabilities.}
#'  \item{\code{AB}} {Heterozygote AB probabilities for. NA if inbred lines}
#'  \item{\code{BB}} {Homozygote BB probabilities. NA if backcross lines}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
decompose_genoprobs = function(genoprobs, contrasts=c(TRUE,TRUE,TRUE)) {
  p <- dim(genoprobs[[1]])[1]
  r <- dim(genoprobs[[1]])[2]

  if(r!=sum(contrasts)) {
    stop("error: allele dimension does not match!")
  }

  geno <- c()
  AA <- c(); AB <- c(); BB <- c()

  if(contrasts[1]) {
    for(chr in 1:length(genoprobs)) {
      q <- dim(genoprobs[[chr]])[3]
      for(i in r*(0:(q-1))) {
        AA <- cbind(AA, genoprobs[[chr]][seq(i*p+1,i*p+p,1)])
      }
    }
  } else { AA <- NA }

  if(contrasts[2]) {
    for(chr in 1:length(genoprobs)) {
      q <- dim(genoprobs[[chr]])[3]
      for(i in 1+r*(0:(q-1))) {
        AB <- cbind(AB, genoprobs[[chr]][seq(i*p+1,i*p+p,1)])
      }
    }

    if(contrasts[3]) {
      for(chr in 1:length(genoprobs)) {
        q <- dim(genoprobs[[chr]])[3]
        for(i in 2+r*(0:(q-1))) {
          BB <- cbind(BB, genoprobs[[chr]][seq(i*p+1,i*p+p,1)])
        }
      }
    } else { BB <- NA }

  } else if(contrasts[3]) {
    AB <- NA
    for(chr in 1:length(genoprobs)) {
      q <- dim(genoprobs[[chr]])[3]
      for(i in 1+r*(0:(q-1))) {
        BB <- cbind(BB, genoprobs[[chr]][seq(i*p+1,i*p+p,1)])
      }
    }
  } else { BB <- NA }

  geno[[1]] <- AA
  geno[[2]] <- AB
  geno[[3]] <- BB

  names(geno) <- c("AA","AB","BB")
  return(geno)
}
