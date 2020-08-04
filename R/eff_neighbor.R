#' Estimation of self and neighbor QTL effects across a genome
#'
#' A function to estimate additive and dominance deviation for self and neighbor QTL effects by a simple regression.
#' @param genoprobs Conditional genotype probabilities as taken from \code{qtl::calc.genoprob()}.
#' @param pheno A vector of individual phenotypes.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial position along an x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param addQTL An optional vector containing marker names that are considered covariates. Namely, this option allows composite interval mapping (Jansen 1993).
#' @param grouping An optional integer vector assigning each individual to a group. This argument can be used when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An optional argument to select trait types. The \code{"quantitative"} or \code{"binary"} calls the \code{"gaussian"} or \code{"binomial"} family in \code{glm()}, respectively.
#' @param fig TRUE/FALSE to plot the effects or not.
#' @param contrasts An optional vector composed of three TRUE/FALSE values, which represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB. If \code{NULL}, it is compiled from \code{genoprobs} automatically.
#' @return A matrix of estimated additive and dominance deviation for self and neighbor effects, with the chromosome numbers and positions. The row names correspond to marker names.
#' \itemize{
#'  \item{\code{chr}} {Chromosome number}
#'  \item{\code{pos}} {Marker position}
#'  \item{\code{a1}} {Additive deviation for self effects}
#'  \item{\code{d1}} {Dominance deviation for self effects}
#'  \item{\code{a2}} {Additive deviation for neighbor effects}
#'  \item{\code{d2}} {Dominance deviation for neighbor effects}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' Similar to Haley-Knott regression (Haley & Knott 1992), the additive and dominance deviations are approximated by a regression of trait values on conditional genotype probabilities.
#' The self QTL effects \code{a1} and \code{d1} are estimated in the same way as the \code{qtl} package performs the Haley-Knott regression.
#' If \code{contrasts = c(TRUE, TRUE, TRUE)}, neighbor QTL effects \code{a1} and \code{d1} are estimated using a quadratic regression; otherwise, the additive neighbor effects are estimated using a linear regression.
#' See also Sato, Takeda & Nagano (2020) for the rationale behind the approximation.
#' @references
#' * Haley CS, Knott SA (1992) A simple regression method for mapping quantitative trait loci in line crosses using flanking markers. Heredity 69:315-324.
#' * Jansen RC (1993) Interval mapping of multiple quantitative trait loci. Genetics 135:205-211.
#' * Sato Y, Takeda K, Nagano AJ (2020) Neighbor QTL: an interval mapping method for quantitative trait loci underlying neighbor effects. bioRxiv \url{https://doi.org/10.1101/2020.05.20.089474}
#' @examples
#' set.seed(1234)
#' test_map <- qtl::sim.map(len=rep(20,5),n.mar=3,include.x=FALSE)
#' test_cross <- qtl::sim.cross(test_map,n.ind=50)
#' test_smap <- cbind(runif(50,1,100),runif(50,1,100))
#' test_genoprobs <- qtl::calc.genoprob(test_cross,step=2)
#'
#' test_eff <- eff_neighbor(genoprobs=test_genoprobs,
#'                          pheno=test_cross$pheno$phenotype,
#'                          smap=test_smap, scale=20, fig=TRUE
#'                          )
#' @export
eff_neighbor = function(genoprobs, pheno, smap, scale, addcovar=NULL, addQTL=NULL, grouping=rep(1,nrow(smap)), response=c("quantitative","binary"), fig=TRUE, contrasts=NULL) {
  response <- match.arg(response)

  switch(response,
         "quantitative" = glm_family <- stats::gaussian(),
         "binary" = glm_family <- stats::binomial()
  )

  makeX = function(addcovar, addQTL, type) {
    X <- c()
    if(is.null(addcovar)==FALSE) {
      X <- cbind(X, addcovar)
    }
    if(is.null(addQTL)==FALSE) {
      if(type=="self") {
        X <- cbind(X, selfprobs[,addQTL])
      } else if(type=="neighbor") {
        X <- cbind(X, neiprobs[,addQTL])
      } else {
        stop("error: type must be 'self' or 'neighbor'")
      }
    }
    return(X)
  }

  self_k = function(k) {
    if(is.null(X)==TRUE) {
      res <- stats::glm.fit(cbind(1,selfprobs[,k]), pheno, family = glm_family)
      coef <- res$coefficients[length(res$coefficients)]
    } else {
      res <- stats::glm.fit(cbind(1,X,selfprobs[,k]), pheno, family = glm_family)
      coef <- res$coefficients[length(res$coefficients)]
      if (is.na(coef)) {
        if(is.null(addcovar)==TRUE) {
          res <- stats::glm.fit(cbind(1,selfprobs[,k]), pheno, family = glm_family)
        } else {
          res <- stats::glm.fit(cbind(1,addcovar,selfprobs[,k]), pheno, family = glm_family)
        }
        coef <- res$coefficients[length(res$coefficients)]
      }
    }
    return(coef)
  }

  selfprobs <- genoprobs2selfprobs(genoprobs=genoprobs, a1=1, d1=0, contrasts=contrasts)
  contrasts <- attr(selfprobs, "contrasts")
  q <- ncol(selfprobs)
  X <- makeX(addcovar=addcovar, addQTL=addQTL, type="self")
  a1 <- as.double(lapply(1:q, self_k))
  
  selfprobs <- genoprobs2selfprobs(genoprobs=genoprobs, a1=0, d1=1, contrasts=contrasts)
  X <- makeX(addcovar=addcovar, addQTL=addQTL, type="self")
  if(contrasts[2]==TRUE) {
    if(contrasts[3]==TRUE) {
      a1 <- a1/2
      d1 <- as.double(lapply(1:q, self_k))
    } else if(contrasts[3]==FALSE) {
      d1 <- a1
    }
  } else {
    a1 <- a1/2
    d1 <- rep(NA, q)
  }

  selfprobs <- genoprobs2selfprobs(genoprobs=genoprobs, a1=a1, d1=d1, contrasts=contrasts)
  nei_k <- function(k) {
    if(is.null(X)==TRUE) {
      res <- stats::glm.fit(cbind(1,selfprobs[,k],neiprobs[,k]), pheno, family = glm_family)
      coef <- res$coefficients[length(res$coefficients)]
    } else {
      res <- stats::glm.fit(cbind(1,X,selfprobs[,k],neiprobs[,k]), pheno, family = glm_family)
      coef <- res$coefficients[length(res$coefficients)]
      if (is.na(coef)) {
        if(is.null(addcovar)==TRUE) {
          res <- stats::glm.fit(cbind(1,selfprobs[,k],neiprobs[,k]), pheno, family = glm_family)
        } else {
          res <- stats::glm.fit(cbind(1,addcovar,selfprobs[,k],neiprobs[,k]), pheno, family = glm_family)
        }
        coef <- res$coefficients[length(res$coefficients)]
      }
    }
    return(coef)
  }

  nei_2k <- function(k) {
    if(is.null(X)==TRUE) {
      res <- stats::glm.fit(cbind(1,selfprobs[,k],stats::poly(neiprobs[,k],2)), pheno, family = glm_family)
      coef <- res$coefficients[-1:0+length(res$coefficients)]
    } else {
      res <- stats::glm.fit(cbind(1,X,selfprobs[,k],stats::poly(neiprobs[,k],2)), pheno, family = glm_family)
      coef <- res$coefficients[-1:0+length(res$coefficients)]
      if (is.na(coef[1])) {
        if(is.null(addcovar)==TRUE) {
          res <- stats::glm.fit(cbind(1,selfprobs[,k],stats::poly(neiprobs[,k],2)), pheno, family = glm_family)
        } else {
          res <- stats::glm.fit(cbind(1,addcovar,selfprobs[,k],stats::poly(neiprobs[,k],2)), pheno, family = glm_family)
        }
        coef <- res$coefficients[-1:0+length(res$coefficients)]
      }
    }
    return(coef)
  }

  if(contrasts[2]==TRUE) {
    if(contrasts[3]==TRUE) {
      neiprobs <- calc_neiprob(genoprobs=genoprobs, a2=1, d2=0.25, contrasts=contrasts, smap=smap, scale=scale, grouping=grouping, d2sq0=TRUE)
      X <- makeX(addcovar=addcovar, addQTL=addQTL, type="neighbor")
      coefs <- do.call(cbind, lapply(1:q, nei_2k))
      a2 <- sign(coefs[1,])*sqrt(abs(coefs[1,]/2))
      d2 <- sqrt(abs(coefs[2,]))
    } else if(contrasts[3]==FALSE) {
      neiprobs <- calc_neiprob(genoprobs=genoprobs, a2=1, d2=0, contrasts=contrasts, smap=smap, scale=scale, grouping=grouping)
      X <- makeX(addcovar=addcovar, addQTL=addQTL, type="neighbor")
      a2 <- as.double(lapply(1:q, nei_k))
      a2 <- sign(a2)*sqrt(abs(a2))
      d2 <- -a2
    }
  } else {
    neiprobs <- calc_neiprob(genoprobs=genoprobs, a2=1, d2=0, contrasts=contrasts, smap=smap, scale=scale, grouping=grouping)
    X <- makeX(addcovar=addcovar, addQTL=addQTL, type="neighbor")
    a2 <- as.double(lapply(1:q, nei_k))
    a2 <- sign(a2)*sqrt(abs(a2/2))
    d2 <- rep(NA, q)
  }

  marker_info <- get_markers(genoprobs=genoprobs)
  coeflist <- data.frame(marker_info, a1, d1, a2, d2)
  colnames(coeflist) <- c("chr","pos","a1","d1","a2","d2")

  if(fig==TRUE) {
    opar <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(opar))

    graphics::par(mfcol=c(2,1))
    plot_eff(res=coeflist, type="self")
    plot_eff(res=coeflist, type="neighbor")
  }
  return(coeflist)
}
