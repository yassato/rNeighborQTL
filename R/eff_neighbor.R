#' Estimation of self and neighbor QTL effects across a genome
#'
#' A function to estimate additive and dominance deviation for self and neighbor QTL effects by a simple regression.
#' @param genoprobs Conditional genotype probabilities as obtained from \code{qtl2::calc_genoprob()}.
#' @param pheno A vector of individual phenotypes.
#' @param gmap Genetic map including observed and pseudomarkers, as obtained from \code{qtl2::insert_pseudomarkers()}.
#' @param contrasts A vector composed of three TRUE/FALSE values. Depending on crossing design, it represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB.
#' @param smap A matrix showing a spatial map. The first and second column include spatial points along a x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param addQTL An optional vector containing marker names that are considered covariates. Namely, this option allows composite interval mapping (Jansen 1993).
#' @param grouping An optional integer vector assigning each individual to a group. This argument can be useful when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An optional argument to select trait types. The \code{"quantitative"} or \code{"binary"} calls the \code{"gaussian"} or \code{"binomial"} family in \code{glm()}, respectively.
#' @param fig TRUE/FALSE to plot the effects or not.
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
#' The self QTL effects \code{a1} and \code{d1} are esimated in the same way as the \code{qtl} package performs the Haley-Knott regression.
#' If \code{contrasts = c(TRUE, TRUE, TRUE)}, neighbor QTL effects \code{a1} and \code{d1} are estimated using a quadratic regression; otherwise, the additive neighbor effects are estimated using a linear regression.
#' @references
#' * Haley CS, Knott SA (1992) A simple regression method for mapping quantitative trait loci in line crosses using flanking markers. Heredity 69:315-324.
#' * Jansen RC (1993) Interval mapping of multiple quantitative trait loci. Genetics 135:205-211.
#' @examples
#' set.seed(1234)
#' data("fake.f2",package="qtl")
#' fake_f2 <- qtl2::convert2cross2(fake.f2)
#' fake_f2 <- subset(fake_f2,chr=c(1:19))
#' smap_f2 <- cbind(runif(qtl2::n_ind(fake_f2),1,100),runif(qtl2::n_ind(fake_f2),1,100))
#' gmap_f2 <- qtl2::insert_pseudomarkers(fake_f2$gmap, step=2)
#' genoprobs_f2 <- qtl2::calc_genoprob(fake_f2,gmap_f2)
#' eff_f2 <- eff_neighbor(genoprobs=genoprobs_f2,
#'                        pheno=fake_f2$pheno[,1],
#'                        gmap=gmap_f2, smap=smap_f2,
#'                        scale=19.37, addcovar=as.matrix(fake_f2$covar)
#'                        )
#' @export
eff_neighbor = function(genoprobs, pheno, gmap, contrasts=c(TRUE,TRUE,TRUE), smap, scale, addcovar=NULL, addQTL=NULL, grouping=rep(1,nrow(smap)), response="quantitative", fig=TRUE) {

  switch(response,
         "quantitative" = glm_family <- "gaussian",
         "binary" = glm_family <- "binomial",
         stop("error: response must be 'quantitative' or 'binary'")
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
      res <- stats::glm(pheno~selfprobs[,k], family=glm_family)
      coef <- res$coef[length(res$coef)]
    } else {
      res <- try(stats::glm(pheno~X+selfprobs[,k], family=glm_family))
      if(is.na(res$coef[length(res$coef)])!=TRUE) {
        coef <- res$coef[length(res$coef)]
      } else {
        if(is.null(addcovar)==TRUE) {
          res <- stats::glm(pheno~selfprobs[,k], family=glm_family)
        } else {
          res <- stats::glm(pheno~addcovar+selfprobs[,k], family=glm_family)
        }
        coef <- res$coef[length(res$coef)]
      }
    }
    return(coef)
  }

  selfprobs <- genoprobs2selfprobs(genoprobs=genoprobs, gmap=gmap, a1=1, d1=0, contrasts=contrasts)
  q <- ncol(selfprobs)
  X <- makeX(addcovar=addcovar, addQTL=addQTL, type="self")
  a1 <- mapply(self_k, 1:q)

  selfprobs <- genoprobs2selfprobs(genoprobs=genoprobs, gmap=gmap, a1=0, d1=1, contrasts=contrasts)
  X <- makeX(addcovar=addcovar, addQTL=addQTL, type="self")
  if(contrasts[2]==TRUE) {
    if(contrasts[3]==TRUE) {
      a1 <- a1/2
      d1 <- mapply(self_k, 1:q)
    } else if(contrasts[3]==FALSE) {
      d1 <- a1
    }
  } else {
    a1 <- a1/2
    d1 <- rep(NA, q)
  }


  nei_k <- function(k) {
    if(is.null(X)==TRUE) {
      res <- stats::glm(pheno~selfprobs[,k]+neiprobs[,k], family=glm_family)
      coef <- res$coef[length(res$coef)]
    } else {
      res <- try(stats::glm(pheno~X+selfprobs[,k]+neiprobs[,k], family=glm_family))
      if(is.na(res$coef[length(res$coef)])!=TRUE) {
        coef <- res$coef[length(res$coef)]
      } else {
        if(is.null(addcovar)==TRUE) {
          res <- stats::glm(pheno~selfprobs[,k]+neiprobs[,k], family=glm_family)
        } else {
          res <- stats::glm(pheno~addcovar+selfprobs[,k]+neiprobs[,k], family=glm_family)
        }
        coef <- res$coef[length(res$coef)]
      }
    }
    return(coef)
  }

  nei_2k <- function(k) {
    if(is.null(X)==TRUE) {
      res <- stats::glm(pheno~selfprobs[,k]+poly(neiprobs[,k],2), family=glm_family)
      coef <- res$coef[c((length(res$coef)-1),length(res$coef))]
    } else {
      res <- try(stats::glm(pheno~X+selfprobs[,k]+poly(neiprobs[,k],2), family=glm_family))
      if(is.na(res$coef[length(res$coef)-1])!=TRUE) {
        coef <- res$coef[c((length(res$coef)-1),length(res$coef))]
      } else {
        if(is.null(addcovar)==TRUE) {
          res <- stats::glm(pheno~selfprobs[,k]+poly(neiprobs[,k],2), family=glm_family)
        } else {
          res <- stats::glm(pheno~addcovar+selfprobs[,k]+poly(neiprobs[,k],2), family=glm_family)
        }
        coef <- res$coef[c((length(res$coef)-1),length(res$coef))]
      }
    }
    return(coef)
  }

  if(contrasts[2]==TRUE) {
    if(contrasts[3]==TRUE) {
      neiprobs <- calc_neiprob(genoprobs=genoprobs, gmap=gmap, a2=1, d2=0.25, contrasts=contrasts, smap=smap, scale=scale, grouping=grouping, d2sq0=TRUE)
      X <- makeX(addcovar=addcovar, addQTL=addQTL, type="neighbor")
      coefs <- mapply(nei_2k, 1:q)
      a2 <- sign(coefs[1,])*sqrt(abs(coefs[1,]/2))
      d2 <- sqrt(abs(coefs[2,]))
    } else if(contrasts[3]==FALSE) {
      neiprobs <- calc_neiprob(genoprobs=genoprobs, gmap=gmap, a2=1, d2=0, contrasts=contrasts, smap=smap, scale=scale, grouping=grouping)
      X <- makeX(addcovar=addcovar, addQTL=addQTL, type="neighbor")
      a2 <- mapply(nei_k, 1:q)
      a2 <- sign(a2)*sqrt(abs(a2))
      d2 <- -a2
    }
  } else {
    neiprobs <- calc_neiprob(genoprobs=genoprobs, gmap=gmap, a2=1, d2=0, contrasts=contrasts, smap=smap, scale=scale, grouping=grouping)
    X <- makeX(addcovar=addcovar, addQTL=addQTL, type="neighbor")
    a2 <- mapply(nei_k, 1:q)
    a2 <- sign(a2)*sqrt(abs(a2/2))
    d2 <- rep(NA, q)
  }

  marker_info <- get_markers(gmap)
  coeflist <- data.frame(marker_info, a1, d1, a2, d2)
  colnames(coeflist) <- c("chr","pos","a1","d1","a2","d2")

  if(fig==TRUE) {
    graphics::par(mfcol=c(2,1))
    plot_eff(res=coeflist, type="self")
    plot_eff(res=coeflist, type="neighbor")
  }
  return(coeflist)
}
