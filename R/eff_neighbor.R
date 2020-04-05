#' Estimated self and neighbor QTL effects across a genome
#'
#' A function to estimate additive and dominance deviation for self and neighbor QTL effects by Haley-Knott regression
#' @param genoprobs Conditional genotype probabilities as obtained from \code{qtl2::calc_genoprob()}.
#' @param pheno A vector of individual phenotypes.
#' @param gmap Genetic map including observed and pseudomarkers, as obtained from \code{qtl2::insert_pseudomarkers()}.
#' @param contrasts  A vector composed of three TRUE/FALSE values. Depending on crossing design, it represents the presence/absence of specific genotypes as c(T/F, T/F, T/F) = AA, AB, BB.
#' @param smap A matrix showing a spatial map. The first and second column include spatial points along a x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param addQTL An optional vector containing marker names that are considered covariates. Namely, this option allows composite interval mapping (Jansen 1993).
#' @param grouping An optional integer vector assigning each individual to a group. This argument can be useful when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An optional argument to select trait types. The \code{"quantitative"}, \code{"binary"}, or \code{"count"} calls the \code{"gaussian"}, \code{"binomial"}, or \code{"poisson"} family in \code{glm()}, respectively.
#' @param fig TRUE/FALSE to plot the effects or not.
#' @return  A matrix of estimated additive and dominance deviation for self and neighbor effects, with the chromosome numbers and positions. The row names correspond to marker names.
#' \itemize{
#'  \item{\code{chr}} {Chromosome number}
#'  \item{\code{pos}} {Marker position}
#'  \item{\code{a1}} {Additive deviation for self effects}
#'  \item{\code{d1}} {Dominance deviation for self effects}
#'  \item{\code{a2}} {Additive deviation for neighbor effects}
#'  \item{\code{d2}} {Dominance deviation for neighbor effects}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @references
#' * Haley CS, Knott SA (1992) A simple regression method for mapping quantitative trait loci in line crosses using flanking markers. Heredity 69:315-324.
#' * Jansen RC (1993) Interval mapping of multiple quantitative trait loci. Genetics 135:205-211.
#' @examples
#' #demo using F2
#' set.seed(1)
#' data("fake.f2",package="qtl")
#' fake_f2 <- qtl2::convert2cross2(fake.f2)
#' fake_f2 <- subset(fake_f2,chr=c(1:19))
#' smap_f2 <- cbind(runif(qtl2::n_ind(fake_f2),1,100),runif(qtl2::n_ind(fake_f2),1,100))
#' gmap_f2 <- qtl2::insert_pseudomarkers(fake_f2$gmap, step=1)
#' genoprobs_f2 <- qtl2::calc_genoprob(fake_f2,gmap_f2,error_prob=0.002)
#' eff_f2_1 <- eff_neighbor(genoprobs=genoprobs_f2,
#'                            pheno=fake_f2$pheno[,1],
#'                            gmap=gmap_f2, smap=smap_f2,
#'                            scale=20, addcovar=as.matrix(fake_f2$covar)
#'                            )
#' @export
eff_neighbor = function(genoprobs, pheno, gmap, contrasts=c(TRUE,TRUE,TRUE), smap, scale, addcovar=NULL, addQTL=NULL, grouping=rep(1,nrow(smap)), response="quantitative", fig=TRUE) {

  switch(response,
         "quantitative" = glm_family <- "gaussian",
         "binary" = glm_family <- "binomial",
         "count" = glm_family <- "poisson",
         print("response must be 'quantitative', 'binary', or 'count'")
  )

  self_k = function(k) {
    if(is.null(X)==TRUE) {
      res = stats::glm(pheno~selfprobs[,k]-1, family=glm_family)
      coef = res$coef[length(res$coef)]
    } else {
      res = try(stats::glm(pheno~X+selfprobs[,k]-1, family=glm_family))
      if(is.na(res$coef[length(res$coef)])!=TRUE) {
        coef = res$coef[length(res$coef)]
      } else {
        if(is.null(addcovar)==TRUE) {
          res = stats::glm(pheno~selfprobs[,k]-1, family=glm_family)
        } else {
          res = stats::glm(pheno~addcovar+selfprobs[,k]-1, family=glm_family)
        }
        coef =res$coef[length(res$coef)]
      }
    }
    return(coef)
  }

  #estimate a1
  selfprobs = genoprobs2selfprobs(genoprobs=genoprobs, gmap=gmap, a1=1, d1=0, contrasts=contrasts)
  q = ncol(selfprobs)

  X = c()
  if(is.null(addcovar)==FALSE) {
    X = cbind(X, addcovar)
  }
  if(is.null(addQTL)==FALSE) {
    X = cbind(X, selfprobs[,addQTL])
  }
  a1 = mapply(self_k, 1:q)

  #estimate d1
  selfprobs = genoprobs2selfprobs(genoprobs=genoprobs, gmap=gmap, a1=0, d1=1, contrasts=contrasts)
  X = c()
  if(is.null(addcovar)==FALSE) {
    X = cbind(X, addcovar)
  }
  if(is.null(addQTL)==FALSE) {
    X = cbind(X, selfprobs[,addQTL])
  }
  if(contrasts[2]==TRUE) {
    d1 = mapply(self_k, 1:q)
  } else {
    d1 = rep(NA, q)
  }

  selfprobs = genoprobs2selfprobs(genoprobs=genoprobs, gmap=gmap, a1=a1, d1=d1, contrasts=contrasts)
  nei_k = function(k) {
    if(is.null(X)==TRUE) {
      res = stats::glm(pheno~selfprobs[,k]+neiprobs[,k]-1, family=glm_family)
      coef = res$coef[length(res$coef)]
    } else {
      res = try(stats::glm(pheno~X+selfprobs[,k]+neiprobs[,k]-1, family=glm_family))
      if(is.na(res$coef[length(res$coef)])!=TRUE) {
        coef = res$coef[length(res$coef)]
      } else {
        if(is.null(addcovar)==TRUE) {
          res = stats::glm(pheno~selfprobs[,k]+neiprobs[,k]-1, family=glm_family)
        } else {
          res = stats::glm(pheno~addcovar+selfprobs[,k]+neiprobs[,k]-1, family=glm_family)
        }
        coef = res$coef[length(res$coef)]
      }
    }
    return(coef)
  }

  #estimate a2
  neiprobs = calc_neiprob(genoprobs=genoprobs, gmap=gmap, a2=1, d2=0, contrasts=contrasts, smap=smap, scale=scale, grouping=grouping)
  X = c()
  if(is.null(addcovar)==FALSE) {
    X = cbind(X, addcovar)
  }
  if(is.null(addQTL)==FALSE) {
    X = cbind(X, neiprobs[,addQTL])
  }
  a2 = mapply(nei_k, 1:q)

  #estimate d2
  neiprobs = calc_neiprob(genoprobs=genoprobs, gmap=gmap, a2=0, d2=1, contrasts=contrasts, smap=smap, scale=scale, grouping=grouping)
  X = c()
  if(is.null(addcovar)==FALSE) {
    X = cbind(X, addcovar)
  }
  if(is.null(addQTL)==FALSE) {
    X = cbind(X, neiprobs[,addQTL])
  }
  if(contrasts[2]==TRUE) {
    d2 = mapply(nei_k, 1:q)
  } else {
    d2 = rep(NA, q)
  }

  marker_info = get_markers(gmap)
  coeflist = data.frame(marker_info, a1, d1, a2, d2)
  coeflist[,5] = sign(coeflist[,5])*sqrt(abs(coeflist[,5]))
  coeflist[,6] = sign(coeflist[,6])*sqrt(abs(coeflist[,6]))
  colnames(coeflist) = c("chr","pos","a1","d1","a2","d2")

  if(fig==TRUE) {
    graphics::par(mfcol=c(2,1))
    plot_eff(res=coeflist, type="self")
    plot_eff(res=coeflist, type="neighbor")
  }
  return(coeflist)
}
