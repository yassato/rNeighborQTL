#' Genome scan for neighbor effects with a QTL model
#'
#' Genome scan with a QTL model by Haley-Knott regression for self and neighbor effects, with possible allowance for additional covariates and non-normal traits.
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
#' @return A matrix of LOD scores for self and neighbor effects, with the chromosome numbers and positions. The row names correspond to marker names.
#' \itemize{
#'  \item{\code{chr}} {Chromosome number}
#'  \item{\code{pos}} {Marker position}
#'  \item{\code{LOD_self}} {LOD score for self effects}
#'  \item{\code{LOD_nei}} {LOD score for neighbor effects}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' This function adopts a stepwise testing from self to neighbor effects; thus, \code{LOD_self} gives the same results as standard QTL mapping.
#' Note that the results return 0 LOD scores for covariate markers when using \code{addQTL} option.
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
#' s_seq <- quantile(dist(smap_f2),c(0.1*(1:10)))
#' scan_f2 <- scan_neighbor(genoprobs=genoprobs_f2,
#'                          pheno=fake_f2$pheno[,1], gmap=gmap_f2,
#'                          contrasts = c(TRUE,TRUE,TRUE), smap=smap_f2,
#'                          scale=19.37, addcovar=as.matrix(fake_f2$covar))
#' plot_nei(scan_f2)
#' @export
scan_neighbor = function(genoprobs, pheno, gmap, contrasts=c(TRUE,TRUE,TRUE), smap, scale, addcovar=NULL, addQTL=NULL, grouping=rep(1,nrow(smap)), response="quantitative") {

  switch(response,
         "quantitative" = glm_family <- "gaussian",
         "binary" = glm_family <- "binomial",
         stop("error: response must be 'quantitative' or 'binary'")
  )

  p <- dim(genoprobs[[1]])[1]
  geno <- decompose_genoprobs(genoprobs=genoprobs,contrasts=contrasts)

  scan_effect <- eff_neighbor(genoprobs=genoprobs, pheno=pheno, gmap=gmap, contrasts=contrasts, smap=smap, scale=scale, addcovar=addcovar, addQTL=addQTL, grouping=grouping, response=response, fig=FALSE)
  q <- nrow(scan_effect)
  p <- dim(genoprobs[[1]])[1]

  y_self_hat <- c()
  for(i in 1:p) y_self_hat <- rbind(y_self_hat, selfprob(i, a1=scan_effect$a1, d1=scan_effect$d1, AA=geno$AA, AB=geno$AB, BB=geno$BB))

  neiprob_i = function(i) {
    id = c(1:p)[grouping == grouping[i]]

    d_i = mapply(function(x) { return(sqrt((smap[x,1]-smap[i,1])^2 + (smap[x,2]-smap[i,2])^2)) },id)
    prob_i = 0
    j_id = id[(d_i>0)&(d_i<=scale)]
    if(length(j_id)==0) {
      return(rep(0,ncol(geno$AA)))
    } else {
      for(j in j_id){
        prob_ij = neiprob(i=i, j=j, a2=scan_effect$a2, d2=scan_effect$d2, AA=geno$AA, AB=geno$AB, BB=geno$BB)
        prob_i = prob_i + prob_ij
      }
      prob_i = prob_i/length(j_id)
      return(prob_i)
    }
  }

  y_nei_hat <- mapply(neiprob_i, 1:p)
  y_nei_hat <- t(y_nei_hat)

  X <- c()
  if(is.null(addcovar)==FALSE) {
    X <- cbind(X, addcovar)
  }
  if(is.null(addQTL)==FALSE) {
    X <- cbind(X, y_self_hat[,match(addQTL, rownames(scan_effect))], y_nei_hat[,match(addQTL, rownames(scan_effect))])
  }

  if(is.null(addcovar)&is.null(addQTL)) {
    LOD_self <- c(); LOD_nei <- c()
    LL_null <- stats::logLik(stats::glm(pheno~1, family=glm_family))
    for(k in 1:q) {
      LL_self <- stats::logLik(stats::glm(pheno~y_self_hat[,k], family=glm_family))
      LL_nei <- stats::logLik(stats::glm(pheno~y_self_hat[,k]+y_nei_hat[,k], family=glm_family))
      LOD_self <- c(LOD_self, log10(exp(LL_self-LL_null)))
      LOD_nei <- c(LOD_nei, log10(exp(LL_nei-LL_self)))
    }
  } else if(is.null(addQTL)==TRUE) {
    LOD_self <- c(); LOD_nei <- c()
    LL_null <- stats::logLik(stats::glm(pheno~X, family=glm_family))
    for(k in 1:q) {
      LL_self <- stats::logLik(stats::glm(pheno~X+y_self_hat[,k], family=glm_family))
      LL_nei <- stats::logLik(stats::glm(pheno~X+y_self_hat[,k]+y_nei_hat[,k], family=glm_family))
      LOD_self <- c(LOD_self, log10(exp(LL_self-LL_null)))
      LOD_nei <- c(LOD_nei, log10(exp(LL_nei-LL_self)))
    }
  } else {
    X <- cbind(y_self_hat[,match(addQTL, rownames(scan_effect))], y_nei_hat[,match(addQTL, rownames(scan_effect))])
    LOD_self <- c(); LOD_nei <- c()
    LL_null <- stats::logLik(stats::glm(pheno~X, family=glm_family))
    for(k in 1:q) {
      LL_self <- stats::logLik(stats::glm(pheno~X+y_self_hat[,k], family=glm_family))
      LL_nei <- stats::logLik(stats::glm(pheno~X+y_self_hat[,k]+y_nei_hat[,k], family=glm_family))
      LOD_self <- c(LOD_self, log10(exp(LL_self-LL_null)))
      LOD_nei <- c(LOD_nei, log10(exp(LL_nei-LL_self)))
    }
  }

  marker_info <- get_markers(gmap)
  LODlist <- data.frame(marker_info, LOD_self, LOD_nei)
  colnames(LODlist) <- c("chr","pos","LOD_self","LOD_nei")

  return(LODlist)
}
