#' Genome scan for neighbor effects with a QTL model
#'
#' Genome scan using a QTL model for self and neighbor effects, with possible allowance for additional covariates and non-normal traits. Theoretical background is described in Sato, Takeda & Nagano (2020).
#' @param genoprobs Conditional genotype probabilities as taken from \code{qtl::calc.genoprob()}.
#' @param pheno A vector of individual phenotypes.
#' @param contrasts A vector composed of three TRUE/FALSE values. Depending on the crossing design, it represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial positions along an x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param addQTL An optional vector containing marker names that are considered covariates. Namely, this option allows composite interval mapping (Jansen 1993).
#' @param grouping An optional integer vector assigning each individual to a group. This argument can be used when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An optional argument to select trait types. The \code{"quantitative"} or \code{"binary"} calls the \code{"gaussian"} or \code{"binomial"} family in \code{glm()}, respectively.
#' @return A matrix of LOD scores for self and neighbor effects, with the chromosome numbers and positions. The row names correspond to marker names.
#' \itemize{
#'  \item{\code{chr}} {Chromosome number}
#'  \item{\code{pos}} {Marker position}
#'  \item{\code{LOD_self}} {LOD score for self effects}
#'  \item{\code{LOD_nei}} {LOD score for neighbor effects}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @seealso \code{\link{eff_neighbor}}
#' @details
#' This function calculates LOD score after the additive and dominance deviation are estimated using \code{eff_neighbor()}.
#' As it adopts a stepwise testing from self to neighbor effects, \code{LOD_self} are the same as standard QTL mapping.
#' Note that the results return 0 LOD scores for covariate markers when using \code{addQTL} option.
#' @references
#' * Jansen RC (1993) Interval mapping of multiple quantitative trait loci. Genetics 135:205-211.
#' * Sato Y, Takeda K, Nagano AJ (2020) Neighbor QTL: an interval mapping method for quantitative trait loci underlying neighbor effects. bioRxiv \url{https://doi.org/10.1101/2020.05.20.089474}
#' @examples
#' set.seed(1234)
#' test_map <- qtl::sim.map(len=rep(20,5),n.mar=3,include.x=FALSE)
#' test_cross <- qtl::sim.cross(test_map,n.ind=50)
#' test_smap <- cbind(runif(50,1,100),runif(50,1,100))
#' test_genoprobs <- qtl::calc.genoprob(test_cross,step=2)
#'
#' test_scan <- scan_neighbor(genoprobs=test_genoprobs,
#'                            pheno=test_cross$pheno$phenotype,
#'                            contrasts=c(TRUE,TRUE,TRUE),
#'                            smap=test_smap, scale=20
#'                            )
#' plot_nei(test_scan)
#' @export
scan_neighbor = function(genoprobs, pheno, contrasts=c(TRUE,TRUE,TRUE), smap, scale, addcovar=NULL, addQTL=NULL, grouping=rep(1,nrow(smap)), response="quantitative") {

  switch(response,
         "quantitative" = glm_family <- "gaussian",
         "binary" = glm_family <- "binomial",
         stop("error: response must be 'quantitative' or 'binary'")
  )

  p <- dim(genoprobs$geno[[1]]$prob)[1]
  geno <- decompose_genoprobs(genoprobs=genoprobs,contrasts=contrasts)

  scan_effect <- eff_neighbor(genoprobs=genoprobs, pheno=pheno, contrasts=contrasts, smap=smap, scale=scale, addcovar=addcovar, addQTL=addQTL, grouping=grouping, response=response, fig=FALSE)
  q <- nrow(scan_effect)
  p <- dim(genoprobs$geno[[1]]$prob)[1]

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

  marker_info <- get_markers(genoprobs=genoprobs)
  LODlist <- data.frame(marker_info, LOD_self, LOD_nei)
  colnames(LODlist) <- c("chr","pos","LOD_self","LOD_nei")

  return(LODlist)
}
