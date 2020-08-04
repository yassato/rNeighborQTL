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
#' @param contrasts An optional vector composed of three TRUE/FALSE values, which represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB. If \code{NULL}, it is compiled from \code{genoprobs} automatically.
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
#'                            smap=test_smap, scale=20
#'                            )
#' plot_nei(test_scan)
#' @export
scan_neighbor = function(genoprobs, pheno, smap, scale, addcovar=NULL, addQTL=NULL, grouping=rep(1,nrow(smap)), response=c("quantitative","binary"), contrasts=NULL) {
  response <- match.arg(response)

  switch(response,
         "quantitative" = glm_family <- stats::gaussian(),
         "binary" = glm_family <- stats::binomial()
  )

  p <- dim(genoprobs$geno[[1]]$prob)[1]
  geno <- decompose_genoprobs(genoprobs=genoprobs,contrasts=contrasts)
  contrasts <- attr(geno, "contrasts")

  scan_effect <- eff_neighbor(genoprobs=genoprobs, pheno=pheno, contrasts=contrasts, smap=smap, scale=scale, addcovar=addcovar, addQTL=addQTL, grouping=grouping, response=response, fig=FALSE)
  q <- nrow(scan_effect)
  p <- dim(genoprobs$geno[[1]]$prob)[1]

  y_self_hat <- genoprobs2selfprobs(geno, a1=scan_effect$a1, d1=scan_effect$d1)
  y_nei_hat <- calc_neiprob(geno, a2=scan_effect$a2, d2=scan_effect$d2, smap=smap, scale=scale, grouping=grouping)

  X <- c()
  if(is.null(addcovar)==FALSE) {
    X <- cbind(X, addcovar)
  }
  if(is.null(addQTL)==FALSE) {
    X <- cbind(X, y_self_hat[,match(addQTL, rownames(scan_effect))], y_nei_hat[,match(addQTL, rownames(scan_effect))])
  }

  if(is.null(addcovar)&is.null(addQTL)) {
    LOD_self <- c(); LOD_nei <- c()
    LL_null <- logLik_glm.fit(rep(1,length(pheno)),pheno,family=glm_family)
    for(k in 1:q) {
      LL_self <- logLik_glm.fit(cbind(1,y_self_hat[,k]),pheno,family=glm_family)
      LL_nei <- logLik_glm.fit(cbind(1,y_self_hat[,k],y_nei_hat[,k]),pheno,family=glm_family)
      LOD_self <- c(LOD_self, log10(exp(LL_self-LL_null)))
      LOD_nei <- c(LOD_nei, log10(exp(LL_nei-LL_self)))
    }
  } else if(is.null(addQTL)==TRUE) {
    LOD_self <- c(); LOD_nei <- c()
    LL_null <- logLik_glm.fit(cbind(1,X),pheno,family=glm_family)
    for(k in 1:q) {
      LL_self <- logLik_glm.fit(cbind(1,X,y_self_hat[,k]),pheno,family=glm_family)
      LL_nei <- logLik_glm.fit(cbind(1,X,y_self_hat[,k],y_nei_hat[,k]),pheno,family=glm_family)
      LOD_self <- c(LOD_self, log10(exp(LL_self-LL_null)))
      LOD_nei <- c(LOD_nei, log10(exp(LL_nei-LL_self)))
    }
  } else {
    X <- cbind(y_self_hat[,match(addQTL, rownames(scan_effect))], y_nei_hat[,match(addQTL, rownames(scan_effect))])
    LOD_self <- c(); LOD_nei <- c()
    LL_null <- logLik_glm.fit(cbind(1,X),pheno,family=glm_family)
    for(k in 1:q) {
      LL_self <- logLik_glm.fit(cbind(1,X,y_self_hat[,k]),pheno,family=glm_family)
      LL_nei <- logLik_glm.fit(cbind(1,X,y_self_hat[,k],y_nei_hat[,k]),pheno,family=glm_family)
      LOD_self <- c(LOD_self, log10(exp(LL_self-LL_null)))
      LOD_nei <- c(LOD_nei, log10(exp(LL_nei-LL_self)))
    }
  }

  marker_info <- get_markers(genoprobs=genoprobs)
  LODlist <- data.frame(marker_info, LOD_self, LOD_nei)
  colnames(LODlist) <- c("chr","pos","LOD_self","LOD_nei")

  return(LODlist)
}
