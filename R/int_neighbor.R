#' Testing marker-by-marker epistasis in neighbor QTL effects
#'
#' A function to test interaction terms between one focal marker and the other markers across a genome.
#' @param genoprobs Conditional genotype probabilities as taken from \code{qtl::calc.genoprob()}.
#' @param pheno A vector of individual phenotypes.
#' @param contrasts A vector composed of three TRUE/FALSE values. Depending on the crossing design, it represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial positions along an x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param addQTL A vector containing marker names that are considered covariates. This argument is necessary for \code{int_neighbor()}, and must match the marker names of \code{gmap}.
#' @param intQTL A name of a focal marker to be tested for its epistasis with the other markers in neighbor effects. The marker name must be included by \code{addQTL}.
#' @param grouping An optional integer vector assigning each individual to a group. This argument can be used when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An optional argument to select trait types. The \code{"quantitative"} or \code{"binary"} calls the \code{"gaussian"} or \code{"binomial"} family in \code{glm()}, respectively.
#' @return A matrix of LOD scores for neighbor epistasis effects, with the chromosome numbers and positions. The row names correspond to marker names.
#' \itemize{
#'  \item{\code{chr}} {Chromosome number}
#'  \item{\code{pos}} {Marker position}
#'  \item{\code{LOD_int}} {LOD score for epistasis in neighbor effects between a focal and the other markers}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' This is an optimal function to test two-way interactions between the main neighbor effect of a focal marker given by \code{intQTL} and the others.
#' All the main neighbor effects are first estimated using \code{eff_neighbor()}, and then a two-way interaction term between the focal marker effect and its counterpart was considered an additional explanatory variable.
#' LOD score was compared between models with or without the two-way interaction.
#' @seealso \code{\link{scan_neighbor}} \code{\link{eff_neighbor}}
#' @examples
#' set.seed(1234)
#' test_map <- qtl::sim.map(len=rep(20,5),n.mar=3,include.x=FALSE)
#' test_cross <- qtl::sim.cross(test_map,n.ind=50)
#' test_smap <- cbind(runif(50,1,100),runif(50,1,100))
#' test_genoprobs <- qtl::calc.genoprob(test_cross,step=2)
#'
#' test_int <- int_neighbor(genoprobs=test_genoprobs,
#'                          pheno=test_cross$pheno$phenotype,
#'                          smap=test_smap,scale=20,
#'                          addQTL=c("c1_D1M1","c1_D1M2"),intQTL="c1_D1M1"
#'                          )
#' plot_nei(test_int, type="int")
#' @export
int_neighbor = function(genoprobs, pheno, contrasts = c(TRUE, TRUE, TRUE), smap, scale, addcovar=NULL, addQTL, intQTL, grouping=rep(1,nrow(smap)), response="quantitative") {

  if(is.na(match(intQTL, addQTL))) {
    warning("A 'intQTL' marker must overlap with 'addQTL'")
    return(NULL)
  }

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

  X <- cbind(y_self_hat[,match(addQTL, rownames(scan_effect))], y_nei_hat[,match(addQTL, rownames(scan_effect))])
  int <- cbind(y_self_hat[,match(intQTL, rownames(scan_effect))], y_nei_hat[,match(intQTL, rownames(scan_effect))])

  if(is.null(addcovar)==FALSE) {
    LOD_int <- c()
    for(k in 1:q) {
      LL_nei <- stats::logLik(stats::glm(pheno~addcovar+X+y_self_hat[,k]+y_self_hat[,k]:int[,1]+y_nei_hat[,k], family=glm_family))
      LL_int <- stats::logLik(stats::glm(pheno~addcovar+X+y_self_hat[,k]+y_self_hat[,k]:int[,1]+y_nei_hat[,k]+y_nei_hat[,k]:int[,2], family=glm_family))
      LOD_int <- c(LOD_int, log10(exp(LL_int-LL_nei)))
    }
  } else if(is.null(addcovar)==TRUE) {
    LOD_int <- c()
    for(k in 1:q) {
      LL_nei <- stats::logLik(stats::glm(pheno~X+y_self_hat[,k]+y_self_hat[,k]:int[,1]+y_nei_hat[,k], family=glm_family))
      LL_int <- stats::logLik(stats::glm(pheno~X+y_self_hat[,k]+y_self_hat[,k]:int[,1]+y_nei_hat[,k]+y_nei_hat[,k]:int[,2], family=glm_family))
      LOD_int <- c(LOD_int, log10(exp(LL_int-LL_nei)))
    }
  }

  marker_info <- get_markers(genoprobs=genoprobs)
  LODlist <- data.frame(marker_info, LOD_int)
  colnames(LODlist) <- c("chr","pos","LOD_int")

  return(LODlist)
}
