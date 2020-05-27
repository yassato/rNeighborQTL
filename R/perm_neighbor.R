#' Permutation tests for neighbor effects with a QTL model
#'
#' A function to calculate a genome-wide LOD threshold using permutation tests for self or neighbor effects.
#' @param genoprobs Conditional genotype probabilities as taken from \code{qtl::calc.genoprob()}.
#' @param pheno A vector of individual phenotypes.
#' @param contrasts A vector composed of three TRUE/FALSE values. Depending on the crossing design, it represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial positions along an x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param addQTL An optional vector containing marker names that are considered covariates. Namely, this option allows composite interval mapping (Jansen 1993).
#' @param intQTL An option when using \code{int_neighbor()}. A name of a focal marker to be tested for its epistasis with the other markers in neighbor effects. The marker name must be included by \code{addQTL}.
#' @param grouping An optional integer vector assigning each individual to a group. This argument can be used when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An optional argument to select trait types. The \code{"quantitative"} or \code{"binary"} calls the \code{"gaussian"} or \code{"binomial"} family in \code{glm()}, respectively.
#' @param times No. of permutation iterations. Default at 99 times
#' @param p_val A vector indicating upper quantiles for permutation LOD scores
#' @param type Select \code{"self"}, \code{"neighbor"}, or \code{"int"} to perform permutation tests for self effects, neighbor effects, or neighbor epistasis, respectively.
#' @param n_core No. of cores for a parallel computation. This does not work for Windows OS. Default is a single-core computation.
#' @return LOD thresholds at given quantiles by \code{p-val}
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @seealso \code{\link{plot_nei}} \code{\link{scan_neighbor}} \code{\link{int_neighbor}}
#' @import parallel
#' @examples
#' set.seed(1234)
#' test_map <- qtl::sim.map(len=rep(20,5),n.mar=3,include.x=FALSE)
#' test_cross <- qtl::sim.cross(test_map,n.ind=50)
#' test_smap <- cbind(runif(50,1,100),runif(50,1,100))
#' test_genoprobs <- qtl::calc.genoprob(test_cross,step=2)
#'
#' test_perm <- perm_neighbor(genoprobs=test_genoprobs,
#'                            pheno=test_cross$pheno$phenotype,
#'                            smap=test_smap,scale=20,
#'                            times=3, p_val=c(1.0,0.5)
#'                            )
#' @export
perm_neighbor = function(genoprobs, pheno, contrasts=c(TRUE,TRUE,TRUE), smap, scale, addcovar=NULL, addQTL=NULL, intQTL=NULL, grouping=rep(1,nrow(smap)), response="quantitative", type="neighbor", times=99, p_val=0.05, n_core=1L) {

  if(type=="self") {
    func = function(x) return(max(scan_neighbor(genoprobs=genoprobs, pheno=sample(pheno), smap=smap, scale=scale, contrasts=contrasts, addcovar=addcovar, addQTL=addQTL, grouping=grouping, response=response)$LOD_self))
  } else if(type=="neighbor") {
    func = function(x) return(max(scan_neighbor(genoprobs=genoprobs, pheno=sample(pheno), smap=smap, scale=scale, contrasts=contrasts, addcovar=addcovar, addQTL=addQTL, grouping=grouping, response=response)$LOD_nei))
  } else if(type=="int") {
    func = function(x) return(max(int_neighbor(genoprobs=genoprobs, pheno=sample(pheno), smap=smap, scale=scale, contrasts=contrasts, addcovar=addcovar, addQTL=addQTL, intQTL=intQTL, grouping=grouping, response=response)$LOD_int))
  } else {
    stop("error: type must be 'self', 'neighbor', or 'int'")
  }

  res = parallel::mcmapply(func, 1:times, mc.cores=getOption("mc.cores",n_core))
  th = stats::quantile(res, 1-p_val, na.rm=TRUE)

  return(th)
}
