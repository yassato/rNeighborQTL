#' Permutation tests for neighbor effects with a QTL model
#'
#' A function to calculate a genome-wide LOD threshold using permutation tests for self or neighbor effects.
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
#' @param times No. of permutation iterations. Default at 99 times
#' @param p_val Upper quantile for permuted LOD scores
#' @param type Select \code{"self"} or \code{"neighbor"} to perform permutation tests for self or neighbor effects. Default is \code{"neighbor"}.
#' @param n_core No. of cores for a parallel computation. This does not work for Windows OS. Default is a single-core computation.
#' @return  A vector of LOD thresholds for all markers, with the chromosome numbers and positions. The row names correspond to marker names.
#' \itemize{
#'  \item{\code{chr}} {Chromosome number}
#'  \item{\code{pos}} {Marker position}
#'  \item{\code{LOD_th}} {A LOD threshold for each marker at \code{p-val}.}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @seealso plot_nei
#' @import parallel
#' @examples
#' #demo using F2
#' set.seed(1)
#' data("fake.f2",package="qtl")
#' fake_f2 <- qtl2::convert2cross2(fake.f2)
#' fake_f2 <- subset(fake_f2,chr=c(1:19))
#' smap_f2 <- cbind(runif(qtl2::n_ind(fake_f2),1,100),runif(qtl2::n_ind(fake_f2),1,100))
#' gmap_f2 <- qtl2::insert_pseudomarkers(fake_f2$gmap, step=1)
#' genoprobs_f2 <- qtl2::calc_genoprob(fake_f2,gmap_f2,error_prob=0.002)
#'
#' scan_f2_1 <- scan_neighbor(genoprobs=genoprobs_f2,
#'                            pheno=fake_f2$pheno[,1],
#'                            gmap=gmap_f2, smap=smap_f2,
#'                            scale=20, addcovar=as.matrix(fake_f2$covar)
#'                            )
#' perm_f2_1 <- perm_neighbor(genoprobs=genoprobs_f2,
#'                            pheno=fake_f2$pheno[,1],
#'                            gmap=gmap_f2, smap=smap_f2,
#'                            scale=20, addcovar=as.matrix(fake_f2$covar),
#'                            times=9, p_val=0.01
#'                            )
#' plot_nei(scan_f2_1,th=perm_f2_1$LOD_th)
#' @export
perm_neighbor = function(genoprobs, pheno, gmap, contrasts=c(TRUE,TRUE,TRUE), smap, scale, addcovar=NULL, addQTL=NULL, grouping=rep(1,nrow(smap)), response="quantitative", type="neighbor", times=99, p_val=0.05, n_core=1L) {

  if(type=="self") {
    func = function(x) return(scan_neighbor(genoprobs=genoprobs, pheno=sample(pheno), gmap=gmap, smap=smap, scale=scale, contrasts=contrasts, addcovar=addcovar, addQTL=addQTL, grouping=grouping, response=response)$LOD_self)
  } else if(type=="neighbor") {
    func = function(x) return(scan_neighbor(genoprobs=genoprobs, pheno=sample(pheno), gmap=gmap, smap=smap, scale=scale, contrasts=contrasts, addcovar=addcovar, addQTL=addQTL, grouping=grouping, response=response)$LOD_nei)
  } else { print("error: type must be 'self' or 'neighbor'") }
  res = parallel::mcmapply(func, 1:times, mc.cores=getOption("mc.cores",n_core))

  th = c()
  for(i in 1:nrow(res)) { th = c(th, stats::quantile(res[i,], 1-p_val, na.rm=TRUE)) }

  marker_info = get_markers(gmap)

  permlist = data.frame(marker_info, th)
  colnames(permlist) = c("chr","pos","LOD_th")

  return(permlist)
}
