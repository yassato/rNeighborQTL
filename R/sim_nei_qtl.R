#' Phenotype simulation for neighbor QTL effects
#'
#' A function to simulate neighbor effects with given QTL effects, distance scale, and causal markers.
#' @param genoprobs Conditional genotype probabilities as obtained from \code{qtl2::calc_genoprob()}.
#' @param gmap Genetic map including observed and pseudomarkers, as obtained from \code{qtl2::insert_pseudomarkers()}.
#' @param contrasts A vector composed of three TRUE/FALSE values, which represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB.
#' @param smap A matrix showing a spatial map. The first and second column include spatial points along a x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param a2 A numeric scalar indicating additive deviation.
#' @param d2 A numeric scalar indicating dominance deviation.
#' @param grouping An integer vector assigning each individual to a group. This argument can be useful when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param n_QTL A positive integer indicating the number of causal markers.
#' @return A numeric matrix containing individuals x marker elements for neighbor QTL effects.
#' \itemize{
#'  \item{\code{true_scale}} {True distance scale of simulated neighbor effects}
#'  \item{\code{true_marker}} {The name(s) of causal markers}
#'  \item{\code{nei_y}} {Simulated neighbor effects standardized to have zero mean and one variance}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' Major genetic effects, \code{a2} and \code{d2}, are allocated to causal loci randomly selected by \code{n_QTL}, while minor polygenic effects (i.e., 1% of \code{a2}) are allocated to the other loci.
#' @examples
#' set.seed(1234)
#' data("fake.f2",package="qtl")
#' fake_f2 <- qtl2::convert2cross2(fake.f2)
#' fake_f2 <- subset(fake_f2,chr=c(1:19))
#' smap_f2 <- cbind(runif(qtl2::n_ind(fake_f2),1,100),runif(qtl2::n_ind(fake_f2),1,100))
#' gmap_f2 <- qtl2::insert_pseudomarkers(fake_f2$gmap, step=2)
#' genoprobs_f2 <- qtl2::calc_genoprob(fake_f2,gmap_f2)
#' s_seq <- quantile(dist(smap_f2),c(0.1*(1:10)))
#'
#' nei_eff <- sim_nei_qtl(genoprobs_f2, gmap_f2, a2=0.5, d2=0.5,
#'                        contrasts=c(TRUE,TRUE,TRUE), smap=smap_f2,
#'                        scale=19, n_QTL=1)
#'
#' scan_f2 <- scan_neighbor(genoprobs=genoprobs_f2,
#'                          pheno=nei_eff$nei_y,
#'                          gmap=gmap_f2, contrasts=c(TRUE,TRUE,TRUE),
#'                          smap=smap_f2, scale=19,
#'                          addcovar=as.matrix(fake_f2$covar))
#' plot_nei(scan_f2)
#' @export
sim_nei_qtl = function(genoprobs, gmap, a2, d2, contrasts=c(TRUE,TRUE,TRUE), smap, scale, grouping=rep(1,nrow(smap)), n_QTL=1) {
  major_eff <- calc_neiprob(genoprobs, gmap, a2=a2, d2=d2, contrasts=contrasts, smap=smap, scale=scale, grouping=grouping)
  minor_eff <- calc_neiprob(genoprobs, gmap, a2=sign(a2)*a2*0.01, d2=sign(d2)*a2*0.01, contrasts=contrasts, smap=smap, scale=scale, grouping=grouping)
  major_loci <- sample(ncol(minor_eff),n_QTL)

  minor_eff[,major_loci] <- major_eff[,major_loci]
  nei_eff <- scale(apply(minor_eff,1,sum))

  true_loci <- colnames(major_eff)[major_loci]

  res <- list(true_scale=scale, true_marker=true_loci, nei_y=nei_eff)
  return(res)
}

