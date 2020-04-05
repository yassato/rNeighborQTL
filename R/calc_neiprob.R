#' Calculating a set of neighbor QTL effects from conditional genotype probabilities
#'
#' A function to calculate self QTL effects for all individuals, with given deviation coefficients and conditional genotype probabilities.
#' @param genoprobs Conditional genotype probabilities as obtained from \code{qtl2::calc_genoprob()}.
#' @param gmap Genetic map including observed and pseudomarkers, as obtained from \code{qtl2::insert_pseudomarkers()}.
#' @param contrasts A vector composed of three TRUE/FALSE values, which represents the presence/absence of specific genotypes as c(T/F, T/F, T/F) = AA, AB, BB.
#' @param smap A matrix showing a spatial map. The first and second column include spatial points along a x-axis and y-axis, respectively.
#' @param scale A numeric scalar indicating the maximum spatial distance between a focal individual and neighbors to define neighbor effects.
#' @param a2 A numeric scalar indicating additive deviation.
#' @param d2 A numeric scalar indicating dominance deviation.
#' @param grouping An integer vector assigning each individual to a group. This argument can be useful when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @return A numeric matrix containing individuals x marker elements for neighbor QTL effects.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @examples
#' set.seed(1)
#' data("fake.f2",package="qtl")
#' fake_f2 <- qtl2::convert2cross2(fake.f2)
#' fake_f2 <- subset(fake_f2,chr=c(1:19))
#' smap_f2 <- cbind(runif(qtl2::n_ind(fake_f2),1,100),runif(qtl2::n_ind(fake_f2),1,100))
#' gmap_f2 <- qtl2::insert_pseudomarkers(fake_f2$gmap, step=1)
#' genoprobs_f2 <- qtl2::calc_genoprob(fake_f2,gmap_f2,error_prob=0.002)
#'
#' neiprobs <- calc_neiprob(genoprobs=genoprobs_f2,
#'                          gmap=gmap_f2, a2=1, d2=0,
#'                          smap=smap_f2, scale=20,
#'                          grouping=fake_f2$covar$pgm
#'                          )
#' @export
calc_neiprob = function(genoprobs, gmap, a2, d2, contrasts=c(TRUE,TRUE,TRUE), smap, scale, grouping=rep(1,nrow(smap))) {
  p = dim(genoprobs[[1]])[1]
  geno = decompose_genoprobs(genoprobs=genoprobs, contrasts=contrasts)

  neiprob_i = function(i) {
    id = c(1:p)[grouping == grouping[i]]

    d_i = mapply(function(x) { return(sqrt((smap[x,1]-smap[i,1])^2 + (smap[x,2]-smap[i,2])^2)) },id)
    prob_i = 0
    j_id = id[(d_i>0)&(d_i<=scale)]
    if(length(j_id)==0) {
      return(rep(0,ncol(geno$AA)))
    } else {
      for(j in j_id){
        prob_ij = neiprob(i=i, j=j, a2=a2, d2=d2, AA=geno$AA, AB=geno$AB, BB=geno$BB)
        prob_i = prob_i + prob_ij
      }
      prob_i = prob_i/length(j_id)
      return(prob_i)
    }
  }
  neiList = mapply(neiprob_i, 1:p)
  neiList = t(neiList)

  marker_info = get_markers(gmap)
  colnames(neiList) = rownames(marker_info)
  rownames(neiList) = rownames(genoprobs[[1]])
  return(neiList)
}
