#' Calculating proportion of phenotypic variation explained (PVE) by neighbor effects
#'
#' A function to calculate PVE by neighbor effects for a series of neighbor distance (\code{s_seq}) using a linear mixed model.
#' @param genoprobs Conditional genotype probabilities as obtained from \code{qtl2::calc_genoprob()}.
#' @param pheno A vector of individual phenotypes.
#' @param gmap Genetic map including observed and pseudomarkers, as obtained from \code{qtl2::insert_pseudomarkers()}.
#' @param contrasts  A vector composed of three TRUE/FALSE values. Depending on crossing design, it represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB.
#' @param smap A matrix showing a spatial map. The first and second column include spatial points along a x-axis and y-axis, respectively.
#' @param s_seq A numeric vector including a set of the maximum spatial distance between a focal individual and neighbors to define neighbor effects. A scalar is also allowed.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param grouping An optional integer vector assigning each individual to a group. This argument can be useful when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An optional argument to select trait types. The \code{"quantitative"} or \code{"binary"} applies the \code{"lmm.aireml()"} or \code{"logistic.mm.aireml()"} for a mixed model, respectively.
#' @param fig TRUE/FALSE to add a figure of PVE or not.
#' @return A matrix containing the maximum neighbor distance, PVE by neighbor effects, and p-value by LRT.
#' \itemize{
#'  \item{\code{scale}} {Maximum neighbor distance given as an argument}
#'  \item{\code{PVE_nei}} {Proportion of phenotypic variation explained (PVE) by neighbor effects}
#'  \item{\code{p-value}} {p-value by a likelihood ratio test between models with or without neighbor effects}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' This function uses mixed models via the \code{gaston} package (Perdry & Dandine-Roulland 2020).
#' If \code{"binary"} is selected, \code{logistic.mm.aireml()} is called via the \code{gaston} package (see Chen et al. 2016 for the theory).
#' In such a case, \code{PVEnei} below is given by the variance component parameter \eqn{\sigma} (i.e., not a proportional value) and p-values are not provided.
#' @references
#' * Perdry H, Dandine-Roulland C (2019) gaston: Genetic Data Handling (QC, GRM, LD, PCA) & Linear Mixed Models. R package version 1.5.5. https://CRAN.R-project.org/package=gaston
#' * Chen H, Wang C, Conomos M. et al. (2016) Control for population structure and relatedness for binary traits in genetic association studies via logistic mixed models. The American Journal of Human Genetics 98: 653-666.
#' @import gaston Matrix
#' @examples
#' set.seed(1234)
#' data("fake.f2",package="qtl")
#' fake_f2 <- qtl2::convert2cross2(fake.f2)
#' fake_f2 <- subset(fake_f2,chr=c(1:19))
#' smap_f2 <- cbind(runif(qtl2::n_ind(fake_f2),1,100),runif(qtl2::n_ind(fake_f2),1,100))
#' gmap_f2 <- qtl2::insert_pseudomarkers(fake_f2$gmap, step=2)
#' genoprobs_f2 <- qtl2::calc_genoprob(fake_f2,gmap_f2)
#' s_seq <- quantile(dist(smap_f2),c(0.1*(1:10)))
#' pve_f2 <- calc_pve(genoprobs=genoprobs_f2,
#'                    pheno=fake_f2$pheno[,1], gmap=gmap_f2,
#'                    smap=smap_f2, s_seq=s_seq,
#'                    contrasts=c(TRUE,TRUE,TRUE),
#'                    addcovar=as.matrix(fake_f2$covar)
#'                    )
#' @export
calc_pve = function(genoprobs, pheno, gmap, contrasts=c(TRUE,TRUE,TRUE), smap, s_seq, addcovar=NULL, grouping=rep(1,nrow(smap)), response="quantitative", fig=TRUE) {

  res <- c()
  for(s in s_seq) {
    if(class(s)=="numeric") { message("scale = ", round(s,3), "\n") }

    selfprobs <- genoprobs2selfprobs(genoprobs=genoprobs, gmap=gmap, a1=1, d1=0, contrasts=contrasts)
    neiprobs <- calc_neiprob(genoprobs=genoprobs, gmap=gmap, contrasts=contrasts, smap=smap, scale=s, a2=1, d2=0, grouping=grouping)

    selfprobs <- (selfprobs-mean(selfprobs))/stats::sd(selfprobs)
    neiprobs <- (neiprobs-mean(neiprobs))/stats::sd(neiprobs)

    K_self <- tcrossprod(selfprobs)/(ncol(selfprobs)-1)
    K_nei <- tcrossprod(neiprobs)/(ncol(neiprobs)-1)
    K_self <- as.matrix(Matrix::nearPD(K_self, maxit=10^6)$mat)
    K_nei <- as.matrix(Matrix::nearPD(K_nei, maxit=10^6)$mat)

    if(is.null(addcovar)) {
      X <- matrix(1, nrow=length(pheno))
    } else {
      X <- stats::model.matrix(~addcovar)
    }

    if(response=="quantitative") {
      aireml2 <- gaston::lmm.aireml(Y=pheno, X=X, K=list(K_self,K_nei), verbose=FALSE)
      aireml1 <- gaston::lmm.aireml(Y=pheno, X=X, K=list(K_self), verbose=FALSE)
      pve <- aireml2$tau[2]/sum(aireml2$tau, aireml2$sigma)
      p_val <- stats::pchisq(-2*(aireml1$logL-aireml2$logL), df=1, lower.tail=FALSE)
    } else if(response=="binary"){
      aireml <- gaston::logistic.mm.aireml(Y=pheno, X=X, K=list(K_self, K_nei), verbose=FALSE)
      pve <- aireml$tau[2]
      p_val <- NA
    } else {
      stop("error: reponse must be 'quantitative' or 'binary'")
    }
    res <- rbind(res, c(s, pve, p_val))
  }
  colnames(res) <- c("scale","PVE_nei", "p-value")

  if(fig==TRUE) {
    PVE <- res[,2]
    deltaPVE <- PVE - c(0, PVE[1:(length(PVE)-1)])
    graphics::plot(s_seq, deltaPVE, type="l", xlab="spatial scale")
    graphics::points(s_seq, deltaPVE)
    max(deltaPVE)
    s_seq[deltaPVE==max(deltaPVE)]
    graphics::points(s_seq[deltaPVE==max(deltaPVE)],deltaPVE[deltaPVE==max(deltaPVE)],pch=16)
  }
  return(res)
}
