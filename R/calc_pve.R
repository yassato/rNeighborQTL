#' Calculating phenotypic variation explained by neighbor effects
#'
#' A function to calculate the proportion or ratio of phenotypic variation explained (PVE or RVE) by neighbor effects for a series of neighbor distance (\code{s_seq}) using mixed models.
#' @param genoprobs Conditional genotype probabilities as taken from \code{qtl::calc.genoprob()}.
#' @param pheno A vector of individual phenotypes.
#' @param contrasts A vector composed of three TRUE/FALSE values. Depending on the crossing design, it represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial positions along an x-axis and y-axis, respectively.
#' @param s_seq A numeric vector including a set of the maximum spatial distance between a focal individual and neighbors to define neighbor effects. A scalar is also allowed.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param grouping An optional integer vector assigning each individual to a group. This argument can be used when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An optional argument to select trait types. The \code{"quantitative"} or \code{"binary"} applies the \code{"lmm.aireml()"} or \code{"logistic.mm.aireml()"} for a mixed model, respectively.
#' @param fig TRUE/FALSE to add a figure of Delta PVE or not.
#' @return A matrix containing the maximum neighbor distance, phenotypic variation explained by neighbor effects, and p-value by a likelihood ratio test.
#' \itemize{
#'  \item{\code{scale}} {Maximum neighbor distance given as an argument}
#'  \item{\code{Var_nei}} {Proportion or ratio of phenotypic variation explained (PVE or RVE) by neighbor effects for linear or logistic mixed models, respectively}
#'  \item{\code{p-value}} {p-value by a likelihood ratio test between models with or without neighbor effects}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' This function calls linear or logistic mixed models via the \code{gaston} package (Perdry & Dandine-Roulland 2020).
#' If \code{"binary"} is selected, \code{Var_nei} in the output is given by the proportion of phenotypic variation explained (PVE) by neighbor effects as PVEnei =\eqn{\sigma^2_2/(\sigma^2_1+\sigma^2_2+\sigma^2_e)}.
#' If \code{"binary"} is selected, \code{Var_nei} is given by the ratio of phenotypic variation explained (RVE) by neighbor effects as RVEnei =\eqn{\sigma^2_2/\sigma^2_1} and p-values are not available.
#' This is because a logistic mixed model \code{logistic.mm.aireml()} called via the \code{gaston} package does not provide \eqn{\sigma^2_e} and log-likelihood (see Chen et al. 2016 for the theory).
#' @references
#' * Perdry H, Dandine-Roulland C (2019) gaston: Genetic Data Handling (QC, GRM, LD, PCA) & Linear Mixed Models. R package version 1.5.5. https://CRAN.R-project.org/package=gaston
#' * Chen H, Wang C, Conomos M. et al. (2016) Control for population structure and relatedness for binary traits in genetic association studies via logistic mixed models. The American Journal of Human Genetics 98: 653-666.
#' @import gaston Matrix
#' @examples
#' set.seed(1234)
#' test_map <- qtl::sim.map(len=rep(20,5),n.mar=3,include.x=FALSE)
#' test_cross <- qtl::sim.cross(test_map,n.ind=50)
#' test_smap <- cbind(runif(50,1,100),runif(50,1,100))
#' test_genoprobs <- qtl::calc.genoprob(test_cross,step=2)
#' s_seq <- quantile(dist(test_smap),c(0.1*(1:10)))
#'
#' test_pve <- calc_pve(genoprobs=test_genoprobs,
#'                      pheno=test_cross$pheno$phenotype,
#'                      smap=test_smap, s_seq=s_seq,
#'                      contrasts=c(TRUE,TRUE,TRUE),
#'                      )
#' @export
calc_pve = function(genoprobs, pheno, contrasts=c(TRUE,TRUE,TRUE), smap, s_seq, addcovar=NULL, grouping=rep(1,nrow(smap)), response="quantitative", fig=TRUE) {

  res <- c()
  for(s in s_seq) {
    if(class(s)=="numeric") { message("scale = ", round(s,3)) }

    selfprobs <- genoprobs2selfprobs(genoprobs=genoprobs, a1=1, d1=0, contrasts=contrasts)

    if((contrasts[2]==TRUE)&(contrasts[3]==FALSE)) {
      neiprobs <- calc_neiprob(genoprobs=genoprobs, contrasts=contrasts, smap=smap, scale=s, a2=1, d2=-1, grouping=grouping)
    } else if(contrasts[2]==TRUE) {
      neiprobs <- calc_neiprob(genoprobs=genoprobs, contrasts=contrasts, smap=smap, scale=s, a2=1, d2=0.5, grouping=grouping, d2sq0=TRUE)
    } else {
      neiprobs <- calc_neiprob(genoprobs=genoprobs, contrasts=contrasts, smap=smap, scale=s, a2=1, d2=0, grouping=grouping)
    }

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
      pve <- aireml$tau[2]/aireml$tau[1]
      p_val <- NA
    } else {
      stop("error: reponse must be 'quantitative' or 'binary'")
    }
    res <- rbind(res, c(s, pve, p_val))
  }
  colnames(res) <- c("scale", "Var_nei", "p-value")

  if(fig==TRUE) {
    PVE <- res[,2]
    deltaPVE <- PVE - c(0, PVE[1:(length(PVE)-1)])

    switch(response,
           "quantitative" = ylab <- "deltaPVE",
           "binary" = ylab <- "deltaRVE"
    )

    graphics::plot(s_seq, deltaPVE, type="l", xlab="spatial scale", ylab=ylab)
    graphics::points(s_seq, deltaPVE)
    graphics::points(s_seq[deltaPVE==max(deltaPVE)[1]],deltaPVE[deltaPVE==max(deltaPVE)[1]],pch=16)
  }
  return(res)
}
