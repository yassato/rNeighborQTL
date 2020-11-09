#' Calculating phenotypic variation explained by neighbor effects
#'
#' A function to calculate the proportion or ratio of phenotypic variation explained (PVE or RVE) by neighbor effects for a series of neighbor distance (\code{s_seq}) using mixed models.
#' @param genoprobs Conditional genotype probabilities as taken from \code{qtl::calc.genoprob()}.
#' @param pheno A vector of individual phenotypes.
#' @param smap A matrix showing a spatial map for individuals. The first and second column include spatial positions along an x-axis and y-axis, respectively.
#' @param s_seq A numeric vector including a set of the maximum spatial distance between a focal individual and neighbors to define neighbor effects. A scalar is also allowed.
#' @param addcovar An optional matrix including additional non-genetic covariates. It contains no. of individuals x no. of covariates.
#' @param grouping An optional integer vector assigning each individual to a group. This argument can be used when \code{smap} contains different experimental replicates. Default setting means that all individuals are belong to a single group.
#' @param response An optional argument to select trait types. The \code{"quantitative"} or \code{"binary"} applies the \code{"lmm.aireml()"} or \code{"logistic.mm.aireml()"} for a mixed model, respectively.
#' @param fig TRUE/FALSE to add a figure of Delta PVE or not.
#' @param contrasts An optional vector composed of three TRUE/FALSE values, which represents the presence/absence of specific genotypes as c(TRUE/FALSE, TRUE/FALSE, TRUE/FALSE) = AA, AB, BB. If \code{NULL}, it is compiled from \code{genoprobs} automatically.
#' @return A matrix containing the maximum neighbor distance, phenotypic variation explained by neighbor effects, and p-value by a likelihood ratio test.
#' \itemize{
#'  \item{\code{scale}} {Maximum neighbor distance given as an argument}
#'  \item{\code{Var_self}} {Proportion or ratio of phenotypic variation explained (PVE or RVE) by self-genotype effects for linear or logistic mixed models, respectively}
#'  \item{\code{Var_nei}} {Proportion or ratio of phenotypic variation explained (PVE or RVE) by neighbor effects for linear or logistic mixed models, respectively}
#'  \item{\code{p-value}} {p-value by a likelihood ratio test between models with or without neighbor effects. Self effects are tested when the scale is zero}
#' }
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' This function calls linear or logistic mixed models via the \code{gaston} package (Perdry & Dandine-Roulland 2020).
#' If \code{"quantitative"} is selected, \code{Var_self} or \code{Var_nei} in the output is given by the proportion of phenotypic variation explained (PVE) by neighbor effects as PVEnei =\eqn{\sigma^2_2/(\sigma^2_1+\sigma^2_2+\sigma^2_e)}.
#' If \code{"binary"} is selected, \code{Var_self} or \code{Var_nei} is given by the ratio of phenotypic variation explained (RVE) by neighbor effects as RVEnei =\eqn{\sigma^2_2/\sigma^2_1} and p-values are not available.
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
#'                      )
#' @export
calc_pve = function(genoprobs, pheno, smap, s_seq, addcovar=NULL, grouping=rep(1,nrow(smap)), response=c("quantitative","binary"), fig=TRUE, contrasts=NULL) {
  
  scaling = function(vec) {
    return((vec-mean(vec))/stats::sd(vec))
  }
  
  response <- match.arg(response)
  
  selfprobs <- genoprobs2selfprobs(genoprobs=genoprobs, a1=1, d1=0, contrasts=contrasts)
  p <- nrow(selfprobs)
  contrasts <- attr(selfprobs, "contrasts")
  selfprobs <- mapply(function(x) { scaling(selfprobs[x,]) }, 1:p)
  selfprobs <- t(selfprobs)
  
  K_self <- tcrossprod(selfprobs)/(ncol(selfprobs)-1)
  K_self <- as.matrix(Matrix::nearPD(K_self, maxit=10^6)$mat)
  
  if(is.null(addcovar)) {
    X <- matrix(1, nrow=length(pheno))
  } else {
    X <- stats::model.matrix(~addcovar)
  }
  
  res <- c()
  if(response=="quantitative") {
    aireml1 <- gaston::lmm.aireml(Y=pheno, X=X, K=list(K_self), verbose=FALSE)
    pve <- aireml1$tau/sum(aireml1$tau, aireml1$sigma)
    p_val <- stats::pchisq(-2*(aireml1$logL0-aireml1$logL), df=1, lower.tail=FALSE)
  } else {
    aireml <- gaston::logistic.mm.aireml(Y=pheno, X=X, K=list(K_self), verbose=FALSE)
    pve <- aireml$tau
    p_val <- NA
  }
  res <- c(0, pve, 0, p_val)
  
  for(s in s_seq) {
    if(class(s)=="numeric") { message("scale = ", round(s,3)) }
    
    if((contrasts[2]==TRUE)&(contrasts[3]==FALSE)) {
      neiprobs <- calc_neiprob(genoprobs=genoprobs, contrasts=contrasts, smap=smap, scale=s, a2=1, d2=-1, grouping=grouping)
    } else if(contrasts[2]==TRUE) {
      neiprobs <- calc_neiprob(genoprobs=genoprobs, contrasts=contrasts, smap=smap, scale=s, a2=1, d2=0.5, grouping=grouping, d2sq0=TRUE)
    } else { #if(response=="binary"){
      neiprobs <- calc_neiprob(genoprobs=genoprobs, contrasts=contrasts, smap=smap, scale=s, a2=1, d2=0, grouping=grouping)
    }
    
    neiprobs <- mapply(function(x) { scaling(neiprobs[x,]) }, 1:p)
    neiprobs <- t(neiprobs)
    
    K_nei <- tcrossprod(neiprobs)/(ncol(neiprobs)-1)
    K_nei <- as.matrix(Matrix::nearPD(K_nei, maxit=10^6)$mat)
    
    if(response=="quantitative") {
      aireml2 <- gaston::lmm.aireml(Y=pheno, X=X, K=list(K_self,K_nei), verbose=FALSE)
      pve_s <- aireml2$tau[1]/sum(aireml2$tau, aireml2$sigma)
      pve_n <- aireml2$tau[2]/sum(aireml2$tau, aireml2$sigma)
      p_val <- stats::pchisq(-2*(aireml1$logL-aireml2$logL), df=1, lower.tail=FALSE)
    } else { #if(response=="binary"){
      aireml <- gaston::logistic.mm.aireml(Y=pheno, X=X, K=list(K_self, K_nei), verbose=FALSE)
      pve_s <- aireml$tau[1]/aireml$tau[2]
      pve_n <- aireml$tau[2]/aireml$tau[1]
      p_val <- NA
    }
    res <- rbind(res, c(s, pve_s, pve_n, p_val))
  }
  colnames(res) <- c("scale", "Var_self", "Var_nei", "p-value")
  rownames(res) <- NULL
  
  if(fig==TRUE) {
    res.sorted <- res[order(res[,1]),]
    PVE <- res.sorted[,3]
    deltaPVE <- diff(PVE)
    
    switch(response,
           "quantitative" = ylab <- "deltaPVE",
           "binary" = ylab <- "deltaRVE"
    )
    
    graphics::plot(res.sorted[-1,1], deltaPVE, type="o",
                   xlab="spatial scale", ylab=ylab,
                   pch=ifelse(deltaPVE==max(deltaPVE)[1],16,1))
  }
  return(res)
}
