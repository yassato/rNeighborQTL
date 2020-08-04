#' Calculating log-likelihood in generalized linear models
#' 
#' An utility function to extract log-likelihood based on AIC of glm.fit()
#' @param ... Arguments to be passed to \code{glm.fit()}.
#' @return Log-likelihood
logLik_glm.fit <- function(...) {
    fit <- stats::glm.fit(...)
    p <- fit$rank
    if (fit$family$family %in% c("gaussian", "Gamma", "inverse.gaussian"))
        p <- p + 1
    p - fit$aic / 2
}
