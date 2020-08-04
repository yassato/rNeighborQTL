#' Plot LOD score for self or neighbor QTL effects
#'
#' Plot LOD curves for a genome scan of self and neighbor QTL effects.
#' @param res Output results of \code{scan_neighbor()}.
#' @param type Plot \code{"self"}, \code{"neighbor"} or \code{"int"} effects. Default is \code{"neighbor"} effects.
#' @param chr An optional vector to select chromosome numbers to be plotted. If \code{NULL}, shown are all chromosomes.
#' @param th Add genome-wide threshold by user-defined vectors or Bonferroni correction. Default is no thresholds added.
#' @param ... Arguments to be passed to \code{plot()}.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @details
#' For the \code{type} argument, \code{"int"} can be selected to draw the results of \code{int_neighbor()}.
#' In this case, the \code{res} object and \code{type} must match, otherwise it returns an error message.
#' @seealso \code{\link{scan_neighbor}} \code{\link{int_neighbor}} \code{\link{perm_neighbor}}
#' @export
plot_nei = function(res, type=c("neighbor","self","int"), chr=NULL, th=NULL, ...) {
  type <- match.arg(type)

  if(is.null(chr)==FALSE) {
    res = res[res$chr==chr,]
  }

  x <- c(1:nrow(res))

  switch(type,
         "self" = y <- res$LOD_self,
         "neighbor" = y <- res$LOD_nei,
         "int" = y <- res$LOD_int,
         stop("error: type must be 'self', 'neighbor', or 'int'")
  )

  if(is.null(y)) {
    stop("error: the output results and type does not match")
    }

  args <- list(...)
  args$type <- "n"; args$main <- type
  args$x <- x;  args$y <- y
  args$xlab <- "";  args$ylab <- "LOD score"
  args$xaxt <- "n"; args$yaxt <- "s"
  do.call(graphics::plot, args)
  for(i in levels(factor(res$chr))) {
    graphics::points(x[res$chr==i], y[res$chr==i], type="l", col=i)
  }
  unobs = grep("_loc", rownames(res))
  graphics::points(x[-unobs], y[-unobs], pch=16, cex=0.75)
  graphics::abline(h=th, col=grDevices::grey(0.5,0.5), lty=2)
}
