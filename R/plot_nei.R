#' Plot a genome scan for self or neighbor QTL effect
#'
#' Plot LOD curves for a genome scan of self and neighbor QTL effects
#' @param res Output results of \code{scan_neighbor()}.
#' @param type Plot \code{"self"} or \code{"neighbor"} effects. Default is \code{"neighbor"} effects.
#' @param chr An optional vector to select chromosome numbers to be plotted. If \code{NULL}, shown are all chromosomes.
#' @param th Add genome-wide threshold by user-defined vectors or Bonferroni correction. Default is no thresholds added.
#' @param ... Arguments to be passed to \code{plot()}.
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
#' @seealso scan_neighbor perm_neighbor
#' @export
plot_nei = function(res, type="neighbor", chr=NULL, th=NULL, ...) {

  if(is.null(chr)==FALSE) {
    res = res[res$chr==chr,]
  }

  x = cumsum(res$pos)

  if(type=="neighbor") {
    y = res$LOD_nei
  } else if(type=="self") {
    y = res$LOD_self
  } else {
    y = NULL
    print("error: type must be 'self' or 'neighbor'")
  }

  args = list(...)
  args$type = "n"; args$main = type
  args$x = x;  args$y = y
  args$xlab = "";  args$ylab = "LOD score"
  args$xaxt = "n"; args$yaxt = "s"
  do.call(graphics::plot, args)
  for(i in levels(factor(res$chr))) {
    graphics::points(x[res$chr==i], y[res$chr==i], type="l", col=i)
  }
  unobs = grep("c.",rownames(res))
  graphics::points(x[-unobs], y[-unobs], pch=16, cex=0.75)

  if(is.null(th)==FALSE) {
    if(th[1]=="bonf") {
      LOD_th = log10(exp(stats::qchisq(0.05/nrow(res), df=2, lower.tail=F)/2))
      graphics::abline(h=LOD_th, col=grDevices::grey(0.5,0.5), lty=2)
    } else {
      graphics::points(x, th, type="l", col=grDevices::grey(0.5,0.5), lty=2)
    }
  }
}
