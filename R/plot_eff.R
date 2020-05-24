#' Plot self and neighbor QTL effects across a genome
#'
#' Plot estimated additive and dominance deviation for self or neighbor effects across a genome
#' @param res Output results of \code{eff_neighbor()}.
#' @param type An option to select \code{"self"} or \code{"neighbor"} effects to be shown. Default is \code{"neighbor"}.
#' @seealso eff_neighbor
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
plot_eff = function(res, type="neighbor") {
  x <- c(1:nrow(res))

  if(type=="neighbor") {
    a <- res$a2
    d <- res$d2
  } else if(type=="self") {
    a <- res$a1
    d <- res$d1
  } else {
    a <- NULL
    d <- NULL
    warning("error: type must be 'self' or 'neighbor'")
    return(NULL)
  }

  graphics::plot(x, (abs(a)+abs(d)), type="n", xlab="", las=1, ylim=c(range(c(stats::na.omit(a),stats::na.omit(d)))),xaxt="n", yaxt="s", main=paste(type,", a:solid; ","d:dashed"), ylab="QTL effect")
  for(i in levels(factor(res$chr))) {
    graphics::points(x[res$chr==i], a[res$chr==i], type="l", col=i, lty=1)
    graphics::points(x[res$chr==i], d[res$chr==i], type="l", col=i, lty=2)
  }
  unobs = grep("_loc", rownames(res))
  graphics::points(x[-unobs], a[-unobs], pch=16, cex=0.75)
  graphics::points(x[-unobs], d[-unobs], pch=1, cex=0.75)
}
