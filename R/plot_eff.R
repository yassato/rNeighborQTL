#' Plot self and neighbor QTL effects across a genome
#'
#' Plot estimated additive and dominance deviation for self or neighbor effects across a genome
#' @param res Output results of \code{eff_neighbor()}.
#' @param type An option to select \code{"self"} or \code{"neighbor"} effects to be shown. Default is \code{"neighbor"}.
#' @seealso eff_neighbor
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
plot_eff = function(res, type=c("neighbor","self")) {
  type <- match.arg(type)

  pos <- res$pos
  chr <- as.factor(res$chr)
  coord <- 0
  M <- 0
  tic <- numeric(nlevels(chr))
  for (i in 1:nlevels(chr)) {
    w <- (chr == levels(chr)[i])
    pos.c <- pos[w]
    coord[w] <- M + pos.c
    mx <- max(pos.c)
    tic[i] <- M + mx/2
    M <- M + mx + max(pos)*0.2
  }
  x <- coord/M
  tic <- tic/M
  
  if(type=="neighbor") {
    a <- res$a2
    d <- res$d2
  } else { #if(type=="self") {
    a <- res$a1
    d <- res$d1
  }

  graphics::plot(x, (abs(a)+abs(d)), type="n", xlab="", las=1, ylim=c(range(c(stats::na.omit(a),stats::na.omit(d)))),xaxt="n", yaxt="s", main=paste(type,", a:solid; ","d:dashed"), ylab="QTL effect")
  for(i in levels(factor(res$chr))) {
    graphics::points(x[res$chr==i], a[res$chr==i], type="l", col=i, lty=1)
    graphics::points(x[res$chr==i], d[res$chr==i], type="l", col=i, lty=2)
  }
  unobs = grep("_loc", rownames(res))
  graphics::points(x[-unobs], a[-unobs], pch=16, cex=0.75)
  graphics::points(x[-unobs], d[-unobs], pch=1, cex=0.75)
  graphics::axis(side=1, at=tic, labels=levels(chr))
}
