#' Reshaping marker information
#'
#' A function to get marker information from a genetic map including observed and pseudo markers
#' @param genoprobs Conditional genotype probabilities as taken from \code{qtl::calc.genoprob()}.
#' @return A matrix showing the chromosome numbers (the first column) and positions (the second column) for all markers (row names).
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
get_markers = function(genoprobs) {
  n_chr <- length(genoprobs$geno)
  marker_info <- c()
  for(chr in 1:n_chr) {
    q <- dim(genoprobs$geno[[chr]]$prob)[2]
    marker_names <- paste0("c",chr,"_",names(attr(genoprobs$geno[[chr]]$prob,"map")))
    chr_info <- rbind(rep(chr, q),attr(genoprobs$geno[[chr]]$prob,"map"))
    colnames(chr_info) <- marker_names
    marker_info <- cbind(marker_info, chr_info)
  }
  marker_info <- t(marker_info)
  colnames(marker_info) <- c("chr","pos")
  return(as.data.frame(marker_info))
}

