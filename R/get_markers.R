#' Reshaping marker information
#'
#' A function to get marker information from a genetic map including observed and pseudo markers
#' @param gmap Genetic map including observed and pseudomarkers, as obtained from \code{qtl2::insert_pseudomarkers()}.
#' @return  A matrix showing the chromosome numbers (the first column) and positions (the second column) for all markers (row names).
#' @author Yasuhiro Sato (\email{sato.yasuhiro.36c@kyoto-u.jp})
get_markers = function(gmap) {
  n_chr = length(gmap)
  marker_info = c()
  for(chr in 1:n_chr) {
    q = length(gmap[[chr]])
    chr_info = rbind(rep(chr, q),gmap[[chr]])
    marker_info = cbind(marker_info, chr_info)
  }
  marker_info = t(marker_info)
  colnames(marker_info) = c("chr","pos")
  return(marker_info)
}
