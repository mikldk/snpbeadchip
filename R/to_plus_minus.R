#' Convert A/B to plus/minus
#' 
#' 
#' @export
to_plus_minus <- function(AB, SNP1, SNP2, ref_strand, sep = "") {
  x <- to_plus_minus_engine(AB, SNP1, SNP2, ref_strand, sep)
  
  return(x)
}

