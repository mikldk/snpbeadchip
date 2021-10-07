#' Compare call alleles
#' 
#' @param d_lst List of `data.frame`'s
#' @param top_n Only take `top_n` rows, if `0` (default), take all rows
#' 
#' @export
call_alleles_ab_compare <- function(d_lst, cutoff, 
                                    nc_int_thres = NULL, 
                                    nc_snr_thres = NULL, 
                                    nc_sd_thres = NULL, 
                                    nc_nbeads_thres = NULL) {
  
  stopifnot(!missing(d_lst))
  stopifnot(!missing(cutoff))
  
  stopifnot(cutoff >= 1)
  
  stopifnot(is.list(d_lst))
  
  for (d in d_lst) {
    stopifnot(is.data.frame(d))
    
    stopifnot("GType" %in% colnames(d))
    #stopifnot("PlusMinus" %in% colnames(d))
    
    stopifnot("Allele_A_base" %in% colnames(d))
    stopifnot("Allele_B_base" %in% colnames(d))
    
    stopifnot("Allele_A_sig_Mean" %in% colnames(d))
    stopifnot("Allele_B_sig_Mean" %in% colnames(d))
    
    stopifnot("Allele_A_sig_SD" %in% colnames(d))
    stopifnot("Allele_A_sig_SD" %in% colnames(d))
    
    stopifnot("Allele_A_sig_N" %in% colnames(d))
    stopifnot("Allele_A_sig_N" %in% colnames(d))
    
    stopifnot("ABmax" %in% colnames(d))
    stopifnot("ABmin" %in% colnames(d))
    stopifnot("AB_frac" %in% colnames(d))
    
  }
  
  if (is.null(nc_int_thres)) {
    nc_int_thres <- -1
  }
  
  if (is.null(nc_snr_thres)) {
    nc_snr_thres <- -1
  }
  
  if (is.null(nc_sd_thres)) {
    nc_sd_thres <- -1
  }
  
  if (is.null(nc_nbeads_thres)) {
    nc_nbeads_thres <- -1
  }
  
  
  ans <- call_alleles_ab_compare_engine(d_lst = d_lst, 
                                        cutoff = cutoff,
                                        nc_int_thres = nc_int_thres, 
                                        nc_snr_thres = nc_snr_thres, 
                                        nc_sd_thres = nc_sd_thres, 
                                        nc_nbeads_thres = nc_nbeads_thres)
  
  return(ans)
}


#' Call alleles
#' 
#' @param d_lst List of `data.frame`'s
#' @param top_n Only take `top_n` rows, if `0` (default), take all rows
#' 
#' @export
call_alleles_ab <- function(d_lst, cutoff, 
                            nc_int_thres = NULL, 
                            nc_snr_thres = NULL, 
                            nc_sd_thres = NULL, 
                            nc_nbeads_thres = NULL) {
  
  stopifnot(!missing(d_lst))
  stopifnot(!missing(cutoff))
  
  stopifnot(cutoff >= 1)
  
  stopifnot(is.list(d_lst))
  
  for (d in d_lst) {
    stopifnot(is.data.frame(d))
    
    #stopifnot("GType" %in% colnames(d))
    #stopifnot("PlusMinus" %in% colnames(d))
    
    stopifnot("Allele_A_base" %in% colnames(d))
    stopifnot("Allele_B_base" %in% colnames(d))
    
    stopifnot("Allele_A_sig_Mean" %in% colnames(d))
    stopifnot("Allele_B_sig_Mean" %in% colnames(d))
    
    stopifnot("Allele_A_sig_SD" %in% colnames(d))
    stopifnot("Allele_A_sig_SD" %in% colnames(d))
    
    stopifnot("Allele_A_sig_N" %in% colnames(d))
    stopifnot("Allele_A_sig_N" %in% colnames(d))
    
    stopifnot("ABmax" %in% colnames(d))
    stopifnot("ABmin" %in% colnames(d))
    stopifnot("AB_frac" %in% colnames(d))
    
  }
  
  if (is.null(nc_int_thres)) {
    nc_int_thres <- -1
  }
  
  if (is.null(nc_snr_thres)) {
    nc_snr_thres <- -1
  }
  
  if (is.null(nc_sd_thres)) {
    nc_sd_thres <- -1
  }
  
  if (is.null(nc_nbeads_thres)) {
    nc_nbeads_thres <- -1
  }
  
  
  ans <- call_alleles_ab_engine(d_lst = d_lst, 
                                cutoff = cutoff,
                                nc_int_thres = nc_int_thres, 
                                nc_snr_thres = nc_snr_thres, 
                                nc_sd_thres = nc_sd_thres, 
                                nc_nbeads_thres = nc_nbeads_thres)
  
  names(ans) <- names(d_lst)
  
  return(ans)
}
