#' Channel probe to alleles
#' 
#' @param d Manifest annotated with idat e.g. from [annotate_manifest_with_idat()]
#' @param top_n Only take `top_n` rows, if `0` (default), take all rows
#' 
#' @export
channel_probe_to_alleles <- function(annotated_manifest, top_n = 0) {
  stopifnot("AddressA_ID" %in% colnames(annotated_manifest))
  stopifnot("AddressB_ID" %in% colnames(annotated_manifest))
  
  stopifnot("SNP" %in% colnames(annotated_manifest))
  stopifnot("IlmnStrand" %in% colnames(annotated_manifest))
  stopifnot("TopGenomicSeqSBE_Left" %in% colnames(annotated_manifest))
  stopifnot("TopGenomicSeqSBE_Right" %in% colnames(annotated_manifest))
  stopifnot("RefStrand" %in% colnames(annotated_manifest))
  stopifnot("ProbeType" %in% colnames(annotated_manifest))
  stopifnot("SNPType" %in% colnames(annotated_manifest))

  stopifnot("A_Red_Mean" %in% colnames(annotated_manifest))
  stopifnot("A_Red_SD" %in% colnames(annotated_manifest))
  stopifnot("A_Red_NBeads" %in% colnames(annotated_manifest))
  stopifnot("A_Grn_Mean" %in% colnames(annotated_manifest))
  stopifnot("A_Grn_SD" %in% colnames(annotated_manifest))
  stopifnot("A_Grn_NBeads" %in% colnames(annotated_manifest))
  
  stopifnot("B_Red_Mean" %in% colnames(annotated_manifest))
  stopifnot("B_Red_SD" %in% colnames(annotated_manifest))
  stopifnot("B_Red_NBeads" %in% colnames(annotated_manifest))
  stopifnot("B_Grn_Mean" %in% colnames(annotated_manifest))
  stopifnot("B_Grn_SD" %in% colnames(annotated_manifest))
  stopifnot("B_Grn_NBeads" %in% colnames(annotated_manifest))
  
  l <- channel_probe_to_alleles_engine(annotated_manifest, top_n)
  d <- as.data.frame(l)
  
  return(d)
}
