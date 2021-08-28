#' Add idat information to manifest
#' 
#' @param manifest Object with manifest, e.g. from [load_manifest()]
#' @param idat Object with idat info, e.g. from [read_idat_single()]
#' 
#' @export
annotate_manifest_with_idat <- function(manifest, idat) {
  stopifnot("AddressA_ID" %in% colnames(manifest))
  stopifnot("AddressB_ID" %in% colnames(manifest))
  stopifnot(!("Name" %in% colnames(manifest)))
  
  
  stopifnot("Name" %in% colnames(idat))
  stopifnot("Red_Mean" %in% colnames(idat))
  stopifnot("Red_SD" %in% colnames(idat))
  stopifnot("Red_NBeads" %in% colnames(idat))
  stopifnot("Grn_Mean" %in% colnames(idat))
  stopifnot("Grn_SD" %in% colnames(idat))
  stopifnot("Grn_NBeads" %in% colnames(idat))

  
  # Avoid work with restoring
  stopifnot(!haskey(idat))
  stopifnot(!haskey(manifest))
  
  #############
  
  idat_intensities <- idat[, c("Name", "Red_Mean", "Red_SD", "Red_NBeads", "Grn_Mean", "Grn_SD", "Grn_NBeads")]
  
  setkey(idat_intensities, Name)
  
  #############
  
  setkey(manifest, AddressA_ID)
  m_idat_A <- idat_intensities[manifest]
  
  setnames(m_idat_A, "Name", "AddressA_ID")
  setnames(m_idat_A, "Red_Mean", "A_Red_Mean")
  setnames(m_idat_A, "Red_SD", "A_Red_SD")
  setnames(m_idat_A, "Red_NBeads", "A_Red_NBeads")
  setnames(m_idat_A, "Grn_Mean", "A_Grn_Mean")
  setnames(m_idat_A, "Grn_SD", "A_Grn_SD")
  setnames(m_idat_A, "Grn_NBeads", "A_Grn_NBeads")
  
  #############
  
  setkey(m_idat_A, AddressB_ID)
  m_idat_AB <- idat_intensities[m_idat_A]
  
  setnames(m_idat_AB, "Name", "AddressB_ID")
  setnames(m_idat_AB, "Red_Mean", "B_Red_Mean")
  setnames(m_idat_AB, "Red_SD", "B_Red_SD")
  setnames(m_idat_AB, "Red_NBeads", "B_Red_NBeads")
  setnames(m_idat_AB, "Grn_Mean", "B_Grn_Mean")
  setnames(m_idat_AB, "Grn_SD", "B_Grn_SD")
  setnames(m_idat_AB, "Grn_NBeads", "B_Grn_NBeads")
  
  stopifnot(nrow(manifest) == nrow(m_idat_AB))
  
  stopifnot(ncol(m_idat_AB) >= 15L)
  setcolorder(m_idat_AB, c(15:ncol(m_idat_AB), 8, 1, 9:14, 2:7))
  
  setkey(manifest, NULL)
  setkey(idat_intensities, NULL)
  setkey(m_idat_A, NULL)
  setkey(m_idat_AB, NULL)
  
  return(m_idat_AB)
}

