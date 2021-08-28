#' @importFrom illuminaio readIDAT
read_idat <- function(path) {
  # d <- read_idat("/home/mikl/gits/mathylation/omni54-sensitivity/data/Sensitivity_firstrunOMNI54/204163550001/204163550001_R01C01_Red.idat")
  stopifnot(length(path) == 1L)
  stopifnot(file.exists(path))
  
  x <- illuminaio::readIDAT(path, what = "all")
  d <- as.data.frame(x$Quants)
  d$Name <- as.integer(rownames(d))
  rownames(d) <- NULL
  #d <- as_tibble(x$Quants, rownames = "Name")
  
  
  stopifnot(isTRUE(all.equal(colnames(d), c("Mean", "SD", "NBeads", "Name"))))
  d <- d[, c("Name", "Mean", "SD", "NBeads")]
  
  return(d)
}

#' Read idat file from red and green channel from sample sheet
#' 
#' @param datapath Folder with `idat` files
#' @param path_red Filename for the `idat` file with the red channel
#' @param path_grn Filename for the `idat` file with the green channel
#' 
#' @return `data.table` with columns `Name`
#' 
#' @importFrom data.table data.table setkey
#' 
#' @export
read_idat_single <- function(datapath, path_red, path_grn) {
  # d <- read_idat_files("/home/mikl/gits/mathylation/omni54-sensitivity/data/Sensitivity_firstrunOMNI54/204163550001", "204163550001_R01C01_Red.idat", "204163550001_R01C01_Grn.idat")
  
  # RED:
  d_R <- read_idat(path = file.path(datapath, path_red))
  stopifnot(isTRUE(all.equal(colnames(d_R), c("Name", "Mean", "SD", "NBeads"))))
  colnames(d_R)[colnames(d_R) == "Mean"] <- "Red_Mean"
  colnames(d_R)[colnames(d_R) == "SD"] <- "Red_SD"
  colnames(d_R)[colnames(d_R) == "NBeads"] <- "Red_NBeads"
  # <- 
  
  # GRN:
  d_G <- read_idat(path = file.path(datapath, path_grn))
  stopifnot(isTRUE(all.equal(colnames(d_G), c("Name", "Mean", "SD", "NBeads"))))
  colnames(d_G)[colnames(d_G) == "Mean"] <- "Grn_Mean"
  colnames(d_G)[colnames(d_G) == "SD"] <- "Grn_SD"
  colnames(d_G)[colnames(d_G) == "NBeads"] <- "Grn_NBeads"
  # <- 
  
  stopifnot(nrow(d_R) == nrow(d_G))
  
  d_R <- data.table::data.table(d_R)
  data.table::setkey(d_R, Name)
  
  d_G <- data.table::data.table(d_G)
  data.table::setkey(d_G, Name)
  
  d <- d_R[d_G]
  stopifnot(nrow(d) == nrow(d_R))
  
  return(d)
}
# 
# 
# read_idat_samplesheet <- function() {

#   
# }