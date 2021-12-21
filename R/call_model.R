#' Convert cluster to A/B genotype
#' 
#' @export
cluster_to_AB <- function(x) {
  y <- x
  y[x == 1L] <- "BB"
  y[x == 2L] <- "AB"
  y[x == 3L] <- "AA"
  y
}


#' Get call model
#' 
#' Note that:
#' 
#' * 1 means BB
#' * 2 means AB
#' * 3 means AA
#' 
#' @export
get_call_model <- function() {
  structure(list(call = NULL, 
                 data = structure(c(NA_real_, NA_real_), .Dim = 1:2), 
                 modelName = "VVV", n = NA_integer_, d = 2L, 
                 G = 3L, BIC = NA_real_, loglik = NA_real_, df = 17, bic = NA_real_, 
                 icl = NA_real_, hypvol = NA_real_, parameters = list(pro = c(0.500479788221025, 
                                                                              0.143228716138262, 0.356291495640446), mean = structure(c(6.96177672254813, 
                                                                                                                                        8.85489676984941, 8.98008865368302, 8.39050131899442, 9.48502667893542, 
                                                                                                                                        6.16325623969781), .Dim = 2:3, .Dimnames = list(c("x", "y"
                                                                                                                                        ), NULL)), variance = list(modelName = "VVV", d = 2L, G = 3L, 
                                                                                                                                                                   sigma = structure(c(0.183408627067186, 0.0587945717557889, 
                                                                                                                                                                                       0.0587945717557889, 0.106750959122843, 0.532778894843797, 
                                                                                                                                                                                       0.455477556606343, 0.455477556606343, 0.443235723673907, 
                                                                                                                                                                                       0.122680156614798, 0.0661433759181621, 0.0661433759181621, 
                                                                                                                                                                                       0.275258041926924), .Dim = c(2L, 2L, 3L), .Dimnames = list(
                                                                                                                                                                                         c("x", "y"), c("x", "y"), NULL)), cholsigma = structure(c(-0.428262334401692, 
                                                                                                                                                                                                                                                   0, -0.137286347719391, 0.296485105616983, -0.729917046549672, 
                                                                                                                                                                                                                                                   0, -0.62401276797054, 0.232042644967799, -0.350257272036996, 
                                                                                                                                                                                                                                                   0, -0.188842263098468, 0.489486099490857), .Dim = c(2L, 
                                                                                                                                                                                                                                                                                                       2L, 3L), .Dimnames = list(c("x", "y"), c("x", "y"), NULL)))), 
                 z = structure(c(NA_real_, NA_real_), .Dim = 1:2), classification = NA_real_, 
                 uncertainty = NA_real_), class = "Mclust")
}

