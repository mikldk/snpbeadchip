#' Get call model
#' 
#' @export
get_call_model <- function() {
  structure(list(call = mclust::Mclust(data = ensemple_trans_data, 
                                       G = 3, modelNames = "VVV"), data = structure(c(NA_real_, 
                                                                                      NA_real_), .Dim = 1:2), modelName = "VVV", n = NA_integer_, d = 2L, 
                 G = 3L, BIC = NA_real_, loglik = NA_real_, df = 17, bic = NA_real_, 
                 icl = NA_real_, hypvol = NA_real_, parameters = list(pro = c(0.498022040140032, 
                                                                              0.146789073617397, 0.355188886242239), mean = structure(c(6.96143507313211, 
                                                                                                                                        8.85414997212262, 8.98014875587864, 8.38963426769399, 9.48389133531632, 
                                                                                                                                        6.16303643323933), .Dim = 2:3, .Dimnames = list(c("x", "y"
                                                                                                                                        ), NULL)), variance = list(modelName = "VVV", d = 2L, G = 3L, 
                                                                                                                                                                   sigma = structure(c(0.183363200566762, 0.0587706663628414, 
                                                                                                                                                                                       0.0587706663628414, 0.106594427810971, 0.522749563933128, 
                                                                                                                                                                                       0.446185116002138, 0.446185116002138, 0.434525934578789, 
                                                                                                                                                                                       0.122676916735384, 0.0661559891154721, 0.0661559891154721, 
                                                                                                                                                                                       0.275099317576411), .Dim = c(2L, 2L, 3L), .Dimnames = list(
                                                                                                                                                                                         c("x", "y"), c("x", "y"), NULL)), cholsigma = structure(c(-0.428209295282998, 
                                                                                                                                                                                                                                                   0, -0.137247525941726, 0.296238998840205, -0.72301422111403, 
                                                                                                                                                                                                                                                   0, -0.617118035817676, -0.231713755412411, -0.350252647006962, 
                                                                                                                                                                                                                                                   0, -0.188880768441864, 0.489309077055823), .Dim = c(2L, 
                                                                                                                                                                                                                                                                                                       2L, 3L), .Dimnames = list(c("x", "y"), c("x", "y"), NULL)))), 
                 z = structure(c(NA_real_, NA_real_), .Dim = 1:2), classification = NA_real_, 
                 uncertainty = NA_real_), class = "Mclust")
}

