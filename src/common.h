#ifndef COMMON_H
#define COMMON_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

inline std::string compl_base(std::string base) {
  if (base == "A") {
    return "T";
  } else if (base == "T") {
    return "A";
  } else if (base == "C") {
    return "G";
  } else if (base == "G") {
    return "C";
  }
  
  Rcpp::print(Rcpp::wrap(base));
  Rcpp::stop("Unknown base ^ to get complementary base of");
} 

#endif /* COMMON_H */