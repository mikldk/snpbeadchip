#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

#include "common.h"

// [[Rcpp::export]]
Rcpp::CharacterVector to_plus_minus_engine(Rcpp::CharacterVector AB,
                                    Rcpp::CharacterVector SNP1,
                                    Rcpp::CharacterVector SNP2,
                                    Rcpp::CharacterVector ref_strand,
                                    Rcpp::CharacterVector sep
                                    ) {
  if (sep.size() != 1) {
    Rcpp::stop("Provide one separator only");
  }
  
  std::string s = as<std::string>(sep[0]);
  
  size_t n = SNP1.size();
  
  if (AB.size() != n) {
    Rcpp::stop("AB.size() != SNP1.size()");
  }
  
  if (SNP2.size() != n) {
    Rcpp::stop("SNP2.size() != SNP1.size()");
  }
  
  if (ref_strand.size() != n) {
    Rcpp::stop("ref_strand.size() != SNP1.size()");
  }
  
  Rcpp::CharacterVector ans(n);
  
  for (size_t i = 0; i < n; ++i) {
    std::string a1 = as<std::string>(SNP1[i]);
    std::string a2 = as<std::string>(SNP2[i]);
    
    if (ref_strand[i] == "-") {
      a1 = compl_base(a1);
      a2 = compl_base(a2);
    }
    
    std::string AB_field = as<std::string>(AB[i]);
    
    if (AB_field == "AA") {
      ans[i] = a1 + s + a1;
    } else if (AB_field == "AB") {
      ans[i] = a1 + s + a2;
    } else if (AB_field == "BB") {
      ans[i] = a2 + s + a2;
    } else {
      ans[i] = "NA";
    }
  }
  
  return ans;
}
