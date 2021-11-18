#include <unordered_map>

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::CharacterVector to_plus_minus(Rcpp::CharacterVector AB,
                                    Rcpp::CharacterVector SNP1,
                                    Rcpp::CharacterVector SNP2,
                                    Rcpp::CharacterVector illumina_strand,
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
  
  if (illumina_strand.size() != n) {
    Rcpp::stop("illumina_strand.size() != SNP1.size()");
  }
  
  
  Rcpp::CharacterVector ans(n);
  
  for (size_t i = 0; i < n; ++i) {
    std::string a1 = as<std::string>(SNP1[i]);
    std::string a2 = as<std::string>(SNP2[i]);
    
    if (illumina_strand[i] == "TOP") {
      std::string tmp = a1;
      a1 = a2;
      a2 = tmp;
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
