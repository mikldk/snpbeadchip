#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List call_alleles_engine(Rcpp::DataFrame d, int top_n = 0) {
  Rcpp::List x;
  
  Rcpp::IntegerVector A_Red_Mean = d["A_Red_Mean"];
  Rcpp::IntegerVector A_Red_SD = d["A_Red_SD"];
  Rcpp::IntegerVector A_Red_NBeads = d["A_Red_NBeads"];
  Rcpp::IntegerVector A_Grn_Mean = d["A_Grn_Mean"];
  Rcpp::IntegerVector A_Grn_SD = d["A_Grn_SD"];
  Rcpp::IntegerVector A_Grn_NBeads = d["A_Grn_NBeads"];
  
  Rcpp::IntegerVector B_Red_Mean = d["B_Red_Mean"];
  Rcpp::IntegerVector B_Red_SD = d["B_Red_SD"];
  Rcpp::IntegerVector B_Red_NBeads = d["B_Red_NBeads"];
  Rcpp::IntegerVector B_Grn_Mean = d["B_Grn_Mean"];
  Rcpp::IntegerVector B_Grn_SD = d["B_Grn_SD"];
  Rcpp::IntegerVector B_Grn_NBeads = d["B_Grn_NBeads"];
  
  Rcpp::CharacterMatrix SNP = d["SNP"];
  Rcpp::CharacterMatrix IlmnStrand = d["IlmnStrand"];
  Rcpp::CharacterMatrix TopGenomicSeqSBE_Left = d["TopGenomicSeqSBE_Left"];
  Rcpp::CharacterMatrix TopGenomicSeqSBE_Right = d["TopGenomicSeqSBE_Right"];
  Rcpp::CharacterMatrix RefStrand = d["RefStrand"];
  Rcpp::CharacterMatrix ProbeType = d["ProbeType"];
  Rcpp::CharacterMatrix SNPType = d["SNPType"];
  
  return x;
}
