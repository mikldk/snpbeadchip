#include <algorithm>
#include <string>

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]

#include "common.h"

inline bool is_base_colour_red(std::string base) {
  if (base == "A" || base == "T") {
    return true;
  } 
  
  return false; // Grn
} 

inline std::string to_upper(std::string base) {
  if (base == "A" || base == "T" || base == "C" || base == "G") {
    return base;
  }
  
  if (base == "a") {
    return "A";
  } else if (base == "t") {
    return "T";
  } else if (base == "c") {
    return "C";
  } else if (base == "g") {
    return "G";
  }
  
  Rcpp::print(Rcpp::wrap(base));
  Rcpp::stop("Unknown lower-case base ^ to convert to upper-case");
} 

// [[Rcpp::export]]
Rcpp::List channel_probe_to_alleles_engine(Rcpp::DataFrame& d, size_t top_n) {
  // All columns below exists
  
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
  
  Rcpp::CharacterVector SNP = d["SNP"];
  Rcpp::CharacterVector IlmnStrand = d["IlmnStrand"];
  Rcpp::CharacterVector TopGenomicSeqSBE_Left = d["TopGenomicSeqSBE_Left"];
  Rcpp::CharacterVector TopGenomicSeqSBE_Right = d["TopGenomicSeqSBE_Right"];
  Rcpp::CharacterVector RefStrand = d["RefStrand"];
  Rcpp::CharacterVector ProbeType = d["ProbeType"];
  Rcpp::CharacterVector SNPType = d["SNPType"];
  
  size_t n = d.nrow();
  
  if (top_n > 0 && top_n < n) {
    n = top_n;
  }
  
  Rcpp::CharacterVector Allele_A_base(n); 
  Rcpp::IntegerVector Allele_A_sig_Mean(n);
  Rcpp::IntegerVector Allele_A_sig_SD(n);
  Rcpp::IntegerVector Allele_A_sig_N(n);
  Rcpp::IntegerVector Allele_A_bg_Mean(n);
  Rcpp::IntegerVector Allele_A_bg_SD(n);
  Rcpp::IntegerVector Allele_A_bg_N(n);
  
  Rcpp::CharacterVector Allele_B_base(n); 
  Rcpp::IntegerVector Allele_B_sig_Mean(n);
  Rcpp::IntegerVector Allele_B_sig_SD(n);
  Rcpp::IntegerVector Allele_B_sig_N(n);
  Rcpp::IntegerVector Allele_B_bg_Mean(n);
  Rcpp::IntegerVector Allele_B_bg_SD(n);
  Rcpp::IntegerVector Allele_B_bg_N(n);
  
  for (size_t i = 0; i < n; ++i) {
    
    std::string ref_strand = as<std::string>(RefStrand[i]);
    std::string probe_type = as<std::string>(ProbeType[i]);
    std::string snp_type = as<std::string>(SNPType[i]);
    
    std::string snp = as<std::string>(SNP[i]);
    std::string allele_A = snp.substr(1, 1);
    std::string allele_B = snp.substr(3, 1);
    
    /*
     * Three types to call:
     * 
     * 1) ProbeType == II & SNPType == UNAMB
     * 2) ProbeType == II & SNPType == INDEL
     * 3) ProbeType == I  & SNPType == AMB
     * 
     */
    
    /////////////////////////////////////////////////////////////////////
    // 1) ProbeType == II & SNPType == UNAMB
    /////////////////////////////////////////////////////////////////////
    if (probe_type == "II" && snp_type == "UNAMB") {
      double red_mean = A_Red_Mean[i];
      double red_sd = A_Red_SD[i];
      double red_nbeads = A_Red_NBeads[i];
      
      double grn_mean = A_Grn_Mean[i];
      double grn_sd = A_Grn_SD[i];
      double grn_nbeads = A_Grn_NBeads[i];
      
      // Red is always allele 'A' (A/B), it comes first in the SNP field
      // Grn is always allele 'B' (A/B), it comes second in the SNP field
      
      Allele_A_sig_Mean[i] = red_mean;
      Allele_A_sig_SD[i] = red_sd;
      Allele_A_sig_N[i] = red_nbeads;
      Allele_A_bg_Mean[i] = NA_INTEGER;
      Allele_A_bg_SD[i] = NA_INTEGER;
      Allele_A_bg_N[i] = NA_INTEGER;
      
      Allele_B_sig_Mean[i] = grn_mean;
      Allele_B_sig_SD[i] = grn_sd;
      Allele_B_sig_N[i] = grn_nbeads;
      Allele_B_bg_Mean[i] = NA_INTEGER;
      Allele_B_bg_SD[i] = NA_INTEGER;
      Allele_B_bg_N[i] = NA_INTEGER;
      
      // FIXME: Common for all?
      Allele_A_base[i] = allele_A;
      Allele_B_base[i] = allele_B;
      if (ref_strand == "-") {
        Allele_A_base[i] = compl_base(as<std::string>(Allele_A_base[i]));
        Allele_B_base[i] = compl_base(as<std::string>(Allele_B_base[i]));
      }
    }
    
    /////////////////////////////////////////////////////////////////////
    // 2) ProbeType == II & SNPType == INDEL
    /////////////////////////////////////////////////////////////////////
    else if (probe_type == "II" && snp_type == "INDEL") {
      Allele_A_base[i] = NA_STRING;
      Allele_A_sig_Mean[i] = NA_INTEGER;
      Allele_A_sig_SD[i] = NA_INTEGER;
      Allele_A_sig_N[i] = NA_INTEGER;
      Allele_A_bg_Mean[i] = NA_INTEGER;
      Allele_A_bg_SD[i] = NA_INTEGER;
      Allele_A_bg_N[i] = NA_INTEGER;
      
      Allele_B_base[i] = NA_STRING; 
      Allele_B_sig_Mean[i] = NA_INTEGER;
      Allele_B_sig_SD[i] = NA_INTEGER;
      Allele_B_sig_N[i] = NA_INTEGER;
      Allele_B_bg_Mean[i] = NA_INTEGER;
      Allele_B_bg_SD[i] = NA_INTEGER;
      Allele_B_bg_N[i] = NA_INTEGER;
    } 
    
    /////////////////////////////////////////////////////////////////////
    // 3) ProbeType == I  & SNPType == AMB
    /////////////////////////////////////////////////////////////////////
    else if (probe_type == "I" && snp_type == "AMB") {
      double a_red_mean = A_Red_Mean[i];
      double a_red_sd = A_Red_SD[i];
      double a_red_nbeads = A_Red_NBeads[i];
      double a_grn_mean = A_Grn_Mean[i];
      double a_grn_sd = A_Grn_SD[i];
      double a_grn_nbeads = A_Grn_NBeads[i];
      
      double b_red_mean = B_Red_Mean[i];
      double b_red_sd = B_Red_SD[i];
      double b_red_nbeads = B_Red_NBeads[i];
      double b_grn_mean = B_Grn_Mean[i];
      double b_grn_sd = B_Grn_SD[i];
      double b_grn_nbeads = B_Grn_NBeads[i];
      
      std::string illstrand = as<std::string>(IlmnStrand[i]);
      std::string tgs_left = as<std::string>(TopGenomicSeqSBE_Left[i]);
      std::string tgs_right = as<std::string>(TopGenomicSeqSBE_Right[i]);
      
      // Look at neighbor base
      std::string neighbor_base = (illstrand == "TOP") ? tgs_right : tgs_left;
      
      // lower case: repetitive elements that are soft-masked
      neighbor_base = to_upper(neighbor_base);
      
      bool take_red = is_base_colour_red(neighbor_base);
      
      if (take_red) {
        // Take Red
        
        // Signal all 'Red'
        Allele_A_sig_Mean[i] = a_red_mean;
        Allele_A_sig_SD[i] = a_red_sd;
        Allele_A_sig_N[i] = a_red_nbeads;
        
        Allele_B_sig_Mean[i] = b_red_mean;
        Allele_B_sig_SD[i] = b_red_sd;
        Allele_B_sig_N[i] = b_red_nbeads;
        
        // Background all 'Grn'
        Allele_A_bg_Mean[i] = a_grn_mean;
        Allele_A_bg_SD[i] = a_grn_sd;
        Allele_A_bg_N[i] = a_grn_nbeads;
        
        Allele_B_bg_Mean[i] = b_grn_mean;
        Allele_B_bg_SD[i] = b_grn_sd;
        Allele_B_bg_N[i] = b_grn_nbeads;
      } else {
        // Take Grn
       
        // Signal all 'Grn'
        Allele_A_sig_Mean[i] = a_grn_mean;
        Allele_A_sig_SD[i] = a_grn_sd;
        Allele_A_sig_N[i] = a_grn_nbeads;
        
        Allele_B_sig_Mean[i] = b_grn_mean;
        Allele_B_sig_SD[i] = b_grn_sd;
        Allele_B_sig_N[i] = b_grn_nbeads; 
        
        // Background all 'Red'
        Allele_A_bg_Mean[i] = a_red_mean;
        Allele_A_bg_SD[i] = a_red_sd;
        Allele_A_bg_N[i] = a_red_nbeads;
        
        Allele_B_bg_Mean[i] = b_red_mean;
        Allele_B_bg_SD[i] = b_red_sd;
        Allele_B_bg_N[i] = b_red_nbeads;
      }
 
      // FIXME: Common for all?
      Allele_A_base[i] = allele_A;
      Allele_B_base[i] = allele_B;
      if (ref_strand == "-") {
        Allele_A_base[i] = compl_base(as<std::string>(Allele_A_base[i]));
        Allele_B_base[i] = compl_base(as<std::string>(Allele_B_base[i]));
      }
    } 
    // Something else:
    else {
      Allele_A_base[i] = NA_STRING;
      Allele_A_sig_Mean[i] = NA_INTEGER;
      Allele_A_sig_SD[i] = NA_INTEGER;
      Allele_A_sig_N[i] = NA_INTEGER;
      Allele_A_bg_Mean[i] = NA_INTEGER;
      Allele_A_bg_SD[i] = NA_INTEGER;
      Allele_A_bg_N[i] = NA_INTEGER;
      
      Allele_B_base[i] = NA_STRING; 
      Allele_B_sig_Mean[i] = NA_INTEGER;
      Allele_B_sig_SD[i] = NA_INTEGER;
      Allele_B_sig_N[i] = NA_INTEGER;
      Allele_B_bg_Mean[i] = NA_INTEGER;
      Allele_B_bg_SD[i] = NA_INTEGER;
      Allele_B_bg_N[i] = NA_INTEGER;
    }
    
    
  }
  
  List ans = List::create(Named("Allele_A_base") = Allele_A_base, 
                          Named("Allele_A_sig_Mean") = Allele_A_sig_Mean, 
                          Named("Allele_A_sig_SD") = Allele_A_sig_SD, 
                          Named("Allele_A_sig_N") = Allele_A_sig_N, 
                          Named("Allele_A_bg_Mean") = Allele_A_bg_Mean, 
                          Named("Allele_A_bg_SD") = Allele_A_bg_SD, 
                          Named("Allele_A_bg_N") = Allele_A_bg_N, 
                          
                          Named("Allele_B_base") = Allele_B_base, 
                          Named("Allele_B_sig_Mean") = Allele_B_sig_Mean, 
                          Named("Allele_B_sig_SD") = Allele_B_sig_SD, 
                          Named("Allele_B_sig_N") = Allele_B_sig_N, 
                          Named("Allele_B_bg_Mean") = Allele_B_bg_Mean, 
                          Named("Allele_B_bg_SD") = Allele_B_bg_SD, 
                          Named("Allele_B_bg_N") = Allele_B_bg_N
  );
  
  return ans;
}
