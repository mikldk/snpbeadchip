#include <algorithm>

#include <Rcpp.h>
using namespace Rcpp;

enum Allele { Allele_NA, Allele_NC, Allele_AA, Allele_BB, Allele_AB };

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
Rcpp::List call_alleles_engine(Rcpp::ListOf<DataFrame>& d_lst, 
                               double cutoff = 3, 
                               int nc_int_thres = -1, 
                               double nc_snr_thres = -1.0, 
                               int nc_sd_thres = -1.0, 
                               int nc_nbeads_thres = -1) {
  // All columns below exists
  
  bool use_int_thres = nc_int_thres > 0;
  bool use_snr_thres = nc_snr_thres > 0;
  bool use_sd_thres = nc_sd_thres > 0;
  bool use_nbeads_thres = nc_nbeads_thres > 0;
  
  size_t n = d_lst.size();
  
  for (size_t i = 0; i < n; ++i) {
    DataFrame d = d_lst[i];
    size_t n_d = d.nrow();
    
    Rcpp::CharacterVector A_base = d["Allele_A_base"];
    Rcpp::CharacterVector B_base = d["Allele_B_base"];
    
    Rcpp::NumericVector AB_frac = d["AB_frac"];
    
    Rcpp::IntegerVector A_mean = d["Allele_A_sig_Mean"];
    Rcpp::IntegerVector B_mean = d["Allele_B_sig_Mean"];
    
    LogicalVector A_NA = is_na(A_mean);
    LogicalVector B_NA = is_na(B_mean);
    
    Rcpp::IntegerVector A_sd;
    Rcpp::IntegerVector B_sd;
    if (use_snr_thres || use_sd_thres) {
      A_sd = d["Allele_A_sig_SD"];
      B_sd = d["Allele_B_sig_SD"];
    }
    
    Rcpp::IntegerVector A_SNR;
    Rcpp::IntegerVector B_SNR;
    if (use_snr_thres) {
      A_SNR = NumericVector(n_d);
      B_SNR = NumericVector(n_d);
      
      for (size_t j = 0; j < n_d; ++j) {
        A_SNR[j] = A_mean[j] / A_sd[j];
        B_SNR[j] = B_mean[j] / B_sd[j];
      }
    }
    
    Rcpp::IntegerVector A_N;
    Rcpp::IntegerVector B_N;
    if (use_nbeads_thres) {
      A_sd = d["Allele_A_sig_N"];
      B_sd = d["Allele_B_sig_N"];
    }
    
    std::vector<Allele> calls(n_d);
    
    for (size_t j = 0; j < n_d; ++j) {
      double ab_frac = AB_frac[j];
      int a = A_mean[j];
      int b = B_mean[j];
      
      int ab_min = 0;
      double snr_min = 0;
      int sd_max = 0;
      int nbeads_min = 0;
      
      Allele call = Allele_NA;
      
      // Then, don't call
      if (A_NA[j] || B_NA[j]) {
        continue;
      }
      
      if (ab_frac < cutoff) {
        call = Allele_AB;
        ab_min = (a < b) ? a : b;
        
        if (use_snr_thres) {
          snr_min = (A_SNR[j] < B_SNR[j]) ? A_SNR[j] : B_SNR[j];
        }
        
        if (use_sd_thres) {
          sd_max = (A_sd[j] > B_sd[j]) ? A_sd[j] : B_sd[j];
        }
        
        if (use_nbeads_thres) {
          nbeads_min = (A_N[j] < B_N[j]) ? A_N[j] : B_N[j];
        }
      } else if (a >= b) {
        call = Allele_AA;
        ab_min = a;
        
        if (use_snr_thres) {
          snr_min = A_SNR[j];
        }
        
        if (use_sd_thres) {
          sd_max = A_sd[j];
        }
        
        if (use_nbeads_thres) {
          nbeads_min = A_N[j];
        }
      } else if (a < b) {
        call = Allele_BB;
        ab_min = b;
        
        if (use_snr_thres) {
          snr_min = B_SNR[j];
        }
        
        if (use_sd_thres) {
          sd_max = B_sd[j];
        }
        
        if (use_nbeads_thres) {
          nbeads_min = B_N[j];
        }
      }
      
      if (use_int_thres && ab_min < nc_int_thres) {
        call = Allele_NC;
      } else if (use_snr_thres && snr_min < nc_snr_thres) { // FIXME: SNR min?
        call = Allele_NC;
      } else if (use_sd_thres && sd_max < nc_sd_thres) {
        call = Allele_NC;
      } else if (use_nbeads_thres && nbeads_min < nc_sd_thres) {
        call = Allele_NC;
      }
      
      calls[j] = call;
    }
    
    
  }
  
  
  List ans = List::create(Named("a") = 0
  );
  
  return ans;
}
