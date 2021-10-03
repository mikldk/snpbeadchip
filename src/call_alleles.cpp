#include <unordered_map>

#include <Rcpp.h>
using namespace Rcpp;

enum class ABAllele : char { 
  NA = 0, 
  NC = 1, 
  AA = 2, 
  BB = 3, 
  AB = 4,
  LAST = 5
};

// Only for pairs of std::hash-able types for simplicity.
// You can of course template this struct to allow other hash functions
struct ABAlleleHash {
  std::size_t operator () (const std::pair<ABAllele, ABAllele> &p) const {
    return static_cast<std::size_t>(p.first)
    * static_cast<std::size_t>(ABAllele::LAST)
    + static_cast<std::size_t>(p.second);
  }
};

DataFrame convert_table_to_df(std::unordered_map<ABAllele, int> table) {
  size_t n = table.size();
  size_t i = 0;
  
  CharacterVector allele(n);
  IntegerVector count(n);
  
  for (auto it = table.begin(); it != table.end(); ++it) {
    ABAllele a = it->first;
    
    if (a == ABAllele::AA) allele[i] = "AA";
    else if (a == ABAllele::BB) allele[i] = "BB";
    else if (a == ABAllele::AB) allele[i] = "AB";
    else if (a == ABAllele::NC) allele[i] = "NC";
    else if (a == ABAllele::NA) allele[i] = "NA";
    else allele[i] = "UNKNOWN";
    
    count[i] = it->second;
    
    ++i;
  }
  
  DataFrame df = DataFrame::create(Named("Call") = allele, 
                                   Named("n") = count);
  return df;
}

DataFrame convert_table_to_df(std::unordered_map<std::pair<ABAllele, ABAllele>, int, ABAlleleHash> table) {
  size_t n = table.size();
  size_t i = 0;
  
  CharacterVector naive_allele(n);
  CharacterVector gs_allele(n);
  IntegerVector count(n);
  
  for (auto it = table.begin(); it != table.end(); ++it) {
    ABAllele naive_a = it->first.first;
    ABAllele gs_a = it->first.second;
    
    if (naive_a == ABAllele::AA) naive_allele[i] = "AA";
    else if (naive_a == ABAllele::BB) naive_allele[i] = "BB";
    else if (naive_a == ABAllele::AB) naive_allele[i] = "AB";
    else if (naive_a == ABAllele::NC) naive_allele[i] = "NC";
    else if (naive_a == ABAllele::NA) naive_allele[i] = "NA";
    else naive_allele[i] = "UNKNOWN";
    
    if (gs_a == ABAllele::AA) gs_allele[i] = "AA";
    else if (gs_a == ABAllele::BB) gs_allele[i] = "BB";
    else if (gs_a == ABAllele::AB) gs_allele[i] = "AB";
    else if (gs_a == ABAllele::NC) gs_allele[i] = "NC";
    else if (gs_a == ABAllele::NA) gs_allele[i] = "NA";
    else gs_allele[i] = "UNKNOWN";
    
    count[i] = it->second;
    
    ++i;
  }
  
  DataFrame df = DataFrame::create(Named("NaiveCall") = naive_allele,
                                   Named("GSCall") = gs_allele, 
                                   Named("n") = count);
  return df;
}

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
  
  std::unordered_map<ABAllele, int> table_naive;
  std::unordered_map<ABAllele, int> table_gs;
  
  std::unordered_map<std::pair<ABAllele, ABAllele>, int, ABAlleleHash> table_compare_naive_gs;

  // FIXME: RcppParallel?
  for (size_t i = 0; i < n; ++i) {
    DataFrame d = d_lst[i];
    size_t n_d = d.nrow();
    
    Rcpp::CharacterVector GS_AB = d["GType"];
    //Rcpp::CharacterVector GS_PM = d["PlusMinus"];
    
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
    
    for (size_t j = 0; j < n_d; ++j) {
      double ab_frac = AB_frac[j];
      int a = A_mean[j];
      int b = B_mean[j];
      
      int ab_min = 0;
      double snr_min = 0;
      int sd_max = 0;
      int nbeads_min = 0;
      
      ABAllele naive_call = ABAllele::NA;
      
      // Then, don't call
      if (A_NA[j] || B_NA[j]) {
        continue;
      }
      
      if (ab_frac < cutoff) {
        naive_call = ABAllele::AB;
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
        naive_call = ABAllele::AA;
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
        naive_call = ABAllele::BB;
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
        naive_call = ABAllele::NC;
      } else if (use_snr_thres && snr_min < nc_snr_thres) { 
        naive_call = ABAllele::NC;
      } else if (use_sd_thres && sd_max < nc_sd_thres) {
        naive_call = ABAllele::NC;
      } else if (use_nbeads_thres && nbeads_min < nc_sd_thres) {
        naive_call = ABAllele::NC;
      }
      
      // naive_call
      ABAllele gs_call = ABAllele::NA;
      
      if (!CharacterVector::is_na(GS_AB[j])) {
        if (GS_AB[j] == "AA") {
          gs_call = ABAllele::AA;
        } else if (GS_AB[j] == "BB") {
          gs_call = ABAllele::BB;
        } else if (GS_AB[j] == "AB") {
          gs_call = ABAllele::AB;
        } else if (GS_AB[j] == "NC") {
          gs_call = ABAllele::NC;
        }
      }
      
      table_naive[naive_call] += 1;
      table_gs[gs_call] += 1;
      
      std::pair<ABAllele, ABAllele> p = std::make_pair(naive_call, gs_call);
      table_compare_naive_gs[p] += 1;
    }
    
    
  }

  List ans = List::create(
    Named("table_naive") = convert_table_to_df(table_naive), 
    Named("table_gs") = convert_table_to_df(table_gs), 
    Named("table_compare_naive_gs") = convert_table_to_df(table_compare_naive_gs)
  );
  
  return ans;
}
