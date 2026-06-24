#include <Rcpp.h>
#include <string.h>
#include "iupac_decode.h"

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

//' @useDynLib MSA2dist, .registration = TRUE
//' @import Rcpp
//' @title rcpp_weightedPi
//' @name rcpp_weightedPi
//' @description calculates weighted pi
//' @return list containing:
//' \itemize{
//'     \item weightedPi
//'     \item model
//'     \item estimator
//' }
//' @param dnavector StringVector [mandatory]
//' @param model string [default: IUPAC]
//' @param estimator string [default: count]
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' rcpp_weightedPi(dnavector=as.character(hiv), model="IUPAC",
//' estimator="count")
//' rcpp_weightedPi(dnavector=as.character(hiv), model="sequence",
//' estimator="count")
//' rcpp_weightedPi(dnavector=as.character(hiv), model="IUPAC",
//' estimator="prob")
//' rcpp_weightedPi(dnavector=as.character(hiv), model="sequence",
//' estimator="prob")
//' @export rcpp_weightedPi
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::List rcpp_weightedPi( Rcpp::StringVector dnavector, std::string model = "IUPAC",
  std::string estimator = "count" ) {
  if(model != "IUPAC" && model != "sequence") {
    Rcpp::stop("Invalid model. Use 'IUPAC' or 'sequence'.");
  }
  if(estimator != "count" && estimator != "prob") {
    Rcpp::stop("Invalid estimator. Use 'count' or 'prob'.");
  }
  int n = dnavector.size();
  CharacterVector dnavectornames = dnavector.attr("names");
  int nsites=dnavector[0].size();
  double pi_sum = 0.0;
  int valid_sites = 0;
  Rcpp::NumericVector out(1, NA_REAL);
  bool useIUPAC = (model == "IUPAC");
  bool useProb = (estimator == "prob");
  for(int s = 0; s < nsites; s++) {
    int A=0, C=0, G=0, T=0;
    int nAlleles = 0;
    for(int i = 0; i < n; i++) {
      std::string seq = Rcpp::as<std::string>(dnavector[i]);
      char base = seq[s];
      if(useIUPAC) {
        decodeIUPAC(seq[s], A, C, G, T, nAlleles);
      } else {
        switch(base) {
          case 'A': A += 1; break;
          case 'C': C += 1; break;
          case 'G': G += 1; break;
          case 'T': T += 1; break;
          default: break; // ignores ambiguous symbols
        }
        nAlleles += 1;
      }
    }
    double pi_site = NA_REAL;
    if(!useProb) {
      if(nAlleles < 2) continue;
      double nDiff = (double)A*C + A*G + A*T + C*G + C*T + G*T;
      double nComp = (double)nAlleles * (nAlleles - 1) / 2.0;
      pi_site = nDiff / nComp;
    } else {
      double total = (double)A + C + G + T;
      if(total == 0) continue;
      double pA = A / total;
      double pC = C / total;
      double pG = G / total;
      double pT = T / total;
      // expected heterozygosity (within-site pi)
      pi_site =
        1.0 -
        (pA*pA + pC*pC + pG*pG + pT*pT);
    }
    pi_sum += pi_site;
    valid_sites++;
  }
  if(valid_sites == 0) {
    out[0] = NA_REAL;
  } else {
    out[0] = pi_sum / valid_sites;
  }
  return Rcpp::List::create(
    Rcpp::Named("weightedPi") = out,
    Rcpp::Named("model") = model,
    Rcpp::Named("estimator") = estimator);
}
