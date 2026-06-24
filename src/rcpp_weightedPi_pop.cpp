#include <Rcpp.h>
#include <string.h>
#include "iupac_decode.h"

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

//' @useDynLib MSA2dist, .registration = TRUE
//' @import Rcpp
//' @title rcpp_weightedPi_pop
//' @name rcpp_weightedPi_pop
//' @description calculates weighted pi per population
//' @return list containing:
//' \itemize{
//'     \item weightedPi
//'     \item populations
//'     \item model
//'     \item estimator
//' }
//' @param dnavector StringVector [mandatory]
//' @param pop_idx IntegerVector [mandatory]
//' @param model string [default: IUPAC]
//' @param estimator string [default: count]
//' @examples
//' ## load example sequence data
//' data("iupac", package="MSA2dist")
//' rcpp_weightedPi_pop(dnavector=as.character(iupac),
//' pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
//' model="IUPAC",
//' estimator="count")
//' rcpp_weightedPi_pop(dnavector=as.character(iupac),
//' pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
//' model="IUPAC",
//' estimator="prob")
//' rcpp_weightedPi_pop(dnavector=as.character(iupac),
//' pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
//' model="sequence",
//' estimator="count")
//' rcpp_weightedPi_pop(dnavector=as.character(iupac),
//' pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
//' model="sequence",
//' estimator="prob")
//' @export rcpp_weightedPi_pop
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::List rcpp_weightedPi_pop(
  Rcpp::StringVector dnavector,
  Rcpp::IntegerVector pop_idx,
  std::string model = "IUPAC",
  std::string estimator = "count" ) {
  if(model != "IUPAC" && model != "sequence") {
    Rcpp::stop("Invalid model. Use 'sequence' or 'IUPAC'.");
  }
  if(estimator != "count" && estimator != "prob") {
    Rcpp::stop("Invalid estimator. Use 'count' or 'prob'.");
  }
  int n = dnavector.size();
  if(pop_idx.size() != n) {
    Rcpp::stop("length(pop_idx) must equal length(dnavector)");
  }
  CharacterVector dnavectornames = dnavector.attr("names");
  std::vector<std::string> seqs(n);
  for(int i = 0; i < n; i++) {
    seqs[i] = Rcpp::as<std::string>(dnavector[i]);
  }
  int nsites = dnavector[0].size();
  int npops = 0;
  for(int i = 0; i < n; i++) {
    if(pop_idx[i] > npops)
      npops = pop_idx[i];
  }
  NumericVector pi_sum(npops, 0.0);
  IntegerVector valid_sites(npops, 0);
  bool useIUPAC = (model == "IUPAC");
  bool useProb = (estimator == "prob");
  for(int s = 0; s < nsites; s++) {
    std::vector<int> A(npops, 0);
    std::vector<int> C(npops, 0);
    std::vector<int> G(npops, 0);
    std::vector<int> T(npops, 0);
    std::vector<int> nAlleles(npops, 0);
    for(int i = 0; i < n; i++) {
      int p = pop_idx[i] - 1; 
      char base = seqs[i][s];
      if(useIUPAC) {
        decodeIUPAC(base, A[p], C[p], G[p], T[p], nAlleles[p]);
      } else {
        switch(base) {
          case 'A': A[p] += 1; nAlleles[p] += 1; break;
          case 'C': C[p] += 1; nAlleles[p] += 1; break;
          case 'G': G[p] += 1; nAlleles[p] += 1; break;
          case 'T': T[p] += 1; nAlleles[p] += 1; break;
          default: break; // ignores ambiguous symbols
        }
      }
    }
    for(int p = 0; p < npops; p++) {
      if(nAlleles[p] == 0) continue;
      double pi_site = NA_REAL;
      if(!useProb) {
        if(nAlleles[p] < 2) continue;
        double nDiff =
          (double)A[p]*C[p] +
          A[p]*G[p] +
          A[p]*T[p] +
          C[p]*G[p] +
          C[p]*T[p] +
          G[p]*T[p]; 
        double nComp =
          (double)nAlleles[p] *
          ((double)nAlleles[p] - 1) / 2.0;
        pi_site = nDiff / nComp;
      } else {
        double total = (double)A[p] + C[p] + G[p] + T[p];
        if(total == 0) continue;
        double pA = A[p] / total;
        double pC = C[p] / total;
        double pG = G[p] / total;
        double pT = T[p] / total;
        // expected heterozygosity (within-site pi)
        pi_site =
          1.0 -
          (pA*pA + pC*pC + pG*pG + pT*pT);
      }
      pi_sum[p] += pi_site;
      valid_sites[p]++;
    }
  }
  NumericVector weightedPi(npops, NA_REAL);
  for(int p = 0; p < npops; p++) {
    if(valid_sites[p] > 0) {
      weightedPi[p] =
        pi_sum[p] /
        (double)valid_sites[p];
    }
  }
  IntegerVector populations(npops);
  for(int p = 0; p < npops; p++) {
    populations[p] = p + 1;
  }
  return Rcpp::List::create(
    Rcpp::Named("weightedPi") = weightedPi,
    Rcpp::Named("populations") = populations,
    Rcpp::Named("model") = model,
    Rcpp::Named("estimator") = estimator);
}
