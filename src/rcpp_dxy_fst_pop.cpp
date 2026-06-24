#include <Rcpp.h>
#include <string.h>
#include "iupac_decode.h"

// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

//' @useDynLib MSA2dist, .registration = TRUE
//' @import Rcpp
//' @title rcpp_dxy_fst_pop
//' @name rcpp_dxy_fst_pop
//' @description calculates dxy and fst per population combination
//' @return list containing:
//' \itemize{
//'     \item dxy
//'     \item fst
//'     \item pi_within
//'     \item pi_between
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
//' rcpp_dxy_fst_pop(dnavector=as.character(iupac),
//' pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
//' model="IUPAC",
//' estimator="count")
//' rcpp_dxy_fst_pop(dnavector=as.character(iupac),
//' pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
//' model="IUPAC",
//' estimator="prob")
//' rcpp_dxy_fst_pop(dnavector=as.character(iupac),
//' pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
//' model="sequence",
//' estimator="count")
//' rcpp_dxy_fst_pop(dnavector=as.character(iupac),
//' pop_idx=c(rep(1,8),rep(2,8),rep(3,8),rep(4,6)),
//' model="sequence",
//' estimator="prob")
//' @export rcpp_dxy_fst_pop
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::List rcpp_dxy_fst_pop(
  Rcpp::StringVector dnavector,
  Rcpp::IntegerVector pop_idx,
  std::string model = "IUPAC",
  std::string estimator = "prob" ) {
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
  NumericMatrix pi_sum(npops, npops);
  IntegerMatrix valid_sites(npops, npops);
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
    for(int p1 = 0; p1 < npops; p1++) {
      if(nAlleles[p1] == 0) continue;
      for(int p2 = p1; p2 < npops; p2++) {
        if(nAlleles[p2] == 0) continue;
        double pi_site = NA_REAL;
        if(p1 == p2) {
          // within-population pi
          if(!useProb) {
            if(nAlleles[p1] < 2) continue;
            double nDiff =
              (double)A[p1]*C[p1] +
              A[p1]*G[p1] +
              A[p1]*T[p1] +
              C[p1]*G[p1] +
              C[p1]*T[p1] +
              G[p1]*T[p1]; 
            double nComp =
              (double)nAlleles[p1] *
              ((double)nAlleles[p1] - 1) / 2.0;
            pi_site = nDiff / nComp;
          } else {
            double total = (double)A[p1] + C[p1] + G[p1] + T[p1];
            if(total == 0) continue;
            double pA = A[p1] / total;
            double pC = C[p1] / total;
            double pG = G[p1] / total;
            double pT = T[p1] / total;
            // expected heterozygosity (within-site pi)
            pi_site =
              1.0 -
              (pA*pA + pC*pC + pG*pG + pT*pT);
          }
        } else {
          // between-population dxy
          if(!useProb) {
            double nComp =
              (double)nAlleles[p1] *
              (double)nAlleles[p2];
            if(nComp == 0) continue;
            double nDiff =
              (double)A[p1] * (C[p2] + G[p2] + T[p2]) +
              (double)C[p1] * (A[p2] + G[p2] + T[p2]) +
              (double)G[p1] * (A[p2] + C[p2] + T[p2]) +
              (double)T[p1] * (A[p2] + C[p2] + G[p2]);
            pi_site = nDiff / nComp;
          } else {
            double total1 = (double)A[p1] + C[p1] + G[p1] + T[p1];
            double total2 = (double)A[p2] + C[p2] + G[p2] + T[p2];
            if(total1 == 0 || total2 == 0) continue;
            double pA1 = A[p1] / total1;
            double pC1 = C[p1] / total1;
            double pG1 = G[p1] / total1;
            double pT1 = T[p1] / total1;
            double pA2 = A[p2] / total2;
            double pC2 = C[p2] / total2;
            double pG2 = G[p2] / total2;
            double pT2 = T[p2] / total2;
            pi_site =
              1.0 -
              (pA1*pA2 + pC1*pC2 + pG1*pG2 + pT1*pT2);
          }
        }
        pi_sum(p1, p2) += pi_site;
        valid_sites(p1, p2)++;
        if(p1 != p2) {
          pi_sum(p2, p1) += pi_site;
          valid_sites(p2, p1)++;
        }
      }
    }
  }
  NumericMatrix pi_between(npops, npops);
  for(int i = 0; i < npops; i++) {
    for(int j = 0; j < npops; j++) {
      if(valid_sites(i, j) > 0) {
        pi_between(i, j) = pi_sum(i, j) /
          (double)valid_sites(i, j);
      } else {
        pi_between(i, j) = NA_REAL;
      }
    }
  }
  NumericMatrix pi_within(npops, npops);
  NumericMatrix dxy(npops, npops);
  NumericMatrix fst(npops, npops);
  for(int i = 0; i < npops; i++) {
    double pi_i = pi_between(i, i);
    for(int j = 0; j < npops; j++) {
      double pi_j = pi_between(j, j);
      pi_within(i, j) =
        (pi_i + pi_j) / 2.0;
      dxy(i, j) = pi_between(i, j);
      if(i == j) {
        fst(i, j) = 0.0;
      } else if(NumericMatrix::is_na(dxy(i, j)) || dxy(i, j) <= 0.0) {
        fst(i, j) = NA_REAL;
      } else {
        fst(i, j) = (dxy(i, j) - pi_within(i, j)) / dxy(i, j);
      }
    }
  }
  IntegerVector populations(npops);
  for(int p = 0; p < npops; p++) {
    populations[p] = p + 1;
  }
  return Rcpp::List::create(
    Rcpp::Named("dxy") = dxy,
    Rcpp::Named("fst") = fst,
    Rcpp::Named("pi_within") = pi_within,
    Rcpp::Named("pi_between") = pi_between,
    Rcpp::Named("valid_sites") = valid_sites,
    Rcpp::Named("populations") = populations,
    Rcpp::Named("model") = model,
    Rcpp::Named("estimator") = estimator);
}
