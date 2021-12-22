#include <Rcpp.h>
#include <string.h>
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
using namespace Rcpp;

//' @useDynLib MSA2dist, .registration = TRUE
//' @import Rcpp
//' @import RcppThread
//' @title rcpp_distSTRING
//' @name rcpp_distSTRING
//' @description calcualtes pairwise distances using a score matrix
//' @return list
//' @param dnavector StringVector
//' @param scoreMatrix NumericMatrix
//' @param ncores number of cores
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' rcpp_distSTRING(dnavector=as.character(hiv), scoreMatrix=iupacMatrix())
//' @export rcpp_distSTRING
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::List rcpp_distSTRING( Rcpp::StringVector dnavector, Rcpp::NumericMatrix scoreMatrix, int ncores = 1 ) {
  std::unordered_map<std::string, double> dist_mat;
  int nCols = scoreMatrix.ncol();
  int nRows = scoreMatrix.nrow();
  CharacterVector Colnames = colnames(scoreMatrix);
  CharacterVector Rownames = rownames(scoreMatrix);
  for( int is = 0; is < nCols; is++ ){
    for( int js = 0; js < nRows; js++){
      std::string isName = "";
      isName = Colnames(is);
      std::string jsName = "";
      jsName = Rownames(js);
      dist_mat[isName+jsName] = scoreMatrix(is,js);
    }
  }
  int n = dnavector.size();
  Rcpp::NumericMatrix distMatrix(n, n);
  CharacterVector dnavectornames = dnavector.attr("names");
  colnames(distMatrix) = dnavectornames;
  rownames(distMatrix) = dnavectornames;
  Rcpp::NumericMatrix sitesMatrix(n, n);
  colnames(sitesMatrix) = dnavectornames;
  rownames(sitesMatrix) = dnavectornames;
  int nsites = dnavector[1].size();
  RcppThread::parallelFor(0, n, [&] (int i) {
    for( int j=i; j < n; j++ ){
      double eqnum = 0;
      int ij_n = nsites;
      for( int s=0; s < nsites; s++){
        std::string is;
        std::string js;
        is = dnavector[i][s];
        js = dnavector[j][s];
        double ij_dist;
        ij_dist = dist_mat[is+js];
        if(ij_dist >= 0.0){
          eqnum = eqnum + ij_dist;
        } else {
          ij_n = ij_n -1;
        };
      }
      distMatrix(i,j) = eqnum / ij_n;
      distMatrix(j,i) = eqnum / ij_n;
      sitesMatrix(i,j) = ij_n;
      sitesMatrix(j,i) = ij_n;
    }
  }, ncores);
  return Rcpp::List::create(Rcpp::Named("distSTRING") = distMatrix, Rcpp::Named("sitesUsed") = sitesMatrix);
}
