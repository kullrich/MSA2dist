#include <Rcpp.h>
#include <string.h>
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
using namespace Rcpp;

//' @useDynLib MSA2dist, .registration = TRUE
//' @import Rcpp
//' @title rcpp_pairwiseDeletionDNA
//' @name rcpp_pairwiseDeletionDNA
//' @description returns number of DNA sites used
//' @return list
//' @param dnavector StringVector
//' @param ncores number of cores
//' @examples
//' ## load example sequence data
//' data("woodmouse", package="ape")
//' w <- woodmouse |> dnabin2dnastring() |> as.character()
//' rcpp_pairwiseDeletionDNA(dnavector=w, ncores=1)
//' @export rcpp_pairwiseDeletionDNA
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::List rcpp_pairwiseDeletionDNA( Rcpp::StringVector dnavector, int ncores = 1 ) {
  std::unordered_map<std::string, double> dist_mat;
  dist_mat["AA"]=0.0;dist_mat["AC"]=1.0;dist_mat["AG"]=1.0;dist_mat["AT"]=1.0;dist_mat["AR"]=-1.0;dist_mat["AY"]=-1.0;dist_mat["AS"]=-1.0;dist_mat["AW"]=-1.0;dist_mat["AK"]=-1.0;dist_mat["AM"]=-1.0;dist_mat["AB"]=-1.0;dist_mat["AD"]=-1.0;dist_mat["AH"]=-1.0;dist_mat["AV"]=-1.0;dist_mat["A."]=-1.0;dist_mat["A-"]=-1.0;dist_mat["AN"]=-1.0;dist_mat["AX"]=-1.0;
  dist_mat["CA"]=1.0;dist_mat["CC"]=0.0;dist_mat["CG"]=1.0;dist_mat["CT"]=1.0;dist_mat["CR"]=-1.0;dist_mat["CY"]=-1.0;dist_mat["CS"]=-1.0;dist_mat["CW"]=-1.0;dist_mat["CK"]=-1.0;dist_mat["CM"]=-1.0;dist_mat["CB"]=-1.0;dist_mat["CD"]=-1.0;dist_mat["CH"]=-1.0;dist_mat["CV"]=-1.0;dist_mat["C."]=-1.0;dist_mat["C-"]=-1.0;dist_mat["CN"]=-1.0;dist_mat["CX"]=-1.0;
  dist_mat["GA"]=1.0;dist_mat["GC"]=1.0;dist_mat["GG"]=0.0;dist_mat["GT"]=1.0;dist_mat["GR"]=-1.0;dist_mat["GY"]=-1.0;dist_mat["GS"]=-1.0;dist_mat["GW"]=-1.0;dist_mat["GK"]=-1.0;dist_mat["GM"]=-1.0;dist_mat["GB"]=-1.0;dist_mat["GD"]=-1.0;dist_mat["GH"]=-1.0;dist_mat["GV"]=-1.0;dist_mat["G."]=-1.0;dist_mat["G-"]=-1.0;dist_mat["GN"]=-1.0;dist_mat["GX"]=-1.0;
  dist_mat["TA"]=1.0;dist_mat["TC"]=1.0;dist_mat["TG"]=1.0;dist_mat["TT"]=0.0;dist_mat["TR"]=-1.0;dist_mat["TY"]=-1.0;dist_mat["TS"]=-1.0;dist_mat["TW"]=-1.0;dist_mat["TK"]=-1.0;dist_mat["TM"]=-1.0;dist_mat["TB"]=-1.0;dist_mat["TD"]=-1.0;dist_mat["TH"]=-1.0;dist_mat["TV"]=-1.0;dist_mat["T."]=-1.0;dist_mat["T-"]=-1.0;dist_mat["TN"]=-1.0;dist_mat["TX"]=-1.0;
  dist_mat["RA"]=-1.0;dist_mat["RC"]=-1.0;dist_mat["RG"]=-1.0;dist_mat["RT"]=-1.0;dist_mat["RR"]=-1.0;dist_mat["RY"]=-1.0;dist_mat["RS"]=-1.0;dist_mat["RW"]=-1.0;dist_mat["RK"]=-1.0;dist_mat["RM"]=-1.0;dist_mat["RB"]=-1.0;dist_mat["RD"]=-1.0;dist_mat["RH"]=-1.0;dist_mat["RV"]=-1.0;dist_mat["R."]=-1.0;dist_mat["R-"]=-1.0;dist_mat["RN"]=-1.0;dist_mat["RX"]=-1.0;
  dist_mat["YA"]=-1.0;dist_mat["YC"]=-1.0;dist_mat["YG"]=-1.0;dist_mat["YT"]=-1.0;dist_mat["YR"]=-1.0;dist_mat["YY"]=-1.0;dist_mat["YS"]=-1.0;dist_mat["YW"]=-1.0;dist_mat["YK"]=-1.0;dist_mat["YM"]=-1.0;dist_mat["YB"]=-1.0;dist_mat["YD"]=-1.0;dist_mat["YH"]=-1.0;dist_mat["YV"]=-1.0;dist_mat["Y."]=-1.0;dist_mat["Y-"]=-1.0;dist_mat["YN"]=-1.0;dist_mat["YX"]=-1.0;
  dist_mat["SA"]=-1.0;dist_mat["SC"]=-1.0;dist_mat["SG"]=-1.0;dist_mat["ST"]=-1.0;dist_mat["SR"]=-1.0;dist_mat["SY"]=-1.0;dist_mat["SS"]=-1.0;dist_mat["SW"]=-1.0;dist_mat["SK"]=-1.0;dist_mat["SM"]=-1.0;dist_mat["SB"]=-1.0;dist_mat["SD"]=-1.0;dist_mat["SH"]=-1.0;dist_mat["SV"]=-1.0;dist_mat["S."]=-1.0;dist_mat["S-"]=-1.0;dist_mat["SN"]=-1.0;dist_mat["SX"]=-1.0;
  dist_mat["WA"]=-1.0;dist_mat["WC"]=-1.0;dist_mat["WG"]=-1.0;dist_mat["WT"]=-1.0;dist_mat["WR"]=-1.0;dist_mat["WY"]=-1.0;dist_mat["WS"]=-1.0;dist_mat["WW"]=-1.0;dist_mat["WK"]=-1.0;dist_mat["WM"]=-1.0;dist_mat["WB"]=-1.0;dist_mat["WD"]=-1.0;dist_mat["WH"]=-1.0;dist_mat["WV"]=-1.0;dist_mat["W."]=-1.0;dist_mat["W-"]=-1.0;dist_mat["WN"]=-1.0;dist_mat["WX"]=-1.0;
  dist_mat["KA"]=-1.0;dist_mat["KC"]=-1.0;dist_mat["KG"]=-1.0;dist_mat["KT"]=-1.0;dist_mat["KR"]=-1.0;dist_mat["KY"]=-1.0;dist_mat["KS"]=-1.0;dist_mat["KW"]=-1.0;dist_mat["KK"]=-1.0;dist_mat["KM"]=-1.0;dist_mat["KB"]=-1.0;dist_mat["KD"]=-1.0;dist_mat["KH"]=-1.0;dist_mat["KV"]=-1.0;dist_mat["K."]=-1.0;dist_mat["K-"]=-1.0;dist_mat["KN"]=-1.0;dist_mat["KX"]=-1.0;
  dist_mat["MA"]=-1.0;dist_mat["MC"]=-1.0;dist_mat["MG"]=-1.0;dist_mat["MT"]=-1.0;dist_mat["MR"]=-1.0;dist_mat["MY"]=-1.0;dist_mat["MS"]=-1.0;dist_mat["MW"]=-1.0;dist_mat["MK"]=-1.0;dist_mat["MM"]=-1.0;dist_mat["MB"]=-1.0;dist_mat["MD"]=-1.0;dist_mat["MH"]=-1.0;dist_mat["MV"]=-1.0;dist_mat["M."]=-1.0;dist_mat["M-"]=-1.0;dist_mat["MN"]=-1.0;dist_mat["MX"]=-1.0;
  dist_mat["BA"]=-1.0;dist_mat["BC"]=-1.0;dist_mat["BG"]=-1.0;dist_mat["BT"]=-1.0;dist_mat["BR"]=-1.0;dist_mat["BY"]=-1.0;dist_mat["BS"]=-1.0;dist_mat["BW"]=-1.0;dist_mat["BK"]=-1.0;dist_mat["BM"]=-1.0;dist_mat["BB"]=-1.0;dist_mat["BD"]=-1.0;dist_mat["BH"]=-1.0;dist_mat["BV"]=-1.0;dist_mat["B."]=-1.0;dist_mat["B-"]=-1.0;dist_mat["BN"]=-1.0;dist_mat["BX"]=-1.0;
  dist_mat["DA"]=-1.0;dist_mat["DC"]=-1.0;dist_mat["DG"]=-1.0;dist_mat["DT"]=-1.0;dist_mat["DR"]=-1.0;dist_mat["DY"]=-1.0;dist_mat["DS"]=-1.0;dist_mat["DW"]=-1.0;dist_mat["DK"]=-1.0;dist_mat["DM"]=-1.0;dist_mat["DB"]=-1.0;dist_mat["DD"]=-1.0;dist_mat["DH"]=-1.0;dist_mat["DV"]=-1.0;dist_mat["D."]=-1.0;dist_mat["D-"]=-1.0;dist_mat["DN"]=-1.0;dist_mat["DX"]=-1.0;
  dist_mat["HA"]=-1.0;dist_mat["HC"]=-1.0;dist_mat["HG"]=-1.0;dist_mat["HT"]=-1.0;dist_mat["HR"]=-1.0;dist_mat["HY"]=-1.0;dist_mat["HS"]=-1.0;dist_mat["HW"]=-1.0;dist_mat["HK"]=-1.0;dist_mat["HM"]=-1.0;dist_mat["HB"]=-1.0;dist_mat["HD"]=-1.0;dist_mat["HH"]=-1.0;dist_mat["HV"]=-1.0;dist_mat["H."]=-1.0;dist_mat["H-"]=-1.0;dist_mat["HN"]=-1.0;dist_mat["HX"]=-1.0;
  dist_mat["VA"]=-1.0;dist_mat["VC"]=-1.0;dist_mat["VG"]=-1.0;dist_mat["VT"]=-1.0;dist_mat["VR"]=-1.0;dist_mat["VY"]=-1.0;dist_mat["VS"]=-1.0;dist_mat["VW"]=-1.0;dist_mat["VK"]=-1.0;dist_mat["VM"]=-1.0;dist_mat["VB"]=-1.0;dist_mat["VD"]=-1.0;dist_mat["VH"]=-1.0;dist_mat["VV"]=-1.0;dist_mat["V."]=-1.0;dist_mat["V-"]=-1.0;dist_mat["VN"]=-1.0;dist_mat["VX"]=-1.0;
  dist_mat[".A"]=-1.0;dist_mat[".C"]=-1.0;dist_mat[".G"]=-1.0;dist_mat[".T"]=-1.0;dist_mat[".R"]=-1.0;dist_mat[".Y"]=-1.0;dist_mat[".S"]=-1.0;dist_mat[".W"]=-1.0;dist_mat[".K"]=-1.0;dist_mat[".M"]=-1.0;dist_mat[".B"]=-1.0;dist_mat[".D"]=-1.0;dist_mat[".H"]=-1.0;dist_mat[".V"]=-1.0;dist_mat[".."]=-1.0;dist_mat[".-"]=-1.0;dist_mat[".N"]=-1.0;dist_mat[".X"]=-1.0;
  dist_mat["-A"]=-1.0;dist_mat["-C"]=-1.0;dist_mat["-G"]=-1.0;dist_mat["-T"]=-1.0;dist_mat["-R"]=-1.0;dist_mat["-Y"]=-1.0;dist_mat["-S"]=-1.0;dist_mat["-W"]=-1.0;dist_mat["-K"]=-1.0;dist_mat["-M"]=-1.0;dist_mat["-B"]=-1.0;dist_mat["-D"]=-1.0;dist_mat["-H"]=-1.0;dist_mat["-V"]=-1.0;dist_mat["-."]=-1.0;dist_mat["--"]=-1.0;dist_mat["-N"]=-1.0;dist_mat["-X"]=-1.0;
  dist_mat["NA"]=-1.0;dist_mat["NC"]=-1.0;dist_mat["NG"]=-1.0;dist_mat["NT"]=-1.0;dist_mat["NR"]=-1.0;dist_mat["NY"]=-1.0;dist_mat["NS"]=-1.0;dist_mat["NW"]=-1.0;dist_mat["NK"]=-1.0;dist_mat["NM"]=-1.0;dist_mat["NB"]=-1.0;dist_mat["ND"]=-1.0;dist_mat["NH"]=-1.0;dist_mat["NV"]=-1.0;dist_mat["N."]=-1.0;dist_mat["N-"]=-1.0;dist_mat["NN"]=-1.0;dist_mat["NX"]=-1.0;
  dist_mat["XA"]=-1.0;dist_mat["XC"]=-1.0;dist_mat["XG"]=-1.0;dist_mat["XT"]=-1.0;dist_mat["XR"]=-1.0;dist_mat["XY"]=-1.0;dist_mat["XS"]=-1.0;dist_mat["XW"]=-1.0;dist_mat["XK"]=-1.0;dist_mat["XM"]=-1.0;dist_mat["XB"]=-1.0;dist_mat["XD"]=-1.0;dist_mat["XH"]=-1.0;dist_mat["XV"]=-1.0;dist_mat["X."]=-1.0;dist_mat["X-"]=-1.0;dist_mat["XN"]=-1.0;dist_mat["XX"]=-1.0;
  int n = dnavector.size();
  Rcpp::NumericMatrix distMatrix(n, n);
  CharacterVector dnavectornames = dnavector.attr("names");
  colnames(distMatrix) = dnavectornames;
  rownames(distMatrix) = dnavectornames;
  Rcpp::NumericMatrix sitesMatrix(n, n);
  colnames(sitesMatrix) = dnavectornames;
  rownames(sitesMatrix) = dnavectornames;
  int nsites = dnavector[1].size();
  RcppThread::ProgressBar bar(n, 1);
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
    };
    bar++;
  }, ncores);
  return Rcpp::List::create(Rcpp::Named("sitesUsed") = sitesMatrix);
}
