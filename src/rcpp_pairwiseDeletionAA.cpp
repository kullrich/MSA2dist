#include <Rcpp.h>
#include <string.h>
#include <RcppThread.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
using namespace Rcpp;

//' @useDynLib MSA2dist, .registration = TRUE
//' @import Rcpp
//' @title rcpp_pairwiseDeletionAA
//' @name rcpp_pairwiseDeletionAA
//' @description returns number of AA sites used
//' @return list
//' @param aavector StringVector [mandatory]
//' @param ncores number of cores [default: 1]
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' h <- hiv |> cds2aa() |> as.character()
//' rcpp_pairwiseDeletionAA(aavector=h, ncores=1)
//' @export rcpp_pairwiseDeletionAA
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::List rcpp_pairwiseDeletionAA( Rcpp::StringVector aavector, int ncores = 1 ) {
  std::unordered_map<std::string, double> dist_mat;
  dist_mat["SS"]=0.0;dist_mat["SR"]=1.0;dist_mat["SL"]=1.0;dist_mat["SP"]=1.0;dist_mat["ST"]=1.0;dist_mat["SA"]=1.0;dist_mat["SV"]=1.0;dist_mat["SG"]=1.0;dist_mat["SI"]=1.0;dist_mat["SF"]=1.0;dist_mat["SY"]=1.0;dist_mat["SC"]=1.0;dist_mat["SH"]=1.0;dist_mat["SQ"]=1.0;dist_mat["SN"]=1.0;dist_mat["SK"]=1.0;dist_mat["SD"]=1.0;dist_mat["SE"]=1.0;dist_mat["SM"]=1.0;dist_mat["SW"]=1.0;dist_mat["S."]=-1.0;;dist_mat["S-"]=-1.0;;dist_mat["SX"]=-1.0;
  dist_mat["RS"]=1.0;dist_mat["RR"]=0.0;dist_mat["RL"]=1.0;dist_mat["RP"]=1.0;dist_mat["RT"]=1.0;dist_mat["RA"]=1.0;dist_mat["RV"]=1.0;dist_mat["RG"]=1.0;dist_mat["RI"]=1.0;dist_mat["RF"]=1.0;dist_mat["RY"]=1.0;dist_mat["RC"]=1.0;dist_mat["RH"]=1.0;dist_mat["RQ"]=1.0;dist_mat["RN"]=1.0;dist_mat["RK"]=1.0;dist_mat["RD"]=1.0;dist_mat["RE"]=1.0;dist_mat["RM"]=1.0;dist_mat["RW"]=1.0;dist_mat["R."]=-1.0;;dist_mat["R-"]=-1.0;;dist_mat["RX"]=-1.0;
  dist_mat["LS"]=1.0;dist_mat["LR"]=1.0;dist_mat["LL"]=0.0;dist_mat["LP"]=1.0;dist_mat["LT"]=1.0;dist_mat["LA"]=1.0;dist_mat["LV"]=1.0;dist_mat["LG"]=1.0;dist_mat["LI"]=1.0;dist_mat["LF"]=1.0;dist_mat["LY"]=1.0;dist_mat["LC"]=1.0;dist_mat["LH"]=1.0;dist_mat["LQ"]=1.0;dist_mat["LN"]=1.0;dist_mat["LK"]=1.0;dist_mat["LD"]=1.0;dist_mat["LE"]=1.0;dist_mat["LM"]=1.0;dist_mat["LW"]=1.0;dist_mat["L."]=-1.0;;dist_mat["L-"]=-1.0;;dist_mat["LX"]=-1.0;
  dist_mat["PS"]=1.0;dist_mat["PR"]=1.0;dist_mat["PL"]=1.0;dist_mat["PP"]=0.0;dist_mat["PT"]=1.0;dist_mat["PA"]=1.0;dist_mat["PV"]=1.0;dist_mat["PG"]=1.0;dist_mat["PI"]=1.0;dist_mat["PF"]=1.0;dist_mat["PY"]=1.0;dist_mat["PC"]=1.0;dist_mat["PH"]=1.0;dist_mat["PQ"]=1.0;dist_mat["PN"]=1.0;dist_mat["PK"]=1.0;dist_mat["PD"]=1.0;dist_mat["PE"]=1.0;dist_mat["PM"]=1.0;dist_mat["PW"]=1.0;dist_mat["P."]=-1.0;;dist_mat["P-"]=-1.0;;dist_mat["PX"]=-1.0;
  dist_mat["TS"]=1.0;dist_mat["TR"]=1.0;dist_mat["TL"]=1.0;dist_mat["TP"]=1.0;dist_mat["TT"]=0.0;dist_mat["TA"]=1.0;dist_mat["TV"]=1.0;dist_mat["TG"]=1.0;dist_mat["TI"]=1.0;dist_mat["TF"]=1.0;dist_mat["TY"]=1.0;dist_mat["TC"]=1.0;dist_mat["TH"]=1.0;dist_mat["TQ"]=1.0;dist_mat["TN"]=1.0;dist_mat["TK"]=1.0;dist_mat["TD"]=1.0;dist_mat["TE"]=1.0;dist_mat["TM"]=1.0;dist_mat["TW"]=1.0;dist_mat["T."]=-1.0;;dist_mat["T-"]=-1.0;;dist_mat["TX"]=-1.0;
  dist_mat["AS"]=1.0;dist_mat["AR"]=1.0;dist_mat["AL"]=1.0;dist_mat["AP"]=1.0;dist_mat["AT"]=1.0;dist_mat["AA"]=0.0;dist_mat["AV"]=1.0;dist_mat["AG"]=1.0;dist_mat["AI"]=1.0;dist_mat["AF"]=1.0;dist_mat["AY"]=1.0;dist_mat["AC"]=1.0;dist_mat["AH"]=1.0;dist_mat["AQ"]=1.0;dist_mat["AN"]=1.0;dist_mat["AK"]=1.0;dist_mat["AD"]=1.0;dist_mat["AE"]=1.0;dist_mat["AM"]=1.0;dist_mat["AW"]=1.0;dist_mat["A."]=-1.0;;dist_mat["A-"]=-1.0;;dist_mat["AX"]=-1.0;
  dist_mat["VS"]=1.0;dist_mat["VR"]=1.0;dist_mat["VL"]=1.0;dist_mat["VP"]=1.0;dist_mat["VT"]=1.0;dist_mat["VA"]=1.0;dist_mat["VV"]=0.0;dist_mat["VG"]=1.0;dist_mat["VI"]=1.0;dist_mat["VF"]=1.0;dist_mat["VY"]=1.0;dist_mat["VC"]=1.0;dist_mat["VH"]=1.0;dist_mat["VQ"]=1.0;dist_mat["VN"]=1.0;dist_mat["VK"]=1.0;dist_mat["VD"]=1.0;dist_mat["VE"]=1.0;dist_mat["VM"]=1.0;dist_mat["VW"]=1.0;dist_mat["V."]=-1.0;;dist_mat["V-"]=-1.0;;dist_mat["VX"]=-1.0;
  dist_mat["GS"]=1.0;dist_mat["GR"]=1.0;dist_mat["GL"]=1.0;dist_mat["GP"]=1.0;dist_mat["GT"]=1.0;dist_mat["GA"]=1.0;dist_mat["GV"]=1.0;dist_mat["GG"]=0.0;dist_mat["GI"]=1.0;dist_mat["GF"]=1.0;dist_mat["GY"]=1.0;dist_mat["GC"]=1.0;dist_mat["GH"]=1.0;dist_mat["GQ"]=1.0;dist_mat["GN"]=1.0;dist_mat["GK"]=1.0;dist_mat["GD"]=1.0;dist_mat["GE"]=1.0;dist_mat["GM"]=1.0;dist_mat["GW"]=1.0;dist_mat["G."]=-1.0;;dist_mat["G-"]=-1.0;;dist_mat["GX"]=-1.0;
  dist_mat["IS"]=1.0;dist_mat["IR"]=1.0;dist_mat["IL"]=1.0;dist_mat["IP"]=1.0;dist_mat["IT"]=1.0;dist_mat["IA"]=1.0;dist_mat["IV"]=1.0;dist_mat["IG"]=1.0;dist_mat["II"]=0.0;dist_mat["IF"]=1.0;dist_mat["IY"]=1.0;dist_mat["IC"]=1.0;dist_mat["IH"]=1.0;dist_mat["IQ"]=1.0;dist_mat["IN"]=1.0;dist_mat["IK"]=1.0;dist_mat["ID"]=1.0;dist_mat["IE"]=1.0;dist_mat["IM"]=1.0;dist_mat["IW"]=1.0;dist_mat["I."]=-1.0;;dist_mat["I-"]=-1.0;;dist_mat["IX"]=-1.0;
  dist_mat["FS"]=1.0;dist_mat["FR"]=1.0;dist_mat["FL"]=1.0;dist_mat["FP"]=1.0;dist_mat["FT"]=1.0;dist_mat["FA"]=1.0;dist_mat["FV"]=1.0;dist_mat["FG"]=1.0;dist_mat["FI"]=1.0;dist_mat["FF"]=0.0;dist_mat["FY"]=1.0;dist_mat["FC"]=1.0;dist_mat["FH"]=1.0;dist_mat["FQ"]=1.0;dist_mat["FN"]=1.0;dist_mat["FK"]=1.0;dist_mat["FD"]=1.0;dist_mat["FE"]=1.0;dist_mat["FM"]=1.0;dist_mat["FW"]=1.0;dist_mat["F."]=-1.0;;dist_mat["F-"]=-1.0;;dist_mat["FX"]=-1.0;
  dist_mat["YS"]=1.0;dist_mat["YR"]=1.0;dist_mat["YL"]=1.0;dist_mat["YP"]=1.0;dist_mat["YT"]=1.0;dist_mat["YA"]=1.0;dist_mat["YV"]=1.0;dist_mat["YS"]=1.0;dist_mat["YI"]=1.0;dist_mat["YF"]=1.0;dist_mat["YY"]=0.0;dist_mat["YC"]=1.0;dist_mat["YH"]=1.0;dist_mat["YQ"]=1.0;dist_mat["YN"]=1.0;dist_mat["YK"]=1.0;dist_mat["YD"]=1.0;dist_mat["YE"]=1.0;dist_mat["YM"]=1.0;dist_mat["YW"]=1.0;dist_mat["Y."]=-1.0;;dist_mat["Y-"]=-1.0;;dist_mat["YX"]=-1.0;
  dist_mat["CS"]=1.0;dist_mat["CR"]=1.0;dist_mat["CL"]=1.0;dist_mat["CP"]=1.0;dist_mat["CT"]=1.0;dist_mat["CA"]=1.0;dist_mat["CV"]=1.0;dist_mat["CG"]=1.0;dist_mat["CI"]=1.0;dist_mat["CF"]=1.0;dist_mat["CY"]=1.0;dist_mat["CC"]=0.0;dist_mat["CH"]=1.0;dist_mat["CQ"]=1.0;dist_mat["CN"]=1.0;dist_mat["CK"]=1.0;dist_mat["CD"]=1.0;dist_mat["CE"]=1.0;dist_mat["CM"]=1.0;dist_mat["CW"]=1.0;dist_mat["C."]=-1.0;;dist_mat["C-"]=-1.0;;dist_mat["CX"]=-1.0;
  dist_mat["HS"]=1.0;dist_mat["HR"]=1.0;dist_mat["HL"]=1.0;dist_mat["HP"]=1.0;dist_mat["HT"]=1.0;dist_mat["HA"]=1.0;dist_mat["HV"]=1.0;dist_mat["HG"]=1.0;dist_mat["HI"]=1.0;dist_mat["HF"]=1.0;dist_mat["HY"]=1.0;dist_mat["HC"]=1.0;dist_mat["HH"]=0.0;dist_mat["HQ"]=1.0;dist_mat["HN"]=1.0;dist_mat["HK"]=1.0;dist_mat["HD"]=1.0;dist_mat["HE"]=1.0;dist_mat["HM"]=1.0;dist_mat["HW"]=1.0;dist_mat["H."]=-1.0;;dist_mat["H-"]=-1.0;;dist_mat["HX"]=-1.0;
  dist_mat["QS"]=1.0;dist_mat["QR"]=1.0;dist_mat["QL"]=1.0;dist_mat["QP"]=1.0;dist_mat["QT"]=1.0;dist_mat["QA"]=1.0;dist_mat["QV"]=1.0;dist_mat["QG"]=1.0;dist_mat["QI"]=1.0;dist_mat["QF"]=1.0;dist_mat["QY"]=1.0;dist_mat["QC"]=1.0;dist_mat["QH"]=1.0;dist_mat["QQ"]=0.0;dist_mat["QN"]=1.0;dist_mat["QK"]=1.0;dist_mat["QD"]=1.0;dist_mat["QE"]=1.0;dist_mat["QM"]=1.0;dist_mat["QW"]=1.0;dist_mat["Q."]=-1.0;;dist_mat["Q-"]=-1.0;;dist_mat["QX"]=-1.0;
  dist_mat["NS"]=1.0;dist_mat["NR"]=1.0;dist_mat["NL"]=1.0;dist_mat["NP"]=1.0;dist_mat["NT"]=1.0;dist_mat["NA"]=1.0;dist_mat["NV"]=1.0;dist_mat["NG"]=1.0;dist_mat["NI"]=1.0;dist_mat["NF"]=1.0;dist_mat["NY"]=1.0;dist_mat["NC"]=1.0;dist_mat["NH"]=1.0;dist_mat["NQ"]=1.0;dist_mat["NN"]=0.0;dist_mat["NK"]=1.0;dist_mat["ND"]=1.0;dist_mat["NE"]=1.0;dist_mat["NM"]=1.0;dist_mat["NW"]=1.0;dist_mat["N."]=-1.0;;dist_mat["N-"]=-1.0;;dist_mat["NX"]=-1.0;
  dist_mat["KS"]=1.0;dist_mat["KR"]=1.0;dist_mat["KL"]=1.0;dist_mat["KP"]=1.0;dist_mat["KT"]=1.0;dist_mat["KA"]=1.0;dist_mat["KV"]=1.0;dist_mat["KG"]=1.0;dist_mat["KI"]=1.0;dist_mat["KF"]=1.0;dist_mat["KY"]=1.0;dist_mat["KC"]=1.0;dist_mat["KH"]=1.0;dist_mat["KQ"]=1.0;dist_mat["KN"]=1.0;dist_mat["KK"]=0.0;dist_mat["KD"]=1.0;dist_mat["KE"]=1.0;dist_mat["KM"]=1.0;dist_mat["KW"]=1.0;dist_mat["K."]=-1.0;;dist_mat["K-"]=-1.0;;dist_mat["KX"]=-1.0;
  dist_mat["DS"]=1.0;dist_mat["DR"]=1.0;dist_mat["DL"]=1.0;dist_mat["DP"]=1.0;dist_mat["DT"]=1.0;dist_mat["DA"]=1.0;dist_mat["DV"]=1.0;dist_mat["DG"]=1.0;dist_mat["DI"]=1.0;dist_mat["DF"]=1.0;dist_mat["DY"]=1.0;dist_mat["DC"]=1.0;dist_mat["DH"]=1.0;dist_mat["DQ"]=1.0;dist_mat["DN"]=1.0;dist_mat["DK"]=1.0;dist_mat["DD"]=0.0;dist_mat["DE"]=1.0;dist_mat["DM"]=1.0;dist_mat["DW"]=1.0;dist_mat["D."]=-1.0;;dist_mat["D-"]=-1.0;;dist_mat["DX"]=-1.0;
  dist_mat["ES"]=1.0;dist_mat["ER"]=1.0;dist_mat["EL"]=1.0;dist_mat["EP"]=1.0;dist_mat["ET"]=1.0;dist_mat["EA"]=1.0;dist_mat["EV"]=1.0;dist_mat["EG"]=1.0;dist_mat["EI"]=1.0;dist_mat["EF"]=1.0;dist_mat["EY"]=1.0;dist_mat["EC"]=1.0;dist_mat["EH"]=1.0;dist_mat["EQ"]=1.0;dist_mat["EN"]=1.0;dist_mat["EK"]=1.0;dist_mat["ED"]=1.0;dist_mat["EE"]=0.0;dist_mat["EM"]=1.0;dist_mat["EW"]=1.0;dist_mat["E."]=-1.0;;dist_mat["E-"]=-1.0;;dist_mat["EX"]=-1.0;
  dist_mat["MS"]=1.0;dist_mat["MR"]=1.0;dist_mat["ML"]=1.0;dist_mat["MP"]=1.0;dist_mat["MT"]=1.0;dist_mat["MA"]=1.0;dist_mat["MV"]=1.0;dist_mat["MG"]=1.0;dist_mat["MI"]=1.0;dist_mat["MF"]=1.0;dist_mat["MY"]=1.0;dist_mat["MC"]=1.0;dist_mat["MH"]=1.0;dist_mat["MQ"]=1.0;dist_mat["MN"]=1.0;dist_mat["MK"]=1.0;dist_mat["MD"]=1.0;dist_mat["ME"]=1.0;dist_mat["MM"]=0.0;dist_mat["MW"]=1.0;dist_mat["M."]=-1.0;;dist_mat["M-"]=-1.0;;dist_mat["MX"]=-1.0;
  dist_mat["WS"]=1.0;dist_mat["WR"]=1.0;dist_mat["WL"]=1.0;dist_mat["WP"]=1.0;dist_mat["WT"]=1.0;dist_mat["WA"]=1.0;dist_mat["WV"]=1.0;dist_mat["WG"]=1.0;dist_mat["WI"]=1.0;dist_mat["WF"]=1.0;dist_mat["WY"]=1.0;dist_mat["WC"]=1.0;dist_mat["WH"]=1.0;dist_mat["WQ"]=1.0;dist_mat["WN"]=1.0;dist_mat["WK"]=1.0;dist_mat["WD"]=1.0;dist_mat["WE"]=1.0;dist_mat["WM"]=1.0;dist_mat["WW"]=0.0;dist_mat["W."]=-1.0;;dist_mat["W-"]=-1.0;;dist_mat["WX"]=-1.0;
  dist_mat[".S"]=-1.0;dist_mat[".R"]=-1.0;dist_mat[".L"]=-1.0;dist_mat[".P"]=-1.0;dist_mat[".T"]=-1.0;dist_mat[".A"]=-1.0;dist_mat[".V"]=-1.0;dist_mat[".G"]=-1.0;dist_mat[".I"]=-1.0;dist_mat[".F"]=-1.0;dist_mat[".Y"]=-1.0;dist_mat[".C"]=-1.0;dist_mat[".H"]=-1.0;dist_mat[".Q"]=-1.0;dist_mat[".N"]=-1.0;dist_mat[".K"]=-1.0;dist_mat[".D"]=-1.0;dist_mat[".E"]=-1.0;dist_mat[".M"]=-1.0;dist_mat[".W"]=-1.0;dist_mat[".."]=-1.0;;dist_mat[".-"]=-1.0;;dist_mat[".X"]=-1.0;
  dist_mat["-S"]=-1.0;dist_mat["-R"]=-1.0;dist_mat["-L"]=-1.0;dist_mat["-P"]=-1.0;dist_mat["-T"]=-1.0;dist_mat["-A"]=-1.0;dist_mat["-V"]=-1.0;dist_mat["-G"]=-1.0;dist_mat["-I"]=-1.0;dist_mat["-F"]=-1.0;dist_mat["-Y"]=-1.0;dist_mat["-C"]=-1.0;dist_mat["-H"]=-1.0;dist_mat["-Q"]=-1.0;dist_mat["-N"]=-1.0;dist_mat["-K"]=-1.0;dist_mat["-D"]=-1.0;dist_mat["-E"]=-1.0;dist_mat["-M"]=-1.0;dist_mat["-W"]=-1.0;dist_mat["-."]=-1.0;;dist_mat["--"]=-1.0;;dist_mat["-X"]=-1.0;
  dist_mat["XS"]=-1.0;dist_mat["XR"]=-1.0;dist_mat["XL"]=-1.0;dist_mat["XP"]=-1.0;dist_mat["XT"]=-1.0;dist_mat["XA"]=-1.0;dist_mat["XV"]=-1.0;dist_mat["XG"]=-1.0;dist_mat["XI"]=-1.0;dist_mat["XF"]=-1.0;dist_mat["XY"]=-1.0;dist_mat["XC"]=-1.0;dist_mat["XH"]=-1.0;dist_mat["XQ"]=-1.0;dist_mat["XN"]=-1.0;dist_mat["XK"]=-1.0;dist_mat["XD"]=-1.0;dist_mat["XE"]=-1.0;dist_mat["XM"]=-1.0;dist_mat["XW"]=-1.0;dist_mat["X."]=-1.0;;dist_mat["X-"]=-1.0;;dist_mat["XX"]=-1.0;
  int n=aavector.size();
  Rcpp::NumericMatrix distMatrix(n, n);
  CharacterVector aavectornames=aavector.attr("names");
  colnames(distMatrix)=aavectornames;
  rownames(distMatrix)=aavectornames;
  Rcpp::NumericMatrix sitesMatrix(n, n);
  colnames(sitesMatrix)=aavectornames;
  rownames(sitesMatrix)=aavectornames;
  int nsites=aavector[1].size();
  RcppThread::ProgressBar bar(n, 1);
  RcppThread::parallelFor(0, n, [&] (int i) {
    for( int j=i; j < n; j++ ) {
      double eqnum=0;
      int ij_n=nsites;
      for( int s=0; s < nsites; s++) {
        std::string is;
        std::string js;
        is=aavector[i][s];
        js=aavector[j][s];
        double ij_dist;
        ij_dist=dist_mat[is+js];
        if(ij_dist >= 0.0) {
          eqnum=eqnum+ij_dist;
        }
        else {
          ij_n=ij_n-1;
        };
      }
      distMatrix(i,j)=eqnum/ij_n;
      distMatrix(j,i)=eqnum/ij_n;
      sitesMatrix(i,j)=ij_n;
      sitesMatrix(j,i)=ij_n;
    };
    bar++;
  }, ncores);
  return Rcpp::List::create(Rcpp::Named("sitesUsed")=sitesMatrix);
}
