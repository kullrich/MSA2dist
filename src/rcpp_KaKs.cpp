/***************************************************************
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.
* Filename: KaKs_Calculator.cpp
* Abstract: including maximum-likelihood and approximate methods.
* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Jan.21, 2005
* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (zhangyb@big.ac.cn)
* Date: Jun.1, 2009
* Modified Version: 2.0.2
* Modified Author: Kristian K Ullrich (ullrich@evolbio.mpg.de)
* Modified Date: July.01, 2022
****************************************************************/
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace std;
#include "KaKs.h"

//' @useDynLib MSA2dist, .registration = TRUE
//' @import Rcpp
//' @title rcpp_KaKs
//' @name rcpp_KaKs
//' @description calculates KaKs as implememted in
//' KaKs Calculator 2.0 \code{MSA2dist} with \code{Rcpp}.
//' @return \code{data.frame}
//' @param cdsstr StringVector [mandatory]
//' @param sgc standard genetic code to use [default: 1]
//' @param method KaKs Calculator 2.0 codon model [default: YN]
//' @param verbose specify if verbose output [default: FALSE]
//' @references Wang et al. (2010) KaKs_Calculator 2.0: a toolkit incorporating
//' gamma-series methods and sliding window strategies.\emph{Genomics,
//' proteomics & bioinformatics.} \bold{8(1)}, 77-80.
//' @examples
//' ## load example sequence data
//' data("hiv", package="MSA2dist")
//' rcpp_KaKs(cdsstr=as.character(hiv[1:3]))
//' @export rcpp_KaKs
//' @author Kristian K Ullrich
// [[Rcpp::export]]
Rcpp::DataFrame rcpp_KaKs( Rcpp::StringVector cdsstr,
    const std::string sgc="1",
    const std::string method="YN",
    bool verbose=false ) {
    Rcpp::DataFrame results_df;
	try {
		KAKS kk;
		if (!kk.Run(cdsstr, sgc, method, verbose)) {
      throw 1;
    }
		results_df=kk.results_df;
	}
	catch (...) {
	}
	return results_df;
}
