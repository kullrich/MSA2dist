/************************************************************
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.
* Filename: MYN.h
* Abstract: Declaration of Modified YN00 (MYN) class.
* Version: 1.0
* Author: Zhang Zhang(zhanghzhang@genomics.org.cn)
* Date: Dec.30, 2005
* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009
* Modified Version: 2.0.2
* Modified Author: Kristian K Ullrich (ullrich@evolbio.mpg.de)
* Modified Date: July.01, 2022
  References: 
	Zhang Zhang, Jun Li, Jun Yu. (2006) Computing Ka and Ks 
	with a consideration of unequal transitional substitutions. 
	BMC Evolutionary Biology, 6:44.
*************************************************************/
#ifndef MYN_H
#define MYN_H
#define CODONFREQ 12
#define square(a) ((a)*(a))
#include "base.h"
#include "YN00.h"

class MYN: public YN00 {
public:
	double GAMMA;
	MYN();
	/* Get the two kappas between purines and between pyrimidines */
	virtual int GetKappa(const std::string seq1, const std::string seq2);	
	/* Calculate the transition probability matrix  */
	int GetPMatCodon(double P[], double kappa, double omega);
	/* Count S and N */
	int CountSites(const std::string z, double &Stot,double &Ntot,double fbS[],double fbN[]);
	/* Correct for multiple substitutions for two kappas */
	int CorrectKappaTN93(double n, double P1, double P2, double Q, double pi4[], double &kappatc_TN93, double &kappaag_TN93);	
	/* Correct for multiple substitutions for Ka and Ks */
	int CorrectKaksTN93(double n, double P1, double P2, double Q, double pi4[], double &kaks, double &SEkaks);
	/* Count Sd and Nd */
	int MYNCountDiffs(const std::string seq1, const std::string seq2, double &Sdts1, double &Sdts2, double &Sdtv,double &Ndts1, double &Ndts2, double &Ndtv,double PMatrix[]);
	/* Main function */
	virtual int DistanceYN00(const std::string seq1, const std::string seq2, double &dS,double &dN, double &SEdS, double &SEdN); 
};

#endif
