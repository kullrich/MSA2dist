/***************************************************************
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.
* Filename: NG86.h
* Abstract: Declaration of NG86 class inherited base class.
* Version: 1.0
* Author: Zhang Zhang  (zhang.zhang@yale.edu)
* Date: Feb.21, 2005
* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009
* Modified Version: 2.0.2
* Modified Author: Kristian K Ullrich (ullrich@evolbio.mpg.de)
* Modified Date: July.01, 2022
  References:
  Nei M, Gojobori T  (1986)  Simple methods for estimating the 
  numbers of synonymous and nonsynonymous nucleotide substitutions.
  Mol Biol Evol 3:418-426.
****************************************************************/
#ifndef NG86_H
#define NG86_H
#include "base.h"
 
/* NG86 class */
class NG86: public Base {
public:
	float GAMMA;//zhangyubin added
	NG86();
	/* Main function of calculating kaks */
	std::string Run(std::string seq1, std::string seq2);
protected:
	/* Count codon's sites */
	void getCondonSite(std::string codon);
	/* Count codon's differences */
	void getCondonDifference(std::string codon1, std::string codon2);
	/* Preprocess */
	void PreProcess(std::string seq1, std::string seq2);
	/* Jukes and Cantor's one-parameter formula */
	double kaks_formula(double p);
public:	
	/* Proportions of sysnonymous(Ps) and nonsysnonymous(Pn): Ps=Sd/S, Pn=Nd/N  */
	double Ps, Pn;
}; 

class NONE: public NG86 {
public:
	NONE();
	/* Main function of calculating kaks */
	std::string Run(std::string seq1, std::string seq2);
};

#endif
