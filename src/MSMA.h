/************************************************************
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.
* Filename: MSMA.h
* Abstract: Declaration of model-selected and model-averged methods' classes.
* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Apr. 2006
* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009
* Modified Version: 2.0.2
* Modified Author: Kristian K Ullrich (ullrich@evolbio.mpg.de)
* Modified Date: July.01, 2022
  References: 
  Goldman N, Yang Z  (1994)  A codon-based model of nucleotide
  substitution for protein-coding DNA sequences. Mol. Biol. 
  Evol. 11:725-736.
  Posada, D. and Buckley, T.R. (2004) Model Selection and Model Averaging
  in Phylogenetics: Advantages of Akaike Information Criterion and Bayesian
  Approaches over Likelihood Ratio Tests, Syst. Biol., 53, 793-808.
  Sullivan, J. and Joyce, P. (2005) Model Selection in Phylogenetics, 
  Annual Review of Ecology, Evolution, and Systematics, 36, 445-466.  
*************************************************************/
#ifndef MA_H
#define MA_H
#include "GY94.h"

class MS: public Base {
public:	
		float GAMMA;  //zhangyubin added
	MS();
	/* Main function */
	std::string Run(const char *seq1, const char *seq2, std::vector<MLResult> &result4MA, std::string &details);
protected:
	/* Calculate Ka and Ks based on a given model */
	void selectModel(const char *seq1, const char *seq2, std::string candidate_model, std::vector<MLResult> &result4MA);
};

class MA: public GY94 {
public:	
	MA();
	/* Main function */
	std::string Run(const char *seq1, const char *seq2, std::vector<MLResult> result4MA);
};

#endif
