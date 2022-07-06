/************************************************************
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.
* Filename: LPB93.h
* Abstract: Declaration of LPB93 and Modified LPB93 class.
* Version: 1.0
* Author: Zhang Zhang  (zhanghzhang@genomics.org.cn)
* Date: Jan.21, 2005
* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009
* Modified Version: 2.0.2
* Modified Author: Kristian K Ullrich (ullrich@evolbio.mpg.de)
* Modified Date: July.01, 2022
  References: 
  Li WH  (1993)  Unbiased estimation of the Rates of synonymous
  and nonsynonymous substitution. J. Mol. Evol. 36:96-99.
  Pamilo P, Bianchi NO  (1993)  Evolution of the Zfx and Zfy 
  genes: rates and interdependence between the genes. Mol. Biol.
  Evol. 10:271-281.
  Tzeng Y-H, Pan R, Li W-H  (2004)  Comparison of Three Methods
  for Estimating Rates of Synonymous and Nonsynonymous Nucleotide
  Substitutions. Mol. Biol. Evol. 21:2290-2298.
*************************************************************/
#ifndef LPB93_H
#define LPB93_H
#include "base.h"
#include "LWL85.h"	//derived from LWL85

/* LPB93 class */
class LPB93: public LWL85 {
public:
	LPB93();
	/* Main function of calculating kaks */
	std::string Run(std::string seq1, std::string seq2);
}; 

class MLPB93: public LPB93 {
public:
	MLPB93();
protected:
	/* Calculate the transition & transversion between two codons at a given position*/
	int TransitionTransversion(std::string codon1, std::string codon2, int pos);
}; 

#endif
