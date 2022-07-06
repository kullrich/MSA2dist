/************************************************************
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.
* Filename: KaKs.h
* Abstract: Declaration of KAKS class including several methods.
* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Jan.21, 2005
* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (zhangyb@big.ac.cn)
* Date: Jun.1, 2009
* Modified Version: 2.0.2
* Modified Author: Kristian K Ullrich (ullrich@evolbio.mpg.de)
* Modified Date: July.01, 2022
*************************************************************/
#ifndef KAKS_H
#define KAKS_H
#include "base.h"
#include "NG86.h"
#include "LWL85.h"
#include "LPB93.h"
#include "GY94.h"
#include "YN00.h"
#include "MYN.h"
#include "MSMA.h"

/* KAKS class */
class KAKS: public Base {
public:	
	/* int tempt; //zhangyubin add for multi-lines data to only gamma */
	KAKS();
	~KAKS();
	/* Main function to call kaks's methods */
	/* bool Run(std::string input, std::string output, std::string sgc, std::string method, bool verbose); */
	bool Run(Rcpp::StringVector cdsstr, std::string sgc, std::string method, bool verbose);
	/* Read and Calculate seq, called in "Run" main function */
	/* bool ReadCalculateSeq(std::string filename); */
	bool ReadCalculateSeq(Rcpp::StringVector cdsstr);
	/* Initialize class, ready for running */
	int Initialize();
	/* Uninitialize class, for unloading */
	int Uninitialize();
	void getGCContent2(std::string str);
protected:		
	/* Use several methods to calculate ka/ks */
	bool calculateKaKs();
    /* NONE: an in-house algorithm in BGI, that is NG86 without correction */
	void start_NONE(float GAMMA);
	/* NG86 */
	void start_NG86(float GAMMA);
	/* LWL85 */
	void start_LWL85(float GAMMA);
	/* Modified LWL85 */
	void start_MLWL85(float GAMMA);
	/* LPB93 */
	void start_LPB93(float GAMMA);
	/* Modified LPB93 */
	void start_MLPB93(float GAMMA);
	/* GY94 */
	void start_GY94(float GAMMA);	
	/* YN00 */
	void start_YN00(float GAMMA);
	/* MYN */
	void start_MYN(float GAMMA);	
	/* Model Selection and Model Averaging */
	void start_MSMA(float GAMMA);
	/* Get GCC of entire sequences and of three codon positions */
	void getGCContent(std::string str);
	/* Check the sequence whether is valid or not */
	//bool checkValid(std::string name, std::string str1, std::string str2);
	bool checkValid(std::string name, std::string name1, std::string name2, std::string comp1, std::string comp2, std::string str1, std::string str2);
	/* Parse the input parameters */
	//bool parseParameter(int argc, const char* argv[]);
	/* bool parseParameter(std::string input_filename, std::string output_filename, std::string select_genetic_code, std::string select_method, bool verbose); */
	bool parseParameter(std::string select_genetic_code, std::string select_method, bool verbose);
	/* Get title information for writing into file */
	std::string getTitleInfo();
public:
	/* Methods' name and reference */
	std::vector<std::string> method_name;
	std::vector<std::string> method_ref;
	/* Parameters' title in outputting file */
	std::vector<std::string> titleInfo;
	/* Results for windows parser that shows results on ListCtrl */
	std::string result4Win;
	/* File name for output */
	std::string output_filename;
	/* Sequence file name */
	std::string seq_filename;
	/* Flag for whether to run NG86, MLWL85, MLPB93, GY94, YN00, MYN, MS/A=model selection/averaging */
	bool none, ng86, gng86, lwl85, glwl85, lpb93,glpb93, yn00, gyn00, mlwl85, gmlwl85, mlpb93, gmlpb93, gy94, myn, ms, ma, gmyn;
	/* verbose */
	bool VERBOSE;
	/* Number of compared pairwise sequences */
	unsigned long number;	//Maybe too many
	/* results data frame */
	Rcpp::DataFrame results_df;
protected:
	/* File name for detailed results for model selection */
	std::string detail_filename;
	/* Detailed results */
	std::string details;
private:
	/* The temporary results for write into file */
	std::string result;
	/* Output stream */
	std::ofstream os;
	/* A pair of sequence */
	std::string seq1, seq2;
}; 

#endif
