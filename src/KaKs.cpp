/************************************************************
* Copyright (C) 2009, Beijing Institute of Genomics of CAS
* All rights reserved.
* Filename: KaKs.cpp
* Abstract: Declaration of KAKS class including several methods.
* Version: 1.0
* Author: Zhang Zhang (zhanghzhang@genomics.org.cn)
* Date: Jan.21, 2005
* Version: 2.0
* Author: Da-Peng Wang(wangdp@big.ac.cn), Yu-Bin Zhang (ybzhang@big.ac.cn)
* Date: Jun.1, 2009
* Modified Version: 2.0.2
* Modified Author: Kristian K Ullrich (ullrich@evolbio.mpg.de)
* Modified Date: July.01, 2022
*************************************************************/
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include "KaKs.h"
#ifdef __APPLE__
#else
    #include <cstring>
#endif

KAKS::KAKS() {
	string items[] = {
	    "Comp1", "Comp2", "seq1", "seq2",
        "Method", "Ka", "Ks", "Ka/Ks",
	    "P-Value(Fisher)", "Length", "S-Sites", "N-Sites",
	    "Fold-Sites(0:2:4)", "Substitutions", "S-Substitutions", "N-Substitutions",
	    "Fold-S-Substitutions(0:2:4)", "Fold-N-Substitutions(0:2:4)", "Divergence-Time", "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)",
	    "GC(1:2:3)", "ML-Score", "AICc", "Akaike-Weight",
	    "Model"}; //25 outputting items
	//Format of outputting parameters
	int i;
	for (i=0; i<sizeof(items)/sizeof(string); i++) titleInfo.push_back(items[i]);
	//Load Methods' Names and References for -h in linux, also for windows' tool tip
	method_name.push_back("NG");
	method_name.push_back("LWL");
	method_name.push_back("LPB");
	method_name.push_back("MLWL");
	method_name.push_back("MLPB");
	method_name.push_back("GY");
	method_name.push_back("YN");
	method_name.push_back("MYN");
	method_name.push_back("MS");
	method_name.push_back("MA");
	method_name.push_back("GNG");
	method_name.push_back("GLWL");
	method_name.push_back("GLPB");
	method_name.push_back("GMLWL");
	method_name.push_back("GMLPB");
	method_name.push_back("GYN");
	method_name.push_back("GMYN");
  method_ref.push_back("Nei, M. and Gojobori, T. (1986) Mol. Biol. Evol., 3, 418-426.");
	method_ref.push_back("Li, W.H., Wu, C.I. and Luo, C.C. (1985) Mol. Biol. Evol., 2, 150-174.");
	method_ref.push_back("Li, W.H. (1993) J. Mol. Evol., 36, 96-99.    Pamilo, P. and Bianchi, N.O. (1993) Mol. Biol. Evol., 10, 271-281.");
	method_ref.push_back("Tzeng, Y.H., Pan, R. and Li, W.H. (2004) Mol. Biol. Evol., 21, 2290-2298.");
	method_ref.push_back("Tzeng, Y.H., Pan, R. and Li, W.H. (2004) Mol. Biol. Evol., 21, 2290-2298.");
	method_ref.push_back("Goldman, N. and Yang, Z. (1994) Mol. Biol. Evol., 11, 725-736.");
	method_ref.push_back("Yang, Z. and Nielsen, R. (2000) Mol. Biol. Evol., 17, 32-43.");
	method_ref.push_back("Zhang, Z., Li, J. and Yu, J. (2006) BMC Evolutionary Biology, 6, 44.");
	method_ref.push_back("Model Selection according to the AICc");
	method_ref.push_back("Model Averaging on a set of candidate models");
	method_ref.push_back("Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.");
	method_ref.push_back("Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.");
	method_ref.push_back("Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.");
	method_ref.push_back("Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.");
	method_ref.push_back("Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.");
	method_ref.push_back("Wang, DP., Zhang, S., He, FH., Zhu, J.,Hu, SN. and Yu, J.(2009) Genomics, Proteomics and Bioinformatics. In press.");
  method_ref.push_back("Wang, DP., Wan, HL., Zhang, S. and Yu, J. (2009) Biology Direct, 4:20 (16 June 2009)");
	Initialize();
}

KAKS::~KAKS() {
	Uninitialize();
	//Free memory
	titleInfo.clear();
	method_name.clear();
	method_ref.clear();
}

int KAKS::Initialize() {
    none=ng86=gng86=lpb93=glpb93=lwl85=glwl85=mlwl85=gmlwl85=mlpb93=gmlpb93=yn00=gyn00=gy94=myn=gmyn=ms=ma=false;
    VERBOSE=false;
	result4Win=result=details="";
	seq_name=seq1_name=seq2_name=seq1=seq2="";
	seq_filename=output_filename=detail_filename="";
	genetic_code=1;
	number=0;
	results_df=Rcpp::DataFrame::create(Named("rownames")=titleInfo);
	return 1;
}

int KAKS::Uninitialize() {
	if (os.is_open()) {
		os.close();
	}
	return 1;
}

/* Get GC of entire sequences GC[0] and of three codon positions GC[1,2,3] */
void KAKS::getGCContent(string str) {
	int i, j;
	initArray(GC, 4);
	for (i=0; i<str.length(); i+=3) {
		string codon=str.substr(i, 3);
		for (j=0; j<3; j++) {
			if (codon[j]=='G' || codon[j]=='C') GC[j+1]++;
		}
	}
	GC[0]=sumArray(GC, 4, 1) / str.length();
	for (i=1; i<4; i++) GC[i] /= (str.length()/3);
}

/****************************************************
* Function: ReadCalculateSeq
* Input Parameter: StringVector
* Output: Read sequences, check sequences' validity
		  and calculate Ka and Ks.
* Return Value: True if succeed, otherwise false.
*****************************************************/
bool KAKS::ReadCalculateSeq(Rcpp::StringVector cdsstr) {
	bool flag=true;
	try	{
	    int i, j;
	    int n=cdsstr.size();
	    Rcpp::StringVector cdsstrnames = cdsstr.attr("names");
		/* result=getTitleInfo();//zhangyubin revised */
		string temp="", name="", name1="", name2="", str="", str1="", str2="", comp1="", comp2="";
		for (i=0; i<n; i++){
		    for (j=i+1; j<n; j++ ){
		        comp1=std::to_string(i+1);
		        comp2=std::to_string(j+1);
		        name1=cdsstrnames[i];
		        name2=cdsstrnames[j];
		        name=name1+"&"+name2;
		        str1=cdsstr[i];
		        str2=cdsstr[j];
		        str=str1+str2;
		        //Get GCC at three codon positions
			    getGCContent(str);
			    //Check str's validility and calculate
			    if (checkValid(name, name1, name2, comp1, comp2, str.substr(0, str.length()/2), str.substr(str.length()/2, str.length()/2))) {
				    if (VERBOSE) {
				        Rcpp::Rcout<<++number<<" "<<name<<endl;
				    }
				    if (!calculateKaKs()) throw 1;	//calculate Ka&Ks using selected method(s)
				    if (VERBOSE) {
				        Rcpp::Rcout<<"[OK]"<<endl;
				    }
			    }
			    name=name1=name2=str=str1=str2=comp1=comp2="";
		    }
		}
	}
	catch (...) {
		flag=false;
	}
	return flag;
}

/**************************************************
* Function: checkValid
* Input Parameter: string, string, string
* Output: Check validity of pairwise sequences
* Return Value: True if succeed, otherwise false.
***************************************************/
bool KAKS::checkValid(string name, string name1, string name2, string comp1, string comp2, string str1, string str2) {
	bool flag=true;
	long i;
	try {
		//Check whether (sequence length)/3==0
		if (str1.length()!=str2.length() || str1.length()%3!=0 || str2.length()%3!=0) {
			Rcpp::Rcout<<endl<<"Error. The size of two sequences in "<<"'"<<name<<"' is not equal."<<endl;
			throw 1;
		}
		//Delete gap and stop codon
		bool found;
		int j;
		for (i=0; i<str1.length(); i+=3) {
			for (found=false, j=0; j<3 && !found; j++) {
				if (str1[j+i]=='-' || str2[j+i]=='-') {
					found=true;
				}
				str1[i+j]=toupper(str1[i+j]);
				str2[i+j]=toupper(str2[i+j]);
				if (convertChar(str1[i+j])==-1 || convertChar(str2[i+j])==-1) {
					found=true;
				}
			}
			if ((getAminoAcid(str1.substr(i,3))=='!') || (getAminoAcid(str2.substr(i,3))=='!')) {
				found=true;
			}
			if (found) {
				str1=str1.replace(i, 3, "");
				str2=str2.replace(i, 3, "");
				i-=3;
			}
		}
		//pass value into private variables
		seq1=str1;
		seq2=str2;
		//pass value into extern variables
		seq_name=name;
		seq1_name=name1;
		seq2_name=name2;
		comp1_name=comp1;
		comp2_name=comp2;
		length=str1.length();
	}
	catch (...) {
		flag=false;
	}
	return flag;
}

/**************************************************
* Function: Run
* Input Parameter: int, const char* []
* Output: Calculate Ka/Ks and output.
* Return Value: void
***************************************************/
bool KAKS::Run(StringVector cdsstr, string sgc, string method, bool verbose) {
	bool flag=true;
	try {
		//Judge whether input parameters are legal
		if (!parseParameter(sgc, method, verbose)) {
			throw 1;
		}
		//Record the start time
		static time_t time_start=time(NULL);
		//Read sequences and calculate Ka & Ks
		if (!ReadCalculateSeq(cdsstr)) {
			throw 1;
		}
		//Time used for running
		time_t t=time(NULL)-time_start;
		int h=t/3600, m=(t%3600)/60, s=t-(t/60)*60;
		//Print on display
		if (VERBOSE) {
		    Rcpp::Rcout<<"Mission accomplished. (Time elapsed: ";
		    if (h) {
		        Rcpp::Rcout<<h<<":"<<m<<":"<<s<<")"<<endl;
		    } else {
		        Rcpp::Rcout<<m<<":"<<s<<")"<<endl;
		    }
		}
	}
	catch (...) {
		flag=false;
	}
	return flag;
}

/**************************************************
* Function: parseParameter
* Input Parameter: int, const char* []
* Output: Parse the input parameters
* Return Value: bool
***************************************************/
bool KAKS::parseParameter(string sgc, string method, bool verbose) {
	bool flag=true;
	string temp;
	VERBOSE=verbose;
	try {
		genetic_code=CONVERT<int>(sgc);
        temp=stringtoUpper(method);
        if (temp=="NONE") none=true;
		else if (temp=="NG") ng86=true;
		else if (temp=="GNG") gng86=true;//added by zhangyubin
		else if (temp=="LWL") lwl85=true;
		else if (temp=="GLWL") glwl85=true;
		else if (temp=="LPB") lpb93=true;
		else if (temp=="GLPB") glpb93=true;
		else if (temp=="MLPB") mlpb93=true;
		else if (temp=="GMLPB") gmlpb93=true;
		else if (temp=="MLWL") mlwl85=true;
		else if (temp=="GMLWL") gmlwl85=true;
		else if (temp=="GY") gy94=true;
		else if (temp=="YN") yn00=true;
		else if (temp=="GYN") gyn00=true;
		else if (temp=="MYN") myn=true;
		else if (temp=="GMYN") gmyn=true;
		else if (temp=="MS") ms=true;
		else if (temp=="MA") ma=true;
		else {
		    Rcpp::Rcout<<"No valid method."<<endl;
		    Rcpp::stop("No valid method.");
		}
		/**
		else if (temp=="ALL") {
			ng86=gng86=lpb93=glpb93=lwl85=glwl85=mlwl85=gmlwl85=mlpb93=gmlpb93=gy94=yn00=gyn00=myn=gmyn=ms=ma=true;  //zhangyubin  added
		}
		**/
	}
	catch (...) {
		flag=false;
	}
	return flag;
}

/*******************************************************
* Function: Calculate
* Input Parameter: void
* Output: Calculate kaks and output results.
* Return Value: void
*
* Note:
********************************************************/
int temp=0;
bool KAKS::calculateKaKs() {
	bool flag=true;
	try {
		//Estimate Ka and Ks
		if (none) start_NONE(0);
		if (ng86||gng86) start_NG86(0);
		if (lwl85||glwl85) start_LWL85(0);
		if (mlwl85||gmlwl85) start_MLWL85(0);
		if (lpb93||glpb93) start_LPB93(0);
		if (mlpb93||gmlpb93) start_MLPB93(0);
		if (gy94) start_GY94(0);
		if (yn00||gyn00) start_YN00(0);
		if (myn||gmyn) start_MYN(0);
		if (ms||ma) start_MSMA(0);
//2009 June added by ZhangYubin for Gamma  started
		if (gmyn||gng86||gyn00||glwl85||gmlwl85||glpb93||gmlpb93) {
 	        float fka,fks,fkaks;
		    string tresult;
		    tresult=result;
		    if (temp==0) {
                tresult=tresult.replace(0, tresult.find('\n')+1, "");
            }
		    temp++;
	        int j,k;
		    j=tresult.find('\n');
		    string linecontent=tresult.substr(0, j+1);
		    linecontent[linecontent.length()-1]='\t';
		    tresult=tresult.replace(0, j+1, "");
		    k=1;
		    while((j=linecontent.find('\t'))>0) {
	            string temp=linecontent.substr(0,j).c_str();
	            if (k==3) {
		            fka=atof(temp.c_str());
	            }
	            if (k==4) {
		            fks=atof(temp.c_str());
	            }
	            if (k==5) {
		            fkaks=atof(temp.c_str());
		            //goto att;
	            }
			    linecontent=linecontent.replace(0, j+1, "");
			    k++;
		    }
            if ((!myn)&&(!ng86)&&(!lwl85)&&(!mlwl85)&&(!lpb93)&&(!mlpb93)&&(!yn00)) {
                //||(gmlpb93&&myn))//&&(!gy94)&&(!ms)&&(!ma))
	            /**
	            if (tempt==0) {
                    result=getTitleInfo();//for Gamma choice??s result delete the first line
                    //Linux version don't need this line
                    tempt++;
                } else {
                    result="";//for Gamma choice??s result delete the first line
                }
	            **/
                result="";
                if (gy94) start_GY94(-1);
                if (ms||ma) start_MSMA(-1);
            }
            if (fkaks<1) {
                //result="";
	            //Estimate Ka and Ks
	            if (none) start_NONE(-1);
	            if (gng86) start_NG86(-1);
	            if (glwl85) start_LWL85(-1);
	            if (gmlwl85) start_MLWL85(4);
	            if (glpb93) start_LPB93(1);
	            if (gmlpb93) start_MLPB93(1);
                //if (gy94) start_GY94(0);
	            if (gyn00) start_YN00(4);
	            if (gmyn) start_MYN(20);
                // if (ms||ma) start_MSMA(0);
            } else if (fkaks>1) {
                //result="";
	            //Estimate Ka and Ks
	            if (none) start_NONE(-1);
	            if (gng86) start_NG86(6);
	            if (glwl85) start_LWL85(0.2);
	            if (gmlwl85) start_MLWL85(0.6);
	            if (glpb93) start_LPB93(1);
	            if (gmlpb93) start_MLPB93(1);
                //if (gy94)	start_GY94(0);
	            if (gyn00) start_YN00(-1);
	            if (gmyn) start_MYN(-1);
                //if (ms||ma) start_MSMA(0);
            } else if (fkaks==1) {
	            if (none) start_NONE(-1);
	            if (gng86) start_NG86(-1);
	            if (glwl85) start_LWL85(-1);
	            if (gmlwl85) start_MLWL85(-1);
	            if (glpb93) start_LPB93(-1);
	            if (gmlpb93) start_MLPB93(-1);
                //if (gy94) start_GY94(-1);
	            if (gyn00) start_YN00(-1);
	            if (gmyn) start_MYN(-1);
                //if (ms||ma) start_MSMA(-1);
            }
        }
        //2009 June added by ZhangYubin for Gamma end
	    result4Win+=result;
		if (VERBOSE) {
		    Rcpp::Rcout<<"result:"<<endl;
		    Rcpp::Rcout<<result<<endl;
		}
		//Add results to data frame
		Rcpp::StringVector tsplit;
		tsplit=splitString(result, "\t");
		if (VERBOSE) {
		    Rcpp::Rcout<<tsplit<<endl;
		}
		string colname_tmp=seq1_name+'_'+seq2_name;
		results_df.push_back(tsplit, colname_tmp);
		result="";
	}
	catch (...) {
		flag=false;
	}
	return flag;
}

//NONE: NG without correction for multiple substitution
void KAKS::start_NONE(float GAMMA) {
	NONE zz;
	zz.GAMMA=GAMMA;
	result+=zz.Run(seq1, seq2);
}

//NG
void KAKS::start_NG86(float GAMMA) {
	NG86 zz;
	zz.GAMMA=GAMMA;
	result+=zz.Run(seq1, seq2);
}

//LWL
void KAKS::start_LWL85(float GAMMA) {
	LWL85 zz;
	zz.GAMMA=GAMMA;
	result+=zz.Run(seq1, seq2);
}

//MLWL
void KAKS::start_MLWL85(float GAMMA) {
	MLWL85 zz;
	zz.GAMMA=GAMMA;
	result+=zz.Run(seq1, seq2);
}

//LPB
void KAKS::start_LPB93(float GAMMA) {
	LPB93 zz;
	zz.GAMMA=GAMMA;
	result+=zz.Run(seq1, seq2);
}

//MLPB
void KAKS::start_MLPB93(float GAMMA) {
	MLPB93 zz;
	zz.GAMMA=GAMMA;
	result+=zz.Run(seq1, seq2);
}

//GY
void KAKS::start_GY94(float GAMMA) {
	GY94 zz("HKY");
	zz.GAMMA=GAMMA;
	result+=zz.Run(seq1.c_str(), seq2.c_str());
}

//YN
void KAKS::start_YN00(float GAMMA) {
	YN00 zz;
	zz.GAMMA=GAMMA;
	result+=zz.Run(seq1, seq2);
}

//MYN
void KAKS::start_MYN(float GAMMA) {
	MYN zz;
	zz.GAMMA=GAMMA;
	result+=zz.Run(seq1, seq2);
}

/************************************************
* Function: start_MSMA
* Input Parameter: void
* Output: Calculate Ka and Ks using the method of
		  model selection or model averaging.
* Return Value: void
*************************************************/
void KAKS::start_MSMA(float GAMMA) {
	vector<MLResult> result4MA;	//generated by MS and used by MA
	//Model Selection
	MS zz1;
	zz1.GAMMA=GAMMA;
	string tmp=zz1.Run(seq1.c_str(), seq2.c_str(), result4MA, details);
	if (ms) {
		result+=tmp;
	}
	//Model Averaging
	if (ma) {
		MA zz2;
		zz2.GAMMA=GAMMA;
		result+=zz2.Run(seq1.c_str(), seq2.c_str(), result4MA);
	}
}

/************************************************
* Function: getTitleInfo
* Input Parameter:
* Output: get title information of outputting file
* Return Value: string
*************************************************/
string KAKS::getTitleInfo() {
	string title="";
	int i=0;
	if (titleInfo.size()>0) {
		//Add "\t" to items except the last one
		for (; i<titleInfo.size()-1; i++) addString(title, titleInfo[i]);
		//The last item is added by "\n"
		addString(title, titleInfo[i], "\n");
	}
	return title;
}
